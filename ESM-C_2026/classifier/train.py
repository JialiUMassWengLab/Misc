"""Training loop for LOF vs GOF variant classifier.

Assumes embedding_diff_lofgof.pt has been produced by classifier.merge_lofgof.
Logs per-epoch train_loss, val_loss, val_auroc to runs/<run_name>/training_log.csv.
Saves best-val-AUROC checkpoint to runs/<run_name>/best.pt.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import torch
from torch.utils.data import DataLoader
from sklearn.metrics import roc_auc_score

from classifier.dataset import (
    SiamesePairDataset,
    VariantEmbeddingDataset,
    VariantPositionDataset,
    VariantPrecomputedPositionDataset,
    collate_fn,
    collate_fn_position,
    collate_fn_siamese,
    load_lofgof_embeddings,
    load_wt_variant_embeddings,
)
from classifier.model import AttentionPoolClassifier
from classifier.model_aa_pos_only import PositionOnlyClassifier
from classifier.model_siamese import SiameseAttentionClassifier


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_METADATA = PROJECT_ROOT / "lof_gof_variant_metadata.tsv"
DEFAULT_CACHE = PROJECT_ROOT / "embedding_diff_lofgof.pt"
DEFAULT_WT_CACHE = PROJECT_ROOT / "wt_embeddings.pt"
DEFAULT_VAR_CACHE = PROJECT_ROOT / "variant_embeddings.pt"


def _parse_hidden_dims(spec: str) -> tuple[int, ...]:
    """Parse a comma-separated int list, e.g. '1024,256' -> (1024, 256)."""
    parts = [p.strip() for p in spec.split(",") if p.strip()]
    if not parts:
        raise argparse.ArgumentTypeError("--hidden-dims must be non-empty")
    try:
        dims = tuple(int(p) for p in parts)
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"--hidden-dims must be ints: {spec!r}") from e
    if any(d <= 0 for d in dims):
        raise argparse.ArgumentTypeError(f"--hidden-dims must be positive: {dims}")
    return dims


def run_epoch(
    model: torch.nn.Module,
    loader: DataLoader,
    criterion: torch.nn.Module,
    optimizer: torch.optim.Optimizer | None,
    device: torch.device,
) -> tuple[float, float | None]:
    is_train = optimizer is not None
    model.train(is_train)

    total_loss = 0.0
    total_n = 0
    all_probs: list[float] = []
    all_labels: list[float] = []

    context = torch.enable_grad() if is_train else torch.no_grad()
    with context:
        for batch in loader:
            labels = batch["label"].to(device, non_blocking=True)

            if "wt_embedding" in batch:
                # Siamese path: pass WT and variant tensors + shared mask.
                wt = batch["wt_embedding"].to(device, non_blocking=True)
                var = batch["var_embedding"].to(device, non_blocking=True)
                mask = batch["attention_mask"].to(device, non_blocking=True)
                logits = model(wt, mask, var, mask).squeeze(-1)
            else:
                # Single-input path (attention / aa_pos_only).
                emb = batch["embedding"].to(device, non_blocking=True)
                raw_mask = batch.get("attention_mask")
                mask = raw_mask.to(device, non_blocking=True) if raw_mask is not None else None
                logits = model(emb, mask).squeeze(-1)

            loss = criterion(logits, labels)

            if is_train:
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

            bs = labels.shape[0]
            total_loss += loss.item() * bs
            total_n += bs

            if not is_train:
                all_probs.extend(torch.sigmoid(logits).detach().cpu().tolist())
                all_labels.extend(labels.detach().cpu().tolist())

    mean_loss = total_loss / max(total_n, 1)
    auroc: float | None = None
    if not is_train and len(set(all_labels)) == 2:
        auroc = float(roc_auc_score(all_labels, all_probs))
    return mean_loss, auroc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE,
                        help="Per-residue or pre-pooled diff cache (used by "
                             "--model attention and --model aa_pos_only).")
    parser.add_argument("--wt-cache", type=Path, default=DEFAULT_WT_CACHE,
                        help="Per-residue WT embedding cache (used by --model siamese).")
    parser.add_argument("--variant-cache", type=Path, default=DEFAULT_VAR_CACHE,
                        help="Per-residue variant embedding cache (used by --model siamese).")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--hidden-dims", type=_parse_hidden_dims, default=(512, 128),
                        help="Comma-separated MLP hidden layer widths, e.g. '1024,256'. "
                             "Applied to both AttentionPoolClassifier and PositionOnlyClassifier. "
                             "Default: 512,128.")
    parser.add_argument("--pool-type", type=str, default="projected",
                        choices=["projected", "direct"],
                        help="Attention pooling variant (only used when --model attention). "
                             "'projected': Linear(x) @ q / sqrt(D). "
                             "'direct': x @ w (learnable weight vector, no projection).")
    parser.add_argument("--model", type=str, default="attention",
                        choices=["attention", "aa_pos_only", "siamese"],
                        help="'attention': attention-pool over the (WT - variant) diff cache. "
                             "'aa_pos_only': MLP on the substitution-position vector alone. "
                             "'siamese': shared-attention pool applied to WT and variant "
                             "separately, then MLP on the difference of the pooled vectors.")
    parser.add_argument("--run-name", type=str, default="default")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    torch.manual_seed(args.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device: {device}")

    if args.model == "siamese":
        wt_embeddings, var_embeddings, rows = load_wt_variant_embeddings(
            args.wt_cache, args.variant_cache, args.metadata
        )
        sample = next(iter(var_embeddings.values()))
        embed_dim = int(sample.shape[-1])
        print(f"siamese cache: per-residue, embed_dim={embed_dim}")

        collate = collate_fn_siamese
        model = SiameseAttentionClassifier(
            embed_dim=embed_dim, hidden_dims=args.hidden_dims,
            dropout=args.dropout, pool_type=args.pool_type,
        ).to(device)
        print(f"model: siamese (pool_type={args.pool_type}, hidden_dims={args.hidden_dims})")

        train_ds = SiamesePairDataset(wt_embeddings, var_embeddings, rows,
                                      split="train", split_seed=args.seed)
        val_ds = SiamesePairDataset(wt_embeddings, var_embeddings, rows,
                                    split="val", split_seed=args.seed)
    else:
        embeddings, rows = load_lofgof_embeddings(args.cache, args.metadata)

        # Peek at one entry to decide whether this cache is per-residue (2-D
        # (seq_len+2, D)) or pre-pooled at the substitution position (1-D (D,)).
        sample_tensor = next(iter(embeddings.values()))
        is_precomputed_position = sample_tensor.ndim == 1
        embed_dim = int(sample_tensor.shape[-1])
        print(f"cache format: {'pre-pooled position (1-D)' if is_precomputed_position else 'per-residue (2-D)'}, "
              f"embed_dim={embed_dim}")

        if args.model == "attention":
            if is_precomputed_position:
                raise ValueError(
                    "--model attention requires a per-residue cache "
                    "(e.g. embedding_diff_lofgof.pt); got a pre-pooled 1-D cache."
                )
            DatasetCls = VariantEmbeddingDataset
            collate = collate_fn
            model = AttentionPoolClassifier(embed_dim=embed_dim, hidden_dims=args.hidden_dims,
                                            dropout=args.dropout, pool_type=args.pool_type).to(device)
            print(f"model: attention (pool_type={args.pool_type}, hidden_dims={args.hidden_dims})")
        else:  # aa_pos_only
            collate = collate_fn_position
            if is_precomputed_position:
                DatasetCls = VariantPrecomputedPositionDataset
            else:
                DatasetCls = VariantPositionDataset
            model = PositionOnlyClassifier(embed_dim=embed_dim, hidden_dims=args.hidden_dims,
                                           dropout=args.dropout).to(device)
            print(f"model: aa_pos_only (embed_dim={embed_dim}, hidden_dims={args.hidden_dims})")

        train_ds = DatasetCls(embeddings, rows, split="train", split_seed=args.seed)
        val_ds = DatasetCls(embeddings, rows, split="val", split_seed=args.seed)

    pin = device.type == "cuda"
    # num_workers=0: workers would fork-copy the 21.6 GB shared dict; single-process
    # is fine since __getitem__ is now O(us).
    train_loader = DataLoader(
        train_ds, batch_size=args.batch_size, shuffle=True,
        collate_fn=collate, num_workers=0, pin_memory=pin,
    )
    val_loader = DataLoader(
        val_ds, batch_size=args.batch_size, shuffle=False,
        collate_fn=collate, num_workers=0, pin_memory=pin,
    )

    pos_weight = torch.tensor([train_ds.n_lof / max(train_ds.n_gof, 1)], device=device)
    print(f"pos_weight (n_lof/n_gof) = {pos_weight.item():.3f}")
    criterion = torch.nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    run_dir = Path(__file__).resolve().parent / "runs" / args.run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "training_log.csv"
    auroc_ckpt_path = run_dir / "best_auroc.pt"
    loss_ckpt_path = run_dir / "best_loss.pt"

    best_auroc = -1.0
    best_loss = float("inf")
    with log_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["epoch", "train_loss", "val_loss", "val_auroc"])
        for epoch in range(1, args.epochs + 1):
            train_loss, _ = run_epoch(model, train_loader, criterion, optimizer, device)
            val_loss, val_auroc = run_epoch(model, val_loader, criterion, None, device)
            writer.writerow([epoch, f"{train_loss:.6f}", f"{val_loss:.6f}",
                             f"{val_auroc:.6f}" if val_auroc is not None else ""])
            f.flush()
            auroc_str = f"{val_auroc:.4f}" if val_auroc is not None else "n/a"
            print(f"epoch {epoch:3d} | train_loss {train_loss:.4f} | val_loss {val_loss:.4f} | val_auroc {auroc_str}")

            if val_auroc is not None and val_auroc > best_auroc:
                best_auroc = val_auroc
                torch.save({"epoch": epoch, "model_state": model.state_dict(),
                            "val_auroc": val_auroc, "val_loss": val_loss}, auroc_ckpt_path)
            if val_loss < best_loss:
                best_loss = val_loss
                torch.save({"epoch": epoch, "model_state": model.state_dict(),
                            "val_auroc": val_auroc, "val_loss": val_loss}, loss_ckpt_path)

    print(f"done. best val_auroc = {best_auroc:.4f} -> {auroc_ckpt_path.name}")
    print(f"      best val_loss  = {best_loss:.4f} -> {loss_ckpt_path.name}")
    print(f"      log at {log_path}")


if __name__ == "__main__":
    main()
