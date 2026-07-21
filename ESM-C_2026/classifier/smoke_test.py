"""Smoke test: load cache -> dataset -> collate -> model forward -> backward."""

from pathlib import Path

import torch
from torch.utils.data import DataLoader

from classifier.dataset import VariantEmbeddingDataset, collate_fn, load_lofgof_embeddings
from classifier.model import AttentionPoolClassifier


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    metadata = PROJECT_ROOT / "combined_variant_metadata.tsv"
    cache = PROJECT_ROOT / "embedding_diff_lofgof.pt"

    embeddings, rows = load_lofgof_embeddings(cache, metadata)
    print(f"embeddings dict: {len(embeddings)} entries")
    print(f"rows: {len(rows)}")

    ds = VariantEmbeddingDataset(embeddings, rows, split="train")
    print(f"train dataset size: {len(ds)}")
    item = ds[0]
    print(f"item[0]: shape={tuple(item['embedding'].shape)}, "
          f"label={item['label'].item()}, name={item['variant_name']!r}")
    assert item["embedding"].ndim == 2
    assert item["embedding"].shape[1] == 1152
    assert item["label"].item() in (0.0, 1.0)

    loader = DataLoader(ds, batch_size=4, shuffle=False, collate_fn=collate_fn, num_workers=0)
    batch = next(iter(loader))
    print(f"batch: emb={tuple(batch['embedding'].shape)}, "
          f"mask={tuple(batch['attention_mask'].shape)} dtype={batch['attention_mask'].dtype}, "
          f"labels={batch['label'].tolist()}")
    assert batch["embedding"].shape[0] == 4
    assert batch["embedding"].shape[2] == 1152
    assert batch["attention_mask"].dtype == torch.bool

    model = AttentionPoolClassifier()
    logits = model(batch["embedding"], batch["attention_mask"])
    print(f"logits: shape={tuple(logits.shape)}, values={logits.detach().flatten().tolist()}")
    assert logits.shape == (4, 1)
    assert torch.isfinite(logits).all()

    criterion = torch.nn.BCEWithLogitsLoss()
    loss = criterion(logits.squeeze(-1), batch["label"])
    loss.backward()
    print(f"loss: {loss.item():.4f}, query grad norm: {model.query.grad.norm().item():.4f}")
    assert model.query.grad is not None
    assert torch.isfinite(model.query.grad).all()

    print("smoke test PASSED")


if __name__ == "__main__":
    main()
