"""Datasets for LOF vs GOF variant classification from ESM-C embedding differences.

Embeddings are loaded into memory once from embedding_diff_lofgof.pt (produced
by classifier/merge_lofgof.py). All splits share the same underlying dict.

Two dataset classes:
- VariantEmbeddingDataset: returns the full (seq_len+2, 1152) tensor per variant.
- VariantPositionDataset:  returns only the (1152,) row at the AA substitution
                           position (BOS+1-indexed, so emb[N]).

`load_lofgof_embeddings` produces a shared `rows` list of
(variant_name, label, sub_pos) tuples. Non-missense and Ter (stop-gain /
stop-loss) variants are dropped there so both datasets operate on identical data.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import torch
from torch.utils.data import Dataset


LABEL_LOF = 0
LABEL_GOF = 1

# Copy of MISSENSE_RE from extract_esm_embeddings.py:81; keep in sync.
# Imported by value (not module import) because that file loads `transformers`
# at import time — heavy for a single regex.
MISSENSE_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)$")


def _parse_missense(variant_name: str) -> tuple[str, int, str] | None:
    """Return (wt_aa3, pos, mut_aa3) for a true missense variant, else None.
    Rejects stop-gain (mut='Ter') and stop-loss/extension (wt='Ter').
    """
    m = MISSENSE_RE.search(variant_name)
    if m is None:
        return None
    wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
    if wt == "Ter" or mut == "Ter":
        return None
    return wt, pos, mut


def load_lofgof_embeddings(
    cache_path: str | Path,
    metadata_tsv: str | Path,
) -> tuple[dict[str, torch.Tensor], list[tuple[str, int, int]]]:
    """Load the merged LOF+GOF embedding cache and the ordered rows list.

    Returns:
        embeddings: variant_name -> float32 Tensor of shape (seq_len+2, 1152)
        rows: (variant_name, label, sub_pos) tuples in metadata-row order.
              label 0=LOF, 1=GOF; sub_pos is the 1-indexed AA position.
              Non-missense and Ter variants have been filtered out.
    """
    cache_path = Path(cache_path)
    if not cache_path.exists():
        raise FileNotFoundError(
            f"{cache_path} not found. Run `python -m classifier.merge_lofgof` first."
        )
    print(f"loading {cache_path} ...")
    embeddings: dict[str, torch.Tensor] = torch.load(
        cache_path, map_location="cpu", weights_only=True
    )
    print(f"  {len(embeddings)} embeddings loaded")

    df = pd.read_csv(metadata_tsv, sep="\t")
    cat = df["category"].astype(str)
    keep = cat.str.startswith("L-") | cat.str.startswith("G-")
    df = df.loc[keep].copy()
    df["label"] = df["category"].str.startswith("G-").astype(int)

    in_cache = df["variant_name"].isin(embeddings)
    missing = int((~in_cache).sum())
    if missing:
        print(f"  dropping {missing} metadata rows without embeddings")
    df = df.loc[in_cache].reset_index(drop=True)

    parsed = df["variant_name"].map(_parse_missense)
    kept = parsed.notna()
    n_dropped = int((~kept).sum())
    if n_dropped:
        print(f"  dropping {n_dropped} non-missense / Ter variants "
              f"(stop-gain, stop-loss, frameshift, splice, etc.)")
    df = df.loc[kept].copy()
    df["sub_pos"] = parsed.loc[kept].map(lambda t: t[1])

    rows = [
        (str(r.variant_name), int(r.label), int(r.sub_pos))
        for r in df.itertuples(index=False)
    ]
    print(f"  {len(rows)} usable variants after filtering")
    return embeddings, rows


class VariantEmbeddingDataset(Dataset):
    """Serves (full-sequence embedding, label) pairs for LOF/GOF variants.

    Labels: L-* -> 0, G-* -> 1. All splits share the same underlying
    `embeddings` dict — pass in the result of load_lofgof_embeddings().
    """

    def __init__(
        self,
        embeddings: dict[str, torch.Tensor],
        rows: list[tuple[str, int, int]],
        split: str = "train",
        split_seed: int = 42,
        split_ratios: tuple[float, float, float] = (0.75, 0.25, 0.0),
    ) -> None:
        if split not in {"train", "val", "test"}:
            raise ValueError(f"split must be train/val/test, got {split!r}")
        if abs(sum(split_ratios) - 1.0) > 1e-6:
            raise ValueError(f"split_ratios must sum to 1, got {split_ratios}")

        self.embeddings = embeddings

        n = len(rows)
        perm = torch.randperm(n, generator=torch.Generator().manual_seed(split_seed)).tolist()
        n_train = int(n * split_ratios[0])
        n_val = int(n * split_ratios[1])
        splits = {
            "train": perm[:n_train],
            "val": perm[n_train : n_train + n_val],
            "test": perm[n_train + n_val :],
        }
        indices = splits[split]
        self.rows_split = [rows[i] for i in indices]

        self.n_lof = sum(1 for _, y, _ in self.rows_split if y == LABEL_LOF)
        self.n_gof = sum(1 for _, y, _ in self.rows_split if y == LABEL_GOF)
        print(f"[{split}] {len(self.rows_split)} variants "
              f"({self.n_lof} LOF, {self.n_gof} GOF)")

    def __len__(self) -> int:
        return len(self.rows_split)

    def __getitem__(self, idx: int) -> dict:
        variant_name, label, _sub_pos = self.rows_split[idx]
        return {
            "embedding": self.embeddings[variant_name],
            "label": torch.tensor(label, dtype=torch.float32),
            "variant_name": variant_name,
        }


class VariantPositionDataset(Dataset):
    """Serves only the (1152,) embedding-diff row at the substitution position.

    Uses the same split logic (seed, ratios) as VariantEmbeddingDataset. Since
    `rows` is pre-filtered by load_lofgof_embeddings() to true missense variants,
    every row here has a valid sub_pos.
    """

    def __init__(
        self,
        embeddings: dict[str, torch.Tensor],
        rows: list[tuple[str, int, int]],
        split: str = "train",
        split_seed: int = 42,
        split_ratios: tuple[float, float, float] = (0.75, 0.25, 0.0),
    ) -> None:
        if split not in {"train", "val", "test"}:
            raise ValueError(f"split must be train/val/test, got {split!r}")
        if abs(sum(split_ratios) - 1.0) > 1e-6:
            raise ValueError(f"split_ratios must sum to 1, got {split_ratios}")

        self.embeddings = embeddings

        n = len(rows)
        perm = torch.randperm(n, generator=torch.Generator().manual_seed(split_seed)).tolist()
        n_train = int(n * split_ratios[0])
        n_val = int(n * split_ratios[1])
        splits = {
            "train": perm[:n_train],
            "val": perm[n_train : n_train + n_val],
            "test": perm[n_train + n_val :],
        }
        indices = splits[split]
        self.rows_split = [rows[i] for i in indices]

        self.n_lof = sum(1 for _, y, _ in self.rows_split if y == LABEL_LOF)
        self.n_gof = sum(1 for _, y, _ in self.rows_split if y == LABEL_GOF)
        print(f"[{split} / pos-only] {len(self.rows_split)} variants "
              f"({self.n_lof} LOF, {self.n_gof} GOF)")

    def __len__(self) -> int:
        return len(self.rows_split)

    def __getitem__(self, idx: int) -> dict:
        variant_name, label, pos = self.rows_split[idx]
        emb = self.embeddings[variant_name]  # (seq_len+2, 1152)
        # BOS at index 0 offsets 1-indexed AA positions, so emb[pos] is the
        # substitution residue's embedding difference.
        assert 1 <= pos <= emb.shape[0] - 2, (
            f"pos {pos} out of range [1, {emb.shape[0] - 2}] for {variant_name}"
        )
        return {
            "embedding": emb[pos],  # (1152,)
            "label": torch.tensor(label, dtype=torch.float32),
            "variant_name": variant_name,
        }


def collate_fn(batch: list[dict]) -> dict:
    """Pad variable-length embeddings to the batch max length; build a bool mask."""
    lengths = [item["embedding"].shape[0] for item in batch]
    max_len = max(lengths)
    embed_dim = batch[0]["embedding"].shape[1]
    b = len(batch)

    padded = torch.zeros(b, max_len, embed_dim, dtype=torch.float32)
    attention_mask = torch.zeros(b, max_len, dtype=torch.bool)
    labels = torch.zeros(b, dtype=torch.float32)
    variant_names: list[str] = []

    for i, item in enumerate(batch):
        seq_len = lengths[i]
        padded[i, :seq_len] = item["embedding"]
        attention_mask[i, :seq_len] = True
        labels[i] = item["label"]
        variant_names.append(item["variant_name"])

    return {
        "embedding": padded,
        "attention_mask": attention_mask,
        "label": labels,
        "variant_name": variant_names,
    }


def collate_fn_position(batch: list[dict]) -> dict:
    """Trivial stack for fixed-size (D,) position-only vectors."""
    return {
        "embedding": torch.stack([b["embedding"] for b in batch]),   # (B, D)
        "label": torch.stack([b["label"] for b in batch]),           # (B,)
        "variant_name": [b["variant_name"] for b in batch],
    }


class VariantPrecomputedPositionDataset(Dataset):
    """Serves pre-pooled (D,) position-only embedding vectors.

    Unlike VariantPositionDataset which slices row `pos` out of a per-residue
    (seq_len+2, D) tensor, this dataset expects each cache entry to already be
    a 1-D vector — someone did the position lookup upstream. Used for caches
    like `6B_embedding_diffs_merged.pt` (D=2560).

    Uses the same split logic (seed, ratios) as the other datasets. Compatible
    with `collate_fn_position`.
    """

    def __init__(
        self,
        embeddings: dict[str, torch.Tensor],
        rows: list[tuple[str, int, int]],
        split: str = "train",
        split_seed: int = 42,
        split_ratios: tuple[float, float, float] = (0.75, 0.25, 0.0),
    ) -> None:
        if split not in {"train", "val", "test"}:
            raise ValueError(f"split must be train/val/test, got {split!r}")
        if abs(sum(split_ratios) - 1.0) > 1e-6:
            raise ValueError(f"split_ratios must sum to 1, got {split_ratios}")

        self.embeddings = embeddings

        n = len(rows)
        perm = torch.randperm(n, generator=torch.Generator().manual_seed(split_seed)).tolist()
        n_train = int(n * split_ratios[0])
        n_val = int(n * split_ratios[1])
        splits = {
            "train": perm[:n_train],
            "val": perm[n_train : n_train + n_val],
            "test": perm[n_train + n_val :],
        }
        indices = splits[split]
        self.rows_split = [rows[i] for i in indices]

        self.n_lof = sum(1 for _, y, _ in self.rows_split if y == LABEL_LOF)
        self.n_gof = sum(1 for _, y, _ in self.rows_split if y == LABEL_GOF)
        print(f"[{split} / precomputed-pos] {len(self.rows_split)} variants "
              f"({self.n_lof} LOF, {self.n_gof} GOF)")

    def __len__(self) -> int:
        return len(self.rows_split)

    def __getitem__(self, idx: int) -> dict:
        variant_name, label, _pos = self.rows_split[idx]
        emb = self.embeddings[variant_name]  # (D,) already the substitution-position vector
        assert emb.ndim == 1, (
            f"expected 1-D pre-pooled vector for {variant_name}, got shape {tuple(emb.shape)}"
        )
        return {
            "embedding": emb.float(),
            "label": torch.tensor(label, dtype=torch.float32),
            "variant_name": variant_name,
        }
