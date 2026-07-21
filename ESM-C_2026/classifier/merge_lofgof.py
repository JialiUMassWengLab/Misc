"""Merge LOF+GOF embedding-diff tensors from all shards into one .pt file.

Mirrors extract_esm_embeddings.merge_shards() but filters to variants whose
category in combined_variant_metadata.tsv starts with L- or G-.

Run once:
    python -m classifier.merge_lofgof

Produces embedding_diff_lofgof.pt at project root (~21.6 GB, float32,
6,464 entries). Subsequent training runs load this single file directly
instead of iterating all 529 shards.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import torch


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def merge_lofgof(metadata_tsv: Path, embedding_dir: Path, out_path: Path) -> None:
    df = pd.read_csv(metadata_tsv, sep="\t")
    cat = df["category"].astype(str)
    keep = cat.str.startswith("L-") | cat.str.startswith("G-")
    wanted: set[str] = set(df.loc[keep, "variant_name"].astype(str))
    print(f"metadata: {keep.sum()} LOF+GOF variants requested", file=sys.stderr)

    shard_paths = sorted(embedding_dir.glob("shard_*.pt"))
    if not shard_paths:
        raise FileNotFoundError(f"No shard_*.pt files under {embedding_dir}")

    merged: dict[str, torch.Tensor] = {}
    for i, p in enumerate(shard_paths):
        d = torch.load(p, map_location="cpu", weights_only=True)
        for k in d.keys() & wanted:
            merged[k] = d[k].float()
        del d
        if (i + 1) % 50 == 0 or i == len(shard_paths) - 1:
            print(f"  shards {i + 1}/{len(shard_paths)}, kept {len(merged)}",
                  file=sys.stderr)

    torch.save(merged, out_path)
    size_gb = out_path.stat().st_size / 1e9
    print(f"wrote {len(merged)} entries to {out_path} ({size_gb:.2f} GB)",
          file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", type=Path,
                        default=PROJECT_ROOT / "combined_variant_metadata.tsv")
    parser.add_argument("--embedding-dir", type=Path,
                        default=PROJECT_ROOT / "embedding_diff")
    parser.add_argument("--out", type=Path,
                        default=PROJECT_ROOT / "embedding_diff_lofgof.pt")
    args = parser.parse_args()
    merge_lofgof(args.metadata, args.embedding_dir, args.out)


if __name__ == "__main__":
    main()
