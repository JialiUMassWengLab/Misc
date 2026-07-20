"""
Extract per-token ESM hidden-state embeddings for wild-type and variant protein
sequences from OMIM LOF/GOF missense variants.

Outputs two shard directories:
  wt_embeddings/        — shard_0000.pt, shard_0001.pt, …
  variant_embeddings/   — shard_0000.pt, shard_0001.pt, …

Each shard is a dict[key -> Tensor(seq_len+2, hidden_dim)].
  - WT shards: key = NP_ accession string
  - Variant shards: key = full variant_name string

Tensor shape is (len(seq) + 2, hidden_dim): BOS token at index 0, amino acids
at indices 1..len(seq), EOS at index len(seq)+1. Token at 1-based position P
is at index P — consistent with how compute_esm_llr.py indexes logits.

Memory-efficient: tensors are freed after each shard is written; only the set
of completed keys is kept in RAM between batches.

Loading all shards into a single dict:
    import torch
    from pathlib import Path
    emb = {}
    for p in sorted(Path("wt_embeddings").glob("shard_*.pt")):
        emb.update(torch.load(p, map_location="cpu"))

Or use merge_shards() provided at the bottom of this file.

Only missense variants matching p.XxxNNNYyy are processed; frameshift, early
termination, indel, splicing, and synonymous variants are skipped.

Run with:
    uv run python extract_esm_embeddings.py
"""

import re
import csv
import sys
import argparse
from pathlib import Path

import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

# ---------------------------------------------------------------------------
# Three-letter to one-letter amino acid map
# ---------------------------------------------------------------------------
AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*",
}

# ---------------------------------------------------------------------------
# Parse FASTA — keyed by NP_ accession
# ---------------------------------------------------------------------------
def load_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current_id = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(chunks)
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current_id is not None:
        sequences[current_id] = "".join(chunks)
    return sequences


# ---------------------------------------------------------------------------
# Parse variant metadata — yield only scorable missense rows
# ---------------------------------------------------------------------------
MISSENSE_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)$")


def parse_variants(path: Path):
    """Yield (variant_name, np_accession, wt_aa1, position_1based, mut_aa1)."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            variant_name = row["variant_name"]
            np_acc       = row["protein_accession"]

            m = MISSENSE_RE.search(variant_name)
            if m is None:
                continue  # frameshift / indel / splicing / synonymous

            wt_aa3, pos_str, mut_aa3 = m.group(1), m.group(2), m.group(3)

            wt_aa1  = AA3TO1.get(wt_aa3)
            mut_aa1 = AA3TO1.get(mut_aa3)
            if wt_aa1 is None or mut_aa1 is None:
                print(f"  [skip] unknown AA code in {variant_name}", file=sys.stderr)
                continue
            if wt_aa1 == mut_aa1:
                continue  # synonymous

            yield variant_name, np_acc, wt_aa1, int(pos_str), mut_aa1


# ---------------------------------------------------------------------------
# Shard-based persistence: only keys are held in RAM between batches
# ---------------------------------------------------------------------------
def load_done_keys(shard_dir: Path) -> set[str]:
    """Scan all existing shard files and collect their keys; tensors are freed immediately."""
    if not shard_dir.exists():
        return set()
    done: set[str] = set()
    for p in sorted(shard_dir.glob("shard_*.pt")):
        d = torch.load(p, map_location="cpu", weights_only=True)
        done.update(d.keys())
        del d
    if done:
        print(f"  resuming — {len(done)} keys already in {shard_dir}", file=sys.stderr)
    return done


def next_shard_index(shard_dir: Path) -> int:
    if not shard_dir.exists():
        return 0
    return len(sorted(shard_dir.glob("shard_*.pt")))


def save_shard(batch_results: dict, shard_dir: Path, shard_idx: int) -> None:
    shard_dir.mkdir(parents=True, exist_ok=True)
    torch.save(batch_results, shard_dir / f"shard_{shard_idx:04d}.pt")


# ---------------------------------------------------------------------------
# Build job lists — single validation pass for both WT and variant jobs
# ---------------------------------------------------------------------------
def build_jobs(
    records,
    sequences: dict[str, str],
    vocab: dict[str, int],
    already_done_wt: set[str],
    already_done_var: set[str],
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """
    Returns (wt_jobs, variant_jobs) where each job is (key, sequence).
    WT jobs are deduplicated by NP_ accession.
    Variant sequence has the mutant amino acid substituted at pos_0.
    """
    wt_jobs: list[tuple[str, str]] = []
    variant_jobs: list[tuple[str, str]] = []
    seen_wt: set[str] = set(already_done_wt)
    skipped = 0

    for variant_name, np_acc, wt_aa1, pos_1based, mut_aa1 in records:
        seq = sequences.get(np_acc)
        if seq is None:
            print(f"  [skip] no sequence for {np_acc} ({variant_name})", file=sys.stderr)
            skipped += 1
            continue

        pos_0 = pos_1based - 1
        if pos_0 >= len(seq):
            print(
                f"  [skip] position {pos_1based} out of range for {np_acc} "
                f"len={len(seq)} ({variant_name})",
                file=sys.stderr,
            )
            skipped += 1
            continue

        if seq[pos_0] != wt_aa1:
            print(
                f"  [skip] sequence mismatch at pos {pos_1based}: "
                f"expected {wt_aa1}, got {seq[pos_0]} ({variant_name})",
                file=sys.stderr,
            )
            skipped += 1
            continue

        if vocab.get(wt_aa1) is None or vocab.get(mut_aa1) is None:
            print(f"  [skip] AA not in vocab: {wt_aa1}/{mut_aa1} ({variant_name})", file=sys.stderr)
            skipped += 1
            continue

        if len(seq) > 10_000:
            print(f"  [skip] protein >10k aa: {np_acc} ({variant_name})", file=sys.stderr)
            skipped += 1
            continue

        if np_acc not in seen_wt:
            seen_wt.add(np_acc)
            wt_jobs.append((np_acc, seq))

        if variant_name not in already_done_var:
            mut_seq = seq[:pos_0] + mut_aa1 + seq[pos_0 + 1:]
            variant_jobs.append((variant_name, mut_seq))

    if skipped:
        print(f"  {skipped} records skipped (see above)", file=sys.stderr)

    return wt_jobs, variant_jobs


# ---------------------------------------------------------------------------
# Core batched embedding extraction
# ---------------------------------------------------------------------------
def run_batch(
    batch: list[tuple[str, str]],
    model,
    tokenizer,
    layer: int,
    device: str,
) -> dict[str, torch.Tensor]:
    """
    Run a forward pass and extract hidden states at the given layer.
    Returns {key: Tensor(actual_len, hidden_dim)} on CPU.
    actual_len = BOS + amino acids + EOS (all non-padding tokens).
    """
    keys = [b[0] for b in batch]
    seqs = [b[1] for b in batch]

    inputs = tokenizer(seqs, return_tensors="pt", padding=True)
    inputs = {k: v.to(device) for k, v in inputs.items()}

    with torch.inference_mode():
        output = model(**inputs, output_hidden_states=True)

    # ESMC hidden_states shape: [n_layers, B, L_padded, hidden_dim]
    layer_hidden = output.hidden_states[layer]  # [B, L_padded, hidden_dim]
    attn_mask = inputs["attention_mask"]         # [B, L_padded]

    result = {}
    for i, key in enumerate(keys):
        actual_len = int(attn_mask[i].sum().item())
        result[key] = layer_hidden[i, :actual_len, :].cpu()

    return result


def extract_embeddings(
    jobs: list[tuple[str, str]],
    model,
    tokenizer,
    layer: int,
    batch_size: int,
    shard_dir: Path,
    already_done: set[str],
    device: str,
    label: str = "items",
) -> None:
    """
    Process jobs in batches, saving each batch as a separate shard file.
    Tensors are freed from memory immediately after each shard is written —
    only the count of completed items is tracked in RAM between batches.
    """
    # Partition by sequence length
    too_long: list[tuple[str, str]] = []
    single:   list[tuple[str, str]] = []
    normal:   list[tuple[str, str]] = []

    for key, seq in jobs:
        n = len(seq)
        if n > 10_000:
            too_long.append((key, seq))
        elif n >= 5_000:
            single.append((key, seq))
        else:
            normal.append((key, seq))

    for key, _ in too_long:
        print(f"  [skip] >10k aa: {key}", file=sys.stderr)

    normal.sort(key=lambda x: len(x[1]))
    single.sort(key=lambda x: len(x[1]))

    total_todo = len(normal) + len(single)
    already_done_count = len(already_done)
    print(
        f"  {len(normal)} {label} at batch_size={batch_size}, "
        f"{len(single)} at batch_size=1 (5k-10k aa), "
        f"{len(too_long)} skipped (>10k aa)",
        file=sys.stderr,
    )

    shard_idx = next_shard_index(shard_dir)
    processed = 0

    def process_batch(batch: list[tuple[str, str]]) -> None:
        nonlocal shard_idx, processed
        batch_results = run_batch(batch, model, tokenizer, layer, device)
        save_shard(batch_results, shard_dir, shard_idx)
        del batch_results           # free tensors immediately
        shard_idx += 1
        processed += len(batch)
        total_so_far = already_done_count + processed
        grand_total  = already_done_count + total_todo
        print(
            f"  {total_so_far}/{grand_total} "
            f"({100 * total_so_far / max(1, grand_total):.1f}%) {label} done",
            file=sys.stderr,
        )

    for start in range(0, len(normal), batch_size):
        process_batch(normal[start : start + batch_size])

    for job in single:
        seq_len = len(job[1])
        process_batch([job])
        print(f"    [singleton: {seq_len} aa]", file=sys.stderr)


# ---------------------------------------------------------------------------
# Optional utility: merge all shards into a single dict (loads everything)
# ---------------------------------------------------------------------------
def merge_shards(shard_dir: Path, output_path: Path) -> None:
    """Concatenate all shard files into a single .pt dict file."""
    merged: dict = {}
    for p in sorted(shard_dir.glob("shard_*.pt")):
        d = torch.load(p, map_location="cpu", weights_only=True)
        merged.update(d)
        del d
    torch.save(merged, output_path)
    print(f"Merged {len(merged)} entries from {shard_dir} → {output_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract per-token ESMC embeddings for WT and variant proteins"
    )
    parser.add_argument("--model",        default="biohub/ESMC-600M")
    parser.add_argument("--layer",        type=int,  default=35)
    parser.add_argument("--batch-size",   type=int,  default=16)
    parser.add_argument(
        "--variants-tsv",
        type=Path,
        default=Path("../OMIM_variant_effects/lof_gof_variant_metadata.tsv"),
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=Path("../OMIM_variant_effects/lof_gof_wt_sequences.fasta"),
    )
    parser.add_argument(
        "--wt-out",
        type=Path,
        default=Path("wt_embeddings"),
        help="Output directory for WT embedding shards (default: wt_embeddings/)",
    )
    parser.add_argument(
        "--variant-out",
        type=Path,
        default=Path("variant_embeddings"),
        help="Output directory for variant embedding shards (default: variant_embeddings/)",
    )
    parser.add_argument(
        "--merge",
        action="store_true",
        help="After embedding, merge shards into wt_out.pt / variant_out.pt (loads all into RAM)",
    )
    args = parser.parse_args()

    print("Loading sequences...", file=sys.stderr)
    sequences = load_fasta(args.fasta)
    print(f"  {len(sequences)} protein sequences loaded", file=sys.stderr)

    print("Parsing variant metadata...", file=sys.stderr)
    records = list(parse_variants(args.variants_tsv))
    print(f"  {len(records)} missense variants parsed", file=sys.stderr)

    print(f"Loading model {args.model}...", file=sys.stderr)
    tokenizer = AutoTokenizer.from_pretrained(args.model)
    model = AutoModelForMaskedLM.from_pretrained(args.model, dtype="auto").eval()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    print(f"  model on {device}", file=sys.stderr)

    vocab = tokenizer.get_vocab()

    print("Scanning existing shards (resume)...", file=sys.stderr)
    wt_done      = load_done_keys(args.wt_out)
    variant_done = load_done_keys(args.variant_out)

    print("Building job lists...", file=sys.stderr)
    wt_jobs, variant_jobs = build_jobs(
        records, sequences, vocab,
        already_done_wt=wt_done,
        already_done_var=variant_done,
    )
    print(f"  {len(wt_jobs)} WT proteins to embed", file=sys.stderr)
    print(f"  {len(variant_jobs)} variant sequences to embed", file=sys.stderr)

    print(f"\n--- Wild-type embeddings (layer {args.layer}) ---", file=sys.stderr)
    extract_embeddings(
        wt_jobs, model, tokenizer, args.layer, args.batch_size,
        args.wt_out, wt_done, device, label="WT proteins",
    )

    print(f"\n--- Variant embeddings (layer {args.layer}) ---", file=sys.stderr)
    extract_embeddings(
        variant_jobs, model, tokenizer, args.layer, args.batch_size,
        args.variant_out, variant_done, device, label="variants",
    )

    wt_total      = len(wt_done)      + len(wt_jobs)
    variant_total = len(variant_done) + len(variant_jobs)
    print(
        f"\nDone. {wt_total} WT embeddings in {args.wt_out}/, "
        f"{variant_total} variant embeddings in {args.variant_out}/",
        file=sys.stderr,
    )

    if args.merge:
        print("\nMerging shards...", file=sys.stderr)
        merge_shards(args.wt_out,      args.wt_out.with_suffix(".pt"))
        merge_shards(args.variant_out, args.variant_out.with_suffix(".pt"))


if __name__ == "__main__":
    main()
