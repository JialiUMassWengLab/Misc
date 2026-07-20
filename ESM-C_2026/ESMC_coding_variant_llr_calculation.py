#!/usr/bin/env python3
"""
ESMC Coding Variant LLR Calculation

Computes ESM masked-marginal log-likelihood ratios (LLR) for missense variants
given as HGVS names. Protein sequences are fetched automatically from NCBI
using the NM accessions embedded in each HGVS name.

    LLR = log P(mutant | masked context) - log P(wildtype | masked context)

Only variants with a parseable missense change p.XxxNNNYyy are scored;
frameshifts, indels, and synonymous variants are skipped automatically.

Usage
-----
    uv run ESMC_coding_variant_llr_calculation.py \\
        --input variants.tsv \\
        [--output results.csv] \\
        [--sequences sequences.fasta] \\
        [--save-sequences sequences.fasta] \\
        [--model biohub/ESMC-600M] \\
        [--batch-size 16] \\
        [--entrez-email you@example.org]
"""

import argparse
import csv
import io
import re
import sys
import time
from pathlib import Path

import torch
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from transformers import AutoModelForMaskedLM, AutoTokenizer

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEFAULT_MODEL      = "biohub/ESMC-600M"
DEFAULT_BATCH_SIZE = 16
EFETCH_CHUNK       = 200  # NM accessions per NCBI efetch call

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Ter": "*",
}

MISSENSE_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
NM_RE       = re.compile(r"(N[MRP]_\d+(?:\.\d+)?)")


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------
def load_variants(path: Path) -> list[str]:
    """Return variant names from a TSV/CSV (with 'variant_name' column) or plain text."""
    text = path.read_text()
    lines = text.splitlines()
    if not lines:
        return []

    first = lines[0]
    if "\t" in first:
        delim = "\t"
    elif "," in first:
        delim = ","
    else:
        # Plain text — one HGVS name per line
        return [ln.strip() for ln in lines if ln.strip()]

    reader = csv.DictReader(io.StringIO(text), delimiter=delim)
    fieldnames = reader.fieldnames or []

    if "variant_name" in fieldnames:
        return [row["variant_name"] for row in reader if row.get("variant_name", "").strip()]

    # Fall back to the first column
    if fieldnames:
        col = fieldnames[0]
        return [row[col] for row in reader if row.get(col, "").strip()]

    return []


def parse_hgvs(variant_name: str) -> tuple[str, str, int, str] | None:
    """Parse an HGVS name → (nm_accession, wt_aa1, pos_1based, mut_aa1), or None.

    Returns None for non-missense variants (frameshift, indel, synonymous).
    """
    nm_m = NM_RE.search(variant_name)
    if nm_m is None:
        return None

    mis_m = MISSENSE_RE.search(variant_name)
    if mis_m is None:
        return None

    wt_aa3, pos_str, mut_aa3 = mis_m.group(1), mis_m.group(2), mis_m.group(3)
    wt_aa1  = AA3TO1.get(wt_aa3)
    mut_aa1 = AA3TO1.get(mut_aa3)
    if wt_aa1 is None or mut_aa1 is None:
        return None
    if wt_aa1 == mut_aa1:
        return None  # synonymous

    return nm_m.group(1), wt_aa1, int(pos_str), mut_aa1


# ---------------------------------------------------------------------------
# Retry wrapper
# ---------------------------------------------------------------------------
def _retry(fn, *args, max_retries=6, base_sleep=2, **kwargs):
    for attempt in range(max_retries):
        try:
            return fn(*args, **kwargs)
        except Exception as exc:
            if attempt == max_retries - 1:
                raise
            sleep = base_sleep * (2 ** attempt)
            print(f"  retry {attempt+1}/{max_retries} in {sleep}s ({exc})", file=sys.stderr)
            time.sleep(sleep)


# ---------------------------------------------------------------------------
# NCBI GenBank sequence fetching
# ---------------------------------------------------------------------------
def _fetch_gb_chunk(ids: list[str]) -> list:
    stream = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="gb", retmode="text")
    records = list(SeqIO.parse(stream, "genbank"))
    stream.close()
    return records


def _match_nm(record_id: str, chunk: list[str]) -> str:
    """Tolerate version differences: NM_000690 matches NM_000690.4."""
    base = record_id.split(".")[0]
    for nm in chunk:
        if nm == record_id or nm.split(".")[0] == base:
            return nm
    return record_id


def _parse_gb_record(record, chunk: list[str], result: dict, needs_protein_fetch: dict) -> None:
    matched_nm = _match_nm(record.id, chunk)
    cds_features = [f for f in record.features if f.type == "CDS"]
    if not cds_features:
        return  # non-coding RNA
    cds = cds_features[0]
    if "pseudo" in cds.qualifiers:
        return  # pseudogene
    translation = cds.qualifiers.get("translation", [""])[0]
    protein_id  = cds.qualifiers.get("protein_id",  [""])[0]
    if translation:
        result[matched_nm] = translation
    elif protein_id:
        needs_protein_fetch[matched_nm] = protein_id


def fetch_sequences(nm_list: list[str], entrez_email: str) -> dict[str, str]:
    """Fetch protein sequences for NM accessions; returns {nm_accession: protein_seq}.

    Three sources are tried in order:
    1. CDS /translation qualifier (bulk GenBank efetch) — covers most NMs.
    2. Solo re-fetch for NMs dropped by NCBI version deduplication.
    3. Protein-db efetch by /protein_id — for records with no inline translation.
    """
    Entrez.email = entrez_email
    result: dict[str, str] = {}
    needs_protein_fetch: dict[str, str] = {}

    chunks = [nm_list[i:i + EFETCH_CHUNK] for i in range(0, len(nm_list), EFETCH_CHUNK)]
    print(f"  Fetching {len(nm_list)} GenBank records in {len(chunks)} chunks...", file=sys.stderr)

    for k, chunk in enumerate(chunks, 1):
        try:
            records = _retry(_fetch_gb_chunk, chunk)
            for record in records:
                _parse_gb_record(record, chunk, result, needs_protein_fetch)

            # Re-fetch NMs silently dropped by NCBI when two versions of the same
            # base accession appear in one chunk (NCBI deduplicates by base accession).
            covered_bases: set[str] = set()
            dropped: list[str] = []
            for nm in chunk:
                base = nm.split(".")[0]
                if nm in result or nm in needs_protein_fetch:
                    covered_bases.add(base)
                elif base in covered_bases:
                    dropped.append(nm)
                else:
                    covered_bases.add(base)

            for nm in dropped:
                try:
                    solo = _retry(_fetch_gb_chunk, [nm])
                    for record in solo:
                        _parse_gb_record(record, [nm], result, needs_protein_fetch)
                    time.sleep(0.3)
                except Exception as e:
                    print(f"    solo re-fetch failed for {nm}: {e}", file=sys.stderr)

        except Exception as e:
            print(f"  chunk {k} failed: {e}", file=sys.stderr)

        time.sleep(0.4)
        if k % 20 == 0:
            print(
                f"  ... {k}/{len(chunks)} chunks done, {len(result)}/{len(nm_list)} resolved",
                file=sys.stderr,
            )

    # Fallback: fetch from protein db for records with protein_id but no inline translation
    if needs_protein_fetch:
        np_unique = list(dict.fromkeys(needs_protein_fetch.values()))
        print(f"  Fallback: fetching {len(np_unique)} NP sequences from protein db...", file=sys.stderr)
        np_to_seq: dict[str, str] = {}
        prot_chunks = [np_unique[i:i + EFETCH_CHUNK] for i in range(0, len(np_unique), EFETCH_CHUNK)]
        for chunk in prot_chunks:
            def _fetch_prot(ids):
                stream = Entrez.efetch(db="protein", id=",".join(ids), rettype="fasta", retmode="text")
                recs = list(SeqIO.parse(stream, "fasta"))
                stream.close()
                return recs
            try:
                prot_recs = _retry(_fetch_prot, chunk)
                for rec in prot_recs:
                    parts = rec.id.split("|")
                    np_acc = rec.id
                    for j, part in enumerate(parts):
                        if part in ("ref", "sp", "gb") and j + 1 < len(parts):
                            np_acc = parts[j + 1]
                            break
                    np_to_seq[np_acc] = str(rec.seq)
            except Exception as e:
                print(f"  protein fallback chunk failed: {e}", file=sys.stderr)
            time.sleep(0.4)

        for nm, np_acc in needs_protein_fetch.items():
            seq = np_to_seq.get(np_acc, "")
            if seq:
                result[nm] = seq
            else:
                print(f"  WARNING: no sequence retrieved for {nm} ({np_acc})", file=sys.stderr)

    return result


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------
def load_fasta(path: Path) -> dict[str, str]:
    """Load sequences from FASTA keyed by first word of the header line."""
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


def save_fasta(sequences: dict[str, str], path: Path) -> None:
    records = [SeqRecord(Seq(seq), id=nm, description="") for nm, seq in sequences.items()]
    with open(path, "w") as fh:
        SeqIO.write(records, fh, "fasta")
    print(f"  Saved {len(records)} sequences to {path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# LLR scoring
# ---------------------------------------------------------------------------
def load_existing_results(path: Path) -> set[str]:
    """Return variant names already present in the output CSV (for resume support)."""
    if not path.exists():
        return set()
    done: set[str] = set()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            done.add(row["variant_name"])
    return done


def score_variants(
    variants: list[str],
    sequences: dict[str, str],
    model,
    tokenizer,
    vocab: dict[str, int],
    batch_size: int,
    output_path: Path,
) -> int:
    """Score all variants with ESMC masked-marginal LLR; writes CSV incrementally.

    Returns total number of variants scored (including previously completed ones).
    Skips variants already in output_path, enabling safe resume after interruption.
    """
    already_done = load_existing_results(output_path)
    if already_done:
        print(f"  Resuming — {len(already_done)} variants already scored", file=sys.stderr)

    skipped = 0
    job_list: list[tuple[str, str, int, int, int]] = []

    for variant_name in variants:
        if variant_name in already_done:
            continue

        parsed = parse_hgvs(variant_name)
        if parsed is None:
            print(f"  [skip] non-missense or unparseable: {variant_name}", file=sys.stderr)
            skipped += 1
            continue

        nm, wt_aa1, pos_1based, mut_aa1 = parsed
        seq = sequences.get(nm)
        if seq is None:
            print(f"  [skip] no sequence for {nm} ({variant_name})", file=sys.stderr)
            skipped += 1
            continue

        pos_0 = pos_1based - 1
        if pos_0 >= len(seq):
            print(
                f"  [skip] position {pos_1based} out of range for {nm} "
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

        wt_idx  = vocab.get(wt_aa1)
        mut_idx = vocab.get(mut_aa1)
        if wt_idx is None or mut_idx is None:
            print(f"  [skip] AA not in model vocab: {wt_aa1}/{mut_aa1} ({variant_name})", file=sys.stderr)
            skipped += 1
            continue

        masked_seq = seq[:pos_0] + tokenizer.mask_token + seq[pos_0 + 1:]
        job_list.append((variant_name, masked_seq, wt_idx, mut_idx, pos_1based))

    # Partition jobs by sequence length to avoid GPU OOM on very long proteins
    too_long, single, normal = [], [], []
    for job in job_list:
        seq_len = len(job[1])
        if seq_len > 10_000:
            too_long.append(job)
        elif seq_len >= 5_000:
            single.append(job)
        else:
            normal.append(job)

    for vname, *_ in too_long:
        print(f"  [skip] protein >10k aa: {vname}", file=sys.stderr)
    skipped += len(too_long)

    print(
        f"  {len(normal)} variants at batch_size={batch_size}, "
        f"{len(single)} at batch_size=1 (5k–10k aa), "
        f"{len(too_long)} skipped (>10k aa), "
        f"{skipped} total skipped",
        file=sys.stderr,
    )

    # Sort ascending by length to minimise padding waste within batches
    normal.sort(key=lambda x: len(x[1]))
    single.sort(key=lambda x: len(x[1]))

    write_header = not output_path.exists() or len(already_done) == 0
    out_fh = open(output_path, "a", newline="")
    writer = csv.writer(out_fh)
    if write_header:
        writer.writerow(["variant_name", "esm_llr"])

    device = model.device
    total_scored = len(already_done)
    total_to_score = len(normal) + len(single) + len(already_done)

    def run_batch(batch):
        vnames    = [b[0] for b in batch]
        seqs      = [b[1] for b in batch]
        wt_idxs   = [b[2] for b in batch]
        mut_idxs  = [b[3] for b in batch]
        positions = [b[4] for b in batch]

        inputs = tokenizer(seqs, return_tensors="pt", padding=True)
        inputs = {k: v.to(device) for k, v in inputs.items()}

        with torch.inference_mode():
            logits = model(**inputs).logits  # (B, L_padded, vocab_size)

        log_probs = torch.log_softmax(logits, dim=-1)

        for i, (vname, pos_1based, wt_idx, mut_idx) in enumerate(
            zip(vnames, positions, wt_idxs, mut_idxs)
        ):
            # pos_1based == token index because ESM prepends a BOS token,
            # shifting 0-indexed AA positions by +1 to match 1-indexed positions.
            lp_wt  = log_probs[i, pos_1based, wt_idx].item()
            lp_mut = log_probs[i, pos_1based, mut_idx].item()
            writer.writerow([vname, f"{lp_mut - lp_wt:.6f}"])

    try:
        for batch_start in range(0, len(normal), batch_size):
            run_batch(normal[batch_start : batch_start + batch_size])
            out_fh.flush()
            total_scored += min(batch_size, len(normal) - batch_start)
            pct = 100 * total_scored / total_to_score if total_to_score else 100
            print(f"  scored {total_scored}/{total_to_score} ({pct:.1f}%)", file=sys.stderr)

        for job in single:
            run_batch([job])
            out_fh.flush()
            total_scored += 1
            pct = 100 * total_scored / total_to_score if total_to_score else 100
            print(
                f"  scored {total_scored}/{total_to_score} ({pct:.1f}%)"
                f"  [singleton: {len(job[1])} aa]",
                file=sys.stderr,
            )
    finally:
        out_fh.close()

    return total_scored


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Compute ESMC masked-marginal LLR scores for HGVS missense variants.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", required=True, type=Path,
        help="Input file: TSV/CSV with 'variant_name' column, or plain text (one HGVS name per line)",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("llr_results.csv"),
        help="Output CSV with variant_name and esm_llr columns",
    )
    parser.add_argument(
        "--sequences", type=Path, default=None,
        help="Pre-fetched sequences FASTA keyed by NM accession (skips NCBI fetch)",
    )
    parser.add_argument(
        "--save-sequences", type=Path, default=None,
        help="Save fetched sequences to this FASTA path for reuse",
    )
    parser.add_argument(
        "--model", default=DEFAULT_MODEL,
        help="HuggingFace masked-LM model identifier",
    )
    parser.add_argument(
        "--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
        help="Sequences per GPU batch (reduce if OOM)",
    )
    parser.add_argument(
        "--entrez-email", default="esmc.scoring@example.org",
        help="Email address for NCBI Entrez API (required when fetching sequences)",
    )
    args = parser.parse_args()

    # 1. Load variant names
    print("Loading variants...", file=sys.stderr)
    variants = load_variants(args.input)
    if not variants:
        sys.exit(f"ERROR: no variants found in {args.input}")
    print(f"  {len(variants)} variants loaded", file=sys.stderr)

    # 2. Get protein sequences
    if args.sequences and args.sequences.exists():
        print(f"\nLoading sequences from {args.sequences}...", file=sys.stderr)
        sequences = load_fasta(args.sequences)
        print(f"  {len(sequences)} sequences loaded", file=sys.stderr)
    else:
        nm_set: set[str] = set()
        for v in variants:
            parsed = parse_hgvs(v)
            if parsed:
                nm_set.add(parsed[0])
        nm_list = sorted(nm_set)
        print(f"\nFetching sequences for {len(nm_list)} unique NM accessions from NCBI...", file=sys.stderr)
        sequences = fetch_sequences(nm_list, args.entrez_email)
        print(f"  {len(sequences)}/{len(nm_list)} sequences fetched", file=sys.stderr)

        if args.save_sequences:
            save_fasta(sequences, args.save_sequences)

    # 3. Load ESMC model
    print(f"\nLoading model {args.model}...", file=sys.stderr)
    tokenizer = AutoTokenizer.from_pretrained(args.model)
    model = AutoModelForMaskedLM.from_pretrained(args.model, dtype="auto").eval()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    print(f"  model on {device}", file=sys.stderr)
    vocab = tokenizer.get_vocab()

    # 4. Score variants and write output
    print(f"\nScoring variants → {args.output}...", file=sys.stderr)
    total = score_variants(
        variants, sequences, model, tokenizer, vocab, args.batch_size, args.output
    )
    print(f"\nDone. {total} LLR scores written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
