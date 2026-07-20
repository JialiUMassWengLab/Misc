# ESMC Coding Variant LLR Calculation

## Background

### What is the ESM masked-marginal LLR?

ESM (Evolutionary Scale Modeling) protein language models are trained using a masked-language-modeling objective: the model learns to predict a masked amino acid from its surrounding sequence context. The **masked-marginal log-likelihood ratio (LLR)** exploits this capability to score the functional impact of a missense substitution:

```
LLR = log P(mutant AA | masked context) − log P(wildtype AA | masked context)
```

A **negative LLR** means the model assigns lower probability to the mutant amino acid than the wildtype given the same sequence context — consistent with the substitution being evolutionarily deleterious. A **positive LLR** indicates the mutant is at least as plausible as the wildtype. The LLR has been validated as a predictor of pathogenicity across several benchmarks (Meier et al. 2021, Brandes et al. 2023).

### The BOS token offset

ESM's tokenizer prepends a beginning-of-sequence (`[CLS]`) token before the first amino acid. This shifts all amino-acid token positions by +1, so the 1-indexed protein position `p` corresponds to token index `p` (not `p−1`). The script exploits this coincidence directly rather than adjusting explicitly.

### Model used

The default model is [`biohub/ESMC-600M`](https://huggingface.co/biohub/ESMC-600M), the 600M-parameter masked protein language model released by EvolutionaryScale. The `--model` flag accepts any HuggingFace `AutoModelForMaskedLM`-compatible model identifier.

---

## Input format

The `--input` file is flexible:

| Format | Description |
|--------|-------------|
| **TSV/CSV with `variant_name` column** | HGVS names are read from the `variant_name` column (e.g., `lof_gof_phase1.tsv`) |
| **Plain text** | One HGVS name per line |

HGVS names must contain an NM accession and a missense protein change. Supported example:

```
NM_000690.4(ALDH2):c.1510G>A (p.Glu504Lys)
```

Variants without a parseable `p.XxxNNNYyy` change (frameshifts, indels, splicing, synonymous) are automatically skipped and reported to stderr.

---

## Output format

The output CSV has two columns:

```
variant_name,esm_llr
NM_000690.4(ALDH2):c.1510G>A (p.Glu504Lys),-1.234567
```

The script writes results **incrementally**: each batch is flushed to disk before the next begins. If the run is interrupted, re-running the same command will resume from where it stopped (already-scored variants are detected from the output file and skipped).

---

## Usage

```bash
uv run ESMC_coding_variant_llr_calculation.py \
    --input variants.tsv \
    --output results.csv \
    --entrez-email your.email@institution.org
```

### All options

| Flag | Default | Description |
|------|---------|-------------|
| `--input` | *(required)* | Input variants file (TSV/CSV with `variant_name` column, or plain text) |
| `--output` | `llr_results.csv` | Output CSV file |
| `--sequences` | — | Pre-fetched sequence FASTA keyed by NM accession; skips NCBI fetch |
| `--save-sequences` | — | Save NCBI-fetched sequences to this FASTA path for later reuse |
| `--model` | `biohub/ESMC-600M` | HuggingFace masked-LM model identifier |
| `--batch-size` | `16` | Sequences per GPU batch (reduce to 1–4 if GPU OOM) |
| `--entrez-email` | `esmc.scoring@example.org` | Email address required by NCBI Entrez API |

### Reusing cached sequences

Fetching sequences from NCBI can take several minutes for large variant sets. Save them once and reuse:

```bash
# First run: fetch and save
uv run ESMC_coding_variant_llr_calculation.py \
    --input variants.tsv \
    --save-sequences sequences.fasta \
    --entrez-email your.email@institution.org

# Subsequent runs: skip NCBI fetch
uv run ESMC_coding_variant_llr_calculation.py \
    --input variants.tsv \
    --sequences sequences.fasta
```

---

## Skipped variants

The script reports each skipped variant to stderr with a reason tag:

| Tag | Reason |
|-----|--------|
| `[skip] non-missense or unparseable` | No NM accession or no `p.XxxNNNYyy` in HGVS name |
| `[skip] no sequence for NM_...` | NM accession not resolved by NCBI fetch |
| `[skip] sequence mismatch at pos N` | Wildtype AA in HGVS does not match the fetched sequence |
| `[skip] position N out of range` | Position exceeds fetched sequence length |
| `[skip] AA not in model vocab` | Amino acid not recognized by the tokenizer |
| `[skip] protein >10k aa` | Sequence too long; would exceed GPU memory |

---

## Memory and runtime notes

- Proteins with sequences 5,000–10,000 aa are scored one at a time (`batch_size=1`) to avoid OOM.
- Proteins longer than 10,000 aa are skipped entirely.
- On a single A10G GPU (24 GB), `--batch-size 16` works for typical human proteins (<2,000 aa). Reduce to 4–8 for larger proteins or smaller GPUs.
- The NCBI fetch respects Entrez rate limits with inter-chunk delays (0.4 s) and exponential-backoff retries (up to 6 attempts).
