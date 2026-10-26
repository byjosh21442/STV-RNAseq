# Raw QC assessment (Phases 2G–2I)

Date: 2026-10-25 · FastQC 0.13.0 · MultiQC 1.35 · run through `workflow/main.nf` (Nextflow 26.10.3, run `kickass_poitras`, 86 s; `-resume` re-run `wise_hodgkin`: all 7 tasks cached).
Inputs: the six original, MD5-validated `data/raw/*.fastq.gz` (not modified, not trimmed).
Outputs: `results/qc/fastqc_raw/<sample>_fastqc.{html,zip}` and `results/qc/multiqc_raw/multiqc_report.html` (all six samples present; git-ignored, regenerate with `nextflow run workflow/main.nf`). Tracked summaries: [`results/qc/raw_qc_summary.tsv`](../results/qc/raw_qc_summary.tsv), [`results/qc/raw_overrepresented.tsv`](../results/qc/raw_overrepresented.tsv) (from `scripts/bash/summarize_fastqc.sh`).

FastQC flags are diagnostic signals with generic thresholds; interpretation below follows the FastQC module documentation (<https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/>, accessed 2026-10-25). No biological interpretation is made.

## FastQC module flags

| Module | free_rep1 | free_rep2 | free_rep3 | inf_rep1 | inf_rep2 | inf_rep3 |
|---|---|---|---|---|---|---|
| Basic statistics | PASS | PASS | PASS | PASS | PASS | PASS |
| Per base sequence quality | PASS | PASS | PASS | PASS | PASS | PASS |
| Per tile sequence quality | PASS | PASS | PASS | PASS | PASS | PASS |
| Per sequence quality scores | PASS | PASS | PASS | PASS | PASS | PASS |
| Per base sequence content | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| Per sequence GC content | PASS | PASS | WARN | WARN | WARN | WARN |
| Per base N content | PASS | PASS | PASS | PASS | PASS | PASS |
| Sequence length distribution | PASS | PASS | PASS | PASS | PASS | PASS |
| Sequence duplication levels | FAIL | FAIL | FAIL | FAIL | FAIL | FAIL |
| Overrepresented sequences | WARN | PASS | WARN | WARN | WARN | WARN |
| Adapter content | PASS | PASS | PASS | PASS | PASS | PASS |

## Metrics

| Metric | free_rep1 | free_rep2 | free_rep3 | inf_rep1 | inf_rep2 | inf_rep3 |
|---|---:|---:|---:|---:|---:|---:|
| Reads (M) | 18.08 | 16.29 | 14.32 | 19.23 | 16.86 | 19.86 |
| Read length (nt) | 75 | 75 | 75 | 75 | 75 | 75 |
| Per-base mean Q, min–max | 33.5–38.8 | 33.5–38.8 | 33.5–38.8 | 33.5–38.8 | 33.5–38.8 | 33.5–38.8 |
| Lowest per-base lower quartile / 10th pct | 34 / 31 | 34 / 31 | 34 / 31 | 34 / 31 | 34 / 31 | 34 / 31 |
| Mode of per-read mean Q | 39 | 39 | 39 | 39 | 39 | 39 |
| Q30 bases (SeqKit, %) | 97 | 97 | 97 | 97 | 97 | 97 |
| Max N content (%) | 0.041 | 0.042 | 0.041 | 0.041 | 0.041 | 0.041 |
| Max adapter content, any type (%) | 0.209 (PolyA) | 0.110 (PolyA) | 0.144 (PolyA) | 0.191 (PolyA) | 0.110 (PolyA) | 0.156 (PolyA) |
| Max Illumina Universal Adapter (%) | 0.063 | 0.056 | 0.076 | 0.057 | 0.046 | 0.056 |
| Duplicate reads (MultiQC, %) | 68.8 | 65.4 | 64.8 | 69.5 | 66.6 | 67.5 |
| Reads in >10× duplication bins (%) | 49.6 | 44.3 | 45.5 | 51.1 | 47.5 | 48.1 |
| Mean GC (%) / mode | 42.8 / 44 | 44.0 / 44 | 44.1 / 45 | 43.8 / 45 | 44.2 / 45 | 43.8 / 45 |
| Overrepresented sequences (n; each % of reads) | 3; 0.10–0.14 | 0 | 3; 0.11–0.16 | 2; 0.14 | 2; 0.12–0.15 | 2; 0.13 |
| Max base-composition bias, pos 1–12 / 13–75 (%) | 23.5 / 5.1 | 25.1 / 4.7 | 26.5 / 5.0 | 24.3 / 5.2 | 25.7 / 4.8 | 24.5 / 4.7 |

## Module-by-module reading

- **Per-base sequence quality — PASS (all).** Mean Phred 33.5–38.8 along the whole read, lowest at position 4, 37.8–37.9 at position 75, so there is no 3′ quality drop-off. The lowest lower quartile is 34 and the lowest 10th percentile is 31. Q30 ≥ 97 %.
- **Per-sequence quality — PASS (all).** Modal read mean quality is Q39; there is no low-quality read sub-population.
- **Read length — PASS.** Every read is 75 nt (also confirmed with SeqKit over all reads).
- **Adapter content — PASS (all).** The highest signal is PolyA, 0.11–0.21 % of reads at the read end; Illumina Universal Adapter is ≤ 0.08 % and PolyG 0.002 %. This is far below the FastQC warning threshold (5 %). Salmon selective alignment uses extension alignment, so short tails need no clipping.
- **Overrepresented sequences — WARN in 5/6, PASS in free_rep2.** There are 0–3 sequences per sample, each 0.10–0.16 % of reads. FastQC labels all of them "No Hit", so they are not adapters or known contaminants. Two of them (`CTTCGGGCTTGTAAGCGATG…`, `GTCCAGTATCTGCCATCGTA…`) appear in 5 samples from both groups. Their identity has **not** been determined (that would need a BLAST lookup; not done, no interpretation). At ≤ 0.16 % each, their effect on library composition is negligible.
- **Duplication — FAIL (all).** About 65–70 % of reads are flagged as duplicates. The FastQC docs explain this: "In RNA-Seq libraries sequences from different transcripts will be present at wildly different levels… this will potentially create large set of duplicates", and duplicates cannot be classed as technical or biological from raw FASTQ. There are no UMIs, so deduplication is neither possible nor appropriate for count-based expression analysis. The levels are similar in both groups (64.8–68.8 % vs 66.6–69.5 %).
- **Per-base sequence content — FAIL (all).** The base-composition bias is confined to positions 1–12 (up to 23.5–26.5 % A/T or G/C imbalance) and falls to ≤ 5.2 % from position 13 onward. This is the random-hexamer priming bias that FastQC says "Nearly all RNA-Seq libraries" show; it "is not a problem which can be fixed by processing, and it doesn't seem to adversely affect the ability to measure expression." Trimming the first bases is not indicated. Salmon `--seqBias` models this kind of bias (option to consider at quantification time, not decided here).
- **Per-sequence GC — PASS in free_rep1/free_rep2; WARN in free_rep3 and all three infected.** Every distribution is unimodal: a single maximum at 44–45 % GC, mean 42.8–44.2 %, SD 6.7–7.4, with no sharp peaks (which would suggest adapter dimers or specific contaminants) and no secondary broad peak (which would suggest another species). My approximate normal-fit deviation is 9.2–10.8 %; FastQC fits its own model against a 15 % warning threshold, so the WARNs are borderline. Trimming or filtering would not change this. **Observation for human review:** WARN occurs in 3/3 infected vs 1/3 free samples. It is recorded, not interpreted.
- **N content — PASS.** Max 0.04 %.
- **Per-tile quality — PASS.**
- **Sequencing batch:** in the sampled read headers, all six libraries come from flowcell HN3JNBCXY lane 1 (`docs/raw_data_validation.md`). Quality profiles are practically identical across samples, consistent with one shared run.

## STV-free vs STV-infected (technical comparison only)

| Aspect | STV-free (rep1–3) | STV-infected (rep1–3) | Comment |
|---|---|---|---|
| Depth (M reads) | 18.08, 16.29, 14.32 (mean 16.23) | 19.23, 16.86, 19.86 (mean 18.65) | Infected ~15 % deeper on average; DESeq2 size factors handle depth |
| Base quality | identical profile | identical profile | same lane |
| Adapters | ≤ 0.21 % | ≤ 0.19 % | negligible in both |
| Duplication | 64.8–68.8 % | 66.6–69.5 % | similar |
| GC mean | 42.8–44.1 % | 43.8–44.2 % | similar; GC WARN 1/3 vs 3/3 (borderline, see above) |
| Overrepresented | 0–3 seqs, ≤ 0.14 % | 2 seqs, ≤ 0.15 % | same sequences recur in both groups |
| Head-bias (pos 1–12) | 23.5–26.5 % | 24.3–25.7 % | same library-prep signature |

No technical metric separates the groups strongly, and no sample is an outlier on quality. free_rep3 has the fewest reads (14.3 M), which is still adequate for gene-level counting.

## Conclusion

**A. No preprocessing appears necessary.**

Evidence:
1. Base quality is high and flat along all 75 nt (no 3′ drop-off) in all samples.
2. Adapter and poly-A/G content is ≤ 0.21 %.
3. N content is negligible.
4. All reads are full length.
5. The FAILs (duplication, per-base content) are the expected RNA-seq patterns that FastQC itself says cannot be fixed by processing.
6. The GC WARNs are borderline, unimodal, and not addressable by trimming.

Following D-003, no read was trimmed, filtered or modified. Two items are recorded for human review only, not as preprocessing triggers:

1. The GC WARN pattern (3/3 infected vs 1/3 free).
2. The identity of the two shared "No Hit" overrepresented sequences.

Confirming the "no trimming" outcome is a human decision; it is logged as D-009 in `docs/decisions.md`.
