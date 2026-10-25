# STV_GSE137303 — tomato leaf RNA-seq re-analysis (Southern tomato virus)

## Objective

Re-analyse public RNA-seq data from *Solanum lycopersicum* cv. M82 leaves, STV-free vs STV-infected, to obtain a reproducible differential-expression result (DESeq2) and plant-specific functional analysis. University project, ~2 weeks; robustness and simplicity over complexity.

## Dataset

| | |
|---|---|
| GEO | [GSE137303](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137303) (PubMed 31620899) |
| BioProject / SRA study | PRJNA565137 / SRP221436 |
| Samples | 6 (GSM4074979–GSM4074984), one SRA run each (SRR10102791–SRR10102796) |
| Library | Illumina HiSeq 2500, **single-end 75 nt**, "TruSeq stranded mRNA" kit (read orientation not documented) |
| Raw data size | 104.6 M reads; 5.7 GB `fastq.gz` at ENA |

Verified sample table: [data/metadata/samples_metadata.tsv](data/metadata/samples_metadata.tsv) · audit: [docs/dataset_audit.md](docs/dataset_audit.md).

## Experimental design

| Condition | Replicates | sample_id |
|---|---|---|
| STV_free (intended reference level) | 3 | STV_free_rep1–3 |
| STV_infected | 3 | STV_infected_rep1–3 |

## Intended workflow (pending approvals — see [docs/decisions.md](docs/decisions.md))

```text
ENA FASTQ (MD5-checked) → FastQC → MultiQC → [trimming only if QC justifies]
→ Salmon (decoy-aware index) → tximport → DESeq2 → plant functional analysis
```

Nextflow will automate FASTQ → QC → quantification; tximport/DESeq2/functional analysis run as R scripts.

## Directory structure

```text
config/            samplesheet.csv (accession-based, generated)
workflow/          Nextflow pipeline (empty: Nextflow not installed yet)
scripts/bash/      audit, metadata and GEO-file scripts
scripts/R/         tximport / DESeq2 / functional analysis (later)
data/metadata/     samples_metadata.tsv + raw source records (sources/)
data/raw/          FASTQ (git-ignored, not downloaded yet)
data/processed/    original_geo/GSE137303_count.txt.gz (authors' table, audit only)
reference/         genome/annotation/indexes (git-ignored, not downloaded yet)
results/           qc/, quantification/, deseq2/, functional/
figures/, envs/    final figures; environment definitions (later)
docs/              audits, design, decisions, run log
```

## Current status (Phase 1, 2026-10-25)

- Done: environment audit, dataset audit (six-sample mapping verified across GEO/SRA/ENA), GEO count-file audit, reference and pipeline research, samplesheet, docs.
- Key findings: data are **single-end 75 nt**; strandedness **unknown**; the GEO "count" table is keyed by **NCBI RefSeq R102 IDs (not ITAG2.4)** and looks **library-size normalized** → not usable for DESeq2.
- Not done (stop conditions): Nextflow/Java/R/tool installation, FASTQ or reference download, any quantification or DESeq2.
- Blocking decisions: reference version (D-001), quantification strategy (D-002).

## Reproduce the completed Phase 1 steps

From the project root, inside WSL (Ubuntu), with `bash`, `curl`, `awk`, `zcat`:

```bash
bash scripts/bash/env_audit.sh                      # read-only environment report
bash scripts/bash/fetch_metadata.sh                 # GEO / NCBI SRA / ENA records -> data/metadata/sources/
bash scripts/bash/build_samples_metadata.sh         # samples_metadata.tsv + config/samplesheet.csv (+ cross-source checks)
bash scripts/bash/fetch_geo_counts.sh               # 1.3 MB GEO supplementary table + structure summary
bash scripts/bash/verify_geo_counts_annotation.sh   # streams NCBI SL2.50 GFF (8 MB, not saved) and matches IDs
```

Remote records can change over time; the copies used for this phase are committed in `data/metadata/sources/`.

## Next steps

1. Approve or amend D-001 … D-008 in [docs/decisions.md](docs/decisions.md).
2. Clone to the WSL filesystem: `git clone "/mnt/c/Users/Cutza/Projects/Proyecto RNASeq Genómica" ~/projects/STV_GSE137303`.
3. Install Java + Nextflow as described in [docs/environment_audit.md](docs/environment_audit.md); then build the trivial Nextflow bootstrap test.
4. After approval: download FASTQ from ENA with MD5 checks, run FastQC/MultiQC, decide on trimming, infer strandedness.
