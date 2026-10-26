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

## Intended workflow (D-002 approved; D-001 reference pending — see [docs/decisions.md](docs/decisions.md))

```text
ENA FASTQ (MD5-checked) → FastQC → MultiQC → [trimming only if QC justifies — not needed so far]
→ Salmon (decoy-aware index) → tximport → DESeq2 → plant functional analysis
```

Nextflow (`workflow/`) automates FASTQ → QC → quantification; tximport/DESeq2/functional analysis will run as R scripts.

## Working location

Active copy: **`~/projects/STV_GSE137303`** on the WSL filesystem (D-006). Activate the tools with `conda activate stv-rnaseq` ([docs/software_versions.md](docs/software_versions.md)).

## Directory structure

```text
config/            samplesheet.csv (generated; accessions + fastq_1)
workflow/          Nextflow DSL2 pipeline: main.nf, modules/, nextflow.config, bootstrap/, tests/stub/
scripts/bash/      audit, metadata, download, sanity-check, FastQC-summary, Nextflow-install scripts
scripts/R/         tximport / DESeq2 / functional analysis (later)
envs/              stv-rnaseq.yml (pinned) + linux-64 lock file
data/metadata/     samples_metadata.tsv + raw source records (sources/)
data/raw/          six ENA FASTQ.gz, MD5-validated (git-ignored)
data/processed/    original_geo/GSE137303_count.txt.gz (authors' table, audit only)
reference/         genome/annotation/indexes (git-ignored; nothing downloaded — D-001)
results/           qc/ (checksum, sanity and QC summary TSVs tracked; FastQC/MultiQC reports ignored),
                   quantification/, deseq2/, functional/, pipeline_info/
figures/           final figures (later)
docs/              audits, design, sources, decisions, validation, QC assessment, run log
```

## Current status (Phase 2, 2026-10-25)

- Phase 1: environment/dataset/GEO-table audits, reference and pipeline research.
- Phase 2 done: WSL-native clone; `stv-rnaseq` env (FastQC 0.13.0, MultiQC 1.35, Salmon 2.7.0, SeqKit, pigz, OpenJDK 21) + Nextflow 26.10.3; bootstrap and `-resume` validated; six FASTQ downloaded from ENA (5.72 GB) with **6/6 MD5 PASS**; sanity checks 6/6 PASS (read counts = ENA, all 75 nt); raw FastQC + MultiQC via the workflow; QC conclusion **A — no preprocessing appears necessary** ([docs/raw_qc_assessment.md](docs/raw_qc_assessment.md)).
- Not done (stop conditions): reference download, Salmon index/quantification, tximport, DESeq2, trimming, interpretation.
- Blocking: reference version **D-001**. To confirm: **D-009** (no trimming).

## Reproduce

Phase 1 (any bash with `curl`, `awk`, `zcat`):

```bash
bash scripts/bash/env_audit.sh                      # read-only environment report
bash scripts/bash/fetch_metadata.sh                 # GEO / NCBI SRA / ENA records -> data/metadata/sources/
bash scripts/bash/build_samples_metadata.sh         # samples_metadata.tsv + config/samplesheet.csv (+ cross-source checks)
bash scripts/bash/fetch_geo_counts.sh               # 1.3 MB GEO supplementary table + structure summary
bash scripts/bash/verify_geo_counts_annotation.sh   # streams NCBI SL2.50 GFF (8 MB, not saved) and matches IDs
```

Remote records can change over time; the copies used for this phase are committed in `data/metadata/sources/`.

Phase 2 (WSL, project root):

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create -f envs/stv-rnaseq.yml   # or: conda create -n stv-rnaseq --file envs/stv-rnaseq.linux-64.lock.txt
conda activate stv-rnaseq
bash scripts/bash/install_nextflow.sh && conda deactivate && conda activate stv-rnaseq
nextflow run workflow/bootstrap/main.nf -w work/bootstrap && nextflow run workflow/bootstrap/main.nf -w work/bootstrap -resume
bash scripts/bash/download_fastq.sh        # ENA download, resume, size+MD5 validation
bash scripts/bash/check_raw_fastq.sh       # gzip/structure/length/read-count checks
nextflow run workflow/main.nf -resume      # FastQC + MultiQC (Salmon off)
bash scripts/bash/summarize_fastqc.sh      # results/qc/raw_qc_summary.tsv
```

## Next steps

1. Decide **D-001** (reference: SL4.0+ITAG4.1 recommended; alternatives SL2.50+ITAG2.4 / SL2.50+NCBI R102) and confirm **D-009** (no trimming).
2. After D-001: download the approved genome + transcript FASTA + annotation into `reference/`, build the decoy-aware Salmon index (`--run_quant true`), run Salmon with `-l A`, inspect `lib_format_counts.json` of all six samples (D-004).
3. Then tximport → DESeq2 in `scripts/R/` (design and thresholds need approval).
