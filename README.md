# STV-RNAseq — reanalysis of GSE137303 (tomato M82 leaves, Southern tomato virus)

Reproducible RNA-seq reanalysis of the transcriptomic response of *Solanum lycopersicum* M82 to Southern tomato virus (STV) infection.

## Biological question

Which leaf genes of tomato cv. M82 differ in expression between plants infected with Southern tomato virus (STV) and STV-free plants? STV is a double-stranded RNA virus of the genus *Amalgavirus*. This repository re-quantifies the public raw reads with a current reference and a current, documented workflow. **Primary result:** no statistically significant differential expression was detected under the prespecified DESeq2 model and FDR threshold in this six-sample dataset (3 vs 3). This is not evidence that STV has no transcriptomic effect, and no biological conclusion is drawn here.

## Dataset

| | |
|---|---|
| GEO series | [GSE137303](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137303) |
| Publication | Fukuhara, Tabara, Koiwa & Takahashi (2020), *Archives of Virology* 165, 11–20, [doi:10.1007/s00705-019-04436-1](https://doi.org/10.1007/s00705-019-04436-1) (PubMed 31620899) |
| BioProject / SRA | PRJNA565137 / SRP221436; runs SRR10102791–SRR10102796 |
| Material | Leaves of *S. lycopersicum* cv. M82 — **3 STV-free vs 3 STV-infected** biological replicates |
| Sequencing | Illumina HiSeq 2500, **single-end 75 nt**, "TruSeq stranded mRNA" kit; 104.6 M reads in total |
| Raw-data source used here | ENA `fastq.gz` (5.7 GB), MD5-verified against ENA (6/6 PASS) |

| Condition | sample_id | Run |
|---|---|---|
| STV_free | STV_free_rep1–3 | SRR10102791–93 |
| STV_infected | STV_infected_rep1–3 | SRR10102794–96 |

Verified metadata: [data/metadata/samples_metadata.tsv](data/metadata/samples_metadata.tsv); audit: [docs/dataset_audit.md](docs/dataset_audit.md).

## Original study vs this reanalysis

| | Original study (GEO record) | This reanalysis |
|---|---|---|
| Genome / annotation | **SL2.50 / ITAG2.4-era analysis** (GEO text); the published count table is keyed by NCBI RefSeq Annotation Release 102 on SL2.50 ([docs/geo_counts_audit.md](docs/geo_counts_audit.md)) | **SGN SL4.0 genome + ITAG4.0 annotation** (D-001) |
| Read processing | TopHat/Bowtie2 alignment; Cuffnorm FPKM / normalized counts | Raw ENA FASTQ → FastQC/MultiQC → **Salmon** (selective alignment, decoy-aware index, `-l SR`, `--gcBias`) → **tximport** → **DESeq2** (Phase 4) |
| Strandedness | Not stated | Measured: all six libraries `SR` (reverse-stranded single-end) |
| Use of the authors' table | — | Audit/cross-check only; never DESeq2 input (normalized values, D-005) |

The reanalysis does not attempt to reproduce the authors' numbers gene by gene. It produces an independent, documented quantification on the current tomato reference.

## Workflow

```text
ENA FASTQ (MD5) → sanity checks (SeqKit) → FastQC → MultiQC        [Phase 2, done]
SGN SL4.0 + ITAG4.0 (SHA-256 manifest) → validation                [Phase 3, done]
→ Salmon 2.7.0 index (ITAG4.0 cDNA + full SL4.0 genome decoys)
→ exploratory Salmon quant (-l A) → orientation SR in 6/6          [Phase 3, done; archived for audit]
→ FINAL Salmon quant (-l SR, --gcBias, SE fragment-length prior 250/25) ×6 → QC   [Phase 3, done]
→ tx2gene (ITAG4.0 GFF) → tximport (gene level, validated)         [Phase 3, done]
→ DESeq2 (~ condition, padj < 0.05, apeglm shrinkage)               [Phase 4, done]
→ functional analysis                                              [not started]
```

- **Architecture:** a small Nextflow DSL2 workflow (`workflow/`, strict syntax, Nextflow 26.10.3) modelled on nf-core/rnaseq 3.27.0 (FastQC, Salmon index/quant, MultiQC modules; same tool versions). It is not a copy of nf-core/rnaseq. Statistics run outside the workflow as R scripts (`scripts/R/`).
- **Software:** pinned conda environments with lock files: `stv-rnaseq` (FastQC 0.13.0, MultiQC 1.35, Salmon 2.7.0, SeqKit 2.14.0, OpenJDK 21) and `stv-r` (R 4.5.3, tximport 1.38.2, DESeq2 1.50.2, apeglm 1.32.0; plotting: ggplot2, ggrepel, pheatmap, RColorBrewer).
- Details: [docs/nextflow_pipeline.md](docs/nextflow_pipeline.md), [docs/software_versions.md](docs/software_versions.md), [docs/methods_references.md](docs/methods_references.md).

## Repository structure

```text
config/            samplesheet.csv (accessions, condition, strandedness, fastq_1)
workflow/          Nextflow pipeline: main.nf, modules/, nextflow.config, bootstrap/, tests/stub/
scripts/bash/      metadata, download, sanity, QC/Salmon summaries, reference download/validation, tx2gene
scripts/R/         01_tximport_import.R (validation + gene-level import), 02_deseq2.R (primary DESeq2 analysis)
envs/              stv-rnaseq.yml, stv-r.yml + linux-64 lock files
data/metadata/     samples_metadata.tsv + raw source records (sources/)
data/raw/          ENA FASTQ (git-ignored)
data/processed/    original_geo/ authors' count table (audit only)
reference/         reference_manifest.tsv, reference_validation.tsv, ITAG4.0_tx2gene*.tsv (tracked);
                   SGN/ downloads (git-ignored)
results/           qc/*.tsv and quantification/*.tsv (tracked summaries); deseq2/*.tsv + plot_data/ (tracked);
                   FastQC/MultiQC reports, Salmon outputs, index, pipeline_info, *.rds (git-ignored)
figures/           deseq2/ PNG figures
docs/              audits, decisions, reference validation, Salmon/tximport reports, methods references, run log
```

## Reproducibility

Everything runs in WSL2/Linux from the project root. Large inputs and outputs are re-created by scripts, never committed.

```bash
# environments
CONDA_CHANNEL_PRIORITY=strict conda env create -f envs/stv-rnaseq.yml   # or: conda create -n stv-rnaseq --file envs/stv-rnaseq.linux-64.lock.txt
CONDA_CHANNEL_PRIORITY=strict conda env create -f envs/stv-r.yml        # or: conda create -n stv-r --file envs/stv-r.linux-64.lock.txt
conda activate stv-rnaseq
bash scripts/bash/install_nextflow.sh && conda deactivate && conda activate stv-rnaseq

# raw data and QC (Phase 2)
bash scripts/bash/download_fastq.sh && bash scripts/bash/check_raw_fastq.sh
nextflow run workflow/main.nf -resume && bash scripts/bash/summarize_fastqc.sh

# reference, index and quantification (Phase 3)
bash scripts/bash/download_reference.sh && bash scripts/bash/validate_reference.sh
# builds the index once and quantifies with the samplesheet strandedness (reverse -> -l SR) and --gcBias (default)
nextflow run workflow/main.nf -resume --run_quant true \
  --transcript_fasta reference/SGN/ITAG4.0/ITAG4.0_cDNA.fasta \
  --genome_fasta reference/SGN/SL4.0/S_lycopersicum_chromosomes.4.00.fa.gz \
  --multiqc_title 'STV_GSE137303 QC + Salmon final (SL4.0 / ITAG4.0, SR, gcBias)' --multiqc_dir qc/multiqc_salmon
# later re-runs: --salmon_index results/genome/index/salmon_index (exploratory settings: --salmon_libtype A --salmon_gc_bias false)
bash scripts/bash/summarize_salmon.sh && bash scripts/bash/build_tx2gene.sh
conda activate stv-r && Rscript scripts/R/01_tximport_import.R

# differential expression (Phase 4)
Rscript scripts/R/02_deseq2.R
```

Provenance:
- SHA-256 for every reference file: `reference/reference_manifest.tsv`.
- ENA MD5 for every FASTQ: `results/qc/raw_fastq_checksums.tsv`.
- SHA-256 of every `quant.sf` (exploratory, check and final runs): `results/quantification/quant_sf_sha256.tsv`.
- Nextflow trace/report: `results/pipeline_info/`.
- Every executed step: [docs/run_log.md](docs/run_log.md).
- Remote metadata can change; the copies used are committed in `data/metadata/sources/`.

## Status (2026-10-25)

| Step | Status |
|---|---|
| Raw FASTQ download + MD5, sanity checks, FastQC/MultiQC | Done; QC conclusion A, no trimming (D-009) — [docs/raw_qc_assessment.md](docs/raw_qc_assessment.md) |
| Reference SL4.0 + ITAG4.0: download, manifest, validation | Done; no material inconsistency — [docs/reference_validation.md](docs/reference_validation.md) |
| Salmon decoy-aware index + quantification (6 samples) | Done. Exploratory `-l A` run measured `SR` in 6/6 (archived). **Final run: explicit `-l SR` + `--gcBias`** (D-010, D-012); mapping 87.95–90.26 % (identical in both runs), no Salmon diagnostics — [docs/salmon_quantification.md](docs/salmon_quantification.md) |
| tx2gene + tximport import check (final run) | Done; 33,976 genes × 6 samples, totals = Salmon — [docs/tximport_preparation.md](docs/tximport_preparation.md) |
| DESeq2 primary analysis (D-013 – D-020) | Done. 18,277 genes tested after the prespecified filter; **0 genes with padj < 0.05; 0 genes with padj < 0.05 and \|log2FC\| ≥ 1**. No statistically significant differential expression was detected under the prespecified DESeq2 model and FDR threshold in this six-sample dataset. Limited power with n = 3 per condition — [docs/results_by_process.md](docs/results_by_process.md) |
| Functional enrichment | **Not started**; DEG thresholds are not relaxed and nominal p-values are not used for gene lists |
