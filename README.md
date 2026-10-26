# STV-RNAseq — reanalysis of GSE137303 (tomato M82 leaves, Southern tomato virus)

Reproducible RNA-seq reanalysis of the transcriptomic response of *Solanum lycopersicum* M82 to Southern tomato virus (STV) infection.

## Biological question

Which leaf genes of tomato cv. M82 differ in expression between plants infected with Southern tomato virus (STV) and STV-free plants? STV is a double-stranded RNA virus of the genus *Amalgavirus*. This repository re-quantifies the public raw reads with a current reference and a current, documented workflow. **Primary result:** no statistically significant differential expression was detected under the prespecified DESeq2 model and FDR threshold in this six-sample dataset (3 vs 3). This is not evidence that STV has no transcriptomic effect, and no biological conclusion is drawn here.

**Pathway level (Phase 6):** this is a separate, pre-specified question, tested with preranked GSEA on the DESeq2 Wald statistic.
- GO gene sets were found enriched toward both ends of the infected-vs-free ranking, at padj < 0.05:
  - positively: e.g. ribosome/translation, photosynthesis, plastid;
  - negatively: e.g. defence/stress-response, endoplasmic reticulum.
- These are coordinated small shifts across many genes, not differentially expressed genes. See [docs/functional_analysis_plan.md](docs/functional_analysis_plan.md).

**Robustness validation (Phase 6B):** pre-specified checks of the Phase 6 signals were run.
- The signals are robust to the Cook's ranking choice and mostly consistent across replicates.
- **None is significant under CAMERA with estimated inter-gene correlation** (same directions).
- ATP-synthase terms depend on organelle-identical genes on the unplaced SL4.0ch00.
- Every theme is therefore classified **EXPLORATORY**: a hypothesis, not a finding.

See [docs/functional_robustness_plan.md](docs/functional_robustness_plan.md) §11 and [docs/biological_comparison_prior_studies.md](docs/biological_comparison_prior_studies.md).

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
→ sample-structure / outlier diagnostics, primary analysis unchanged [Phase 5, done]
→ Wald-statistic ranking → GO (PLAZA 5.0, ITAG4.0) → fgsea preranked GSEA [Phase 6, done]
→ robustness validation: CAMERA, Cook's, replicate, leave-one-out, ch00 checks [Phase 6B, done]
→ closure: decision audit; organelle-aware quantification ASSESSED (not run)  [Phase 6B closure, done]
→ biological integration: evidence matrix, ethylene case study, GO redundancy, literature [Phase 7, done]
→ MapMan/Mercator4 inputs (manual web step pending); agriGO assessed, not used
```

- **Architecture:** a small Nextflow DSL2 workflow (`workflow/`, strict syntax, Nextflow 26.10.3) modelled on nf-core/rnaseq 3.27.0 (FastQC, Salmon index/quant, MultiQC modules; same tool versions). It is not a copy of nf-core/rnaseq. Statistics run outside the workflow as R scripts (`scripts/R/`).
- **Software:** pinned conda environments with lock files: `stv-rnaseq` (FastQC 0.13.0, MultiQC 1.35, Salmon 2.7.0, SeqKit 2.14.0, OpenJDK 21) and `stv-r` (R 4.5.3, tximport 1.38.2, DESeq2 1.50.2, apeglm 1.32.0; plotting: ggplot2, ggrepel, pheatmap, RColorBrewer).
- Phase 5 added `rrcov` 1.7-7 to `stv-r` for the robust-PCA diagnostic.
- Phase 6 added `fgsea` 1.36.2 and `ontologyIndex` 2.12 to `stv-r`.
- Phase 6B added `limma` 3.66.0 and `edgeR` 4.8.2 (+ `statmod`) to `stv-r`. BLAST+ 2.17.0+ was used read-only from an existing environment.
- Details: [docs/nextflow_pipeline.md](docs/nextflow_pipeline.md), [docs/software_versions.md](docs/software_versions.md), [docs/methods_references.md](docs/methods_references.md).

## Repository structure

```text
config/            samplesheet.csv (accessions, condition, strandedness, fastq_1)
workflow/          Nextflow pipeline: main.nf, modules/, nextflow.config, bootstrap/, tests/stub/
scripts/bash/      metadata, download, sanity, QC/Salmon summaries, reference download/validation, tx2gene
scripts/R/         01_tximport_import.R (validation + gene-level import), 02_deseq2.R (primary DESeq2 analysis)
                   03_sample_structure_diagnostics.R, 04_geo_historical_structure_check.R (Phase 5 diagnostics)
                   05_functional_ranking_annotation.R, 06_functional_gsea.R, 07_mapman_mercator_prep.R,
                   08_functional_posthoc_descriptives.R (Phase 6)
                   09 – 16 functional robustness scripts + functional_robustness_helpers.R (Phase 6B);
                   17_prior_study_locus_lookup.R, 18_functional_ch00_posthoc_identity_diagnostic.R (Phase 6B, post-hoc)
                   19_organelle_scope_diagnostic.R (Phase 6B closure, post-hoc; D-044)
                   20_biological_integration.R, 21_go_redundancy_modules.R, 22_phase7_figures.R (Phase 7)
envs/              stv-rnaseq.yml, stv-r.yml + linux-64 lock files
data/metadata/     samples_metadata.tsv + raw source records (sources/)
data/raw/          ENA FASTQ (git-ignored)
data/processed/    original_geo/ authors' count table (audit only)
reference/         reference_manifest.tsv, reference_validation.tsv, ITAG4.0_tx2gene*.tsv (tracked);
                   SGN/ and functional/ downloads (git-ignored)
results/           qc/*.tsv and quantification/*.tsv (tracked summaries); deseq2/*.tsv + plot_data/ (tracked);
                   FastQC/MultiQC reports, Salmon outputs, index, pipeline_info, *.rds (git-ignored)
                   qc/sample_structure/, qc/stv_signal/*.tsv (tracked Phase 5 diagnostics)
                   functional/ (Phase 6: ranking, annotation, go_gsea, plant_reactome_gsea, mapman, agrigo,
                   prespecification); functional/annotation/go_gene_sets.rds git-ignored
                   functional_robustness/ (Phase 6B: prespecification, leading_edge, cooks_sensitivity, camera,
                   comparison, replicate_consistency, loo, ch00, themes, prior_studies)
                   organelle_assessment/ (Phase 6B closure, D-044: prespecification, blast, provenance, scope tables)
                   biological_integration/ (Phase 7: theme matrix, ethylene case study, GO redundancy, literature search)
figures/           deseq2/ PNG figures
                   sample_structure/ PNG figures (Phase 5)
                   functional/ PNG figures and enrichment_curves/ (Phase 6)
                   functional_robustness/ PNG figures F1 – F5 (Phase 6B)
                   biological_integration/ PNG figures F1 – F5b (Phase 7)
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
# sample-structure diagnostics (Phase 5)
conda activate stv-rnaseq && bash scripts/bash/fastq_header_audit.sh && bash scripts/bash/stv_signal_check.sh  # Phase 5
conda activate stv-r && Rscript scripts/R/03_sample_structure_diagnostics.R && Rscript scripts/R/04_geo_historical_structure_check.R  # Phase 5
# functional analysis (Phase 6); 06 requires the frozen pre-specification in results/functional/prespecification/
bash scripts/bash/download_functional_annotation.sh
Rscript scripts/R/05_functional_ranking_annotation.R && Rscript scripts/R/06_functional_gsea.R
Rscript scripts/R/07_mapman_mercator_prep.R && Rscript scripts/R/08_functional_posthoc_descriptives.R
# robustness validation (Phase 6B); every script checks the frozen pre-specification in results/functional_robustness/prespecification/
for s in 09_functional_leading_edge_summary 10_functional_cooks_sensitivity 11_functional_camera 12_functional_method_comparison \
         13_functional_replicate_consistency; do Rscript scripts/R/$s.R; done
bash scripts/bash/download_organelle_reference.sh && bash scripts/bash/ch00_organelle_blast.sh   # BLAST+ from an existing env
for s in 15_functional_ch00_investigation 14_functional_loo_sensitivity 16_functional_robustness_figures_themes \
         17_prior_study_locus_lookup 18_functional_ch00_posthoc_identity_diagnostic; do Rscript scripts/R/$s.R; done
# Phase 6B closure (D-044; post-hoc, optional): organelle assessment; 19 checks the frozen criteria in results/organelle_assessment/prespecification/
bash scripts/bash/organelle_multihsp_blast.sh && Rscript scripts/R/19_organelle_scope_diagnostic.R
bash scripts/bash/ch00_ncbi_provenance.sh   # network: NCBI Datasets + E-utilities; system python3
# Phase 7 (integration of frozen results; no new test)
bash scripts/bash/literature_search_europepmc.sh   # network: Europe PMC; metadata only
Rscript scripts/R/20_biological_integration.R && Rscript scripts/R/21_go_redundancy_modules.R && Rscript scripts/R/22_phase7_figures.R
```

BLAST+ 2.17.0 (bioconda) is expected in a conda env named `rachis-qiime2-2026.4`, and conda is sourced from
`$HOME/miniconda3`. See [docs/software_versions.md](docs/software_versions.md). The full stage map, with inputs, outputs,
environments and manual steps, is in [docs/pipeline_execution_map.md](docs/pipeline_execution_map.md). Known limitations
are in [docs/reproducibility_audit.md](docs/reproducibility_audit.md).

Provenance:
- SHA-256 for every reference file: `reference/reference_manifest.tsv`.
- ENA MD5 for every FASTQ: `results/qc/raw_fastq_checksums.tsv`.
- SHA-256 of every `quant.sf` (exploratory, check and final runs): `results/quantification/quant_sf_sha256.tsv`.
- Nextflow trace/report: `results/pipeline_info/`.
- Every executed step: [docs/run_log.md](docs/run_log.md).
- Remote metadata can change; the copies used are committed in `data/metadata/sources/`.

## Status (2026-10-26)

| Step | Status |
|---|---|
| Raw FASTQ download + MD5, sanity checks, FastQC/MultiQC | Done; QC conclusion A, no trimming (D-009) — [docs/raw_qc_assessment.md](docs/raw_qc_assessment.md) |
| Reference SL4.0 + ITAG4.0: download, manifest, validation | Done; no material inconsistency — [docs/reference_validation.md](docs/reference_validation.md) |
| Salmon decoy-aware index + quantification (6 samples) | Done. Exploratory `-l A` run measured `SR` in 6/6 (archived). **Final run: explicit `-l SR` + `--gcBias`** (D-010, D-012); mapping 87.95–90.26 % (identical in both runs), no Salmon diagnostics — [docs/salmon_quantification.md](docs/salmon_quantification.md) |
| tx2gene + tximport import check (final run) | Done; 33,976 genes × 6 samples, totals = Salmon — [docs/tximport_preparation.md](docs/tximport_preparation.md) |
| DESeq2 primary analysis (D-013 – D-020) | Done. 18,277 genes tested after the prespecified filter; **0 genes with padj < 0.05; 0 genes with padj < 0.05 and \|log2FC\| ≥ 1**. No statistically significant differential expression was detected under the prespecified DESeq2 model and FDR threshold in this six-sample dataset. Limited power with n = 3 per condition — [docs/results_by_process.md](docs/results_by_process.md) |
| Sample-structure diagnostics (Phase 5) | Done. STV_free_rep1 is consistently separated in expression-space diagnostics, but no independent technical failure or documented experimental covariate justifies its removal; STV labels supported by STV reads; all six samples kept (D-022) — [docs/sample_structure_investigation.md](docs/sample_structure_investigation.md) |
| Functional analysis (Phase 6) | Done. Pre-specified preranked GSEA (fgsea 1.36.2) on the DESeq2 Wald ranking (18,238 genes) with PLAZA 5.0 GO for ITAG4.0: significant sets (padj < 0.05) BP 120, MF 51, CC 44. Plant Reactome not run (projected; 3.4 % coverage). MapMan/Mercator4 inputs prepared (manual web step pending). agriGO not suitable (conclusion C). DEG thresholds are not relaxed and nominal p-values are not used for gene lists — [docs/functional_analysis_plan.md](docs/functional_analysis_plan.md) |
| Functional robustness validation (Phase 6B) | Done; no STOP condition. CAMERA (estimated correlation) 0 significant sets → Tier 1 = Tier 2 = 0, Tier 3 = 210, Tier 4 = 5. Cook's-stable 210/215. Replicate classes: 144 broadly / 69 partially / 2 ambiguous; 0 driven by STV_free_rep1. Leave-one-out maximum sample dependence 47.7 %. ch00-sensitive ≤ 12.5 %. **All themes EXPLORATORY**. Phase 6 outputs byte-identical — [docs/functional_robustness_plan.md](docs/functional_robustness_plan.md) §11, [docs/ch00_functional_investigation.md](docs/ch00_functional_investigation.md), [docs/biological_comparison_prior_studies.md](docs/biological_comparison_prior_studies.md) |
| Phase 6B closure (decision audit, organelle assessment) | Done; **no result changed**, nothing re-quantified. Decision audit D-033 – D-044: 3 APPROVE, 7 APPROVE WITH CAVEAT; D-042 resolved by the owner as REGENERATE ON DEMAND. SL4.0ch00 (a Salmon decoy) already contains the complete plastome and mitogenome; **Option A approved (D-044): no organelle-aware re-quantification**; no gene-exclusion fgsea. Mercator4 still under maintenance. Committed `13e87be` — [docs/phase6b_decision_audit.md](docs/phase6b_decision_audit.md), [docs/organelle_requantification_assessment.md](docs/organelle_requantification_assessment.md), [docs/reproducibility_audit.md](docs/reproducibility_audit.md), [docs/pipeline_execution_map.md](docs/pipeline_execution_map.md) |
| Biological integration (Phase 7) | Done; frozen results only, no new test. Gene level unchanged (**0 genes padj < 0.05**). Ten themes: **7 EXPLORATORY** (photosynthesis/plastid, defence, hormone/ethylene, heme, glutathione, water/abiotic stress, cell wall) and **3 UNSUPPORTED** (ribosome/translation: single-sample dependence; ER: STV_free_rep1 dependence; ATP synthase: organelle-identical features). Ethylene: lower-rank directional tendency (same direction as Fukuhara et al. and Badeli et al.), not significant at gene level (ACO1 padj 0.297) and not CAMERA-supported. The 215 GO sets form about 73 redundancy modules. Mercator4 still under maintenance — [docs/phase7_biological_interpretation.md](docs/phase7_biological_interpretation.md), [docs/biological_integration.md](docs/biological_integration.md), [docs/literature_biological_context.md](docs/literature_biological_context.md) |
