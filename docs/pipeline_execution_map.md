# Pipeline execution map

Each stage lists the script, its inputs and outputs, the environment, whether it is automatic or manual, and whether it is
required to reproduce the primary results. Run everything from the project root in WSL2/Linux. Status: Phase 6B closure,
2026-10-26.

**Required** means needed for the primary chain:
- FASTQ → `quant.sf` → DESeq2 (Phase 4 result);
- → ranking → fgsea (Phase 6 result);
- → robustness classes (Phase 6B).

**Optional** means diagnostic, audit or post-hoc. Nothing downstream depends on it, unless stated.

| # | Stage | Script / command | Main inputs | Main outputs (tracked unless marked *ign.*) | Env | Mode | Required? | Decisions |
|---|---|---|---|---|---|---|---|---|
| 0 | Environments | `conda env create -f envs/*.yml` (or lock files); `scripts/bash/install_nextflow.sh` | `envs/` | conda envs `stv-rnaseq`, `stv-r` | — | auto | Required | D-006, D-008 |
| 0b | Environment audit | `scripts/bash/env_audit.sh` | — | stdout | any | auto | Optional | — |
| 1 | Metadata | `fetch_metadata.sh` → `build_samples_metadata.sh` | GEO / SRA / ENA | `data/metadata/`, `config/samplesheet.csv` | stv-rnaseq | auto (network) | Required | — |
| 1b | GEO table audit | `fetch_geo_counts.sh`, `verify_geo_counts_annotation.sh` | GEO supplementary | `data/processed/original_geo/` | stv-rnaseq | auto (network) | Optional | D-005 |
| 2 | Raw data | `download_fastq.sh` → `check_raw_fastq.sh` | `samples_metadata.tsv` | `data/raw/` (*ign.*), `results/qc/raw_fastq_checksums.tsv` | stv-rnaseq | auto (network) | Required | D-007 |
| 3 | Read QC | `nextflow run workflow/main.nf -resume` → `summarize_fastqc.sh` | FASTQ, samplesheet | `results/qc/fastqc_raw/`, `multiqc_raw/` (*ign.*); `results/qc/*.tsv` | stv-rnaseq | auto | Required (QC gate) | D-003, D-009 |
| 4 | Reference | `download_reference.sh` → `validate_reference.sh` | SGN URLs | `reference/SGN/` (*ign.*), `reference/reference_manifest.tsv`, `reference_validation.tsv` | stv-rnaseq | auto (network) | Required | D-001 |
| 5 | Index + quantification | `nextflow run workflow/main.nf -resume --run_quant true --transcript_fasta … --genome_fasta …` → `summarize_salmon.sh` | ITAG4.0 cDNA, SL4.0 genome (decoys), FASTQ | `results/genome/index/` and `results/quantification/salmon/` (*ign.*); `quant_sf_sha256.tsv`, Salmon summaries | stv-rnaseq | auto | Required | D-001, D-002, D-004, D-010 – D-012 |
| 5b | Run comparison | `compare_salmon_runs.sh` | two Salmon result sets | comparison TSV | stv-rnaseq | auto | Optional | — |
| 6 | tx2gene | `build_tx2gene.sh` | ITAG4.0 GFF | `reference/ITAG4.0_tx2gene*.tsv` | stv-rnaseq | auto | Required | — |
| 7 | Gene-level import | `Rscript scripts/R/01_tximport_import.R` | `quant.sf`, tx2gene | `txi_final.rds` (*ign.*), import checks | stv-r | auto | Required | D-002 |
| 8 | DESeq2 | `02_deseq2.R` | `txi_final.rds` | `results/deseq2/*.tsv`, `deseq2_objects.rds` (*ign.*), `figures/deseq2/` | stv-r | auto | **Required (primary result)** | D-013 – D-020 |
| 9 | Sample-structure diagnostics | `fastq_header_audit.sh`, `stv_signal_check.sh` (stv-rnaseq); `03_*.R`, `04_*.R` (stv-r) | FASTQ, DESeq2 objects, GEO table | `results/qc/sample_structure/`, `results/qc/stv_signal/` | both | auto | Optional (diagnostic; D-022 keeps all samples) | D-022 |
| 10 | Functional resources | `download_functional_annotation.sh` | PLAZA 5.0, GO OBO, Plant Reactome | `reference/functional/` (*ign.*), `annotation_sources_manifest.tsv` | stv-rnaseq | auto (network) | Required | D-024, D-025 |
| 11 | Ranking + GO sets | `05_functional_ranking_annotation.R` | DESeq2 results / objects, PLAZA, GO | `results/functional/ranking/`, `annotation/`, `go_gene_sets.rds` (*ign.*) | stv-r | auto | Required | D-023 – D-026 |
| 12 | Preranked GSEA | `06_functional_gsea.R` (checks the Phase 6 freeze) | ranking, GO sets | `results/functional/go_gsea/`, `figures/functional/` | stv-r | auto | **Required (Phase 6 result)** | D-027, D-031 |
| 13 | MapMan inputs | `07_mapman_mercator_prep.R` | ITAG4.0 proteins, DESeq2 | `results/functional/mapman/` | stv-r | auto | Optional | D-029 |
| 13b | Mercator4 annotation | Web submission (`MERCATOR4_SUBMISSION.md`) | `mapman/` FASTA | — | web | **manual; pending** (site under maintenance 2026-10-26) | Optional | D-029 |
| 14 | Post-hoc descriptives | `08_functional_posthoc_descriptives.R` | leading edges | `go_gsea/posthoc_*` | stv-r | auto | Optional | — |
| 15 | Leading-edge summary | `09_functional_leading_edge_summary.R` (checks the 6B freeze; same for 10 – 16) | Phase 6 leading edges | `functional_robustness/leading_edge/` | stv-r | auto | Required (6B) | D-041 |
| 16 | Cook's sensitivity | `10_functional_cooks_sensitivity.R` | DESeq2 objects, GO sets | `cooks_sensitivity/` | stv-r | auto | Required (6B) | D-033 |
| 17 | CAMERA | `11_functional_camera.R` | `quant.sf` (SHA-256 checked), tx2gene, GO sets | `camera/` (voom table *ign.*) | stv-r | auto | Required (6B) | D-034, D-042 |
| 18 | Method comparison | `12_functional_method_comparison.R` | fgsea + CAMERA tables | `comparison/` | stv-r | auto | Required (6B) | D-035 |
| 19 | Replicate consistency | `13_functional_replicate_consistency.R` | DESeq2 objects, GO sets, comparison | `replicate_consistency/` | stv-r | auto | Required (6B) | D-036 |
| 20 | Organelle reference + BLAST | `download_organelle_reference.sh` → `ch00_organelle_blast.sh` | NCBI RefSeq, ITAG4.0 | `reference/functional/organelle/` (*ign.*), `ch00/` BLAST tables | stv-rnaseq + external BLAST env | auto (network) | Required (6B) | D-038 |
| 21 | ch00 investigation | `15_functional_ch00_investigation.R` | BLAST tables, GFF, ranking, GO sets | `ch00/` | stv-r | auto | Required (6B) | D-038 |
| 22 | Leave-one-out | `14_functional_loo_sensitivity.R` | `txi_final.rds` | `loo/` | stv-r | auto | Required (6B) | D-037 |
| 23 | Themes + figures | `16_functional_robustness_figures_themes.R` | rows 15 – 22 | `themes/`, `figures/functional_robustness/` | stv-r | auto | Required (6B) | D-039, D-040 |
| 24 | Prior-study lookup | `17_prior_study_locus_lookup.R` | DESeq2, ranking, themes | `prior_studies/` | stv-r | auto | Optional (post-hoc) | D-043 |
| 25 | ch00 post-hoc exclusion | `18_functional_ch00_posthoc_identity_diagnostic.R` | `ch00/`, ranking, GO sets | `ch00/posthoc_identity98_*` | stv-r | auto | Optional (post-hoc) | D-043 |
| 26 | Organelle assessment freeze | `cp` + `sha256sum` (see `run_log.md` entry 115) | `docs/organelle_requantification_assessment.md` §1 | `results/organelle_assessment/prespecification/` | — | manual (one-off) | Optional (D-044) | D-044 |
| 27 | Multi-HSP BLAST | `organelle_multihsp_blast.sh` | organelle FASTA (SHA-256 checked), ITAG4.0 cDNA, SL4.0 genome | `organelle_assessment/blast/` | external BLAST env | auto | Optional (post-hoc) | D-044 |
| 28 | Organelle scope diagnostic | `19_organelle_scope_diagnostic.R` (checks the D-044 freeze) | BLAST tables, DESeq2, leading edges, 6B tables, `quant.sf` (read-only), `duplicate_clusters.tsv` | `organelle_assessment/*.tsv` | stv-r | auto | Optional (post-hoc) | D-044 |
| 29 | GenBank provenance | `ch00_ncbi_provenance.sh` | NCBI Datasets / E-utilities, SL4.0 genome, row 28 output | `organelle_assessment/provenance/` | system python3, curl | auto (network) | Optional (post-hoc) | D-044 |
| 30 | Literature search (Phase 7) | `literature_search_europepmc.sh` | Europe PMC REST | `biological_integration/literature/` (metadata) | system python3 | auto (network) + manual selection | Optional (interpretation only) | D-046 |
| 31 | Evidence matrix and deep dives (Phase 7) | `20_biological_integration.R` | Phase 6 / 6B / D-044 tables, `go_gene_sets.rds`, `deseq2_objects.rds`, `curated_theme_literature.tsv` (manual) | `biological_integration/*.tsv` | stv-r | auto | Required for Phase 7 interpretation | D-045 |
| 32 | GO redundancy (Phase 7) | `21_go_redundancy_modules.R` | row 31 output, `go_gene_sets.rds`, GO OBO | `go_redundancy_*.tsv`, F5, F5b | stv-r | auto | Optional (descriptive) | D-047 |
| 33 | Phase 7 figures | `22_phase7_figures.R` | rows 31 – 32 | `figures/biological_integration/F1 – F4` | stv-r | auto | Optional | D-048 |

## Dependency outline

```text
1 → 2 → 3 ─┐
4 ─────────┴→ 5 → 6 → 7 → 8 ──→ 11 → 12 ──→ 15, 16, 17 → 18 → 19
                     │     │          │       20 → 21 ─────────────┐
                     │     └→ 9       └→ 13 (→ 13b manual), 14     22 ──→ 23 → 24, 25
                     └──────────────────────────────── 17 (quant.sf), 28 (quant.sf)
                                             26 → 27 → 28 → 29
          (18 – 25, 28) + 30 (manual literature summary) → 31 → 32 → 33        [Phase 7]
```

Not part of any stage: the Windows/WSL editing helpers (`.stv_tools`, staging). They are not needed to reproduce any
result; see [reproducibility_audit.md](reproducibility_audit.md) R11.
