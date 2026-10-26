# Results by process — STV_GSE137303

One section per process, with the key verified numbers and links to the detailed documents. Phases 1 – 4 are summarized retrospectively from their committed outputs. Phase 5 is a **diagnostic checkpoint** that did not change the primary analysis. Phase 6 adds a separate **pathway-level** analysis (preranked GSEA); it does not change the gene-level result.

## Primary result (Phase 4; D-013 – D-019)

- **18,277 genes tested** after the prespecified filter (count ≥ 10 in ≥ 3 samples).
- **0 genes with padj < 0.05** (Benjamini–Hochberg, `alpha = 0.05`).
- **0 genes with padj < 0.05 and |log2FC| ≥ 1** (MLE log2FC).
- No statistically significant differential expression was detected under the prespecified DESeq2 model and FDR threshold in this six-sample dataset. This is **not** proof that STV has no transcriptomic effect.

## QC caveat (Phase 5; owner decision D-022)

STV_free_rep1 is consistently separated in expression-space diagnostics, but no independent technical failure or documented experimental covariate justifies its removal. All six samples are kept; the six-sample analysis remains the primary analysis.

## Limitations kept for later discussion

- **n = 3 biological replicates per condition.** Power to detect moderate changes is limited: with 3 replicates, Schurch et al. (2016) report that only 20 – 40 % of significantly DE genes are detected. Only the 0-gene outcome at the prespecified thresholds is reported; no relaxed or nominal-p gene list is used.
- **Observed sample structure.** STV_free_rep1 inflates within-group variance in the STV-free group. Each STV-free sample is closest to the STV-infected sample with the same replicate number, but this is undocumented and is not modelled.
- **Undocumented design factors.** Plant identity, growth batch, extraction date and library batch are unknown. The STV-free and STV-infected plants may come from different M82 sub-lines (a secondary, unverified statement in the 2026 reanalysis).
- **No leave-one-out run.** A sensitivity analysis without STV_free_rep1 was not run (D-022). If ever run, it would have n = 2 STV-free samples and could not replace the primary analysis.

## Results by process

| # | Process | Key result | Detail |
|---|---|---|---|
| 1 | Dataset and metadata (GEO GSE137303, PRJNA565137) | 6 single-end 75-nt leaf libraries, cv. M82, 1-month-old, 3 STV-free + 3 STV-infected; HiSeq 2500. **No documented pairing/blocking factor** | [dataset_audit.md](dataset_audit.md), [sample_structure_investigation.md](sample_structure_investigation.md) §2 |
| 2 | Raw-data acquisition (ENA) | 6/6 files MD5 = ENA; reads 14,318,362 – 19,862,460 (SeqKit = ENA = SRA) | [raw_data_validation.md](raw_data_validation.md) |
| 3 | Read QC (SeqKit, FastQC, MultiQC) | Q30 97 % in all; mean per-base Q 33.5 – 38.8; adapters ≤ 0.21 %; per-base content and duplication FAIL in all six (typical of RNA-seq). No trimming (D-009) | [raw_qc_assessment.md](raw_qc_assessment.md) |
| 4 | Sequencing layout (read names, Phase 5) | All six libraries on instrument SN679, run 271, flowcell HN3JNBCXY, lanes 1 + 2 (lane-1 fraction 0.496 – 0.498): **one sequencing batch** | `results/qc/sample_structure/fastq_header_lanes.tsv` |
| 5 | Reference (SGN SL4.0 + ITAG4.0) | 34,075 genes/mRNAs; decoy-aware Salmon index (k = 31) | [reference_validation.md](reference_validation.md), [reference_decision.md](reference_decision.md) |
| 6 | Quantification (Salmon 2.7.0, `SR`, `--gcBias`) | Mapping 87.95 – 90.26 %; decoy 8.45 – 10.57 %; SR-compatible ratio 0.9924 – 0.9943; no diagnostics | [salmon_quantification.md](salmon_quantification.md) |
| 7 | Gene-level import (tximport 1.38.2) | 33,976 genes × 6; gene totals = Salmon totals; 9,808 genes with zero counts in all samples | [tximport_preparation.md](tximport_preparation.md) |
| 8 | Pre-filter (D-014) | count ≥ 10 in ≥ 3 samples → 18,277 genes retained (53.79 %) | `results/deseq2/deseq2_summary.tsv` |
| 9 | DESeq2 model (D-013, D-015) | `~ condition`, reference STV_free; size factors (counts/avgTxLength) 0.848 – 1.212; parametric dispersion trend (asymptDisp 0.0835, extraPois 5.16); median final dispersion 0.088; 0 non-converged genes | [decisions.md](decisions.md) D-013 – D-015 |
| 10 | Testing (D-016) | 18,238 of the 18,277 tested genes have p-values (39 Cook's-filtered, cutoff 18.0); independent-filtering threshold at quantile 0 %; **padj < 0.05: 0 genes; padj < 0.05 and \|log2FC\| ≥ 1: 0 genes**; raw p < 0.05: 192 (descriptive only, never used as a gene list) | `results/deseq2/deseq2_results_all.tsv` |
| 11 | Shrinkage (apeglm, D-017) | Completed; all genes converged | `results/deseq2/deseq2_results_shrunken.tsv` |
| 12 | Sample QC (D-018) | VST PCA top 500: PC1 49.28 %, PC2 27.17 %; no separation by condition; hclust k = 2 = {STV_free_rep1} vs rest. Figures drawn with pheatmap/ggrepel/ggplot2 (visualization only, D-020) | `figures/deseq2/` |
| 13 | STV label support (Phase 5) | 0 STV reads in free libraries; 304 / 386 / 580 in infected (15.3 – 34.4 per million; 97 – 98 % of LC429302.1 covered). Labels supported | [sample_structure_investigation.md](sample_structure_investigation.md) §3 |
| 14 | Correlation / distance / PCA sensitivity (Phase 5) | STV_free_rep1 most remote (median distance 93.4 vs 68.9 – 76.0; median within-condition Pearson 0.928 vs 0.941 – 0.961); the same configuration in top 100/500/1000/all-gene PCA; every sample's nearest neighbour is from the other condition; same-number proximity observed but not documented | §4 |
| 15 | Cook's diagnostics (Phase 5) | STV_free_rep1 has the largest Cook's distance in 29.9 % of genes (expected 16.7 %); above the cutoff in 33 of the 39 filtered genes (top sample in 19). Filtering unchanged | §5 |
| 16 | Robust PCA (Phase 5) | `rrcov` defaults are degenerate at n = 6 (`PcaGrid` k = n; `PcaHubert` h = n); subtask stopped; no flag used | §6 |
| 17 | Upstream QC synthesis (Phase 5) | **No unique technical failure** of STV_free_rep1: small edge values (mapping −1.08 points, decoy +0.93 points vs next sample); same run/flowcell/lanes; all QC passed | §7 |
| 18 | GEO processed table — historical structure check (Phase 5) | Not used for DE inference. STV_free_rep1 most remote in both FPKM and Cuffnorm "count" columns; same nearest-neighbour pattern | §8 |
| 19 | External reanalysis (Phase 5) | Badeli et al. 2026 (preprint) independently saw STV−_rep1 separated on PC1; no outlier handling | [external_reanalysis_comparison.md](external_reanalysis_comparison.md) |
| 20 | Decision matrix and owner decision (Phase 5) | **Option A approved (D-022): all six samples kept**; the six-sample DESeq2 analysis remains primary. B not applicable (no documented covariate). C (leave-one-out) not run at this checkpoint: n_free would be 2 | §10 |

## Phase 6 — functional analysis at the pathway level (primary DESeq2 result unchanged)

The gene-level result above is **unchanged**: 18,277 genes tested; 0 with padj < 0.05; 0 with padj < 0.05 and \|log2FC\| ≥ 1; six samples kept.

Phase 6 asks a different, competitive question: are GO gene sets concentrated at either end of the full DESeq2 Wald ranking? It uses preranked GSEA without any gene threshold. Pre-specification, frozen before testing: [functional_analysis_plan.md](functional_analysis_plan.md) §1 – §7; decisions D-023 – D-032.

| # | Process | Input | Tool / version | Parameters | Reference pipeline / study | Quantitative result | Technical interpretation | Decision | Evidence / output files |
|---|---|---|---|---|---|---|---|---|---|
| 21 | Ranking (Task 1) | `results/deseq2/deseq2_results_all.tsv` (primary fit) | R 4.5.3, data.table; DESeq2 1.50.2 object for the direction check | Wald `stat`, STV_infected vs STV_free; universe = finite stat **and** DESeq2 p-value | Love et al. 2014; Kelsey et al. 2025 (TE-Seq); Tomaszewski et al. 2026; Li et al. 2026 | **18,238 genes ranked**; 39 Cook's-filtered excluded (stat −3.16 to 1.34); 0 duplicates, 0 ties, 0 NA/Inf; stat −4.405 to 2.864 (sd 0.668); 9,382 positive / 8,856 negative | Direction verified (DESeq2 description "condition STV infected vs STV free"; sign(stat) = sign(log2FC) for all genes) | D-023 | `results/functional/ranking/`; `figures/functional/ranking_wald_statistic.png` |
| 22 | GO annotation (Task 2) | SGN ITAG4.0 / ITAG4.1 GO files; PLAZA 5.0 `go.sly.csv.gz` | ontologyIndex 2.12; GO release 2026-08-05 | PLAZA primary sources, direct annotations; exact IDs; alt_id/obsolete mapping | Van Bel et al. 2022; GO docs | PLAZA primary: **82.88 %** of ranked genes; SGN ITAG4.0 exact: 2.52 %; SGN ITAG4.1 exact tx + identical protein: 45.57 %; PLAZA incl. projections 89.14 % (not used) | Identifiers exact (25,750/25,750 PLAZA genes are ITAG4.0 IDs); no version stripping; 192 obsolete terms (1,477 pairs) dropped | D-024, D-025 | `results/functional/annotation/`; `results/functional/annotation_sources_manifest.tsv` |
| 23 | GO gene sets (Task 3) | 103,133 direct gene–term pairs | ontologyIndex 2.12 | is_a + part_of propagation; roots excluded; ranked universe; **15 – 500 genes** | GSEA User Guide; fgsea vignette; Li et al. 2026 | Tested sets **BP 952, MF 459, CC 237**; ranked genes in tested sets: BP 49.7 %, MF 50.9 %, CC 22.4 % | Size range fixed before testing; identical-membership terms reported (BP 187, MF 82, CC 57) | D-025, D-026 | `annotation/go_gene_set_catalog.tsv`; `figures/functional/go_gene_set_sizes.png` |
| 24 | Plant Reactome (Task 4) | `gene_ids_by_pathway_and_species.tab` (current release) | R | Exact ITAG4.0 IDs only | Plant Reactome docs; Gupta et al. 2024 | 222 tomato pathways, 1,242 IDs; 733 exact ITAG4.0 IDs; **622 ranked genes (3.41 %)**; **5 pathways** with 15 – 500 ranked genes; 464 IDs differ only by version (not re-mapped) | Tomato pathways are orthology **projections** from curated rice pathways | **D-028: collection D not run (STOP)** | `results/functional/plant_reactome_gsea/` |
| 25 | Primary preranked GSEA (Task 5) | Ranking + GO sets | **fgsea 1.36.2** `fgseaMultilevel` | eps 0; nPermSimple 10,000; scoreType std; gseaParam 1; seed 20261025; SerialParam; BH within collection; **padj < 0.05** | Subramanian et al. 2005; Korotkevich et al.; GSEA User Guide | **BP 120** (75 positive / 45 negative NES); **MF 51** (21 / 30); **CC 44** (40 / 4) significant sets | Positive NES: ribosome/translation, photosynthesis, plastid/thylakoid. Negative NES: defence/stress-response, oxidoreductase-binding and ER terms (see plan §8) | D-027 | `results/functional/go_gsea/gsea_go_*_all.tsv`, `*_significant_fdr05.tsv` |
| 26 | Diagnostics (Task 6) | GSEA output | fgsea, R | Warnings classified in advance; rerun identity; seed 20261026 diagnostic; leading-edge composition | fgsea vignette/source; GSEA User Guide | 0 warnings, 0 NA, rerun identical; seed diagnostic: 15 BP calls differ (all primary padj 0.031 – 0.070), MF/CC 0; leading edges: 2,858 genes, **0 with gene padj < 0.05**, median \|stat\| 0.70, max single-gene share 0.23 | Set-level signals reflect many small coordinated shifts; the gene-permutation null ignores inter-gene correlation; borderline BP sets are seed-sensitive | D-027, D-031 | `gsea_go_summary.tsv`, `gsea_warnings.tsv`, `gsea_go_*_seed_diagnostic.tsv`, `leading_edge_*.tsv`; `figures/functional/gsea_go_*.png`, `enrichment_curves/` (29 curves by the fixed rule) |
| 27 | Post-hoc descriptive (not pre-specified) | Leading edges + ITAG4.0 GFF | R | Share of leading-edge genes on SL4.0ch00 | — | 15 of 2,858 leading-edge genes on ch00; ≥ 10 % of the leading edge in 27 of 215 significant sets (ATP-synthase / proton-transport terms) | Interpretive caveat only | — | `go_gsea/posthoc_leading_edge_ch00_share.tsv` |
| 28 | MapMan / Mercator4 preparation (Task 7) | ITAG4.0 proteins (unmodified), apeglm/MLE LFC, Wald stat | Mercator4 **v8.0** (web; site under maintenance on 2026-10-26) | Manual web submission documented; no ORA | Schwacke et al. 2019; Bolger et al. 2021 | FASTA 34,075 records (sha256 `84c15e…`); 32,522 terminal and 334 internal stops; expression table 18,277 rows (apeglm LFC −1.556 to 1.092) | Descriptive visualization only; colours are not significance | D-029 | `results/functional/mapman/` (incl. `MERCATOR4_SUBMISSION.md`) |
| 29 | agriGO feasibility (Task 8) | agriGO v2 manual/download pages | — | SEA, PAGE, custom annotation assessed | Tian et al. 2017; Kim & Volsky 2005 | ITAG4.0 background "user defined", GO 2016; custom upload ≤ 4 MB placed on agriGO's 2016 ontology; web-only, results kept for days | Not reproducible; redundant with fgsea | **D-030: conclusion C (not used)** | [agrigo_feasibility.md](agrigo_feasibility.md); `results/functional/agrigo/` |

**Phase 6 statement for reports:**
- No individual gene is differentially expressed at FDR < 0.05.
- At the pathway level, pre-specified preranked GSEA on the DESeq2 Wald statistic found GO gene sets positively enriched in the infected-vs-free ranking, notably ribosome/translation, photosynthesis and plastid terms.
- It found other sets negatively enriched, notably defence/stress-response, oxidoreductase-binding and endoplasmic-reticulum terms.
- These are coordinated, small shifts across many genes, assessed with a gene-permutation null. They are hypotheses for follow-up, not evidence of activation or inhibition. See the limits in [functional_analysis_plan.md](functional_analysis_plan.md) §8.4.
