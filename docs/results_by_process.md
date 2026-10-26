# Results by process — STV_GSE137303

One section per process, with the key verified numbers and links to the detailed documents. Phases 1 – 4 are summarized retrospectively from their committed outputs. Phase 5 is a **diagnostic checkpoint** that did not change the primary analysis. Phase 6 adds a separate **pathway-level** analysis (preranked GSEA); it does not change the gene-level result. Phase 6B validates the Phase 6 pathway signals with pre-specified robustness checks; it changes neither result.

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
- **No gene-level leave-one-out run.** A DE analysis without STV_free_rep1 was not run (D-022). It would have n = 2 STV-free samples and could not replace the primary analysis. Phase 6B ran a **pathway-level** leave-one-sample-out as a labelled sensitivity diagnostic only (D-037, row 36).

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

## Phase 6B — functional robustness validation (Phase 4 and Phase 6 primary results unchanged)

Phase 6B asks which Phase 6 signals hold up under:
- a correlation-adjusted test (CAMERA);
- the Cook's ranking choice;
- replicate-level checks;
- single-sample removal;
- the SL4.0ch00 organelle-like genes.

Rules were frozen before any result: [functional_robustness_plan.md](functional_robustness_plan.md) §1 – §10 (frozen 2026-10-26T07:36:28Z). The outcome is in §11. Decisions D-033 – D-043. The 76 Phase 6 primary outputs are byte-identical after 6B (`sha256sum -c`). **No STOP condition was triggered.**

| # | Process | Input | Tool / version | Parameters | Reference pipeline / study | Quantitative result | Technical interpretation | Decision | Evidence / output files |
|---|---|---|---|---|---|---|---|---|---|
| 30 | Leading-edge summary, tracked reproduction (Task 1) | Phase 6 leading-edge tables | R 4.5.3, data.table | Descriptive | — | Exact reproduction: 2,858 genes (1,639 + / 1,219 −), median \|stat\| 0.704, max share 0.228, median LE fraction 0.417 / 0.387 / 0.581 | The numbers no longer depend on an untracked helper | D-041 | `scripts/R/09_*`; `results/functional_robustness/leading_edge/` |
| 31 | Cook's-filtered genes: inclusion sensitivity (D-023 final) | 18,277 finite Wald statistics | fgsea 1.36.2 | Phase 6 call and seed; sets rebuilt and identity-checked | Love et al. 2014; Cook 1977; Kelsey et al. 2025 | 210 / 215 sets Cook's-stable, 5 supported (padj < 0.25), 0 lost; 10 sensitivity-only; ρ(NES) ≥ 0.994. STV_free_rep1 has the largest Cook's distance in 19 of 39 genes | The ranking choice does not drive the Phase 6 calls | D-033 | `cooks_sensitivity/`; F5 `cooks_sensitivity_nes.png` |
| 32 | Independent expression model (Task 4) | Same six `quant.sf` (SHA-256 checked) | tximport 1.38.2 (`lengthScaledTPM`), edgeR 4.8.2 (TMM), limma 3.66.0 (voom, eBayes) | 18,238-gene universe; `~ condition`; coefficient 2 | Law et al. 2014, 2016; Soneson et al. 2015; Robinson & Oshlack 2010 | TMM 0.920 – 1.061; prior df 3.56; 0 genes at adj. P < 0.05; ρ(limma t, Wald) **0.963**, sign agreement 0.933 | Equivalent gene-level ranking; the model passes the equivalence check | D-034 | `camera/camera_summary.tsv`, `voom_mean_variance.png` |
| 33 | CAMERA, correlation-adjusted competitive test (Task 5) | voom object + Phase 6 GO sets | limma `camera` | Primary `inter.gene.cor = NA`; secondary 0.01 (descriptive); BH within collection | Wu & Smyth 2012; Goeman & Bühlmann 2007 | **0 sets at FDR < 0.05** in BP (952), MF (459), CC (237); minimum FDR 0.997 – 0.999. Median estimated ρ̄ in fgsea-significant sets 0.056 / 0.065 / 0.091 (median VIF 5.6 / 5.4 / 14.5). Secondary 0.01: 40 / 7 / 26 | With estimated correlation, no set-level shift is significant | D-034 | `camera/camera_go_*_all.tsv` |
| 34 | fgsea vs CAMERA; tiers (Task 6) | Rows 25 and 33 | R | Pre-specified tiers | GSEA User Guide (FDR 0.25) | Tier 1 **0**, Tier 2 **0**, Tier 3 **210**, Tier 3C 0, Tier 4 **5**. Direction agreement 100 % for fgsea-significant sets (89 % for all); ρ(signed significance) 0.923 / 0.914 / 0.893 | Same direction and similar ranking, but different significance calibration: every Phase 6 call is **method-dependent** | D-035 | `comparison/`; F1 `method_agreement_signed_significance.png` |
| 35 | Replicate consistency (Task 7) | Blind VST; leading edges | R | Pre-specified classes | — | 144 broadly consistent, 69 partially, 2 ambiguous, **0 strongly influenced by STV_free_rep1** | Most fgsea signals are carried by more than one replicate pair | D-036 | `replicate_consistency/`; F2 `replicate_consistency_*.png` |
| 36 | Leave-one-sample-out (Task 8; sensitivity only) | `txi_final.rds`; 6 five-sample fits | DESeq2 1.50.2, fgsea 1.36.2 | Genes fixed to 18,277; Phase 6 calls | Cook 1977; Schurch et al. 2016 | Retained + supported: 34 % (CC, drop STV_free_rep2) to 99 % (BP, drop STV_free_rep3). Maximum sample-dependent share **47.7 %** (CC, STV_infected_rep1); STV_free_rep1: 13.3 / 37.3 / 11.4 % | No single sample drives > 50 % of a collection, but several sets depend on one sample | D-037 | `loo/`; F3 `loo_retention.png` |
| 37 | SL4.0ch00 investigation (Task 9) | ITAG4.0 cDNA/proteins; NC_007898.3, NC_035963.1 | BLAST+ 2.17.0+ | Pre-specified classes; exclusion of 22 near-identical ranked ch00 genes | Camacho et al. 2009; Kahlau et al. 2006 | 186 / 513 ch00 genes near-identical to an organelle genome. ch00-sensitive sets 12.5 / 9.8 / 9.1 %; each of the 9 lost sets (ATP-synthase / proton-transport and related terms) contains 5 excluded genes. Post-hoc broader exclusion: 15.8 / 5.9 / 9.1 % | ATP-synthase / proton-transport signals depend on organelle-identical gene models (technical caveat) | D-038, D-043 | [ch00_functional_investigation.md](ch00_functional_investigation.md); `ch00/`; F4 |
| 38 | Theme classification (Task 14) | Rows 34 – 37 | R | Pre-specified regexes and rules | — | **All 10 themes EXPLORATORY** (no Tier 1/2 set anywhere) | No theme meets the ROBUST or SUPPORTED WITH CAUTION criteria | D-039 | `themes/theme_evidence_classification.tsv` |

### Phase 6B evidence table by theme (Task 10)

Columns from tracked tables:
- sets: significant in ≥ 1 method; + / − = NES sign; "+" means higher ranked in STV-infected leaves;
- replicate classes: broadly / partially / ambiguous;
- LOO dep.: sets dependent on STV_free_rep1;
- ch00: ch00-sensitive sets;
- Cook's: sets not retained.

Literature columns: what each source reports, with its directness (details in [biological_comparison_prior_studies.md](biological_comparison_prior_studies.md)). Literature is used for plausibility only and does not affect any class.

| Theme | Sets (+ / −) | Tiers | Min fgsea padj | Min CAMERA FDR | Replicate classes | LOO dep. | ch00 | Cook's | Class | STV literature (same dataset) | Other literature (indirect) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Ribosome / translation / RNA processing | 34 (34 / 0) | T3 33, T4 1 | 5.0 × 10⁻³⁰ | 0.997 | 0 / 32 / 2 | 6 | 0 | 0 | EXPLORATORY | No claim | — |
| Photosynthesis / plastid / thylakoid | 25 (25 / 0) | T3 25 | 4.2 × 10⁻¹¹ | 0.997 | 15 / 10 / 0 | 1 | 0 | 1 | EXPLORATORY | No claim | — |
| ATP synthase / proton transport | 5 (5 / 0) | T3 5 | 3.3 × 10⁻³ | 0.997 | 5 / 0 / 0 | 0 | **4** | 0 | EXPLORATORY | No claim | — (technical ch00 caveat) |
| Defence / biotic response | 8 (0 / 8) | T3 8 | 1.2 × 10⁻⁸ | 0.999 | 7 / 1 / 0 | 1 | 0 | 0 | EXPLORATORY | Badeli 2026: PR-protein locus and "response to fungus" (ORA p ≈ 0.1) lower in STV+, not significant | Whitham et al. 2003: defence/stress gene changes after acute RNA viruses (Arabidopsis); not comparable |
| Hormone signalling | 5 (0 / 5) | T3 4, T4 1 | 4.0 × 10⁻³ | 0.999 | 4 / 1 / 0 | 1 | 1 | 0 | EXPLORATORY | Fukuhara 2020 (abstract): ethylene biosynthesis/signalling genes down. Badeli 2026: ACO1 trend down, not significant | — |
| Endoplasmic reticulum / protein folding | 4 (0 / 4) | T3 4 | 1.9 × 10⁻⁴ | 0.997 | 1 / 3 / 0 | 4 | 0 | 0 | EXPLORATORY | No claim | — |
| Glutathione / detoxification | 3 (0 / 3) | T3 3 | 1.4 × 10⁻³ | 0.998 | 3 / 0 / 0 | 1 | 0 | 1 | EXPLORATORY | No claim | Khankhum et al. 2016: oxidation-reduction dominates the GO distribution of endornavirus-associated DEGs (bean; not an enrichment test) |
| Heme / oxidoreductase | 11 (5 / 6) | T3 11 | 3.0 × 10⁻⁴ | 0.997 | 9 / 2 / 0 | 4 | 1 | 0 | EXPLORATORY | No claim | Khankhum et al. 2016 (as above) |
| Water deprivation / abiotic stress | 3 (0 / 3) | T3 3 | 1.6 × 10⁻⁶ | 0.999 | 3 / 0 / 0 | 1 | 0 | 0 | EXPLORATORY | No claim | — |
| Carbohydrate binding / cell wall | 4 (2 / 2) | T3 4 | 1.7 × 10⁻⁴ | 0.998 | 4 / 0 / 0 | 0 | 1 | 1 | EXPLORATORY | Badeli 2026: β-galactosidase and PGIP loci in Table 2 (n.s.; ITAG4.0 names differ) | — |
| Other (no theme regex) | 113 (65 / 48) | T3 110, T4 3 | 3.6 × 10⁻²⁷ | 0.997 | 93 / 20 / 0 | 21 | 17 | 2 | — | — | — |

Persistent-virus context, applying to all themes and indirect: Roossinck 2010 and Takahashi et al. 2019 are reviews describing persistent/latent plant viruses as generally asymptomatic. Sabanadzovic et al. 2009 describe STV as seed-transmitted. No tomato–virus transcriptome source was cited (see the comparison document).

**Phase 6B statement for reports:**
- The Phase 6 pathway-level signals are robust to the Cook's ranking choice. They are mostly consistent across replicate pairs.
- They are partly sensitive to the removal of single samples. ATP-synthase / proton-transport terms depend on organelle-identical SL4.0ch00 gene models.
- **None** of the 215 fgsea-significant GO sets is significant under CAMERA with estimated inter-gene correlation. Directions agree.
- Under the pre-specified rules, every theme is therefore **EXPLORATORY**: a hypothesis for an adequately replicated experiment, not a finding.
- No gene is differentially expressed at FDR < 0.05.
