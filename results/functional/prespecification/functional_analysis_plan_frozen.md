# Functional analysis plan — Phase 6 (pre-specified)

**Status.** Sections 1 – 7 were written on 2026-10-26 **before any gene-set enrichment test was run**. They were frozen and checksummed in `docs/run_log.md` (Phase 6) before `scripts/R/06_functional_gsea.R` was executed. Section 8 (outcome) was appended after the run and does not change sections 1 – 7.

Decisions: `docs/decisions.md` D-023 – D-032. Methods and references: `docs/methods_references.md` §25 – §33. agriGO assessment: [`agrigo_feasibility.md`](agrigo_feasibility.md).

---

## 1. Two statistical questions, kept separate

| | Gene level (Phase 4, **unchanged**) | Pathway level (this phase) |
|---|---|---|
| Question | Is gene *g* differentially expressed between STV_infected and STV_free? | Are the genes of a set concentrated at one end of the **full** infected-vs-free ranking more than random sets of the same size? (competitive, preranked GSEA) |
| Method | DESeq2 Wald test, BH, `alpha = 0.05` (D-013 – D-019) | fgsea multilevel preranked GSEA on the DESeq2 Wald statistic, BH within each collection (D-023 – D-027) |
| Result | 18,277 genes tested; **0** genes with padj < 0.05; **0** with padj < 0.05 and \|log2FC\| ≥ 1; all six samples kept (D-022) | see §8 |

A significant gene set does **not** make any of its genes differentially expressed, and it does not change the Phase 4 result. A non-significant set is not evidence of no biological effect (n = 3 per group).

## 2. Why classical DEG-list enrichment is not used

Over-representation analysis (ORA; agriGO "SEA", Mercator4/MapMan bin ORA, clusterProfiler `enrichGO`) needs a list of differentially expressed genes. The primary DESeq2 analysis produced **zero FDR-significant DEGs**. Relaxing the individual-gene threshold only to obtain an enrichment list would violate the prespecified analysis. The 192 genes with raw p < 0.05 are therefore never used as an input list: they have no multiplicity control, and at 18,238 tests about 912 raw p < 0.05 are expected by chance alone. Preranked GSEA uses every ranked gene and needs no threshold, so it answers the pathway-level question without redefining differential expression.

## 3. Original project plan vs final evidence-based roles

| Tool / resource | Original plan | Final role (evidence-based) | Phase 6 status | Evidence |
|---|---|---|---|---|
| GSEA (fgsea) | — | **Primary threshold-free statistical pathway analysis** on GO BP, MF, CC | Run (§4) | D-023 – D-027 |
| GO annotation | SGN ITAG4.0 `goterms` file | GO gene sets from **PLAZA 5.0** (ITAG4.0 gene models, exact IDs, primary sources only). The SGN files fail on version or coverage (§6) | Built | D-024, D-025 |
| Plant Reactome | Pathway enrichment | Plant/tomato pathway gene-set source and pathway interpretation, **only if validated**. Tomato pathways are orthology **projections** from curated rice pathways. Exact mapping covers 3.41 % of ranked genes; only 5 pathways have 15 – 500 ranked genes | **Not run as GSEA** (STOP for collection D); audit and exact mapping tables kept | D-028 |
| MapMan / Mercator4 | Functional bins + DEG enrichment | Plant-specific annotation and **descriptive** expression visualization | Inputs prepared; Mercator4 web submission is a documented manual step; no ORA | D-029 |
| agriGO v2 | GO enrichment (SEA) | Complementary GO analysis **only if** PAGE/ID compatibility is validated | Conclusion **C**: not suitable for this dataset/version; not run | D-030; [`agrigo_feasibility.md`](agrigo_feasibility.md) |
| clusterProfiler | — | Cited as a GSEA/GO framework (Yu 2012; Wu 2021; Xu 2024). Not used: fgsea is called directly so that every parameter is explicit | Not used | D-027 |

## 4. Pre-specified analysis (frozen before testing)

### 4.1 Ranking (D-023)

- Metric: DESeq2 Wald statistic `stat` of the primary six-sample fit (`results/deseq2/deseq2_results_all.tsv`), contrast **STV_infected vs STV_free**. `stat = log2FC_MLE / lfcSE`, so **positive = higher in STV_infected**.
  - Verified against the stored DESeq2 object: description "Wald statistic: condition STV infected vs STV free"; identical values; sign(stat) = sign(MLE log2FC) for all 18,277 genes.
- Universe: the **18,238** genes with a finite statistic **and** a DESeq2 p-value.
  - **39 genes excluded**: DESeq2 set their p-value to NA because a sample exceeded the Cook's-distance cutoff (count outlier). DESeq2 withholds inference for them, so their statistic is not a valid basis for a test. Their statistics (−3.16 to 1.34) and the positions they would have taken are listed in `results/functional/ranking/ranking_excluded_genes.tsv`.
  - Excluded: none for NA/Inf statistics (0).
- Diagnostics (no enrichment involved):
  - duplicate IDs 0; ties 0 (fgsea definition); statistic range −4.405 to 2.864; sd 0.668.
  - 9,382 positive and 8,856 negative statistics.
  - \|stat\| > 3 in 19 genes vs 49.2 expected under N(0,1); \|stat\| > 4 in 3 genes vs 1.2 expected.
  - Not used: padj ranking, significant-only ranking, apeglm LFC (no reference establishes a reason to replace the Wald statistic). No jitter or tie-breaking is needed (0 ties).

### 4.2 Gene sets (D-024, D-025, D-026)

- **Source**: PLAZA 5.0 dicots, `go.sly.csv.gz` (instance `dicots_05`, generated 2021-09-20). Assembly/annotation `ITAG4.0`; provider solgenomics.net.
  - All 25,750 PLAZA gene IDs are exact ITAG4.0 gene IDs. **No identifier version was stripped.**
- **Annotations used**: `go_source = primary` rows with `is_propagated_from_child_term = 0` (direct annotations): UniProt-GOA (IEA, IBA and experimental codes), InterProScan IEA, and GO-Consortium imports.
  - Not used: PLAZA orthology / gene-family / tree projections (`from_iortho`, `from_family_enrichment`, `from_tree`). They transfer GO from other species, and cross-species transfer is not silently accepted.
  - No `NOT` qualifiers are present.
- **Ontology**: GO `go-basic.obo`, dated release 2026-08-05 (`data-version: releases/2026-07-26`). PLAZA term IDs were mapped to that release:
  - `alt_id` → primary term: 41 terms;
  - obsolete term with a single `replaced_by` → its replacement: 220 terms;
  - obsolete without replacement → dropped: 192 terms, 1,477 gene–term pairs; listed in `go_term_version_mapping.tsv`.
- **Propagation**: to all ancestors over `is_a` and `part_of` only (GO docs: safe for grouping; `regulates` is not). No namespace is crossed (checked).
- **Collections**: BP, MF and CC separately.
  - The three root terms are excluded.
  - Sets are restricted to the 18,238 ranked genes.
  - Genes without GO stay in the ranked list, as in GSEA.
- **Redundancy**: all terms are kept, as in standard GO GSEA. Terms with identical gene membership are reported but not merged (BP 187, MF 82, CC 57 of the tested terms).
- **Size range (D-026): 15 – 500 genes after restriction to the ranked universe**, fixed before any test.
  - GSEA User Guide: "by default, GSEA ignores gene sets that contain fewer than 15 genes or more than 500 genes; defaults that are appropriate for datasets with 10,000 to 20,000 features" (here 18,238).
  - fgsea vignette: `minSize = 15, maxSize = 500`.
  - Plant precedent: Li et al. 2026 (potato, DESeq2 Wald ranking) used 10 – 500.
- **Sets tested**: BP **952**, MF **459**, CC **237**.

### 4.3 Test (D-027)

| Parameter | Value | Why |
|---|---|---|
| Implementation | **fgsea 1.36.2** (Bioconductor 3.22), `fgseaMultilevel()` (what `fgsea()` dispatches to) | fgsea vignette: adaptive multi-level split Monte-Carlo p-values, arbitrarily small p-values; `fgsea()` warns users of `fgseaSimple` that multilevel is recommended |
| `minSize`, `maxSize` | 15, 500 | §4.2 |
| `scoreType` | `"std"` | Two-sided; the statistic has both signs |
| `gseaParam` | 1 | Weighted running sum (GSEA default, p = 1); the Wald statistic magnitude is meaningful |
| `eps` | 0 | Exact p-value estimation (vignette: set `eps = 0` "to estimate P-value more accurately") |
| `nPermSimple` | 10,000 | Value suggested by fgsea's own warning for unbalanced statistics; fixed in advance, not tuned |
| `sampleSize` | 101 (default) | — |
| Seed / parallelism | `set.seed(20261025)` before each collection; `BiocParallel::SerialParam()` | Deterministic; checked by re-running each collection and requiring identical output |
| Multiple testing | BH within each collection (fgsea `padj`) | Collections run separately |
| **Criterion** | **padj < 0.05** | Prespecified; GSEA User Guide: with few samples and gene-set permutation, "use a more stringent FDR cutoff, such as 5%". **Not changed after results** |

**Wording.**
- Positive NES = "positively enriched in the infected-vs-free ranking" (set members tend toward the STV_infected end).
- Negative NES = "negatively enriched in the infected-vs-free ranking" (toward the STV_free end).
- No set is called "activated" or "inhibited".

### 4.4 Warnings (classified in advance)

| fgsea warning | Classification |
|---|---|
| "There are ties in the preranked stats" | Not expected (0 ties). If raised: material → STOP and report |
| "…P-values were not calculated properly due to unbalanced … gene-level statistic values … set to NA" | **Material** → STOP and report (results unusable for those sets) |
| "For some of the pathways the P-values were likely overestimated…" | Informational: the affected p-values are conservative (log2err = NA). Counted and reported |
| Any other warning or error | Material → STOP and report |

### 4.5 Diagnostics (all collections, complete results)

- Ranking distribution (`figures/functional/ranking_wald_statistic.png`).
- Set-size distributions (`figures/functional/go_gene_set_sizes.png`).
- Nominal p-value histograms.
- NES distribution and positive/negative counts.
- padj distribution; `log2err` NA count.
- Run-twice reproducibility check.
- **Seed diagnostic**: one extra run with seed 20261026, counting how many significance calls change. It is diagnostic only and never replaces the primary seed.
- Leading-edge composition of significant sets:
  - size, and fraction of the set;
  - median and maximum \|stat\|;
  - share of the largest gene;
  - pairwise Jaccard overlap.

### 4.6 Figure-selection criteria (fixed before any pathway identity is seen; D-031)

1. Always produced:
   - ranking distribution;
   - set-size histograms;
   - per-collection NES vs −log10(padj) for **all** tested sets (no text labels);
   - per-collection nominal p-value histograms.
2. Enrichment curves (`fgsea::plotEnrichment`) **only for FDR-significant sets** (padj < 0.05).
   - At most 5 per collection and direction, ordered by padj ascending, then \|NES\| descending, then GO ID.
   - If no set is significant, no enrichment curve is drawn and this is reported.
3. No figure is made for a set chosen by name or by expected biology. Non-significant sets appear only in the complete tables and the overview plots.

### 4.7 STOP conditions applied

- GO coverage of the ranked universe unexpectedly poor.
- Ambiguous or version-mismatched identifiers.
- Plant Reactome poor or undocumented.
- Ambiguous duplicates.
- Material fgsea warnings.
- Too few genes or sets.
- Ranking direction inconsistent with the contrast.

At the time of freezing:
- GO passes: 82.88 % of ranked genes have ≥ 1 usable direct GO annotation; exact IDs; no duplicates.
- Plant Reactome **fails** (projected, 3.41 % coverage, 5 testable pathways), so collection D is stopped.
- Ranking direction passes.

## 5. Reference supporting each analytical choice

| Choice | Reference(s) |
|---|---|
| Preranked, threshold-free gene-set test; ES, NES, leading edge | Subramanian et al. 2005 (10.1073/pnas.0506580102); GSEA User Guide |
| DESeq2 Wald statistic as the gene-level metric | Love et al. 2014 (10.1186/s13059-014-0550-8). Workflow: Kelsey, Kalekar & Sedivy 2025, TE-Seq, *Mobile DNA* (10.1186/s13100-025-00381-w), "ranked by their Wald statistic (produced by DESeq2)" |
| Plant precedent for signed DESeq2 Wald ranking + GSEA | Tomaszewski et al. 2026, *Front Plant Sci* (10.3389/fpls.2026.1733427; fgsea, Solanaceae); Li et al. 2026, *Plants* (10.3390/plants15182842; PSTVd-infected potato, 10 – 500 genes, BH q < 0.05) |
| fgsea multilevel, `eps = 0`, 15 – 500 | Korotkevich et al., bioRxiv 10.1101/060012; fgsea 1.36/1.38 vignette and source |
| Size range 15 – 500; FDR 5 % for gene-set permutation; no duplicate ranking values; avoid near-duplicate sets | GSEA User Guide (docs.gsea-msigdb.org) |
| BH adjustment | Benjamini & Hochberg 1995 (10.1111/j.2517-6161.1995.tb02031.x) |
| Competitive tests ignore inter-gene correlation (interpretation caveat) | Tamayo et al. 2016 (10.1177/0962280212460441); Wu & Smyth 2012 (10.1093/nar/gks461) |
| GO resource, true-path propagation over `is_a` / `part_of`, go-basic | Ashburner et al. 2000 (10.1038/75556); GO Consortium 2023 (10.1093/genetics/iyad031); GO docs "Relations in the Gene Ontology", "Download the ontology" |
| OBO parsing / ancestors | ontologyIndex 2.12; Greene, Richardson & Turro 2017 (10.1093/bioinformatics/btw763) |
| Tomato GO annotation on ITAG4.0 | PLAZA 5.0: Van Bel et al. 2022 (10.1093/nar/gkab1024); ITAG4.0/SL4.0: Hosmani et al. 2019 preprint (10.1101/767764) |
| Plant Reactome projection status | Plant Reactome "Computational Inference" and Release Summary V24; Gupta et al. 2024 (10.1093/nar/gkad1052); Naithani et al. 2020 (10.1093/nar/gkz996) |
| MapMan / Mercator4 preparation | Schwacke et al. 2019 (10.1016/j.molp.2019.01.003); Bolger, Schwacke & Usadel 2021 (10.1007/978-1-0716-1609-3_9) |
| agriGO / PAGE assessment | Tian et al. 2017 (10.1093/nar/gkx382); Kim & Volsky 2005 (10.1186/1471-2105-6-144) |
| clusterProfiler (framework, not used) | Yu et al. 2012 (10.1089/omi.2011.0118); Wu et al. 2021 (10.1016/j.xinn.2021.100141); Xu et al. 2024 (10.1038/s41596-024-01020-z) |

## 6. Annotation coverage known before testing

Ranked universe: 18,238 genes. Source: `results/functional/annotation/go_annotation_coverage.tsv`.

| Source | Ranked genes covered | % | Used? |
|---|---|---|---|
| PLAZA 5.0 primary direct GO (after GO-version mapping) | 15,115 | 82.88 | **yes** |
| — BP / MF / CC after propagation (non-root) | 10,686 / 12,202 / 10,272 | 58.59 / 66.90 / 56.32 | yes |
| — in tested sets (15 – 500): BP / MF / CC | 9,058 / 9,287 / 4,089 | 49.67 / 50.92 / 22.42 | yes |
| PLAZA 5.0 any source incl. orthology projections | 16,258 | 89.14 | no (cross-species transfer) |
| SGN `ITAG4.0_goterms.txt`, exact gene ID | 459 | 2.52 | no (21,311 of its IDs carry other versions; no stripping) |
| SGN `ITAG4.1_goterms.txt`, exact transcript ID + identical protein | 8,311 | 45.57 | no (other release; lower coverage) |
| Plant Reactome tomato, exact ID | 622 | 3.41 | no (projected; 5 testable pathways) |

## 7. Not done in this phase (by design)

- ORA on nominal-p genes, and any DEG list.
- Relaxed FDR.
- KEGG.
- Sensitivity runs with other sizes, metrics or permutation settings.
- PLAZA projected GO.
- Plant Reactome GSEA.
- agriGO.
- Mercator4 bin ORA.
- Any change to the Phase 4 DESeq2 analysis.

## 8. Outcome (appended after the run)

_To be completed after `scripts/R/06_functional_gsea.R`._
