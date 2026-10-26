# Comparison with prior analyses of GSE137303 (Phase 6B, Task 11)

**Scope.** This document compares the Phase 4 – 6B results of this project with the two published analyses of the same RNA-seq
reads:
- **Fukuhara et al. (2020)**, the study that generated GSE137303;
- **Badeli et al. (2026)**, a reanalysis available as a preprint.

It records agreement and disagreement claim by claim. It does not interpret biology, and it uses no gene list, relaxed
threshold or nominal-p over-representation analysis.

**Independence of this pipeline.**
- The project's ranking, gene sets, tests, thresholds, tiers and theme rules were defined in Phase 6 / 6B and frozen before
  CAMERA, the sensitivity analyses and this comparison:
  - pre-specification in `results/functional/prespecification/`;
  - Phase 6B pre-specification in `results/functional_robustness/prespecification/`.
- None of these rules was derived from the gene or pathway results of either study. In Phase 5, Badeli et al. 2026 was read
  for methods and sample structure only ([external_reanalysis_comparison.md](external_reanalysis_comparison.md)).
- Disclosure: the Phase 6B theme regex "Hormone signalling" contains the word *ethylene*. The term came from the Phase 6
  results, and it also appears in the Fukuhara et al. (2020) abstract, which the project has known since Phase 1.

**Not independent data.** All three analyses use the same six libraries. Agreement between them is agreement between
**analysis pipelines**. It is not independent biological replication.

## 1. The three analyses

| | Fukuhara et al. 2020 (original) | Badeli et al. 2026 (preprint, not peer reviewed) | This project |
|---|---|---|---|
| Access | Abstract (PubMed 31620899); full text subscription-only, **methods not read** | HTML full text, read 2026-10-25 and 2026-10-26 | — |
| Reference | SL2.50 (GEO record) | Ensembl SL3.0 (ITAG3.0 models) + STV KX949574 | SGN SL4.0 + ITAG4.0 (D-001) |
| Quantification | TopHat/Bowtie2 → Cuffnorm FPKM (GEO record) | Trimmomatic → Rsubread align → featureCounts | Salmon 2.7.0 (`-l SR`, `--gcBias`, decoys) → tximport |
| Gene-level test | Not stated in the abstract. According to Badeli et al., FPKM differences without genome-wide FDR (**secondary, unverified**) | DESeq2 1.42.0, BH FDR < 0.05 and \|log2FC\| ≥ 1; exploratory FDR < 0.10 / 0.20 | DESeq2 1.50.2, padj < 0.05 (D-013 – D-016); no relaxed list |
| Functional analysis | Not stated in the abstract | clusterProfiler hypergeometric GO-BP ORA on relaxed lists (top 100 genes or nominal p < 0.05) | Pre-specified preranked GSEA (fgsea) on the full Wald ranking; CAMERA validation; no ORA |

## 2. Claim-by-claim comparison

Gene values are from `results/functional_robustness/prior_studies/prior_study_locus_lookup.tsv`
(`scripts/R/17_prior_study_locus_lookup.R`). Loci are matched **by Solyc locus prefix only**: the source uses ITAG3.0 versions,
this project ITAG4.0. Set values are from `results/functional_robustness/comparison/method_comparison_BP.tsv`.

| # | Claim in the prior analysis | Source | This project (primary model; Phase 6B status) | Agreement |
|---|---|---|---|---|
| 1 | STV infection affects host gene expression | Fukuhara 2020 abstract | Gene level: **0 genes at padj < 0.05** (18,277 tested). Pathway level: 215 GO sets at fgsea padj < 0.05, but **0 at CAMERA FDR < 0.05**; every theme EXPLORATORY | Not supported at gene level; pathway-level shifts are method-dependent |
| 2 | Genes of ethylene biosynthesis and signalling are down-regulated in STV-infected plants | Fukuhara 2020 abstract | ACO1 locus (Solyc07g049530.3): MLE log2FC −1.81, apeglm −1.56, Wald −4.06 (rank 18,236 of 18,238), padj 0.297. GO sets, all lower ranked in STV-infected leaves: "ethylene-activated signaling pathway" (NES −2.01, padj 0.004), "cellular response to ethylene stimulus" (−2.02, 0.004), "response to ethylene" (−1.69, 0.041); all **Tier 3** (CAMERA FDR 0.999). "ethylene biosynthetic process" (GO:0009693) has 4 genes in the universe and was **not tested** (< 15). Theme "Hormone signalling": EXPLORATORY | Same direction; not significant at gene level; set-level signal not robust across methods |
| 3 | No gene meets FDR < 0.05 and \|log2FC\| ≥ 1 | Badeli 2026 | 0 genes at padj < 0.05 (with or without the fold-change criterion) | Agree (same data; different pipelines) |
| 4 | ACO1 log2FC −1.47, padj 0.225 | Badeli 2026 Table 2 | Solyc07g049530.3: apeglm −1.56, padj 0.297; in the leading edge of 8 fgsea-significant sets (Tier 3) | Agree in direction and non-significance |
| 5 | "Pathogenesis-related (PR)-protein" Solyc12g009240, log2FC −1.12, FDR 0.083 | Badeli 2026 Table 2 | Solyc12g009240.1. Its **ITAG4.0 description is "Ethylene-responsive transcription factor"**. apeglm −1.10, padj 0.193, the most negative Wald statistic of the ranking (rank 18,238); in 1 significant leading edge (Defence theme, Tier 3) | Same direction, not significant; **the gene name differs between annotations** (not resolved here) |
| 6 | NAC1 (Solyc01g095080) log2FC ≈ −0.4, FDR > 0.5; ERF017 (Solyc07g007770) +3.97, FDR 0.157 | Badeli 2026 | Both loci exist in ITAG4.0 (Solyc01g095080.3, Solyc07g007770.2) but are **not among the 18,277 tested genes**: they fail the count ≥ 10 in ≥ 3 samples prefilter (D-014) | Not comparable |
| 7 | Other Table 2 loci: DUF247 Solyc08g067160, beta-galactosidase Solyc04g008330, PGIP Solyc02g086990, glyoxalase II Solyc02g031950, GMP synthase Solyc01g109010 | Badeli 2026 Table 2 | Not significant (padj 0.83 – 1.00). Same signs. ITAG4.0 descriptions differ for all five: Solyc08g067160.3 "Alpha/beta-Hydrolases superfamily protein", Solyc04g008330.1 "Glycosyltransferase", Solyc02g086990.4 "Cyclic nucleotide-gated ion channel-like", Solyc02g031950.3 "pathogenesis-related family protein", Solyc01g109010.3 "phox domain-containing protein" (rank 1 of the ranking, padj 1.00) | Same signs; annotation names differ |
| 8 | GO-BP ORA on relaxed lists: no term at FDR < 0.05; "ethylene biosynthetic process" and "response to fungus" at p ≈ 0.1 | Badeli 2026 | Different test and no gene list. "response to fungus" (GO:0009620): NES −1.89, fgsea padj 0.013, **Tier 3** (CAMERA FDR 0.999). "ethylene biosynthetic process": not tested (4 genes) | Not directly comparable; neither analysis supports these terms under a correlation-aware or FDR-controlled test |
| 9 | No induction of SA/JA defence markers; "STV does not trigger broad defense responses" | Badeli 2026 Discussion | "Defence / biotic response" theme: 8 sets, all lower ranked in STV-infected leaves, fgsea-only (Tier 3), EXPLORATORY. "regulation of jasmonic acid mediated signaling pathway": Tier 4 (seed-unstable) | No contradiction at a robust level; neither analysis shows a robust defence signal |
| 10 | STV−_rep1 separated on PC1 | Badeli 2026 | Same sample is the most remote in all Phase 5 diagnostics (D-022). In Phase 6B: 0 sets "strongly influenced by STV_free_rep1"; leave-one-out dependence on STV_free_rep1 ≤ 37.3 % per collection | Agree |
| 11 | Plant phenotype (taller plants, more fruit, higher germination; ~86 % seed transmission) | Fukuhara 2020 abstract | Not assessable from leaf RNA-seq | — |

## 3. Results without a counterpart in the prior analyses

Neither prior analysis reports a pathway-level test on the full ranking. The Phase 6 positive enrichment of
ribosome/translation and photosynthesis/plastid GO sets has **no counterpart** in either study:
- Ribosome / translation / RNA processing: 34 sets.
- Photosynthesis / plastid / thylakoid: 25 sets.

Badeli et al. state that "none of the highly expressed genes or major metabolic enzymes were differently expressed". That is a
gene-level statement and does not test set-level shifts. In Phase 6B, both themes are **EXPLORATORY**:
- all their sets are Tier 3;
- 0 / 34 ribosome sets are "broadly condition-consistent" (32 partially consistent, 2 ambiguous);
- 15 / 25 photosynthesis sets are broadly consistent.

## 4. Wider literature (plausibility only; directness stated)

Priority order as fixed in plan §8. None of these sources tested the hypotheses examined here.

| Tier | Source | What it actually reports | Directness for this dataset |
|---|---|---|---|
| STV | Sabanadzovic et al. (2009) | Description of STV: a 3.5 kb dsRNA virus from tomato. Seed-transmitted, but not transmitted mechanically or by grafting | Virus biology only; no host transcriptome |
| STV | Fukuhara et al. (2020) | Phenotypes and an ethylene-related expression claim (abstract) | Same data (rows 1, 2, 11) |
| STV | Badeli et al. (2026) | Reanalysis; 0 DEGs; trends only | Same data (rows 3 – 10) |
| Tomato–virus transcriptomics | — | No source cited. Two targeted Europe PMC queries (2026-10-26) returned no open-access study with a tomato host, a virus-only design and a leaf transcriptome: the hits were melon/ToLCNDV and tomato with whitefly + ToCV | — |
| Persistent / endornavirus-like viruses | Roossinck (2010) | Review: persistent plant viruses (Partitiviridae, *Endornavirus*) "do not generally cause disease"; contrasts persistent, chronic and acute lifestyles | Conceptual; not about STV transcriptomes |
| Persistent viruses | Takahashi et al. (2019) | Review of latent/persistent plant viruses and possible benefits to hosts (co-authored by the GSE137303 authors) | Conceptual |
| Persistent viruses | Khankhum et al. (2016) | Data paper. Common bean leaves, 3 endornavirus-infected vs 3 endornavirus-free libraries, edgeR. 132 DEGs (84 down, 48 up). The Blast2GO distribution (not an enrichment test) is dominated by oxidation-reduction processes | Different host, virus and statistics; indirect |
| Broad plant–virus | Whitham et al. (2003) | Arabidopsis microarrays after inoculation with five acute RNA viruses: virus-general and virus-specific changes in defence/stress-associated genes | Acute viruses, other host; indirect |

## 5. Summary (technical)

- **Gene level.** Badeli et al. and this project agree that no gene reaches a genome-wide FDR threshold. The gene-level test
  of the original study is not documented in the accessible record. The ethylene-related trend (ACO1 locus) has the same
  sign in Badeli et al. and here; it is not significant.
- **Pathway level.** This project adds a threshold-free test on the full ranking. Its positive signals (ribosome/translation,
  photosynthesis/plastid) and negative signals (defence, ethylene response, ER, oxidoreductase) are fgsea-only. None
  passes CAMERA with estimated inter-gene correlation. All are classified **EXPLORATORY**
  ([functional_robustness_plan.md](functional_robustness_plan.md) §11).
- **Annotation.** Several loci carry different functional names in ITAG3.0/Ensembl (Badeli et al.) and ITAG4.0 (this project).
  Gene names from the prior analyses should not be transferred to ITAG4.0 without checking.

## References

Full entries with DOIs and access dates are in [methods_references.md](methods_references.md):
- §1 (Fukuhara et al. 2020);
- §24 (Badeli et al. 2026);
- §41 (Sabanadzovic et al. 2009; Roossinck 2010; Takahashi et al. 2019; Khankhum et al. 2016; Whitham et al. 2003).
