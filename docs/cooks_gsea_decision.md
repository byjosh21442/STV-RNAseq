# Cook's-filtered genes in the GSEA ranking — decision D-023 (final)

Written before the Cook's-inclusion sensitivity run. Sources accessed 2026-10-26. The rule below was **not** chosen by
comparing the number of significant pathways under each option.

## 1. What Cook's filtering is in DESeq2

- Cook's distance "is a measure of how much a single sample is influencing the fitted coefficients for a gene"
  (DESeq2 vignette, Bioconductor 3.22). Love et al. (2014) define it as "the scaled distance that the coefficient vector
  ... would move if the sample were removed and the model refit", i.e. a case-deletion (leave-one-out) influence measure
  (Cook 1977).
- DESeq2 flags a gene when, for a sample, Cook's distance exceeds the 0.99 quantile of F(p, m − p) (p = model
  parameters, m = samples). Only samples in groups with ≥ 3 replicates are considered.
- With ≤ 6 replicates per condition the whole gene is flagged; with ≥ 7 the outlier count is replaced and the gene
  re-fitted (`minReplicatesForReplace = 7`, used unchanged in D-015).
- In this dataset: cutoff 18.0, 39 of 18,277 tested genes flagged (`padj_NA_reason = cooks_outlier`).

## 2. Is the Wald statistic retained? Why are p and padj NA?

- **Yes, `stat` is retained.** The vignette's Theory section plots `res$stat` against the maximum Cook's distance
  (`idx <- !is.na(W)`), and our 39 flagged genes all have finite `stat` (range −3.16 to 1.34;
  `results/functional/ranking/ranking_excluded_genes.tsv`).
- **p and padj are set to NA by design:** "If a row contains a sample with an extreme count outlier then the p value and
  adjusted p value will be set to NA" (vignette). Love et al. (2014): in small experiments "merely ignoring the outlier
  may even be considered data cherry-picking – and therefore, it is more prudent to exclude the whole gene from
  downstream analysis"; such genes are "removed from subsequent analysis, including P value adjustment".
- The vignette FAQ shows how to obtain unfiltered values (`results(dds, cooksCutoff = FALSE, independentFiltering =
  FALSE)`), i.e. the filter is an analysis choice, not a computational failure.

## 3. What fgsea does with them

`fgsea:::preparePathwaysAndStats` (fgsea 1.36.2, installed source) stops with "Not all stats values are finite numbers"
when any statistic is NA/Inf. fgsea therefore never decides: the analyst chooses which genes enter `stats`. Cook's-flagged
genes have finite statistics, so they enter **only if** the analyst builds the vector from `stat` alone.

## 4. Published DESeq2-Wald → GSEA workflows

| Source | What it does with genes whose p-value is NA | Evidence level |
|---|---|---|
| TE-Seq pipeline (Kelsey et al. 2025, *Mobile DNA*; public code `maxfieldk/TE-Seq`, `workflow/srna/scripts/ea.R`) | Builds the ranking as `setNames(res[[stat]], res$gene_id) %>% na.omit()` after default `results()`: genes with NA **stat** are dropped, Cook's-flagged genes (finite stat, NA p) are **retained implicitly**; no discussion of Cook's filtering | Code inspected (GitHub code search, 2026-10-26) |
| Phase 6 plant precedents (Tomaszewski et al. 2026; Li et al. 2026) | Handling of NA-p genes not stated in the parts verified for Phase 6 | Not verified; not used as evidence |
| DESeq2 vignette / Love et al. 2014 | Flagged genes are "set aside"/excluded "from downstream analysis" (statement about the DE analysis, not specifically GSEA) | Official documentation / primary paper |

No workflow was found that discusses the question explicitly.

## 5. Methodological literature on influential genes in ranked gene-set analysis

A Europe PMC search (`"Cook's distance" AND "GSEA" AND "DESeq2"`, 2026-10-26) returned 10 applied studies and no
methodological paper. No published evaluation of including vs excluding Cook's-flagged genes in preranked GSEA was
found. The evidence is therefore **ambiguous**.

## 6. Risks

| Option | Risks |
|---|---|
| Exclude the 39 genes (primary since Phase 6) | Non-random removal: genes with an influential sample are removed, which could remove part of a real signal if that sample is biologically variable (the maximum-Cook's sample of each gene is reported in the sensitivity output). Set sizes shrink slightly (39 genes = 0.21 % of the ranking). |
| Include the 39 genes | Their Wald statistic is, by definition, strongly influenced by one sample (DESeq2's own case-deletion diagnostic); placing them in the ranking lets single-sample artefacts enter leading edges — the situation Love et al. (2014) call "data cherry-picking" if ignored. Ranks of all other genes shift slightly. |

Because fgsea is rank-based and the flagged statistics lie inside the range of the ranked statistics (−3.16 to 1.34 vs
−4.41 to 2.86), the practical impact is expected to be small in either direction; this is checked, not assumed.

## 7. Recommended rule for this project (D-023 final; D-033)

1. **Primary: keep the Phase 6 ranking unchanged** — 18,238 genes with a finite Wald statistic **and** a DESeq2
   p-value (Cook's-flagged genes excluded). This follows the only official statement available (DESeq2 sets these
   genes aside; Love et al. 2014 recommend excluding the whole gene in small experiments) and keeps Phase 6 outputs
   byte-identical.
2. **Sensitivity (labelled, never primary):** fgsea with all 18,277 finite-statistic genes (Cook's-flagged genes
   included), gene sets rebuilt with the identical GO processing, same parameters and seed. Interpretation rules are fixed
   in [functional_robustness_plan.md](functional_robustness_plan.md) §2.
3. The choice is not revisited based on which option yields more significant pathways.

## References

- DESeq2 vignette, Bioconductor 3.22: <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html> (accessed 2026-10-26).
- Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology, 15*, 550. https://doi.org/10.1186/s13059-014-0550-8 (full text PMC4302049 read 2026-10-26).
- Cook, R. D. (1977). Detection of influential observation in linear regression. *Technometrics, 19*(1), 15–18. https://doi.org/10.2307/1268249
- Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., & Sergushichev, A. (2021). Fast gene set enrichment analysis [Preprint]. *bioRxiv*. https://doi.org/10.1101/060012 ; fgsea 1.36.2 source `preparePathwaysAndStats` (installed package, inspected 2026-10-26).
- Kelsey, M. M. G., Kalekar, R. L., & Sedivy, J. M. (2025). TE-Seq: A transposable element annotation and RNA-Seq pipeline. *Mobile DNA, 16*(1), 44. https://doi.org/10.1186/s13100-025-00381-w ; code <https://github.com/maxfieldk/TE-Seq> (`workflow/srna/scripts/ea.R`, inspected 2026-10-26).
