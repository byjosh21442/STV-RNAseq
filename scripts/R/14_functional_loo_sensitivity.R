#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6B, Task 8 (D-037): leave-one-sample-out influence at the pathway level. LABELLED SENSITIVITY
# ANALYSIS: six reduced DESeq2 fits (five samples each), never a replacement for the primary analysis and never a reason
# to remove a sample. Criteria fixed in docs/functional_robustness_plan.md §6.
# Usage (project root, env stv-r): Rscript scripts/R/14_functional_loo_sensitivity.R
# Inputs : results/quantification/txi_final.rds (git-ignored, local), config/samplesheet.csv, results/deseq2/deseq2_results_all.tsv,
#          tracked GO mapping table + GO release (via helpers), Phase 6 fgsea tables
# Outputs: results/functional_robustness/loo/*
source("scripts/R/functional_robustness_helpers.R")
suppressPackageStartupMessages(library(DESeq2))
run_utc <- fr_check_frozen()
out <- file.path(fr_out, "loo")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cd <- fr_coldata()
txi <- readRDS("results/quantification/txi_final.rds")
fr_check(identical(colnames(txi$counts), fr_samples) && identical(txi$countsFromAbundance, "no"), "txi_final.rds changed")
tested <- fread("results/deseq2/deseq2_results_all.tsv", sep = "\t", quote = "")$gene_id
fr_check(length(tested) == 18277L, "tested gene set changed")
pairs <- fr_go_pairs()

runs <- list(); cmp <- list()
for (s in fr_samples) {
  keep <- setdiff(fr_samples, s)
  tx <- txi
  for (a in c("abundance", "counts", "length")) tx[[a]] <- tx[[a]][, keep]
  dds <- suppressMessages(DESeqDataSetFromTximport(tx, colData = cd[keep, , drop = FALSE], design = ~ condition))
  dds <- dds[tested, ]
  dds <- DESeq(dds, test = "Wald", fitType = "parametric", minReplicatesForReplace = 7, parallel = FALSE, quiet = TRUE)
  fr_check(identical(resultsNames(dds), c("Intercept", fr_coef)), "unexpected resultsNames in a reduced fit")
  r <- results(dds, name = fr_coef, alpha = fr_fdr)
  stats <- fr_ranking(rownames(r), r$stat, r$pvalue)
  sets <- fr_sets(pairs, names(stats))
  run_lab <- paste0("drop_", s)
  info <- data.table(run = run_lab, dropped_sample = s, samples_used = length(keep), genes_tested = nrow(r),
                     genes_ranked = length(stats), genes_pvalue_NA = sum(is.na(r$pvalue)),
                     genes_padj_lt_0.05 = sum(r$padj < fr_fdr, na.rm = TRUE),
                     spearman_wald_vs_primary = NA_real_)
  prim <- fread("results/functional/ranking/deseq2_wald_ranking.tsv", sep = "\t", quote = "")
  common <- intersect(prim$gene_id, names(stats))
  info$spearman_wald_vs_primary <- round(cor(prim[match(common, gene_id), stat], stats[common], method = "spearman"), 4)
  for (k in fr_collections) {
    rk <- fr_fgsea(sets[[k]], stats)
    cmp[[paste(s, k)]] <- fr_compare_primary(k, rk, run_lab)[, dropped_sample := s]
    info[, (paste0("material_warnings_", k)) := fr_material_warnings(rk)]
    info[, (paste0("sets_tested_", k)) := nrow(rk)]
  }
  runs[[s]] <- info
  cat(run_lab, "done\n")
}
runs <- rbindlist(runs)
cmp <- rbindlist(cmp)
fr_wt(runs, file.path(out, "loo_runs.tsv"))
fr_wt(cmp[significant_fdr05 %in% TRUE | new_significant, .(collection, run, dropped_sample, go_id, term_name, NES, padj,
                                                            NES_s, padj_s, status, new_significant)][order(collection, go_id, run)],
      file.path(out, "loo_primary_significant_sets_status.tsv"))
summ <- cmp[, .(primary_significant = sum(significant_fdr05 %in% TRUE), retained = sum(status == "retained", na.rm = TRUE),
                supported = sum(status == "supported", na.rm = TRUE), lost = sum(status == "lost", na.rm = TRUE),
                new_significant = sum(new_significant),
                spearman_NES_common_sets = round(cor(NES, NES_s, method = "spearman", use = "complete.obs"), 4)),
            by = .(collection, dropped_sample)]
fr_wt(summ[order(collection, match(dropped_sample, fr_samples))], file.path(out, "loo_summary.tsv"))

# sample-s dependence (plan §6): supported in both runs dropping the other two samples of the same condition, lost when s dropped
w <- dcast(cmp[significant_fdr05 %in% TRUE], collection + go_id + term_name ~ dropped_sample, value.var = "status")
dep <- list()
for (s in fr_samples) {
  grp <- if (startsWith(s, "STV_free")) fr_samples[1:3] else fr_samples[4:6]
  ctrl <- setdiff(grp, s)
  ok_ctrl <- w[[ctrl[1]]] %in% c("retained", "supported") & w[[ctrl[2]]] %in% c("retained", "supported")
  w[, (paste0("dependent_on_", s)) := ok_ctrl & w[[s]] == "lost"]
  dep[[s]] <- w[, .(sample = s, primary_significant = .N, dependent_sets = sum(get(paste0("dependent_on_", s)))), by = collection]
}
dep <- rbindlist(dep)
dep[, fraction_dependent := round(dependent_sets / primary_significant, 4)]
fr_wt(w, file.path(out, "loo_dependence_by_set.tsv"))
fr_wt(dep[order(collection, match(sample, fr_samples))], file.path(out, "loo_dependence_summary.tsv"))
print(runs); print(summ); print(dep)
# fgsea warnings of the reduced runs are reported in loo_runs.tsv (not a STOP condition in plan §6)
fr_check(all(dep$fraction_dependent <= 0.5), "pathways predominantly driven by one sample (LOO dependence > 50%)")
writeLines(c(paste("run_start_utc:", run_utc), capture.output(sessionInfo())), file.path(out, "sessionInfo_14.txt"))
cat("14_functional_loo_sensitivity.R finished\n")
