#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6B, Task 1 (D-041): tracked reproduction of the Phase 6 leading-edge summary that was first
# computed by an untracked helper (p6_le_v2.R). Reads tracked Phase 6 outputs only; writes nothing under
# results/functional/. Descriptive only.
# Usage (project root, env stv-r): Rscript scripts/R/09_functional_leading_edge_summary.R
# Inputs : results/functional/go_gsea/{leading_edge_diagnostics_significant_sets.tsv,leading_edge_significant_sets_long.tsv,
#          gsea_go_{BP,MF,CC}_all.tsv}
# Outputs: results/functional_robustness/leading_edge/leading_edge_summary_{by_collection,overall}.tsv,
#          leading_edge_gene_recurrence_top20.tsv
source("scripts/R/functional_robustness_helpers.R")
out <- file.path(fr_out, "leading_edge")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- fread("results/functional/go_gsea/leading_edge_diagnostics_significant_sets.tsv", sep = "\t", quote = "")
d[, leading_edge_size := as.numeric(leading_edge_size)]   # integer in the file; median() of integers needs a double column
by_col <- d[, .(significant_sets = .N, leading_edge_size_median = median(leading_edge_size),
                leading_edge_fraction_median = round(median(leading_edge_fraction), 3),
                median_abs_stat_median = round(median(median_abs_stat), 3), max_abs_stat_max = round(max(max_abs_stat), 3),
                top_gene_share_max = round(max(top_gene_share_of_sum_abs_stat), 3)), by = collection]
l <- fread("results/functional/go_gsea/leading_edge_significant_sets_long.tsv", sep = "\t", quote = "")
by_col <- merge(by_col, l[, .(unique_leading_edge_genes = uniqueN(gene_id)), by = collection], by = "collection")
for (k in fr_collections) {
  a <- fr_phase6(k)
  by_col[collection == k, `:=`(sets_padj_lt_0.01 = a[padj < 0.01, .N],
                               median_size_significant_sets = median(as.numeric(a[significant_fdr05 %in% TRUE, size])))]
}
by_col <- by_col[match(fr_collections, collection)]
u <- unique(l[, .(gene_id, stat, padj_gene)])
overall <- data.table(metric = c("unique_leading_edge_genes", "positive_stat", "negative_stat", "gene_padj_lt_0.05",
                                 "median_abs_stat_unique_genes", "top_gene_share_max_all_sets"),
                      value = c(nrow(u), sum(u$stat > 0), sum(u$stat < 0), sum(u$padj_gene < 0.05, na.rm = TRUE),
                                round(median(abs(u$stat)), 3), round(max(d$top_gene_share_of_sum_abs_stat), 3)))
rec <- l[, .(n_significant_sets = .N), by = .(gene_id, description)][order(-n_significant_sets, gene_id)][1:20]
fr_wt(by_col, file.path(out, "leading_edge_summary_by_collection.tsv"))
fr_wt(overall, file.path(out, "leading_edge_summary_overall.tsv"))
fr_wt(rec, file.path(out, "leading_edge_gene_recurrence_top20.tsv"))
print(by_col); print(overall)
cat("09_functional_leading_edge_summary.R finished\n")
