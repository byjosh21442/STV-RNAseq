#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6B, D-023/D-033 sensitivity: preranked GSEA with the 39 Cook's-filtered genes INCLUDED
# (all 18,277 genes with a finite Wald statistic). Labelled sensitivity; the Phase 6 primary result is not modified.
# Also verifies (a) that the shared helper reproduces the Phase 6 primary fgsea result exactly and (b) that GO sets
# rebuilt from the tracked mapping table equal the Phase 6 sets on the primary universe (STOP otherwise).
# Usage (project root, env stv-r): Rscript scripts/R/10_functional_cooks_sensitivity.R
# Inputs : results/deseq2/deseq2_results_all.tsv, results/deseq2/deseq2_objects.rds (Cook's matrix, read-only),
#          results/functional/ranking/deseq2_wald_ranking.tsv, results/functional/annotation/go_gene_sets.rds,
#          results/functional/annotation/gene_annotation_mapping_tested_genes.tsv, GO go-basic 2026-08-05
# Outputs: results/functional_robustness/cooks_sensitivity/*, results/functional_robustness/helper_validation_phase6.tsv
source("scripts/R/functional_robustness_helpers.R")
suppressPackageStartupMessages(library(DESeq2))
run_utc <- fr_check_frozen()
out <- file.path(fr_out, "cooks_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

res <- fread("results/deseq2/deseq2_results_all.tsv", sep = "\t", quote = "", na.strings = "NA")
fr_check(nrow(res) == 18277L, "DESeq2 result rows changed")
rk <- fread("results/functional/ranking/deseq2_wald_ranking.tsv", sep = "\t", quote = "")
stats_primary <- fr_ranking(res$gene_id, res$stat, res$pvalue)
fr_check(identical(names(stats_primary), rk$gene_id) && isTRUE(all.equal(unname(stats_primary), rk$stat, tolerance = 0)),
         "primary ranking rule does not reproduce the Phase 6 ranking")
stats_incl <- fr_ranking(res$gene_id, res$stat, res$pvalue, include_na_p = TRUE)
fr_check(length(stats_incl) == 18277L, "inclusion ranking is not 18,277 genes")
added <- setdiff(names(stats_incl), names(stats_primary))
fr_check(length(added) == 39L && all(res[gene_id %in% added, padj_NA_reason] == "cooks_outlier"),
         "genes added by the inclusion rule are not exactly the 39 Cook's-flagged genes")

# --- GO sets: rebuilt == Phase 6 on the primary universe (else STOP: GO mapping changed) -------------------------------
pairs <- fr_go_pairs()
gs <- readRDS("results/functional/annotation/go_gene_sets.rds")
sets_primary_rebuilt <- fr_sets(pairs, names(stats_primary))
same <- vapply(fr_collections, function(k) fr_same_sets(sets_primary_rebuilt[[k]], gs$sets[[k]]), TRUE)
fr_check(all(same), paste("rebuilt GO sets differ from Phase 6 in:", paste(fr_collections[!same], collapse = ",")))
sets_incl <- fr_sets(pairs, names(stats_incl))

# --- helper validation: Phase 6 primary reproduced exactly -------------------------------------------------------------
val <- list(); comp <- list(); allres <- list()
for (k in fr_collections) {
  r0 <- fr_fgsea(gs$sets[[k]], stats_primary)
  p6 <- fr_phase6(k)
  m <- match(p6$go_id, r0$pathway)
  le0 <- vapply(r0$leadingEdge[m], paste, "", collapse = ";")
  val[[k]] <- data.table(collection = k, sets = nrow(r0), sets_phase6 = nrow(p6),
                         max_abs_diff_pval = max(abs(r0$pval[m] - p6$pval)), max_abs_diff_padj = max(abs(r0$padj[m] - p6$padj)),
                         max_abs_diff_NES = max(abs(r0$NES[m] - p6$NES)), leading_edges_identical = identical(le0, p6$leading_edge_genes),
                         material_warnings = fr_material_warnings(r0))
  # inclusion run
  r1 <- fr_fgsea(sets_incl[[k]], stats_incl)
  allres[[k]] <- r1[, .(collection = k, go_id = pathway, size, ES, NES, pval, padj, log2err,
                        leading_edge_size = lengths(leadingEdge), leading_edge_genes = vapply(leadingEdge, paste, "", collapse = ";"))]
  allres[[k]][, n_cooks_genes_in_set := vapply(sets_incl[[k]][go_id], function(g) sum(g %in% added), 0L)]
  allres[[k]][, n_cooks_genes_in_leading_edge := vapply(strsplit(leading_edge_genes, ";", fixed = TRUE), function(g) sum(g %in% added), 0L)]
  comp[[k]] <- fr_compare_primary(k, r1, "cooks_inclusion")
  attr(comp[[k]], "warn") <- fr_material_warnings(r1)
}
val <- rbindlist(val)
fr_wt(val, file.path(fr_out, "helper_validation_phase6.tsv"))
fr_check(all(val$sets == val$sets_phase6) && all(val$max_abs_diff_pval < 1e-12) && all(val$max_abs_diff_NES < 1e-12) &&
         all(val$leading_edges_identical), "helper does not reproduce the Phase 6 primary fgsea result")

cmp <- rbindlist(comp)
fr_wt(rbindlist(allres), file.path(out, "gsea_go_cooks_inclusion_all.tsv"))
fr_wt(cmp[order(collection, padj, go_id)], file.path(out, "cooks_inclusion_vs_primary.tsv"))
summ <- cmp[, .(sets_tested_primary = sum(!is.na(padj)), sets_tested_inclusion = sum(!is.na(padj_s)),
                primary_significant = sum(significant_fdr05 %in% TRUE),
                retained = sum(status == "retained", na.rm = TRUE), supported = sum(status == "supported", na.rm = TRUE),
                lost = sum(status == "lost", na.rm = TRUE), new_significant_inclusion_only = sum(new_significant),
                spearman_NES_common_sets = round(cor(NES, NES_s, method = "spearman", use = "complete.obs"), 4)),
            by = collection]
summ[, material_warnings := vapply(collection, function(k) attr(comp[[k]], "warn"), 0L)]
fr_wt(summ, file.path(out, "cooks_inclusion_summary.tsv"))

# --- descriptive: which sample carries the maximum Cook's distance for each of the 39 genes ----------------------------
obj <- readRDS("results/deseq2/deseq2_objects.rds")
ck <- assays(obj$dds)[["cooks"]][added, , drop = FALSE]
fr_check(identical(colnames(ck), fr_samples), "Cook's matrix sample order changed")
cs <- data.table(gene_id = added, max_cooks = apply(ck, 1, max), sample_max_cooks = colnames(ck)[apply(ck, 1, which.max)])
cs <- merge(cs, res[, .(gene_id, stat, baseMean, maxCooks, description)], by = "gene_id")
cs[, rank_in_inclusion_ranking := match(gene_id, names(stats_incl))]
fr_wt(cs[order(sample_max_cooks, -max_cooks)], file.path(out, "cooks_genes_max_sample.tsv"))
fr_wt(cs[, .N, by = sample_max_cooks][order(-N)], file.path(out, "cooks_genes_max_sample_counts.tsv"))
writeLines(c(paste("run_start_utc:", run_utc), capture.output(sessionInfo())), file.path(out, "sessionInfo_10.txt"))
print(val); print(summ); print(cs[, .N, by = sample_max_cooks])
cat("10_functional_cooks_sensitivity.R finished\n")
