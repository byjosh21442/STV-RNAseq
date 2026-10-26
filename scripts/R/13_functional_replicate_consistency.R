#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6B, Task 7 (D-036): descriptive replicate-consistency diagnostic for every GO set significant in
# fgsea or CAMERA. Blind VST (D-018), all six samples, STV_free_rep1 shown explicitly, no inferential p-values.
# Rules fixed in docs/functional_robustness_plan.md §5; figure selection in §9 (F2).
# Usage (project root, env stv-r): Rscript scripts/R/13_functional_replicate_consistency.R
# Inputs : results/deseq2/deseq2_objects.rds (vsd_blind; git-ignored, local), results/functional_robustness/comparison/*,
#          results/functional/go_gsea/gsea_go_*_all.tsv (leading edges), results/functional/annotation/go_gene_sets.rds
# Outputs: results/functional_robustness/replicate_consistency/*, figures/functional_robustness/replicate_consistency_*.png
source("scripts/R/functional_robustness_helpers.R")
suppressPackageStartupMessages({ library(SummarizedExperiment); library(ggplot2) })
run_utc <- fr_check_frozen()
out <- file.path(fr_out, "replicate_consistency")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fr_fig, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS("results/deseq2/deseq2_objects.rds")
vsd <- assay(obj$vsd_blind)
fr_check(identical(colnames(vsd), fr_samples), "VST sample order/labels changed")
fr_check(identical(as.character(obj$vsd_blind$condition), rep(c("STV_free", "STV_infected"), each = 3)), "VST condition labels changed")
rm(obj); invisible(gc())
sdv <- apply(vsd, 1, sd)
z <- (vsd[sdv > 0, ] - rowMeans(vsd[sdv > 0, ])) / sdv[sdv > 0]
gs <- readRDS("results/functional/annotation/go_gene_sets.rds")
free <- 1:3; inf <- 4:6

rows <- list()
for (k in fr_collections) {
  cmp <- fread(file.path(fr_out, "comparison", sprintf("method_comparison_%s.tsv", k)), sep = "\t", quote = "")
  p6 <- fr_phase6(k)
  ev <- cmp[!is.na(tier)]
  for (i in seq_len(nrow(ev))) {
    id <- ev$go_id[i]
    if (ev$sig_f[i]) {
      genes <- strsplit(p6[go_id == id, leading_edge_genes], ";", fixed = TRUE)[[1]]; gene_basis <- "fgsea_leading_edge"
      s <- sign(ev$NES[i])
    } else {
      genes <- gs$sets[[k]][[id]]; gene_basis <- "full_set_camera_only"
      s <- ev$dir_c[i]
    }
    g <- intersect(genes, rownames(z))
    sc <- colMeans(z[g, , drop = FALSE])
    d_all <- mean(sc[inf]) - mean(sc[free])
    P <- mean(outer(sc[inf], sc[free], "-") * s > 0)
    r_drop <- vapply(seq_along(sc), function(j) {
      ii <- setdiff(inf, j); ff <- setdiff(free, j)
      (mean(sc[ii]) - mean(sc[ff])) / d_all
    }, 0)
    names(r_drop) <- paste0("r_drop_", fr_samples)
    r1 <- r_drop[["r_drop_STV_free_rep1"]]
    cls <- if (s * d_all <= 0) "ambiguous" else if (r1 < 0.25) "strongly influenced by STV_free_rep1" else
      if (P >= 8/9 && r1 >= 0.5) "broadly condition-consistent" else if (P >= 6/9) "partially consistent" else "ambiguous"
    row <- data.table(collection = k, go_id = id, term_name = ev$term_name[i], tier = ev$tier[i],
      sig_fgsea = ev$sig_f[i], sig_camera = ev$sig_c[i], NES = ev$NES[i], padj_fgsea = ev$padj[i], FDR_camera = ev$FDR[i],
      expected_sign = s, gene_basis = gene_basis, n_genes_scored = length(g))
    row[, (paste0("score_", fr_samples)) := as.list(round(unname(sc), 4))]
    row[, `:=`(delta_all = round(d_all, 4), pairs_in_expected_direction = round(P, 4))]
    row[, (names(r_drop)) := as.list(round(unname(r_drop), 4))]
    row[, classification := cls]
    rows[[length(rows) + 1]] <- row
  }
}
rc <- rbindlist(rows)
fr_wt(rc, file.path(out, "replicate_consistency_all_significant_sets.tsv"))
cls_levels <- c("broadly condition-consistent", "partially consistent", "strongly influenced by STV_free_rep1", "ambiguous")
summ <- rc[, .N, by = .(collection, tier, classification)]
summ[, classification := factor(classification, levels = cls_levels)]
fr_wt(summ[order(collection, tier, classification)], file.path(out, "replicate_consistency_summary.tsv"))
t12 <- rc[tier %in% c("Tier 1", "Tier 2"), .(tier12_sets = .N, strongly_free_rep1 = sum(classification == "strongly influenced by STV_free_rep1"),
                                           broadly = sum(classification == "broadly condition-consistent")), by = collection]
fr_wt(t12, file.path(out, "replicate_consistency_tier12_check.tsv"))
print(summ[order(collection, tier, classification)]); print(t12)

# --- F2: pre-specified selection (plan §9) --------------------------------------------------------------------------
sel <- list()
for (k in fr_collections) {
  x <- rc[collection == k]
  t1 <- head(x[tier == "Tier 1"][order(FDR_camera, padj_fgsea, go_id)], 6)
  fo <- head(x[sig_fgsea == TRUE & sig_camera == FALSE & tier %in% c("Tier 3", "Tier 4")][order(padj_fgsea, go_id)], 2)
  if (nrow(t1)) sel[[length(sel) + 1]] <- cbind(t1, selection = "Tier 1 (ordered by CAMERA FDR, fgsea padj, GO ID)")
  if (nrow(fo)) sel[[length(sel) + 1]] <- cbind(fo, selection = "diagnostic example: strongest fgsea-only set")
}
sel <- rbindlist(sel)
fr_wt(sel[, .(collection, go_id, term_name, tier, selection, classification)], file.path(out, "figure_F2_selection.tsv"))
for (k in fr_collections) {
  x <- sel[collection == k]
  if (!nrow(x)) next
  long <- melt(x, id.vars = c("go_id", "term_name", "tier", "selection", "classification"),
               measure.vars = paste0("score_", fr_samples), variable.name = "sample", value.name = "score")
  long[, sample := sub("^score_", "", sample)]
  long[, condition := factor(sub("_rep[0-9]$", "", sample), levels = c("STV_free", "STV_infected"))]
  long[, marker := fifelse(sample == "STV_free_rep1", "STV_free_rep1", "other samples")]
  long[, panel := sprintf("%s %s\n%s | %s\n%s", go_id, substr(term_name, 1, 38), tier,
                          fifelse(grepl("diagnostic", selection), "diagnostic example", "Tier 1 selection"), classification)]
  means <- long[, .(score = mean(score)), by = .(panel, condition)]
  p <- ggplot(long, aes(condition, score)) +
    geom_crossbar(data = means, aes(ymin = score, ymax = score), width = 0.5, colour = "grey40", linewidth = 0.3) +
    geom_point(aes(shape = marker, colour = condition), size = 2.4, position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    geom_text(aes(label = sub("STV_(free|infected)_", "", sample)), size = 2.2, nudge_x = 0.22, colour = "grey30") +
    scale_colour_manual(values = c(STV_free = "#1B9E77", STV_infected = "#D95F02")) +
    scale_shape_manual(values = c("STV_free_rep1" = 17, "other samples" = 16)) +
    facet_wrap(~ panel, ncol = 4) +
    labs(title = sprintf("GO %s: replicate consistency of selected sets (all six samples)", k),
         subtitle = "Sample score = mean per-gene z-score (blind VST) over the fgsea leading edge; bar = condition mean; triangle = STV_free_rep1",
         x = NULL, y = "Sample score", colour = NULL, shape = NULL) +
    theme_bw(base_size = 8) + theme(legend.position = "bottom", strip.text = element_text(size = 6.5))
  ggsave(file.path(fr_fig, sprintf("replicate_consistency_%s.png", k)), p, width = 10, height = 2.6 * ceiling(nrow(x) / 4) + 1.2, dpi = 200)
}

fr_check(all(t12$strongly_free_rep1 <= 0.5 * t12$tier12_sets),
         "pathways predominantly driven by one sample: > 50% of Tier 1+2 sets strongly influenced by STV_free_rep1")
writeLines(c(paste("run_start_utc:", run_utc), capture.output(sessionInfo())), file.path(out, "sessionInfo_13.txt"))
cat("13_functional_replicate_consistency.R finished\n")
