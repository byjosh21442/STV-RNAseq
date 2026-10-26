#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6, Tasks 5-6: primary preranked GSEA (fgsea multilevel) on GO BP, MF and CC, with diagnostics.
# All parameters were pre-specified in docs/functional_analysis_plan.md sections 1-7 (frozen copy and sha256 in
# results/functional/prespecification/) and docs/decisions.md D-023 - D-031. Nothing here is tuned on results.
# Usage (project root, env stv-r): Rscript scripts/R/06_functional_gsea.R
# Inputs : results/functional/ranking/deseq2_wald_ranking.tsv, results/functional/annotation/go_gene_sets.rds and
#          go_gene_set_catalog.tsv (scripts/R/05_functional_ranking_annotation.R)
# Outputs: results/functional/go_gsea/*, figures/functional/gsea_*.png, figures/functional/enrichment_curves/*.png
# Material fgsea warnings or NA results: all tables are written, then the script stops with an error (STOP condition).
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(BiocParallel)
  library(ggplot2)
})
options(warn = 1)

# --- pre-specified parameters (D-026, D-027, D-031) -------------------------------------------------------------
seed       <- 20261025L
seed_diag  <- 20261026L   # diagnostic only; never replaces the primary result
fdr        <- 0.05
max_curves <- 5L
fg_args <- list(minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000, sampleSize = 101,
                scoreType = "std", gseaParam = 1, BPPARAM = SerialParam())
collections <- c("BP", "MF", "CC")

rank_tsv  <- "results/functional/ranking/deseq2_wald_ranking.tsv"
sets_rds  <- "results/functional/annotation/go_gene_sets.rds"
cat_tsv   <- "results/functional/annotation/go_gene_set_catalog.tsv"
frozen    <- "results/functional/prespecification/functional_analysis_plan_frozen.md"
frozen_ck <- "results/functional/prespecification/prespecification_checksum.tsv"
out_dir   <- "results/functional/go_gsea"
fig_dir   <- "figures/functional"
curve_dir <- file.path(fig_dir, "enrichment_curves")
for (d in c(out_dir, fig_dir, curve_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

check <- function(ok, msg) if (!isTRUE(ok)) stop("STOP CONDITION: ", msg, call. = FALSE)
wt <- function(x, file) fwrite(x, file, sep = "\t", quote = FALSE, na = "NA")

# --- the plan must have been frozen before this run ---------------------------------------------------------------
check(file.exists(frozen) && file.exists(frozen_ck), "pre-specification snapshot missing")
ck <- fread(frozen_ck, colClasses = "character")
sha_now <- sub(" .*$", "", system2("sha256sum", frozen, stdout = TRUE))
check(identical(sha_now, ck[file == frozen, sha256]), "frozen plan checksum mismatch")
run_start_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
check(run_start_utc >= ck[file == frozen, frozen_utc], "GSEA run started before the pre-specification was frozen")

# --- inputs ---------------------------------------------------------------------------------------------------
rk <- fread(rank_tsv, sep = "\t", quote = "")
check(!anyDuplicated(rk$gene_id) && !is.unsorted(rev(rk$stat)), "ranking not unique/sorted")
stats <- setNames(rk$stat, rk$gene_id)
gs <- readRDS(sets_rds)
check(identical(sort(gs$ranked_genes), sort(rk$gene_id)), "gene-set universe differs from the ranking")
check(gs$min_size == fg_args$minSize && gs$max_size == fg_args$maxSize, "size range differs from 05 script")
catalog <- fread(cat_tsv, sep = "\t", quote = "")
term_name <- setNames(catalog$term_name, catalog$term)

# --- run ----------------------------------------------------------------------------------------------------------
run_fgsea <- function(pw, s) {
  warns <- character()
  set.seed(s)
  res <- withCallingHandlers(do.call(fgseaMultilevel, c(list(pathways = pw, stats = stats), fg_args)),
                             warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(res = res, warnings = warns)
}
classify <- function(w) fifelse(grepl("likely overestimated", w, fixed = TRUE), "informational_conservative_pvalues", "material")
key_cols <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size")

summ <- list(); warn_tab <- list(); all_tabs <- list(); le_long <- list(); le_diag <- list(); sel <- list()
for (k in collections) {
  pw <- gs$sets[[k]]
  r1 <- run_fgsea(pw, seed)
  r2 <- run_fgsea(pw, seed)
  rd <- run_fgsea(pw, seed_diag)
  res <- r1$res
  expected_n <- catalog[namespace == k & in_size_range_15_500 == TRUE, .N]
  identical_rerun <- isTRUE(all.equal(as.data.frame(r1$res[, ..key_cols]), as.data.frame(r2$res[, ..key_cols]), tolerance = 0)) &&
    identical(r1$res$leadingEdge, r2$res$leadingEdge) && identical(r1$warnings, r2$warnings)
  m <- match(res$pathway, rd$res$pathway)
  sig1 <- res$padj < fdr; sigd <- rd$res$padj[m] < fdr
  sd_tab <- data.table(collection = k, go_id = res$pathway, term_name = unname(term_name[res$pathway]), NES = res$NES,
                       pval_primary = res$pval, padj_primary = res$padj, pval_seed_diag = rd$res$pval[m],
                       padj_seed_diag = rd$res$padj[m], significant_primary = sig1, significant_seed_diag = sigd,
                       call_differs = xor(sig1, sigd))
  wt(sd_tab[order(-call_differs, padj_primary, go_id)], file.path(out_dir, sprintf("gsea_go_%s_seed_diagnostic.tsv", k)))
  if (length(r1$warnings)) warn_tab[[k]] <- data.table(collection = k, warning = r1$warnings, class = classify(r1$warnings))

  tab <- res[, .(collection = k, go_id = pathway, term_name = unname(term_name[pathway]), size, ES, NES, pval, padj, log2err,
                 direction = fifelse(NES > 0, "positive (toward STV_infected end)", "negative (toward STV_free end)"),
                 significant_fdr05 = padj < fdr,
                 leading_edge_size = lengths(leadingEdge),
                 leading_edge_genes = vapply(leadingEdge, paste, "", collapse = ";"))]
  setorder(tab, pval, go_id)
  all_tabs[[k]] <- tab
  wt(tab, file.path(out_dir, sprintf("gsea_go_%s_all.tsv", k)))
  sig <- tab[significant_fdr05 %in% TRUE][order(padj, -abs(NES), go_id)]
  wt(sig, file.path(out_dir, sprintf("gsea_go_%s_significant_fdr05.tsv", k)))

  # leading-edge composition (significant sets only)
  if (nrow(sig)) {
    sub_res <- res[pathway %in% sig$go_id]
    le <- data.table(go_id = rep(sub_res$pathway, lengths(sub_res$leadingEdge)),
                     gene_id = unlist(sub_res$leadingEdge, use.names = FALSE))
    le <- merge(le, rk[, .(gene_id, rank, stat, log2FoldChange_MLE, padj_gene = padj, description)], by = "gene_id", sort = FALSE)
    le <- merge(sig[, .(go_id, term_name, NES, padj)], le, by = "go_id", sort = FALSE)
    setorder(le, padj, go_id, rank)
    le_long[[k]] <- cbind(collection = k, le)
    d <- le[, .(leading_edge_size = .N, median_abs_stat = median(abs(stat)), max_abs_stat = max(abs(stat)),
                top_gene = gene_id[which.max(abs(stat))], top_gene_share_of_sum_abs_stat = max(abs(stat)) / sum(abs(stat))),
            by = go_id]
    d <- merge(sig[, .(go_id, term_name, size, NES, padj)], d, by = "go_id")
    d[, leading_edge_fraction := leading_edge_size / size]
    le_diag[[k]] <- cbind(collection = k, d)
  }
  jac_pairs <- 0L; jac_ge05 <- 0L; jac_max <- NA_real_
  if (nrow(sig) > 1) {
    les <- setNames(res[match(sig$go_id, pathway), leadingEdge], sig$go_id)
    cmb <- combn(length(les), 2)
    j <- apply(cmb, 2, function(ix) { a <- les[[ix[1]]]; b <- les[[ix[2]]]; length(intersect(a, b)) / length(union(a, b)) })
    jac_pairs <- length(j); jac_ge05 <- sum(j >= 0.5); jac_max <- max(j)
  }

  # enrichment curves: FDR-significant only, <= 5 per direction (D-031)
  for (dir in c("positive", "negative")) {
    cand <- if (dir == "positive") sig[NES > 0] else sig[NES < 0]
    cand <- head(cand[order(padj, -abs(NES), go_id)], max_curves)
    if (!nrow(cand)) {
      sel[[paste(k, dir)]] <- data.table(collection = k, direction = dir, order = NA_integer_, go_id = NA_character_,
                                         term_name = "none: no FDR-significant set in this direction", NES = NA_real_,
                                         padj = NA_real_, figure_file = NA_character_)
      next
    }
    for (i in seq_len(nrow(cand))) {
      id <- cand$go_id[i]
      f <- file.path(curve_dir, sprintf("%s_%s_%d_%s.png", k, dir, i, sub(":", "_", id)))
      p <- plotEnrichment(pw[[id]], stats, gseaParam = 1) +
        labs(title = paste(strwrap(sprintf("GO %s %s: %s", k, id, cand$term_name[i]), 75), collapse = "\n"),
             subtitle = sprintf("NES = %.2f, padj = %.3g, size = %d; DESeq2 Wald ranking, infected vs free",
                                cand$NES[i], cand$padj[i], cand$size[i]),
             x = "Rank (1 = most positive statistic)", y = "Enrichment score")
      ggsave(f, p, width = 7.5, height = 4.5, dpi = 200)
      sel[[paste(k, dir, i)]] <- data.table(collection = k, direction = dir, order = i, go_id = id, term_name = cand$term_name[i],
                                            NES = cand$NES[i], padj = cand$padj[i], figure_file = f)
    }
  }

  summ[[k]] <- data.table(
    collection = k, sets_with_ranked_genes = length(pw), sets_tested_15_500 = nrow(res), sets_expected_from_catalog = expected_n,
    significant_fdr05 = sum(sig1, na.rm = TRUE), significant_positive_NES = sum(sig1 & res$NES > 0, na.rm = TRUE),
    significant_negative_NES = sum(sig1 & res$NES < 0, na.rm = TRUE),
    all_positive_NES = sum(res$NES > 0, na.rm = TRUE), all_negative_NES = sum(res$NES < 0, na.rm = TRUE),
    min_padj = min(res$padj, na.rm = TRUE), min_pval = min(res$pval, na.rm = TRUE),
    NES_range = paste(signif(range(res$NES, na.rm = TRUE), 4), collapse = " to "),
    na_pval_or_NES = sum(is.na(res$pval) | is.na(res$NES)), log2err_NA = sum(is.na(res$log2err)),
    warnings_total = length(r1$warnings), warnings_material = sum(classify(r1$warnings) == "material"),
    rerun_same_seed_identical = identical_rerun,
    seed_diag_significance_calls_changed = sum(xor(sig1, sigd), na.rm = TRUE),
    seed_diag_max_abs_log10_pval_diff = max(abs(log10(res$pval) - log10(rd$res$pval[m])), na.rm = TRUE),
    le_jaccard_pairs_significant = jac_pairs, le_jaccard_pairs_ge_0.5 = jac_ge05, le_jaccard_max = jac_max)
}

summary_tab <- rbindlist(summ)
wt(summary_tab, file.path(out_dir, "gsea_go_summary.tsv"))
wt(if (length(warn_tab)) rbindlist(warn_tab) else data.table(collection = character(), warning = character(), class = character()),
   file.path(out_dir, "gsea_warnings.tsv"))
wt(rbindlist(all_tabs)[significant_fdr05 %in% TRUE], file.path(out_dir, "gsea_go_all_collections_significant_fdr05.tsv"))
wt(if (length(le_long)) rbindlist(le_long) else data.table(collection = character(), go_id = character(), gene_id = character()),
   file.path(out_dir, "leading_edge_significant_sets_long.tsv"))
wt(if (length(le_diag)) rbindlist(le_diag) else data.table(collection = character(), go_id = character()),
   file.path(out_dir, "leading_edge_diagnostics_significant_sets.tsv"))
wt(rbindlist(sel, fill = TRUE), file.path(out_dir, "figure_selection.tsv"))

# --- overview figures (all tested sets, no labels) -----------------------------------------------------------------
allt <- rbindlist(all_tabs)
allt[, collection := factor(collection, levels = collections)]
allt[, status := fifelse(significant_fdr05 %in% TRUE, "padj < 0.05", "not significant")]
p <- ggplot(allt, aes(NES, -log10(padj), colour = status)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = -log10(fdr), linetype = 2, colour = "grey40") +
  geom_vline(xintercept = 0, colour = "grey70") +
  facet_wrap(~ collection, nrow = 1) +
  scale_colour_manual(values = c("padj < 0.05" = "#D95F02", "not significant" = "grey55")) +
  labs(title = "Preranked GSEA (fgsea multilevel) on GO, all tested gene sets",
       subtitle = "Ranking: DESeq2 Wald statistic, STV_infected vs STV_free; NES > 0 = toward STV_infected end; dashed line: padj = 0.05",
       x = "Normalized enrichment score (NES)", y = "-log10(BH-adjusted p)", colour = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "gsea_go_nes_vs_padj.png"), p, width = 11, height = 4.5, dpi = 200)

png(file.path(fig_dir, "gsea_go_pvalue_histograms.png"), width = 2700, height = 900, res = 200)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
for (k in collections) {
  hist(all_tabs[[k]]$pval, breaks = seq(0, 1, 0.05), col = "grey80", border = "grey60",
       main = sprintf("GO %s: nominal p-values (%d sets)", k, nrow(all_tabs[[k]])), xlab = "fgsea nominal p-value")
}
dev.off()

# --- versions --------------------------------------------------------------------------------------------------
pk <- c("fgsea", "BiocParallel", "data.table", "fastmatch", "cowplot", "ggplot2", "ontologyIndex", "DESeq2")
wt(rbind(data.table(software = "R", version = paste(R.version$major, R.version$minor, sep = ".")),
         data.table(software = pk, version = vapply(pk, function(x) as.character(packageVersion(x)), ""))),
   "results/functional/software_versions_phase6.tsv")
writeLines(c(paste("run_start_utc:", run_start_utc), capture.output(sessionInfo())),
           file.path(out_dir, "sessionInfo_06_functional_gsea.txt"))

print(summary_tab)
material <- summary_tab$warnings_material > 0 | summary_tab$na_pval_or_NES > 0 |
  summary_tab$sets_tested_15_500 != summary_tab$sets_expected_from_catalog | !summary_tab$rerun_same_seed_identical
check(!any(material), paste("material GSEA problem in:", paste(summary_tab$collection[material], collapse = ",")))
cat("06_functional_gsea.R finished\n")
