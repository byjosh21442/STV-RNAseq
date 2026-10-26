#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 4: DESeq2 differential expression, STV_infected vs STV_free (docs/decisions.md D-013..D-019).
# Usage (project root, env stv-r): Rscript scripts/R/02_deseq2.R
# Inputs : results/quantification/txi_final.rds (scripts/R/01_tximport_import.R), config/samplesheet.csv,
#          results/quantification/tximport_validation.tsv, reference/SGN/ITAG4.0/ITAG4.0_descriptions.txt
# Outputs: results/deseq2/*.tsv|txt, results/deseq2/plot_data/*.tsv, figures/deseq2/*.png,
#          results/deseq2/deseq2_objects.rds (git-ignored)
# Any STOP CONDITION error halts the script before results are written; no sample or gene is ever removed
# except by the approved pre-filter.
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})
options(warn = 1)

# --- approved parameters (D-013..D-019) ----------------------------------------------------------------------
levels_condition <- c("STV_free", "STV_infected")
coef_name        <- "condition_STV_infected_vs_STV_free"
min_count        <- 10
min_samples      <- 3
fdr              <- 0.05
min_abs_lfc      <- 1
pca_ntop         <- 500
ma_ylim          <- c(-5, 5)
heatmap_min_rows <- 2
seed             <- 20261025L   # vst() draws a random gene subset (nsub); fixed for reproducibility

txi_rds   <- "results/quantification/txi_final.rds"
txi_check <- "results/quantification/tximport_validation.tsv"
sheet_csv <- "config/samplesheet.csv"
desc_txt  <- "reference/SGN/ITAG4.0/ITAG4.0_descriptions.txt"
out_dir   <- "results/deseq2"
plot_dir  <- file.path(out_dir, "plot_data")
fig_dir   <- "figures/deseq2"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

summ <- data.frame(section = character(), metric = character(), value = character())
add <- function(section, metric, value) {
  summ[nrow(summ) + 1, ] <<- list(section, metric, paste(value, collapse = ","))
}
check <- function(ok, msg) if (!isTRUE(ok)) stop("STOP CONDITION: ", msg, call. = FALSE)
wt <- function(x, file) write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cond_cols <- c(STV_free = "#1B9E77", STV_infected = "#D95F02")

# --- inputs and metadata checks --------------------------------------------------------------------------------
txi <- readRDS(txi_rds)
check(identical(dim(txi$counts), c(33976L, 6L)),
      paste("tximport dimensions changed:", paste(dim(txi$counts), collapse = " x ")))
check(identical(txi$countsFromAbundance, "no"), "tximport countsFromAbundance is not 'no'")
val <- read.delim(txi_check, stringsAsFactors = FALSE)
expected_totals <- as.numeric(strsplit(val$value[val$check == "gene_count_totals"], ",")[[1]])
check(isTRUE(all.equal(unname(round(colSums(txi$counts))), expected_totals)),
      "tximport column totals differ from results/quantification/tximport_validation.tsv")

sheet <- read.csv(sheet_csv, stringsAsFactors = FALSE)
check(identical(colnames(txi$counts), sheet$sample), "txi sample order differs from config/samplesheet.csv")
coldata <- data.frame(row.names = sheet$sample, sample = sheet$sample,
                      condition = factor(sheet$condition, levels = levels_condition),
                      replicate = sheet$replicate, run_accession = sheet$run_accession)
check(!anyNA(coldata$condition), "condition values outside the approved levels")
check(all(table(coldata$condition) == 3), "expected 3 samples per condition")

add("input", "txi_dim", dim(txi$counts))
add("input", "txi_countsFromAbundance", txi$countsFromAbundance)
add("input", "samples", rownames(coldata))
add("input", "condition", as.character(coldata$condition))

# --- DESeqDataSet, design and factor levels ---------------------------------------------------------------------
dds_all <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
check(identical(levels(dds_all$condition), levels_condition), "condition factor levels / reference level")
add("design", "design", deparse(design(dds_all)))
add("design", "condition_levels", levels(dds_all$condition))
add("design", "reference_level", levels(dds_all$condition)[1])
add("design", "log2FC_definition", "log2(STV_infected / STV_free); positive = higher in infected")

# --- approved pre-filter (D-014) --------------------------------------------------------------------------------
cts_all <- counts(dds_all)
keep <- rowSums(cts_all >= min_count) >= min_samples
all_zero_rounded <- sum(rowSums(cts_all) == 0)
add("prefilter", "rule", sprintf("count >= %d in >= %d samples", min_count, min_samples))
add("prefilter", "genes_before", nrow(cts_all))
add("prefilter", "genes_all_zero_estimated_counts", sum(rowSums(txi$counts) == 0))
add("prefilter", "genes_all_zero_rounded_counts", all_zero_rounded)
add("prefilter", "genes_removed_by_rule", sum(!keep))
add("prefilter", "genes_removed_not_all_zero", sum(!keep) - all_zero_rounded)
add("prefilter", "genes_retained", sum(keep))
add("prefilter", "percent_genes_retained", round(100 * mean(keep), 2))
add("prefilter", "percent_counts_retained_per_sample",
    sprintf("%.3f", 100 * colSums(cts_all[keep, ]) / colSums(cts_all)))
dds <- dds_all[keep, ]

# --- model fit (D-015) ------------------------------------------------------------------------------------------
dds <- DESeq(dds, test = "Wald", fitType = "parametric", minReplicatesForReplace = 7, parallel = FALSE)
check(identical(resultsNames(dds), c("Intercept", coef_name)),
      paste("unexpected resultsNames:", paste(resultsNames(dds), collapse = ",")))
check(!is.null(normalizationFactors(dds)), "tximport normalization factors missing")
add("model", "resultsNames", resultsNames(dds))
add("model", "test", "Wald")
add("model", "normalization", "tximport avgTxLength offsets -> normalizationFactors (sizeFactors not used)")
add("model", "minReplicatesForReplace", 7)

nf <- normalizationFactors(dds)
sf_len <- estimateSizeFactorsForMatrix(counts(dds) / assays(dds)[["avgTxLength"]])
add("normalization", "size_factor_counts_over_avgTxLength", sprintf("%.4f", sf_len))
add("normalization", "normalizationFactor_column_median", sprintf("%.4f", apply(nf, 2, median)))
add("normalization", "normalizationFactor_range", sprintf("%.4f", range(nf)))

dfun <- dispersionFunction(dds)
mc <- mcols(dds)
fit_type <- attr(dfun, "fitType")
disp_coef <- attr(dfun, "coefficients")
add("dispersion", "fitType_used", fit_type)
add("dispersion", "trend_coefficients", if (is.null(disp_coef)) "NA" else
    paste(names(disp_coef), signif(disp_coef, 5), sep = "=", collapse = ";"))
add("dispersion", "dispPriorVar", signif(attr(dfun, "dispPriorVar"), 5))
add("dispersion", "varLogDispEsts", signif(attr(dfun, "varLogDispEsts"), 5))
add("dispersion", "median_gene_wise", signif(median(mc$dispGeneEst), 5))
add("dispersion", "median_fitted", signif(median(mc$dispFit), 5))
add("dispersion", "median_final", signif(median(mc$dispersion), 5))
add("dispersion", "dispersion_outliers_kept_gene_wise", sum(mc$dispOutlier, na.rm = TRUE))
add("dispersion", "beta_not_converged", sum(!mc$betaConv, na.rm = TRUE))

# --- results (D-016) --------------------------------------------------------------------------------------------
res <- results(dds, name = coef_name, alpha = fdr)
res_contrast <- results(dds, contrast = c("condition", "STV_infected", "STV_free"), alpha = fdr)
check(isTRUE(all.equal(res$log2FoldChange, res_contrast$log2FoldChange)) &&
        identical(is.na(res$padj), is.na(res_contrast$padj)),
      "results by coefficient name and by explicit contrast disagree")
lfc_description <- mcols(res_contrast)$description[colnames(res_contrast) == "log2FoldChange"]
check(grepl("condition STV_infected vs STV_free", lfc_description, fixed = TRUE),
      paste("unexpected contrast description:", lfc_description))

nc <- counts(dds, normalized = TRUE)
mean_inf <- rowMeans(nc[, dds$condition == "STV_infected"])
mean_free <- rowMeans(nc[, dds$condition == "STV_free"])
both_pos <- mean_inf > 0 & mean_free > 0
rho <- cor(res$log2FoldChange[both_pos], log2(mean_inf / mean_free)[both_pos], method = "spearman")
big <- both_pos & abs(res$log2FoldChange) >= min_abs_lfc
sign_agree <- mean(sign(res$log2FoldChange[big]) == sign(log2(mean_inf / mean_free)[big]))
check(rho > 0, sprintf("MLE log2FC anti-correlated with infected/free normalized means (rho = %.3f)", rho))
add("direction_check", "lfc_description", lfc_description)
add("direction_check", "spearman_MLE_lfc_vs_log2_group_mean_ratio", sprintf("%.4f", rho))
add("direction_check", "sign_agreement_genes_abs_MLE_lfc_ge_1", sprintf("%.4f (n = %d)", sign_agree, sum(big)))

cooks <- assays(dds)[["cooks"]]
p <- length(resultsNames(dds))
m <- ncol(dds)
cooks_cutoff <- qf(0.99, p, m - p)
check(isTRUE(all.equal(unname(apply(cooks, 1, max)), unname(mc$maxCooks))), "maxCooks recomputation mismatch")
na_reason <- ifelse(!is.na(res$padj), "",
             ifelse(is.na(res$pvalue), ifelse(res$baseMean == 0, "all_zero", "cooks_outlier"),
                    "independent_filtering"))
check(all(mc$maxCooks[na_reason == "cooks_outlier"] > cooks_cutoff), "NA p-value not explained by Cook's distance")
add("cooks", "cutoff_qf_0.99_p_m-p", sprintf("%.4f (p = %d, m = %d)", cooks_cutoff, p, m))
add("cooks", "genes_maxCooks_above_cutoff", sum(mc$maxCooks > cooks_cutoff, na.rm = TRUE))
add("cooks", "genes_pvalue_set_NA_by_cooks", sum(na_reason == "cooks_outlier"))
add("cooks", "genes_above_cutoff_per_sample",
    paste(colnames(dds), colSums(cooks > cooks_cutoff, na.rm = TRUE), sep = "=", collapse = ";"))

sig  <- !is.na(res$padj) & res$padj < fdr
strg <- sig & abs(res$log2FoldChange) >= min_abs_lfc
category <- ifelse(is.na(res$padj), "padj_NA",
            ifelse(strg, "stringent", ifelse(sig, "fdr05_only", "not_significant")))
direction <- ifelse(sig, ifelse(res$log2FoldChange > 0, "up", "down"), "")
md <- metadata(res)
add("results", "alpha", md$alpha)
add("results", "pAdjustMethod", "BH")
add("results", "independent_filtering_threshold_baseMean", sprintf("%.4f", md$filterThreshold))
add("results", "independent_filtering_threshold_quantile", names(md$filterThreshold))
add("results", "genes_in_results", nrow(res))
add("results", "genes_with_pvalue", sum(!is.na(res$pvalue)))
add("results", "genes_with_padj", sum(!is.na(res$padj)))
add("results", "raw_pvalue_lt_0.05_descriptive_only", sum(res$pvalue < 0.05, na.rm = TRUE))
add("results", "padj_lt_0.05", sum(sig))
add("results", "padj_lt_0.05_up", sum(sig & res$log2FoldChange > 0))
add("results", "padj_lt_0.05_down", sum(sig & res$log2FoldChange < 0))
add("results", "stringent_padj_lt_0.05_abs_MLE_lfc_ge_1", sum(strg))
add("results", "stringent_up", sum(strg & res$log2FoldChange > 0))
add("results", "stringent_down", sum(strg & res$log2FoldChange < 0))
add("results", "padj_NA_total", sum(is.na(res$padj)))
add("results", "padj_NA_all_zero", sum(na_reason == "all_zero"))
add("results", "padj_NA_cooks_outlier", sum(na_reason == "cooks_outlier"))
add("results", "padj_NA_independent_filtering", sum(na_reason == "independent_filtering"))
writeLines(capture.output(summary(res, alpha = fdr)), file.path(out_dir, "deseq2_results_summary_DESeq2.txt"))
wt(data.frame(theta = md$filterNumRej$theta, numRej = md$filterNumRej$numRej),
   file.path(out_dir, "deseq2_independent_filtering.tsv"))

# --- apeglm shrinkage (D-017) -----------------------------------------------------------------------------------
shr <- tryCatch(
  lfcShrink(dds, coef = coef_name, type = "apeglm", res = res, returnList = TRUE),
  error = function(e) stop("STOP CONDITION: apeglm shrinkage failed: ", conditionMessage(e), call. = FALSE))
res_shr <- shr$res
check(identical(res_shr$pvalue, res$pvalue) && identical(res_shr$padj, res$padj), "lfcShrink altered p-values")
check(identical(rownames(res_shr), rownames(res)), "lfcShrink changed gene order")
conv <- if ("conv" %in% colnames(shr$fit$diag)) shr$fit$diag[, "conv"] else rep(NA_real_, nrow(res))
check(length(conv) == nrow(res), "apeglm diagnostics do not match the number of genes")
add("apeglm", "status", "completed")
add("apeglm", "coef", coef_name)
add("apeglm", "conv_code_counts", paste(names(table(conv, useNA = "ifany")), table(conv, useNA = "ifany"),
                                        sep = "=", collapse = ";"))
add("apeglm", "genes_NA_shrunken_lfc", sum(is.na(res_shr$log2FoldChange)))
add("apeglm", "significant_genes_sign_change_vs_MLE",
    sum(sig & sign(res_shr$log2FoldChange) != sign(res$log2FoldChange), na.rm = TRUE))

# --- annotation and tables --------------------------------------------------------------------------------------
ann <- read.delim(desc_txt, header = FALSE, quote = "", col.names = c("gene_id", "description"),
                  colClasses = "character")
check(!anyDuplicated(ann$gene_id), "duplicated gene IDs in ITAG4.0_descriptions.txt")
gene_id <- rownames(res)
description <- ann$description[match(gene_id, ann$gene_id)]
add("annotation", "genes_in_results_without_description", sum(is.na(description)))

all_df <- data.frame(gene_id, description, baseMean = res$baseMean,
                     log2FoldChange_MLE = res$log2FoldChange, lfcSE_MLE = res$lfcSE, stat = res$stat,
                     pvalue = res$pvalue, padj = res$padj, maxCooks = mc$maxCooks,
                     padj_NA_reason = na_reason, de_category = category, direction = direction)
shr_df <- data.frame(gene_id, description, baseMean = res_shr$baseMean,
                     log2FoldChange_apeglm = res_shr$log2FoldChange, lfcSE_apeglm = res_shr$lfcSE,
                     pvalue = res_shr$pvalue, padj = res_shr$padj, log2FoldChange_MLE = res$log2FoldChange,
                     apeglm_conv = conv, de_category = category, direction = direction)
sig_df <- data.frame(gene_id, description, baseMean = res$baseMean,
                     log2FoldChange_MLE = res$log2FoldChange, lfcSE_MLE = res$lfcSE,
                     log2FoldChange_apeglm = res_shr$log2FoldChange, lfcSE_apeglm = res_shr$lfcSE,
                     stat = res$stat, pvalue = res$pvalue, padj = res$padj,
                     direction = direction, de_category = category)[sig, ]
by_padj <- function(d) d[order(d$padj, d$pvalue, d$gene_id, na.last = TRUE), ]
by_shrunken <- function(d) d[order(-abs(d$log2FoldChange_apeglm), d$padj, d$gene_id), ]
wt(by_padj(all_df), file.path(out_dir, "deseq2_results_all.tsv"))
wt(by_padj(shr_df), file.path(out_dir, "deseq2_results_shrunken.tsv"))
wt(by_shrunken(sig_df), file.path(out_dir, "deseq2_significant_fdr05.tsv"))
wt(by_shrunken(sig_df[sig_df$de_category == "stringent", ]), file.path(out_dir, "deseq2_stringent_degs.tsv"))

# --- sample QC (D-018) ------------------------------------------------------------------------------------------
set.seed(seed)
vsd <- vst(dds, blind = TRUE)
pca_df <- plotPCA(vsd, intgroup = "condition", ntop = pca_ntop, returnData = TRUE)
pv <- attr(pca_df, "percentVar")
rv <- MatrixGenerics::rowVars(assay(vsd))
top <- order(rv, decreasing = TRUE)[seq_len(min(pca_ntop, length(rv)))]
pc_full <- prcomp(t(assay(vsd)[top, ]))
pve <- pc_full$sdev^2 / sum(pc_full$sdev^2)
check(isTRUE(all.equal(unname(pve[1:2]), unname(pv))), "PCA variance does not match plotPCA")
add("sample_qc", "vst", sprintf("vst(blind = TRUE, nsub = 1000), seed = %d", seed))
add("sample_qc", "pca_ntop", length(top))
add("sample_qc", "pca_percent_variance_PC1_to_PC6", sprintf("%.2f", 100 * pve))
wt(data.frame(sample = pca_df$name, condition = pca_df$condition, PC1 = pca_df$PC1, PC2 = pca_df$PC2),
   file.path(plot_dir, "pca_vst_top500.tsv"))
wt(data.frame(PC = paste0("PC", seq_along(pve)), percent_variance = 100 * pve),
   file.path(plot_dir, "pca_vst_top500_variance.tsv"))

g <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 3, show.legend = FALSE, seed = seed) +
  scale_color_manual(values = cond_cols) +
  labs(x = sprintf("PC1: %.1f%% variance", 100 * pv[1]), y = sprintf("PC2: %.1f%% variance", 100 * pv[2]),
       title = "PCA, VST (blind) - top 500 most variable genes") +
  coord_fixed() + theme_bw()
ggsave(file.path(fig_dir, "pca_vst_top500.png"), g, width = 7, height = 6, dpi = 300)

sample_dists <- dist(t(assay(vsd)))
dm <- as.matrix(sample_dists)
hc <- hclust(sample_dists, method = "complete")
anno <- data.frame(condition = coldata$condition, row.names = rownames(coldata))
pheatmap(dm, clustering_distance_rows = sample_dists, clustering_distance_cols = sample_dists,
         clustering_method = "complete", color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         annotation_col = anno, annotation_row = anno, annotation_colors = list(condition = cond_cols),
         display_numbers = TRUE, number_format = "%.1f",
         main = "Euclidean distances, VST (blind), all retained genes",
         filename = file.path(fig_dir, "sample_distance_heatmap.png"), width = 8, height = 7)
wt(data.frame(sample = rownames(dm), dm, check.names = FALSE), file.path(plot_dir, "sample_distances_vst.tsv"))
png(file.path(fig_dir, "sample_dendrogram.png"), width = 7, height = 5, units = "in", res = 300)
plot(hc, main = "Hierarchical clustering of samples",
     sub = "complete linkage; Euclidean distance on VST (blind)", xlab = "")
dev.off()
wt(data.frame(step = seq_len(nrow(hc$merge)), merge1 = hc$merge[, 1], merge2 = hc$merge[, 2], height = hc$height),
   file.path(plot_dir, "sample_hclust_merges.tsv"))

k2 <- cutree(hc, k = 2)
diag(dm) <- NA
same <- outer(coldata$condition, coldata$condition, "==")
sample_qc <- data.frame(
  sample = rownames(coldata), condition = coldata$condition,
  total_estimated_counts = round(colSums(txi$counts)),
  detected_genes_count_gt_0 = colSums(cts_all > 0),
  counts_retained_after_prefilter = colSums(counts(dds)),
  size_factor_counts_over_avgTxLength = sf_len,
  normalizationFactor_median = apply(nf, 2, median),
  normalized_count_total = colSums(nc),
  PC1 = pca_df$PC1, PC2 = pca_df$PC2,
  nearest_sample = colnames(dm)[apply(dm, 1, which.min)],
  nearest_distance = apply(dm, 1, min, na.rm = TRUE),
  mean_distance_same_condition = sapply(seq_len(nrow(dm)), function(i) mean(dm[i, same[i, ]], na.rm = TRUE)),
  mean_distance_other_condition = sapply(seq_len(nrow(dm)), function(i) mean(dm[i, !same[i, ]], na.rm = TRUE)),
  hclust_k2_cluster = k2[rownames(coldata)])
wt(sample_qc, file.path(out_dir, "deseq2_sample_qc.tsv"))
add("sample_qc", "hclust_k2_clusters_match_condition",
    length(unique(paste(sample_qc$condition, sample_qc$hclust_k2_cluster))) == 2)

png(file.path(fig_dir, "library_size_normalization.png"), width = 11, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(10, 5, 3, 1))
barplot(sample_qc$total_estimated_counts / 1e6, names.arg = sample_qc$sample, las = 2,
        col = cond_cols[as.character(sample_qc$condition)],
        ylab = "Estimated counts (millions)", main = "Library size (tximport, all genes)")
boxplot(log2(nc + 1), las = 2, col = cond_cols[as.character(coldata$condition)], outline = FALSE,
        ylab = "log2(normalized count + 1)", main = "Normalized counts (retained genes)")
dev.off()
qs <- apply(log2(nc + 1), 2, quantile, probs = c(0, 0.25, 0.5, 0.75, 1))
wt(data.frame(sample = colnames(nc), t(qs), check.names = FALSE),
   file.path(plot_dir, "normalized_log2_count_quantiles.tsv"))

png(file.path(fig_dir, "dispersion_estimates.png"), width = 7, height = 6, units = "in", res = 300)
plotDispEsts(dds, main = sprintf("Dispersion estimates (fitType = %s)", fit_type))
dev.off()
wt(data.frame(gene_id, baseMean = mc$baseMean, dispGeneEst = mc$dispGeneEst, dispFit = mc$dispFit,
              dispMAP = mc$dispMAP, dispersion = mc$dispersion, dispOutlier = mc$dispOutlier),
   file.path(plot_dir, "dispersion_estimates.tsv"))

# --- result figures (D-019) -------------------------------------------------------------------------------------
png(file.path(fig_dir, "ma_plot_mle.png"), width = 7, height = 6, units = "in", res = 300)
plotMA(res, alpha = fdr, ylim = ma_ylim, main = "MA plot - MLE log2FC (STV_infected / STV_free)")
dev.off()
png(file.path(fig_dir, "ma_plot_apeglm.png"), width = 7, height = 6, units = "in", res = 300)
plotMA(res_shr, alpha = fdr, ylim = ma_ylim, main = "MA plot - apeglm log2FC (STV_infected / STV_free)")
dev.off()

cat_levels <- c("not significant", "padj < 0.05 only", "stringent: padj < 0.05 & |MLE log2FC| >= 1")
vol <- data.frame(gene_id, baseMean = res$baseMean, log2FoldChange_MLE = res$log2FoldChange,
                  log2FoldChange_apeglm = res_shr$log2FoldChange, padj = res$padj,
                  neg_log10_padj = -log10(res$padj), de_category = category)
wt(vol, file.path(plot_dir, "ma_volcano_data.tsv"))
vol_plot <- vol[!is.na(vol$padj), ]
vol_plot$category <- factor(c(not_significant = cat_levels[1], fdr05_only = cat_levels[2],
                              stringent = cat_levels[3])[vol_plot$de_category], levels = cat_levels)
add("figures", "volcano_genes_not_drawn_padj_NA", sum(is.na(vol$padj)))
g <- ggplot(vol_plot, aes(log2FoldChange_apeglm, neg_log10_padj, color = category)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = -log10(fdr), linetype = 2, color = "grey40") +
  scale_color_manual(values = setNames(c("grey70", "#3182BD", "#E6550D"), cat_levels), drop = FALSE) +
  labs(x = "apeglm log2 fold change (STV_infected / STV_free)", y = "-log10(BH-adjusted p)",
       color = NULL, title = "Volcano plot",
       subtitle = sprintf("%d genes with padj = NA not drawn", sum(is.na(vol$padj)))) +
  theme_bw() + theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(fig_dir, "volcano_apeglm.png"), g, width = 7, height = 7, dpi = 300)

stringent_ids <- gene_id[strg]
if (length(stringent_ids) >= heatmap_min_rows) {
  set.seed(seed)
  vsd_design <- vst(dds, blind = FALSE)
  mat <- assay(vsd_design)[stringent_ids, , drop = FALSE]
  mat <- mat - rowMeans(mat)
  pheatmap(mat, annotation_col = anno, annotation_colors = list(condition = cond_cols),
           show_rownames = nrow(mat) <= 50, clustering_method = "complete",
           main = sprintf("Stringent DEGs (n = %d), VST (blind = FALSE), row-centred", nrow(mat)),
           filename = file.path(fig_dir, "deg_heatmap_stringent.png"), width = 7, height = 8)
  wt(data.frame(gene_id = rownames(mat), mat, check.names = FALSE), file.path(plot_dir, "deg_heatmap_stringent.tsv"))
  add("figures", "deg_heatmap", sprintf("drawn (%d stringent genes)", nrow(mat)))
} else {
  add("figures", "deg_heatmap", sprintf("not drawn (%d stringent genes < %d)", length(stringent_ids), heatmap_min_rows))
}

# --- provenance -------------------------------------------------------------------------------------------------
for (pkg in c("DESeq2", "tximport", "apeglm", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer")) {
  add("versions", pkg, as.character(packageVersion(pkg)))
}
add("versions", "R", R.version.string)
add("versions", "salmon", val$value[val$check == "salmon_version"])
summ <- rbind(summ[summ$section == "versions", ], summ[summ$section != "versions", ])
wt(summ, file.path(out_dir, "deseq2_summary.tsv"))
saveRDS(list(dds = dds, res = res, res_apeglm = res_shr, vsd_blind = vsd),
        file.path(out_dir, "deseq2_objects.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "deseq2_sessionInfo.txt"))
print(summ, right = FALSE, row.names = FALSE)
cat("DESeq2 analysis complete ->", out_dir, "and", fig_dir, "\n")
