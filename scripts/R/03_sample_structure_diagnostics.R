#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 5 diagnostic checkpoint: sample-structure diagnostics (Tasks 3-6).
# READ-ONLY with respect to the primary analysis: reads the existing Phase 4 objects, never refits, filters or
# removes samples/genes, and writes nothing under results/deseq2/ or figures/deseq2/.
# Usage (project root, env stv-r): Rscript scripts/R/03_sample_structure_diagnostics.R
# Inputs : results/deseq2/deseq2_objects.rds (git-ignored; written by scripts/R/02_deseq2.R),
#          results/deseq2/plot_data/{pca_vst_top500,sample_distances_vst}.tsv, results/deseq2/deseq2_sample_qc.tsv,
#          upstream QC tables listed in the "upstream QC" section, config/samplesheet.csv
# Outputs: results/qc/sample_structure/*.tsv|txt, figures/sample_structure/*.png
# Method sources (docs/methods_references.md): correlation/distance/PCA - DESeq2 vignette + rnaseqGene workflow;
# median replicate correlation - Gierlinski et al. 2015; Cook's boxplot - DESeq2 vignette "Approach to count
# outliers"; robust PCA (rrcov PcaGrid/PcaHubert, default parameters, rlog input) - Chen et al. 2020.
suppressPackageStartupMessages({
  library(DESeq2)
  library(rrcov)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})
options(warn = 1)

seed                    <- 20261025L           # same seed as Phase 4 (scripts/R/02_deseq2.R)
pca_ntop_sets           <- c(100L, 500L, 1000L) # plus all retained genes; sensitivity diagnostics only
expected_genes          <- 18277L
expected_cooks_filtered <- 39L
rpca_crit               <- 0.975                # rrcov default crit.pca.distances (Chen et al. 2020: defaults)

obj_rds   <- "results/deseq2/deseq2_objects.rds"
sheet_csv <- "config/samplesheet.csv"
out_dir   <- "results/qc/sample_structure"
fig_dir   <- "figures/sample_structure"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

summ <- data.frame(section = character(), metric = character(), value = character())
add <- function(section, metric, value) {
  summ[nrow(summ) + 1, ] <<- list(section, metric, paste(value, collapse = ","))
}
check <- function(ok, msg) if (!isTRUE(ok)) stop("STOP CONDITION: ", msg, call. = FALSE)
wt <- function(x, file) write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
rd <- function(file) read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)
cond_cols <- c(STV_free = "#1B9E77", STV_infected = "#D95F02")

# --- inputs: existing Phase 4 objects (read-only) ----------------------------------------------------------------
obj <- readRDS(obj_rds)
dds <- obj$dds
res <- obj$res
vsd <- obj$vsd_blind
sheet <- read.csv(sheet_csv, stringsAsFactors = FALSE)
samples <- sheet$sample
check(identical(colnames(dds), samples) && identical(colnames(vsd), samples), "sample order differs from samplesheet")
check(nrow(dds) == expected_genes && nrow(vsd) == expected_genes, "retained gene count differs from Phase 4")
check(identical(deparse(design(dds)), "~condition"), "design is not ~condition")
cond <- setNames(as.character(dds$condition), samples)
rep_label <- setNames(sheet$replicate, samples)
v <- assay(vsd)
add("input", "objects", obj_rds)
add("input", "genes_retained", nrow(v))
add("input", "transform", "vst(blind = TRUE) stored by Phase 4 (not recomputed)")

# --- Task 3a: distances; reproduce Phase 4 ------------------------------------------------------------------------
dm <- as.matrix(dist(t(v)))
p4_dist <- rd("results/deseq2/plot_data/sample_distances_vst.tsv")
check(identical(p4_dist$sample, samples), "Phase 4 distance table sample order")
check(isTRUE(all.equal(unname(dm), unname(as.matrix(p4_dist[, samples])), tolerance = 1e-8)),
      "VST Euclidean distances do not reproduce Phase 4")
add("reproducibility", "phase4_vst_distances_reproduced", TRUE)

# --- Task 3b: Pearson / Spearman correlation on VST ----------------------------------------------------------------
cor_p <- cor(v, method = "pearson")
cor_s <- cor(v, method = "spearman")
wt(data.frame(sample = samples, cor_p, check.names = FALSE), file.path(out_dir, "correlation_pearson_vst.tsv"))
wt(data.frame(sample = samples, cor_s, check.names = FALSE), file.path(out_dir, "correlation_spearman_vst.tsv"))
anno <- data.frame(condition = factor(cond, levels = names(cond_cols)), row.names = samples)
cor_heatmap <- function(mat, label, file) {
  off <- mat[upper.tri(mat)]
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
           breaks = seq(min(off), 1, length.out = 101),
           color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
           display_numbers = TRUE, number_format = "%.4f",
           annotation_col = anno, annotation_row = anno, annotation_colors = list(condition = cond_cols),
           main = sprintf("%s correlation, VST (blind)\n%d retained genes; samplesheet order, not clustered",
                          label, nrow(v)),
           fontsize = 9, filename = file, width = 9, height = 7.5)
}
cor_heatmap(cor_p, "Pearson", file.path(fig_dir, "correlation_heatmap_pearson_vst.png"))
cor_heatmap(cor_s, "Spearman", file.path(fig_dir, "correlation_heatmap_spearman_vst.png"))

# --- Task 3c: pairwise table and per-sample medians -----------------------------------------------------------------
pr <- t(combn(samples, 2))
pw <- data.frame(sample_a = pr[, 1], sample_b = pr[, 2],
                 condition_a = cond[pr[, 1]], condition_b = cond[pr[, 2]],
                 same_condition = cond[pr[, 1]] == cond[pr[, 2]],
                 replicate_label_a = rep_label[pr[, 1]], replicate_label_b = rep_label[pr[, 2]],
                 euclidean_vst = dm[pr], pearson_vst = cor_p[pr], spearman_vst = cor_s[pr])
pw$euclidean_rank_1_closest <- rank(pw$euclidean_vst, ties.method = "min")
pw <- pw[order(pw$euclidean_vst), ]
wt(pw, file.path(out_dir, "pairwise_sample_distances_vst.tsv"))

off_diag <- function(m) { diag(m) <- NA; m }
dmo <- off_diag(dm); cpo <- off_diag(cor_p); cso <- off_diag(cor_s)
same <- outer(cond, cond, "==")
every <- matrix(TRUE, length(samples), length(samples))
row_median <- function(m, mask) sapply(seq_len(nrow(m)), function(i) median(m[i, mask[i, ]], na.rm = TRUE))
per_sample <- data.frame(
  sample = samples, condition = cond,
  median_euclidean_all = row_median(dmo, every),
  median_euclidean_same_condition = row_median(dmo, same),
  median_euclidean_other_condition = row_median(dmo, !same),
  nearest_sample = samples[apply(dmo, 1, which.min)],
  nearest_euclidean = apply(dmo, 1, min, na.rm = TRUE),
  median_pearson_all = row_median(cpo, every),
  median_pearson_same_condition = row_median(cpo, same),
  median_spearman_all = row_median(cso, every),
  median_spearman_same_condition = row_median(cso, same))
per_sample$rank_median_euclidean_all_1_most_remote <- rank(-per_sample$median_euclidean_all, ties.method = "min")
wt(per_sample, file.path(out_dir, "per_sample_distance_correlation_summary.tsv"))
add("structure", "nearest_neighbour", paste(samples, per_sample$nearest_sample, sep = "->", collapse = ";"))
add("structure", "most_remote_sample_median_euclidean_all", samples[which.max(per_sample$median_euclidean_all)])

# --- Task 3d: PCA gene-set sensitivity (diagnostic only; primary QC PCA stays top 500, D-018) -----------------------
rv <- MatrixGenerics::rowVars(v)
ord <- order(rv, decreasing = TRUE)
sets <- c(setNames(pca_ntop_sets, paste0("top", pca_ntop_sets)), all_retained = nrow(v))
pca_var <- list()
pca_coord <- list()
for (nm in names(sets)) {
  sel <- ord[seq_len(sets[[nm]])]
  pc <- prcomp(t(v[sel, ]), center = TRUE, scale. = FALSE)
  pve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  pca_var[[nm]] <- data.frame(gene_set = nm, n_genes = length(sel), PC = paste0("PC", seq_along(pve)),
                              percent_variance = pve)
  pca_coord[[nm]] <- data.frame(gene_set = nm, n_genes = length(sel), sample = samples, condition = cond,
                                PC1 = pc$x[, 1], PC2 = pc$x[, 2], PC3 = pc$x[, 3],
                                distance_to_centroid_all_PCs = sqrt(rowSums(pc$x^2)))
}
p4_pca <- rd("results/deseq2/plot_data/pca_vst_top500.tsv")
t5 <- pca_coord[["top500"]]
check(isTRUE(all.equal(t5$PC1, p4_pca$PC1[match(samples, p4_pca$sample)], tolerance = 1e-8)) &&
        isTRUE(all.equal(t5$PC2, p4_pca$PC2[match(samples, p4_pca$sample)], tolerance = 1e-8)),
      "top-500 PCA does not reproduce Phase 4 coordinates")
add("reproducibility", "phase4_pca_top500_reproduced", TRUE)
pca_var_df <- do.call(rbind, pca_var)
pca_coord_df <- do.call(rbind, pca_coord)
wt(pca_var_df, file.path(out_dir, "pca_sensitivity_variance.tsv"))
wt(pca_coord_df, file.path(out_dir, "pca_sensitivity_coordinates.tsv"))
for (nm in names(sets)) {
  x <- pca_var[[nm]]$percent_variance
  add("pca_sensitivity", paste0(nm, "_percent_variance_PC1_to_PC5"), sprintf("%.2f", x[1:5]))
  add("pca_sensitivity", paste0(nm, "_most_remote_sample_from_centroid"),
      samples[which.max(pca_coord[[nm]]$distance_to_centroid_all_PCs)])
}
panel_lab <- sapply(names(sets), function(nm) {
  x <- pca_var[[nm]]$percent_variance
  sprintf("%s (n = %d genes): PC1 %.1f%%, PC2 %.1f%%", nm, sets[[nm]], x[1], x[2])
})
pca_coord_df$panel <- factor(panel_lab[pca_coord_df$gene_set], levels = panel_lab)
g <- ggplot(pca_coord_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = sample), size = 2.6, show.legend = FALSE, seed = seed) +
  facet_wrap(~panel, scales = "free") +
  scale_color_manual(values = cond_cols) +
  labs(title = "SENSITIVITY DIAGNOSTIC - PCA by gene set, VST (blind), centred, not scaled",
       subtitle = "Not used to choose a PCA; the primary QC PCA remains top 500 (D-018). PC signs are arbitrary per panel.") +
  theme_bw()
ggsave(file.path(fig_dir, "pca_sensitivity_vst.png"), g, width = 11, height = 9, dpi = 300)

# --- Task 4: Cook's distances from the existing fit (filtering untouched) ---------------------------------------------
cooks <- assays(dds)[["cooks"]]
check(identical(colnames(cooks), samples), "Cook's matrix sample order")
p <- length(resultsNames(dds))
m <- ncol(dds)
cutoff <- qf(0.99, p, m - p)
n_zero <- colSums(cooks == 0, na.rm = TRUE)
qs <- t(apply(cooks, 2, quantile, probs = c(0, 0.25, 0.5, 0.75, 0.95, 0.99, 1), na.rm = TRUE))
max_idx <- apply(cooks, 1, function(r) if (all(is.na(r))) NA_integer_ else which.max(r))
n_max <- tabulate(max_idx[!is.na(max_idx)], nbins = m)
above <- colSums(cooks > cutoff, na.rm = TRUE)
ck <- data.frame(sample = samples, condition = cond, n_genes = nrow(cooks), n_NA = colSums(is.na(cooks)),
                 n_zero = n_zero, min = qs[, 1], q25 = qs[, 2], median = qs[, 3],
                 mean = colMeans(cooks, na.rm = TRUE), q75 = qs[, 4], q95 = qs[, 5], q99 = qs[, 6], max = qs[, 7],
                 cutoff_qf_0.99 = cutoff, n_genes_above_cutoff = above,
                 n_genes_sample_has_max_cooks = n_max,
                 pct_genes_sample_has_max_cooks = 100 * n_max / sum(n_max))
check(sum(above) > 0 && all(n_max >= 0), "Cook's summary failed")
add("cooks", "cutoff_qf_0.99_p_m-p", sprintf("%.4f (p = %d, m = %d)", cutoff, p, m))
add("cooks", "genes_above_cutoff_per_sample", paste(samples, above, sep = "=", collapse = ";"))
add("cooks", "genes_sample_has_max_cooks", paste(samples, n_max, sep = "=", collapse = ";"))
add("cooks", "zero_cooks_values_per_sample", paste(samples, n_zero, sep = "=", collapse = ";"))

lc <- log10(cooks)
lc[is.infinite(lc)] <- NA   # only relevant if a Cook's distance is exactly 0 (count recorded above)
png(file.path(fig_dir, "cooks_log10_boxplot.png"), width = 8, height = 7, units = "in", res = 300)
par(mar = c(11, 5, 3, 2))  # vignette uses mar = c(8,5,2,2): boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)
boxplot(lc, range = 0, las = 2, col = cond_cols[cond], ylab = "log10(Cook's distance)",
        main = "Cook's distance per sample (all retained genes)")
abline(h = log10(cutoff), lty = 2, col = "red")
legend("bottomleft", lty = 2, col = "red", bg = "white", cex = 0.8,
       legend = sprintf("DESeq2 cutoff qf(0.99, %d, %d) = %.2f", p, m - p, cutoff))
dev.off()

na_cooks <- is.na(res$pvalue) & res$baseMean > 0
check(sum(na_cooks) == expected_cooks_filtered, sprintf("Cook's-filtered genes = %d, expected %d",
                                                        sum(na_cooks), expected_cooks_filtered))
check(all(mcols(dds)$maxCooks[na_cooks] > cutoff), "a filtered gene is not explained by Cook's cutoff")
g39 <- rownames(dds)[na_cooks]
cm <- cooks[g39, , drop = FALSE]
infl <- samples[apply(cm, 1, which.max)]
cnt <- counts(dds)[g39, , drop = FALSE]
ncnt <- counts(dds, normalized = TRUE)[g39, , drop = FALSE]
t39 <- data.frame(gene_id = g39, baseMean = res$baseMean[na_cooks], maxCooks = apply(cm, 1, max),
                  influential_sample = infl, influential_condition = cond[infl],
                  n_samples_above_cutoff = rowSums(cm > cutoff),
                  samples_above_cutoff = apply(cm > cutoff, 1, function(z) paste(samples[z], collapse = ";")),
                  setNames(as.data.frame(signif(cm, 5)), paste0("cooks_", samples)),
                  setNames(as.data.frame(cnt), paste0("count_", samples)),
                  setNames(as.data.frame(round(ncnt, 1)), paste0("normcount_", samples)),
                  check.names = FALSE)
t39 <- t39[order(t39$influential_sample, -t39$maxCooks), ]
wt(t39, file.path(out_dir, "cooks_filtered_genes_influential_sample.tsv"))
infl_tab <- table(factor(infl, levels = samples))
ck$n_cooks_filtered_genes_influential <- as.integer(infl_tab[samples])
wt(ck, file.path(out_dir, "cooks_per_sample_summary.tsv"))
add("cooks", "filtered_genes_influential_sample", paste(samples, infl_tab[samples], sep = "=", collapse = ";"))

# --- Task 5: robust PCA (rrcov; Chen et al. 2020) - classification only, never removal --------------------------------
set.seed(seed)
rld <- rlog(dds, blind = TRUE)
check(identical(colnames(rld), samples), "rlog sample order")
rpca_inputs <- list(rlog_blind = t(assay(rld)), vst_blind = t(v))
rpca_methods <- list(PcaGrid = function(X) PcaGrid(X, crit.pca.distances = rpca_crit),
                     PcaHubert = function(X) PcaHubert(X, crit.pca.distances = rpca_crit))
rpca_rows <- list()
rpca_fits <- list()
for (inp in names(rpca_inputs)) for (meth in names(rpca_methods)) {
  key <- paste(meth, inp, sep = "_on_")
  set.seed(seed)
  fit <- tryCatch(rpca_methods[[meth]](rpca_inputs[[inp]]), error = function(e) e)
  if (inherits(fit, "error")) {
    add("robust_pca", paste0(key, "_error"), conditionMessage(fit))
    next
  }
  rpca_fits[[key]] <- fit
  sdv <- fit@sd; odv <- fit@od; csd <- fit@cutoff.sd; cod <- fit@cutoff.od
  cls <- ifelse(sdv <= csd & odv <= cod, "regular",
         ifelse(sdv > csd & odv <= cod, "good_leverage",
         ifelse(sdv <= csd & odv > cod, "orthogonal_outlier", "bad_leverage")))
  quan <- if ("quan" %in% slotNames(fit)) fit@quan else NA
  # Degeneracy checks (n = 6): k equal to the number of samples leaves no orthogonal complement (OD ~ 0), and
  # h = n means the "robust" covariance step uses every sample, i.e. no resistance to an outlier.
  degenerate <- c(if (fit@k >= nrow(rpca_inputs[[inp]])) "k_equals_n_no_orthogonal_space",
                  if (isTRUE(quan >= nrow(rpca_inputs[[inp]]))) "h_equals_n_no_robustness")
  assessment <- if (length(degenerate)) paste(degenerate, collapse = ";") else "not_degenerate"
  rpca_rows[[key]] <- data.frame(method = meth, input = inp, n_samples = nrow(rpca_inputs[[inp]]),
                                 n_genes = ncol(rpca_inputs[[inp]]), k_selected = fit@k, h_quan = quan,
                                 crit_pca_distances = rpca_crit, degeneracy = assessment,
                                 sample = samples, condition = cond,
                                 score_distance = sdv[samples], orthogonal_distance = odv[samples],
                                 cutoff_sd = csd, cutoff_od = cod,
                                 rrcov_flag = ifelse(fit@flag[samples], "regular", "outlier"),
                                 outlier_map_class = cls[samples])
  add("robust_pca", paste0(key, "_k"), fit@k)
  add("robust_pca", paste0(key, "_h_quan"), quan)
  add("robust_pca", paste0(key, "_degeneracy"), assessment)
  add("robust_pca", paste0(key, "_flagged_outliers"),
      if (any(!fit@flag)) paste(samples[!fit@flag[samples]], collapse = ";") else "none")
}
if (length(rpca_rows)) {
  rpca_df <- do.call(rbind, rpca_rows)
  wt(rpca_df, file.path(out_dir, "robust_pca_rrcov_classification.tsv"))
  png(file.path(fig_dir, "robust_pca_outlier_maps.png"), width = 10, height = 9, units = "in", res = 300)
  par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
  for (key in names(rpca_rows)) {
    d <- rpca_rows[[key]]
    plot(d$score_distance, d$orthogonal_distance, pch = 19, col = cond_cols[d$condition],
         xlim = c(0, 1.2 * max(d$score_distance, d$cutoff_sd[1])),
         ylim = c(0, 1.2 * max(d$orthogonal_distance, d$cutoff_od[1])),
         xlab = "Robust score distance", ylab = "Orthogonal distance",
         main = sprintf("%s on %s (k = %d, h = %s)\n%s", d$method[1], d$input[1], d$k_selected[1],
                        d$h_quan[1], d$degeneracy[1]), cex.main = 0.85)
    abline(v = d$cutoff_sd[1], h = d$cutoff_od[1], lty = 2, col = "grey40")
    text(d$score_distance, d$orthogonal_distance, labels = d$sample, pos = 3, cex = 0.65)
  }
  dev.off()
}

# --- Task 6: six-sample upstream QC table ------------------------------------------------------------------------------
sq <- rd("results/quantification/salmon_qc_summary.tsv")
rf <- rd("results/qc/raw_fastq_sanity.tsv")
rq <- rd("results/qc/raw_qc_summary.tsv")
lf <- rd("results/quantification/salmon_libformat_summary.tsv")
ln <- rd("results/qc/sample_structure/fastq_header_lanes.tsv")
dq <- rd("results/deseq2/deseq2_sample_qc.tsv")
st <- rd("results/qc/stv_signal/stv_signal_summary.tsv")
pick <- function(d, col) {
  x <- d[[col]][match(samples, d$sample)]
  check(length(x) == length(samples) && !anyNA(x), paste("missing QC column", col))
  setNames(x, samples)
}
lane1_frac <- sapply(samples, function(s) { z <- ln[ln$sample == s, ]; z$reads[z$lane == 1] / sum(z$reads) })
run_id <- sapply(samples, function(s) {
  z <- ln[ln$sample == s, ]
  paste(unique(paste(z$instrument, z$run_number, z$flowcell, sep = ":")), collapse = ";")
})
lanes <- sapply(samples, function(s) paste(sort(unique(ln$lane[ln$sample == s])), collapse = "+"))
qc_rows <- list(
  list("raw_reads", "SeqKit/ENA (salmon_qc_summary.tsv)", pick(sq, "raw_reads")),
  list("q30_bases_pct", "SeqKit stats -a (raw_fastq_sanity.tsv)", pick(rf, "q30_pct")),
  list("gc_pct_mean", "SeqKit (raw_fastq_sanity.tsv)", pick(rf, "gc_pct")),
  list("fastqc_per_sequence_gc_flag", "FastQC (raw_qc_summary.tsv)", pick(rq, "flag:Per sequence GC content")),
  list("max_adapter_pct_any_type", "FastQC (raw_qc_summary.tsv)", pick(rq, "max_adapter_pct")),
  list("duplicate_reads_pct", "FastQC via MultiQC, 100 - dedup_pct_remaining", 100 - pick(rq, "dedup_pct_remaining")),
  list("fastqc_overrepresented_sequences_n", "FastQC (raw_qc_summary.tsv)", pick(rq, "n_overrepresented")),
  list("instrument_run_flowcell", "FASTQ read names (fastq_header_lanes.tsv)", run_id),
  list("lanes", "FASTQ read names (fastq_header_lanes.tsv)", lanes),
  list("lane1_read_fraction", "FASTQ read names (fastq_header_lanes.tsv)", lane1_frac),
  list("salmon_percent_mapped", "Salmon meta_info (salmon_qc_summary.tsv)", pick(sq, "percent_mapped")),
  list("salmon_pct_decoy_of_processed", "Salmon meta_info (salmon_qc_summary.tsv)", pick(sq, "pct_decoy_of_processed")),
  list("salmon_SR_compatible_fragment_ratio", "Salmon lib_format_counts (salmon_libformat_summary.tsv)",
       pick(lf, "compatible_fragment_ratio")),
  list("salmon_strand_mapping_bias", "Salmon lib_format_counts (salmon_libformat_summary.tsv)",
       pick(lf, "strand_mapping_bias")),
  list("library_size_estimated_counts", "tximport (deseq2_sample_qc.tsv)", pick(dq, "total_estimated_counts")),
  list("detected_genes_count_gt_0", "tximport (deseq2_sample_qc.tsv)", pick(dq, "detected_genes_count_gt_0")),
  list("size_factor_counts_over_avgTxLength", "DESeq2 (deseq2_sample_qc.tsv)",
       pick(dq, "size_factor_counts_over_avgTxLength")),
  list("PC1_top500", "VST PCA (deseq2_sample_qc.tsv)", pick(dq, "PC1")),
  list("PC2_top500", "VST PCA (deseq2_sample_qc.tsv)", pick(dq, "PC2")),
  list("median_euclidean_all", "this script", setNames(per_sample$median_euclidean_all, samples)),
  list("median_euclidean_same_condition", "this script", setNames(per_sample$median_euclidean_same_condition, samples)),
  list("median_pearson_same_condition", "this script", setNames(per_sample$median_pearson_same_condition, samples)),
  list("cooks_genes_above_cutoff", "this script", setNames(ck$n_genes_above_cutoff, samples)),
  list("cooks_genes_sample_has_max", "this script", setNames(ck$n_genes_sample_has_max_cooks, samples)),
  list("cooks_filtered_genes_influential_of_39", "this script",
       setNames(ck$n_cooks_filtered_genes_influential, samples)),
  list("stv_reads_salmon", "STV signal check (stv_signal_summary.tsv)", pick(st, "stv_numreads_total")),
  list("stv_reads_per_million", "STV signal check (stv_signal_summary.tsv)", pick(st, "stv_reads_per_million_processed")))
fmt <- function(x) if (is.numeric(x)) formatC(x, digits = 6, format = "g") else as.character(x)
qc_tab <- do.call(rbind, lapply(qc_rows, function(r) {
  x <- r[[3]][samples]
  num <- is.numeric(x)
  data.frame(metric = r[[1]], source = r[[2]], t(setNames(fmt(x), samples)),
             min_sample = if (num) samples[which.min(x)] else "", max_sample = if (num) samples[which.max(x)] else "",
             check.names = FALSE)
}))
wt(qc_tab, file.path(out_dir, "upstream_qc_six_samples.tsv"))

# --- provenance -------------------------------------------------------------------------------------------------------
for (pkg in c("DESeq2", "rrcov", "robustbase", "pcaPP", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer")) {
  add("versions", pkg, as.character(packageVersion(pkg)))
}
add("versions", "R", R.version.string)
add("parameters", "seed", seed)
add("parameters", "pca_gene_sets", paste(names(sets), sets, sep = "=", collapse = ";"))
add("parameters", "rrcov", sprintf("PcaGrid(X) and PcaHubert(X), defaults, crit.pca.distances = %.3f; X = samples x genes", rpca_crit))
summ <- rbind(summ[summ$section == "versions", ], summ[summ$section != "versions", ])
wt(summ, file.path(out_dir, "sample_structure_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sample_structure_sessionInfo.txt"))
print(summ, right = FALSE, row.names = FALSE)
cat("Sample-structure diagnostics complete ->", out_dir, "and", fig_dir, "\n")
