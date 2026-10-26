#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 5 diagnostic checkpoint, Task 7.
# HISTORICAL PROCESSED-DATA STRUCTURE CHECK - NOT USED FOR DIFFERENTIAL-EXPRESSION INFERENCE.
# Describes sample-to-sample structure in the authors' GEO processed table (Cuffnorm 2.2.1 FPKM, SL2.50 + NCBI
# RefSeq annotation; docs/geo_counts_audit.md). Nothing here feeds DESeq2 or any decision on the primary analysis.
# Usage (project root, env stv-r): Rscript scripts/R/04_geo_historical_structure_check.R
# Input : data/processed/original_geo/GSE137303_count.txt.gz (tracked; MD5 sidecar)
# Output: results/qc/sample_structure/geo_historical/*.tsv, figures/sample_structure/geo_historical_*.png
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})
options(warn = 1)
label <- "Historical processed-data structure check - not used for differential-expression inference"
seed <- 20261025L
geo_file <- "data/processed/original_geo/GSE137303_count.txt.gz"
out_dir <- "results/qc/sample_structure/geo_historical"
fig_dir <- "figures/sample_structure"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
check <- function(ok, msg) if (!isTRUE(ok)) stop("STOP CONDITION: ", msg, call. = FALSE)
wt <- function(x, file) {
  writeLines(paste("#", label), file)
  suppressWarnings(write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", append = TRUE))
}
cond_cols <- c(STV_free = "#1B9E77", STV_infected = "#D95F02")

# Column label -> sample mapping VERIFIED in docs/geo_counts_audit.md (GEO Sample_description + SRA file names).
label_map <- c("12" = "STV_free_rep1", "432" = "STV_free_rep2", "433" = "STV_free_rep3",
               "22" = "STV_infected_rep1", "412" = "STV_infected_rep2", "422" = "STV_infected_rep3")
samples <- unname(label_map)
cond <- setNames(sub("_rep[0-9]+$", "", samples), samples)
md5_expected <- trimws(strsplit(readLines(paste0(geo_file, ".md5"))[1], "\\s+")[[1]][1])
check(identical(unname(tools::md5sum(geo_file)), md5_expected), "GEO file MD5 differs from sidecar")

geo <- read.delim(gzfile(geo_file), check.names = FALSE, quote = "", stringsAsFactors = FALSE)
check(nrow(geo) == 30389L && ncol(geo) == 22L, "GEO table structure differs from docs/geo_counts_audit.md")
fpkm_cols <- paste0(names(label_map), "-FPKM")
count_cols <- sapply(names(label_map), function(l) if (paste0(l, "-count") %in% names(geo)) paste0(l, "-count")
                     else paste0(l, "count"))   # GEO header has "412count" without hyphen
check(all(c(fpkm_cols, count_cols) %in% names(geo)), "expected FPKM/count columns not found")

summ <- data.frame(section = character(), metric = character(), value = character())
add <- function(section, metric, value) summ[nrow(summ) + 1, ] <<- list(section, metric, paste(value, collapse = ","))
add("label", "scope", label)
add("input", "file", geo_file)
add("input", "md5", md5_expected)

analyse <- function(mat, tag) {
  colnames(mat) <- samples
  keep <- rowSums(mat) > 0
  lx <- log2(mat[keep, , drop = FALSE] + 1)
  add(tag, "rows_total", nrow(mat))
  add(tag, "rows_nonzero_any_sample_used", sum(keep))
  add(tag, "transform", "log2(x + 1)")
  cp <- cor(lx, method = "pearson")
  cs <- cor(lx, method = "spearman")
  dm <- as.matrix(dist(t(lx)))
  wt(data.frame(sample = samples, cp, check.names = FALSE), file.path(out_dir, paste0(tag, "_pearson.tsv")))
  wt(data.frame(sample = samples, cs, check.names = FALSE), file.path(out_dir, paste0(tag, "_spearman.tsv")))
  pr <- t(combn(samples, 2))
  pw <- data.frame(sample_a = pr[, 1], sample_b = pr[, 2], same_condition = cond[pr[, 1]] == cond[pr[, 2]],
                   euclidean = dm[pr], pearson = cp[pr], spearman = cs[pr])
  pw$euclidean_rank_1_closest <- rank(pw$euclidean, ties.method = "min")
  wt(pw[order(pw$euclidean), ], file.path(out_dir, paste0(tag, "_pairwise.tsv")))
  dmo <- dm; diag(dmo) <- NA
  cso <- cs; diag(cso) <- NA
  same <- outer(cond, cond, "==")
  ps <- data.frame(sample = samples, condition = cond,
                   median_euclidean_all = apply(dmo, 1, median, na.rm = TRUE),
                   median_euclidean_same_condition = sapply(seq_along(samples), function(i)
                     median(dmo[i, same[i, ]], na.rm = TRUE)),
                   nearest_sample = samples[apply(dmo, 1, which.min)],
                   median_spearman_same_condition = sapply(seq_along(samples), function(i)
                     median(cso[i, same[i, ]], na.rm = TRUE)))
  ps$rank_median_euclidean_all_1_most_remote <- rank(-ps$median_euclidean_all, ties.method = "min")
  wt(ps, file.path(out_dir, paste0(tag, "_per_sample.tsv")))
  add(tag, "most_remote_sample_median_euclidean", samples[which.max(ps$median_euclidean_all)])
  add(tag, "nearest_neighbours", paste(samples, ps$nearest_sample, sep = "->", collapse = ";"))
  rv <- apply(lx, 1, var)
  top <- order(rv, decreasing = TRUE)[seq_len(min(500L, nrow(lx)))]
  pc <- prcomp(t(lx[top, ]), center = TRUE, scale. = FALSE)
  pve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  add(tag, "pca_top500_percent_variance_PC1_to_PC5", sprintf("%.2f", pve[1:5]))
  pcd <- data.frame(sample = samples, condition = cond, PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                    distance_to_centroid_all_PCs = sqrt(rowSums(pc$x^2)))
  wt(pcd, file.path(out_dir, paste0(tag, "_pca_top500.tsv")))
  list(pcd = pcd, pve = pve, cs = cs)
}

fp <- analyse(as.matrix(geo[, fpkm_cols]), "geo_fpkm")
ct <- analyse(as.matrix(geo[, count_cols]), "geo_cuffnorm_count")

pcd <- rbind(cbind(fp$pcd, panel = sprintf("Cuffnorm FPKM: PC1 %.1f%%, PC2 %.1f%%", fp$pve[1], fp$pve[2])),
             cbind(ct$pcd, panel = sprintf("Cuffnorm 'count': PC1 %.1f%%, PC2 %.1f%%", ct$pve[1], ct$pve[2])))
g <- ggplot(pcd, aes(PC1, PC2, color = condition)) +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = sample), size = 2.6, show.legend = FALSE, seed = seed) +
  facet_wrap(~panel, scales = "free") + scale_color_manual(values = cond_cols) +
  labs(title = "GEO processed table (GSE137303_count.txt.gz): PCA, log2(x + 1), top 500 variable genes",
       subtitle = label) +
  theme_bw()
ggsave(file.path(fig_dir, "geo_historical_pca.png"), g, width = 11, height = 5.5, dpi = 300)
anno <- data.frame(condition = factor(cond, levels = names(cond_cols)), row.names = samples)
off <- fp$cs[upper.tri(fp$cs)]
pheatmap(fp$cs, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(min(off), 1, length.out = 101),
         color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100), display_numbers = TRUE, number_format = "%.4f",
         annotation_col = anno, annotation_row = anno, annotation_colors = list(condition = cond_cols),
         main = "GEO Cuffnorm FPKM: Spearman correlation\nHistorical check - not used for DE inference",
         fontsize = 9, filename = file.path(fig_dir, "geo_historical_spearman_fpkm.png"), width = 9, height = 7.5)

wt(summ, file.path(out_dir, "geo_historical_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "geo_historical_sessionInfo.txt"))
print(summ, right = FALSE, row.names = FALSE)
