#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6B, Tasks 4-5 (D-034): CAMERA (limma) competitive gene-set test on voom-transformed counts,
# GO BP/MF/CC, identical gene sets and universe as the Phase 6 fgsea analysis. Validation only: the primary DESeq2
# model and the Phase 6 fgsea results are not modified. All choices are fixed in docs/functional_robustness_plan.md §3.
# Usage (project root, env stv-r): Rscript scripts/R/11_functional_camera.R
# Inputs : results/quantification/salmon/<sample>/quant.sf (git-ignored, local), reference/ITAG4.0_tx2gene.tsv,
#          config/samplesheet.csv, results/functional/annotation/go_gene_sets.rds, results/functional/ranking/deseq2_wald_ranking.tsv
# Outputs: results/functional_robustness/camera/*
source("scripts/R/functional_robustness_helpers.R")
suppressPackageStartupMessages({ library(tximport); library(limma); library(edgeR) })
run_utc <- fr_check_frozen()
out <- file.path(fr_out, "camera")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
summ <- data.table(section = character(), metric = character(), value = character())
add <- function(s, m, v) summ <<- rbind(summ, data.table(section = s, metric = m, value = paste(v, collapse = ",")))

# --- input matrix: tximport lengthScaledTPM (tximport vignette, limma-voom section) ------------------------------------
cd <- fr_coldata()
files <- setNames(file.path("results/quantification/salmon", cd$sample, "quant.sf"), cd$sample)
fr_check(all(file.exists(files)), "quant.sf files missing")
qsha <- fread("results/quantification/quant_sf_sha256.tsv", colClasses = "character")
add("input", "quant_sf_sha256_reference_file", "results/quantification/quant_sf_sha256.tsv")
sha_now <- vapply(files, function(f) sub(" .*$", "", system2("sha256sum", f, stdout = TRUE)), "")
add("input", "quant_sf_sha256_now", paste(names(sha_now), substr(sha_now, 1, 12), sep = ":"))
fr_check(all(sha_now %in% unlist(qsha)), "quant.sf checksums differ from results/quantification/quant_sf_sha256.tsv")
t2g <- fread("reference/ITAG4.0_tx2gene.tsv", sep = "\t")
txi <- tximport(files, type = "salmon", tx2gene = as.data.frame(t2g), countsFromAbundance = "lengthScaledTPM")
fr_check(identical(colnames(txi$counts), fr_samples), "tximport sample order changed")
add("input", "tximport_countsFromAbundance", txi$countsFromAbundance)
add("input", "tximport_dim", dim(txi$counts))

rk <- fread("results/functional/ranking/deseq2_wald_ranking.tsv", sep = "\t", quote = "")
gs <- readRDS("results/functional/annotation/go_gene_sets.rds")
universe <- rk$gene_id
fr_check(identical(sort(universe), sort(gs$ranked_genes)), "Phase 6 universe and gene sets disagree")
fr_check(all(universe %in% rownames(txi$counts)), "ranked genes missing from the tximport matrix")
cnt <- txi$counts[universe, ]
add("input", "genes_universe", nrow(cnt))
add("input", "lengthScaledTPM_column_sums_universe", round(colSums(cnt)))
add("input", "lengthScaledTPM_column_sums_all_genes", round(colSums(txi$counts)))

# --- DGEList, TMM, voom, linear model ---------------------------------------------------------------------------------
design <- model.matrix(~ condition, data = cd)
fr_check(identical(colnames(design), c("(Intercept)", "conditionSTV_infected")), "unexpected design columns")
fr_check(identical(rownames(design), fr_samples), "design row order")
y <- DGEList(counts = cnt, samples = cd)
y <- calcNormFactors(y, method = "TMM")
add("normalization", "TMM_norm_factors", sprintf("%s=%.4f", rownames(y$samples), y$samples$norm.factors))
add("normalization", "effective_library_size", sprintf("%s=%.0f", rownames(y$samples), y$samples$lib.size * y$samples$norm.factors))
fe <- filterByExpr(y, design)
add("diagnostic", "filterByExpr_would_keep_of_universe", sum(fe))
v <- voom(y, design, save.plot = TRUE)
add("voom", "span", v$span)
add("voom", "weights_range", signif(range(v$weights), 4))
fit <- eBayes(lmFit(v, design))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
add("limma", "df_prior", signif(fit$df.prior, 5))
add("limma", "s2_prior", signif(fit$s2.prior, 5))
add("limma", "genes_adj_P_lt_0.05_descriptive_only", sum(tt$adj.P.Val < 0.05))

# equivalence diagnostic with the DESeq2 Wald statistic (plan §3): rho <= 0 -> STOP
rho <- cor(tt[rk$gene_id, "t"], rk$stat, method = "spearman")
add("equivalence", "spearman_limma_t_vs_DESeq2_wald", round(rho, 4))
add("equivalence", "sign_agreement_limma_t_vs_DESeq2_wald", round(mean(sign(tt[rk$gene_id, "t"]) == sign(rk$stat)), 4))
fr_check(rho > 0, "CAMERA model not equivalent: limma t and DESeq2 Wald statistic are not positively correlated")
fr_wt(data.table(gene_id = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t, P.Value = tt$P.Value,
                 adj.P.Val = tt$adj.P.Val), file.path(out, "voom_limma_gene_statistics.tsv"))

# --- CAMERA ------------------------------------------------------------------------------------------------------------
cat_tab <- fread("results/functional/annotation/go_gene_set_catalog.tsv", sep = "\t", quote = "")
term_name <- setNames(cat_tab$term_name, cat_tab$term)
tabs <- list()
for (k in fr_collections) {
  tested <- cat_tab[namespace == k & in_size_range_15_500 == TRUE, term]
  sets_k <- gs$sets[[k]][tested]
  fr_check(all(lengths(sets_k) >= 15 & lengths(sets_k) <= 500), "set size outside 15-500")
  idx <- ids2indices(sets_k, rownames(v))
  fr_check(identical(lengths(idx), lengths(sets_k)), "ids2indices lost genes (identifier mismatch)")
  c_na <- camera(v, idx, design, contrast = 2, inter.gene.cor = NA, allow.neg.cor = FALSE, use.ranks = FALSE,
                 trend.var = FALSE, sort = FALSE)
  c_01 <- camera(v, idx, design, contrast = 2, inter.gene.cor = 0.01, allow.neg.cor = FALSE, use.ranks = FALSE,
                 trend.var = FALSE, sort = FALSE)
  fr_check(identical(rownames(c_na), names(sets_k)) && identical(rownames(c_01), names(sets_k)), "camera row order")
  tab <- data.table(collection = k, go_id = rownames(c_na), term_name = unname(term_name[rownames(c_na)]),
                    NGenes = c_na$NGenes, Correlation = c_na$Correlation,
                    VIF = 1 + (c_na$NGenes - 1) * pmax(c_na$Correlation, 0),
                    Direction = c_na$Direction, PValue = c_na$PValue, FDR = c_na$FDR,
                    significant_fdr05 = c_na$FDR < fr_fdr,
                    Direction_cor001 = c_01$Direction, PValue_cor001 = c_01$PValue, FDR_cor001 = c_01$FDR,
                    significant_fdr05_cor001 = c_01$FDR < fr_fdr)
  setorder(tab, PValue, go_id)
  tabs[[k]] <- tab
  fr_wt(tab, file.path(out, sprintf("camera_go_%s_all.tsv", k)))
  fr_wt(tab[significant_fdr05 == TRUE], file.path(out, sprintf("camera_go_%s_significant_fdr05.tsv", k)))
  q <- function(x) paste(signif(quantile(x, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), 3), collapse = "/")
  add(paste0("camera_", k), "sets_tested", nrow(tab))
  add(paste0("camera_", k), "significant_fdr05_primary_cor_NA", sum(tab$significant_fdr05))
  add(paste0("camera_", k), "significant_up", sum(tab$significant_fdr05 & tab$Direction == "Up"))
  add(paste0("camera_", k), "significant_down", sum(tab$significant_fdr05 & tab$Direction == "Down"))
  add(paste0("camera_", k), "min_FDR", signif(min(tab$FDR), 4))
  add(paste0("camera_", k), "sets_FDR_lt_0.25", sum(tab$FDR < 0.25))
  add(paste0("camera_", k), "correlation_quantiles_all_sets", q(tab$Correlation))
  add(paste0("camera_", k), "correlation_negative_sets", sum(tab$Correlation < 0))
  add(paste0("camera_", k), "VIF_quantiles_all_sets", q(tab$VIF))
  add(paste0("camera_", k), "significant_fdr05_secondary_cor_0.01", sum(tab$significant_fdr05_cor001))
  add(paste0("camera_", k), "direction_agreement_NA_vs_0.01", round(mean(tab$Direction == tab$Direction_cor001), 4))
}
fr_wt(summ, file.path(out, "camera_summary.tsv"))
png(file.path(out, "voom_mean_variance.png"), width = 1400, height = 1100, res = 200)
plot(v$voom.xy, pch = 16, cex = 0.25, xlab = v$voom.xy$xlab, ylab = v$voom.xy$ylab, main = "voom: mean-variance trend (18,238 genes)")
lines(v$voom.line, col = "red", lwd = 2)
dev.off()
writeLines(c(paste("run_start_utc:", run_utc), capture.output(sessionInfo())), file.path(out, "sessionInfo_11.txt"))
print(summ, nrows = 200)
cat("11_functional_camera.R finished\n")
