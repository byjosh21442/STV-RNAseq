#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6B, Task 6 (D-035): fgsea (Phase 6 primary) vs CAMERA (primary, inter.gene.cor = NA) per GO
# collection: overlap, direction agreement, rank correlation and pre-specified Tiers 1-4 (docs/functional_robustness_plan.md §4).
# Usage (project root, env stv-r): Rscript scripts/R/12_functional_method_comparison.R
# Inputs : results/functional/go_gsea/gsea_go_{BP,MF,CC}_{all,seed_diagnostic}.tsv; results/functional_robustness/camera/camera_go_*_all.tsv
# Outputs: results/functional_robustness/comparison/*
source("scripts/R/functional_robustness_helpers.R")
run_utc <- fr_check_frozen()
out <- file.path(fr_out, "comparison")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

tier_of <- function(sig_f, sig_c, conc, padj, FDR, unstable) {
  supp_c <- FDR < fr_support_fdr & conc;  opp_c <- FDR < fr_support_fdr & !conc
  supp_f <- padj < fr_support_fdr & conc; opp_f <- padj < fr_support_fdr & !conc
  fifelse(!sig_f & !sig_c, NA_character_,
  fifelse((sig_f & opp_c) | (sig_c & opp_f) | (sig_f & sig_c & !conc) | (sig_f & !sig_c & !supp_c & unstable), "Tier 4",
  fifelse(sig_f & sig_c & conc & !unstable, "Tier 1",
  fifelse((sig_f & sig_c & conc & unstable) | (sig_f & !sig_c & supp_c) | (sig_c & !sig_f & supp_f), "Tier 2",
  fifelse(sig_f & !sig_c, "Tier 3", "Tier 3C")))))
}

tabs <- list(); summ <- list()
for (k in fr_collections) {
  f <- fr_phase6(k)[, .(collection, go_id, term_name, size, NES, pval, padj, sig_f = significant_fdr05 %in% TRUE,
                        leading_edge_size)]
  sd <- fr_phase6_seed(k)[, .(go_id, fgsea_seed_unstable = call_differs %in% TRUE)]
  cm <- fread(file.path(fr_out, "camera", sprintf("camera_go_%s_all.tsv", k)), sep = "\t", quote = "")
  fr_check(setequal(f$go_id, cm$go_id) && nrow(f) == nrow(cm), paste("fgsea and CAMERA tested different sets in", k))
  m <- merge(f, cm[, .(go_id, NGenes, Correlation, VIF, Direction, PValue, FDR, sig_c = significant_fdr05,
                       Direction_cor001, PValue_cor001, FDR_cor001, sig_c01 = significant_fdr05_cor001)], by = "go_id")
  m <- merge(m, sd, by = "go_id")
  fr_check(all(m$size == m$NGenes), paste("set sizes differ between methods in", k))
  m[, dir_f := sign(NES)][, dir_c := fifelse(Direction == "Up", 1, -1)][, concordant := dir_f == dir_c]
  m[, tier := tier_of(sig_f, sig_c, concordant, padj, FDR, fgsea_seed_unstable)]
  m[, signed_sig_fgsea := dir_f * -log10(pval)][, signed_sig_camera := dir_c * -log10(PValue)]
  m[, signed_sig_camera_cor001 := fifelse(Direction_cor001 == "Up", 1, -1) * -log10(PValue_cor001)]
  setorder(m, tier, FDR, padj, go_id, na.last = TRUE)
  tabs[[k]] <- m
  fr_wt(m, file.path(out, sprintf("method_comparison_%s.tsv", k)))
  both <- m[sig_f & sig_c]
  s <- data.table(
    collection = k, sets_tested = nrow(m), significant_fgsea = sum(m$sig_f), significant_camera = sum(m$sig_c),
    significant_both = nrow(both), fgsea_only = sum(m$sig_f & !m$sig_c), camera_only = sum(m$sig_c & !m$sig_f),
    both_concordant = sum(both$concordant), both_discordant = sum(!both$concordant),
    overlap_jaccard = round(nrow(both) / max(1, sum(m$sig_f | m$sig_c)), 4),
    overlap_fraction_of_fgsea = round(nrow(both) / max(1, sum(m$sig_f)), 4),
    direction_agreement_all_sets = round(mean(m$concordant), 4),
    direction_agreement_fgsea_significant = round(mean(m[sig_f == TRUE, concordant]), 4),
    spearman_signed_significance = round(cor(m$signed_sig_fgsea, m$signed_sig_camera, method = "spearman"), 4),
    spearman_signed_significance_cor001 = round(cor(m$signed_sig_fgsea, m$signed_sig_camera_cor001, method = "spearman"), 4),
    camera_cor001_significant = sum(m$sig_c01), camera_cor001_significant_and_fgsea = sum(m$sig_c01 & m$sig_f),
    tier1 = sum(m$tier == "Tier 1", na.rm = TRUE), tier2 = sum(m$tier == "Tier 2", na.rm = TRUE),
    tier3 = sum(m$tier == "Tier 3", na.rm = TRUE), tier3C = sum(m$tier == "Tier 3C", na.rm = TRUE),
    tier4 = sum(m$tier == "Tier 4", na.rm = TRUE),
    median_correlation_fgsea_significant = round(median(m[sig_f == TRUE, Correlation]), 4),
    median_VIF_fgsea_significant = round(median(m[sig_f == TRUE, VIF]), 3),
    median_correlation_all = round(median(m$Correlation), 4))
  summ[[k]] <- s
}
summ <- rbindlist(summ)
fr_wt(summ, file.path(out, "method_comparison_summary.tsv"))
tier_dir <- rbindlist(tabs)[!is.na(tier), .N, by = .(collection, tier, fgsea_direction = fifelse(dir_f > 0, "positive", "negative"))]
fr_wt(tier_dir[order(collection, tier, fgsea_direction)], file.path(out, "tier_counts_by_direction.tsv"))
print(summ)

# STOP (plan §4): major contradiction
contra <- summ$significant_both > 0 & summ$both_discordant / pmax(1, summ$significant_both) >= 0.10
fr_check(!any(contra), paste("major CAMERA contradiction (>= 10% of jointly significant sets discordant) in:",
                             paste(summ$collection[contra], collapse = ",")))
fr_check(all(summ$spearman_signed_significance >= 0), "negative between-method correlation of signed significance")
writeLines(c(paste("run_start_utc:", run_utc), capture.output(sessionInfo())), file.path(out, "sessionInfo_12.txt"))
cat("12_functional_method_comparison.R finished\n")
