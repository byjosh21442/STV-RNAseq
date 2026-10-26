#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6B, Tasks 10, 13, 14 (D-039, D-040): pre-specified robustness figures (F1, F3, F4, F5; F2 is
# drawn by script 13), theme assignment by fixed regular expressions and the frozen ROBUST / SUPPORTED WITH CAUTION /
# EXPLORATORY rule (docs/functional_robustness_plan.md §8 - §9). Literature columns of the evidence table are added in
# docs/results_by_process.md, not here.
# Usage (project root, env stv-r): Rscript scripts/R/16_functional_robustness_figures_themes.R
# Inputs : results/functional_robustness/{comparison,replicate_consistency,loo,ch00,cooks_sensitivity}/*, Phase 6 tables
# Outputs: results/functional_robustness/themes/*, figures/functional_robustness/*.png
source("scripts/R/functional_robustness_helpers.R")
suppressPackageStartupMessages(library(ggplot2))
run_utc <- fr_check_frozen()
out <- file.path(fr_out, "themes")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fr_fig, recursive = TRUE, showWarnings = FALSE)

cmp <- rbindlist(lapply(fr_collections, function(k) fread(file.path(fr_out, "comparison", sprintf("method_comparison_%s.tsv", k)),
                                                          sep = "\t", quote = "")))
cmp[, collection := factor(collection, levels = fr_collections)]
rc <- fread(file.path(fr_out, "replicate_consistency", "replicate_consistency_all_significant_sets.tsv"), sep = "\t", quote = "")
loo <- fread(file.path(fr_out, "loo", "loo_dependence_by_set.tsv"), sep = "\t", quote = "")
loo_sum <- fread(file.path(fr_out, "loo", "loo_summary.tsv"), sep = "\t", quote = "")
chs <- fread(file.path(fr_out, "ch00", "ch00_sensitivity_vs_primary.tsv"), sep = "\t", quote = "")
cks <- fread(file.path(fr_out, "cooks_sensitivity", "cooks_inclusion_vs_primary.tsv"), sep = "\t", quote = "")
le <- fread("results/functional/go_gsea/leading_edge_significant_sets_long.tsv", sep = "\t", quote = "")

# --- themes (fixed regular expressions, plan §8) ---------------------------------------------------------------------
themes <- c(
  "Ribosome / translation / RNA processing" = "ribosom|translation|rRNA|ribonucleoprotein|preribosome|nucleolus|tRNA",
  "Photosynthesis / plastid / thylakoid" = "photosynth|plastid|chloroplast|thylakoid|photosystem|chlorophyll|light.harvesting|stroma",
  "ATP synthase / proton transport" = "ATP synthase|proton-transporting|proton transmembrane|ATP biosynthetic|ATP metabolic",
  "Defence / biotic response" = "defense|immune|response to (other organism|external biotic|biotic|bacterium|fungus|virus|oomycete|insect|wounding)|chitin",
  "Hormone signalling" = "ethylene|jasmonic|salicylic|abscisic|hormone|auxin|gibberellin|cytokinin",
  "Water deprivation / abiotic stress" = "water deprivation|osmotic|salt|desiccation|drought|cold|heat|hypoxia|oxidative stress|abiotic",
  "Endoplasmic reticulum / protein folding" = "endoplasmic reticulum|Golgi|unfolded protein|protein folding",
  "Glutathione / detoxification" = "glutathione|detoxification|toxin|xenobiotic",
  "Heme / oxidoreductase" = "heme|tetrapyrrole|monooxygenase|oxidoreductase",
  "Carbohydrate binding / cell wall" = "carbohydrate binding|cell wall|polysaccharide|glucan|pectin|xyloglucan")
sig <- cmp[!is.na(tier)]
memb <- rbindlist(lapply(names(themes), function(t) sig[grepl(themes[[t]], term_name, ignore.case = TRUE), .(theme = t, collection, go_id)]))
other <- sig[!paste(collection, go_id) %in% paste(memb$collection, memb$go_id), .(theme = "Other (no theme regex matched)", collection, go_id)]
memb <- rbind(memb, other)
x <- merge(memb, sig[, .(collection, go_id, term_name, tier, sig_f, sig_c, NES, padj, Direction, FDR)], by = c("collection", "go_id"))
x <- merge(x, rc[, .(collection, go_id, classification)], by = c("collection", "go_id"), all.x = TRUE)
x <- merge(x, loo[, .(collection, go_id, free1_dependent = dependent_on_STV_free_rep1)], by = c("collection", "go_id"), all.x = TRUE)
x <- merge(x, chs[, .(collection, go_id, ch00_sensitive)], by = c("collection", "go_id"), all.x = TRUE)
x <- merge(x, cks[, .(collection, go_id, cooks_status = status)], by = c("collection", "go_id"), all.x = TRUE)
x[, direction := fifelse(sig_f, fifelse(NES > 0, "positive", "negative"), fifelse(Direction == "Up", "positive", "negative"))]
fr_wt(x[order(theme, tier, collection, go_id)], file.path(out, "theme_membership_significant_sets.tsv"))

classify <- function(d) {
  t12 <- d[tier %in% c("Tier 1", "Tier 2")]; t1 <- d[tier == "Tier 1"]; fsig <- d[sig_f == TRUE]
  crit <- c(has_tier1 = nrow(t1) > 0,
            broadly_ge_50pct_tier12 = nrow(t12) > 0 && mean(t12$classification == "broadly condition-consistent") >= 0.5,
            no_tier12_strongly_free_rep1 = !any(t12$classification == "strongly influenced by STV_free_rep1"),
            loo_free_rep1_dependent_le_25pct = nrow(fsig) == 0 || mean(fsig$free1_dependent %in% TRUE) <= 0.25,
            tier1_not_ch00_sensitive = any(!(t1$ch00_sensitive %in% TRUE)))
  cls <- if (all(crit)) "ROBUST" else if (nrow(t12) > 0) "SUPPORTED WITH CAUTION" else "EXPLORATORY"
  c(as.list(crit), theme_classification = cls)
}
ev <- x[, c(list(sets_significant_any = .N, tier1 = sum(tier == "Tier 1"), tier2 = sum(tier == "Tier 2"),
                 tier3 = sum(tier == "Tier 3"), tier3C = sum(tier == "Tier 3C"), tier4 = sum(tier == "Tier 4"),
                 positive_sets = sum(direction == "positive"), negative_sets = sum(direction == "negative"),
                 fgsea_significant = sum(sig_f), fgsea_min_padj = signif(min(padj[sig_f], Inf), 3),
                 camera_significant = sum(sig_c), camera_min_FDR = signif(min(FDR), 3),
                 tier12_broadly = sum(tier %in% c("Tier 1", "Tier 2") & classification == "broadly condition-consistent"),
                 tier12_partially = sum(tier %in% c("Tier 1", "Tier 2") & classification == "partially consistent"),
                 tier12_free_rep1 = sum(tier %in% c("Tier 1", "Tier 2") & classification == "strongly influenced by STV_free_rep1"),
                 tier12_ambiguous = sum(tier %in% c("Tier 1", "Tier 2") & classification == "ambiguous"),
                 loo_free_rep1_dependent = sum(free1_dependent %in% TRUE), ch00_sensitive = sum(ch00_sensitive %in% TRUE),
                 cooks_not_retained = sum(sig_f & !(cooks_status %in% "retained")),
                 tier12_go_ids = paste(go_id[tier %in% c("Tier 1", "Tier 2")], collapse = ";")),
            classify(.SD)), by = theme]
# leading-edge genes recurring most often within each theme's fgsea-significant sets (top 5)
lex <- merge(x[sig_f == TRUE, .(theme, collection, go_id)], le[, .(collection, go_id, gene_id, description)], by = c("collection", "go_id"))
top_le <- lex[, .N, by = .(theme, gene_id, description)][order(theme, -N, gene_id)][, head(.SD, 5), by = theme]
top_le <- top_le[, .(top_leading_edge_genes = paste(sprintf("%s (%s; %d sets)", gene_id, substr(sub(" \\(AHRD.*$", "", description), 1, 40), N),
                                                    collapse = " | ")), by = theme]
ev <- merge(ev, top_le, by = "theme", all.x = TRUE)
setorder(ev, -tier1, -tier2, theme)
fr_wt(ev, file.path(out, "theme_evidence_classification.tsv"))
# descriptive breakdown by direction (NOT pre-specified; reported for transparency only)
evd <- x[, c(list(sets = .N, tier1 = sum(tier == "Tier 1"), tier2 = sum(tier == "Tier 2")), classify(.SD)), by = .(theme, direction)]
fr_wt(evd[order(theme, direction)], file.path(out, "theme_by_direction_descriptive.tsv"))
print(ev[, .(theme, sets_significant_any, tier1, tier2, tier3, tier3C, tier4, positive_sets, negative_sets, theme_classification)])

# --- F1: method agreement --------------------------------------------------------------------------------------------
cmp[, tier_plot := fifelse(is.na(tier), "not significant in either", tier)]
cmp[, tier_plot := factor(tier_plot, levels = c("Tier 1", "Tier 2", "Tier 3", "Tier 3C", "Tier 4", "not significant in either"))]
p1 <- ggplot(cmp[order(tier_plot, decreasing = TRUE)], aes(signed_sig_fgsea, signed_sig_camera, colour = tier_plot)) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 1.1, alpha = 0.75) + facet_wrap(~ collection, nrow = 1) +
  scale_colour_manual(values = c("Tier 1" = "#1B7837", "Tier 2" = "#7FBC41", "Tier 3" = "#E08214", "Tier 3C" = "#8073AC",
                                 "Tier 4" = "#B2182B", "not significant in either" = "grey75"), drop = FALSE) +
  labs(title = "Q1 - fgsea (DESeq2 Wald ranking) vs CAMERA (voom, estimated inter-gene correlation), all tested GO sets",
       subtitle = "Signed significance = sign(direction) x -log10(nominal p); positive = toward STV_infected; tiers pre-specified (D-035)",
       x = "fgsea signed significance", y = "CAMERA signed significance", colour = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(fr_fig, "method_agreement_signed_significance.png"), p1, width = 11, height = 4.6, dpi = 200)

# --- F3: leave-one-out retention -------------------------------------------------------------------------------------
ls_long <- melt(loo_sum, id.vars = c("collection", "dropped_sample"), measure.vars = c("retained", "supported", "lost"),
                variable.name = "status", value.name = "sets")
ls_long[, dropped_sample := factor(dropped_sample, levels = fr_samples)][, collection := factor(collection, levels = fr_collections)]
p3 <- ggplot(ls_long, aes(dropped_sample, sets, fill = status)) + geom_col() + facet_wrap(~ collection, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(retained = "#1B7837", supported = "#A6DBA0", lost = "#B2182B")) +
  labs(title = "Q3 - Leave-one-sample-out sensitivity: fate of the Phase 6 FDR-significant GO sets",
       subtitle = "retained: padj < 0.05 same sign; supported: padj < 0.25 same sign; lost: otherwise. Five-sample refits; sensitivity only",
       x = "Sample left out", y = "Phase 6 significant sets", fill = NULL) +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")
ggsave(file.path(fr_fig, "loo_retention.png"), p3, width = 11, height = 4.6, dpi = 200)

# --- F4 / F5: ch00 and Cook's sensitivity ----------------------------------------------------------------------------
plast <- unique(x[theme %in% c("Photosynthesis / plastid / thylakoid", "ATP synthase / proton transport"), paste(collection, go_id)])
nes_plot <- function(d, title, subtitle, file) {
  d <- d[!is.na(NES) & !is.na(NES_s)]
  d[, group := fifelse(paste(collection, go_id) %in% plast, "photosynthesis/plastid or ATP-synthase theme (significant set)",
                       fifelse(significant_fdr05 %in% TRUE, "other Phase 6 significant set", "not significant in Phase 6"))]
  p <- ggplot(d[order(group != "not significant in Phase 6")], aes(NES, NES_s, colour = group)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) + geom_point(size = 1, alpha = 0.75) +
    facet_wrap(~ factor(collection, levels = fr_collections), nrow = 1) +
    scale_colour_manual(values = c("photosynthesis/plastid or ATP-synthase theme (significant set)" = "#B35806",
                                   "other Phase 6 significant set" = "#542788", "not significant in Phase 6" = "grey75")) +
    labs(title = title, subtitle = subtitle, x = "NES, Phase 6 primary", y = "NES, sensitivity run", colour = NULL) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file.path(fr_fig, file), p, width = 11, height = 4.6, dpi = 200)
}
nes_plot(chs, "Q4 - Ranking without near-identical organelle-type SL4.0ch00 genes vs Phase 6 primary",
         "Same fgsea call and seed; gene sets rebuilt on the reduced universe (D-038)", "ch00_sensitivity_nes.png")
nes_plot(cks, "Q1 - Cook's-flagged genes included (18,277 genes) vs Phase 6 primary (18,238 genes)",
         "Same fgsea call and seed; gene sets rebuilt on the inclusion universe (D-033)", "cooks_sensitivity_nes.png")
writeLines(c(paste("run_start_utc:", run_utc), capture.output(sessionInfo())), file.path(out, "sessionInfo_16.txt"))
cat("16_functional_robustness_figures_themes.R finished\n")
