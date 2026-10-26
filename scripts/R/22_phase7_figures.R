#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 7 (Task 15): exploratory figures F1 - F4 from the Phase 7 tables (no new statistics).
# F1 evidence overview (themes x evidence layers); F2 fgsea vs CAMERA on identical axes; F3 theme summary; F4 ethylene case.
# Colour rules: fgsea significance is never drawn in a colour that suggests validation; CAMERA cells are grey (not significant).
# Usage (project root, env stv-r): Rscript scripts/R/22_phase7_figures.R   (after scripts 20 and 21)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(cowplot) })
options(warn = 1)
out <- "results/biological_integration"; fig <- "figures/biological_integration"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) fread(file.path(out, f), sep = "\t", quote = "")
m <- rd("theme_evidence_matrix.tsv"); s <- rd("theme_set_level_evidence.tsv")
eth <- rd("ethylene_case_study.tsv"); es <- rd("ethylene_go_sets.tsv")
themes10 <- m$theme
short <- c("Ribosome / translation / RNA processing" = "Ribosome / translation", "Photosynthesis / plastid / thylakoid" = "Photosynthesis / plastid",
           "ATP synthase / proton transport" = "ATP synthase / proton transport", "Defence / biotic response" = "Defence response",
           "Hormone signalling" = "Hormone / ethylene", "Endoplasmic reticulum / protein folding" = "Endoplasmic reticulum",
           "Heme / oxidoreductase" = "Heme / oxidoreductase", "Glutathione / detoxification" = "Glutathione",
           "Water deprivation / abiotic stress" = "Water / abiotic stress", "Carbohydrate binding / cell wall" = "Cell wall / carbohydrate")
theme_cols <- c("Ribosome / translation / RNA processing" = "#1f78b4", "Photosynthesis / plastid / thylakoid" = "#33a02c",
                "ATP synthase / proton transport" = "#b2df8a", "Defence / biotic response" = "#e31a1c", "Hormone signalling" = "#ff7f00",
                "Endoplasmic reticulum / protein folding" = "#6a3d9a", "Heme / oxidoreductase" = "#b15928",
                "Glutathione / detoxification" = "#fb9a99", "Water deprivation / abiotic stress" = "#a6cee3",
                "Carbohydrate binding / cell wall" = "#cab2d6", "Other (no theme regex matched)" = "#d9d9d9")

# ---------------------------------------------------------------- F1 integrated evidence overview
pct <- function(x) paste0(round(100 * x), " %")
f1 <- rbindlist(lapply(seq_len(nrow(m)), function(i) {
  r <- m[i]; n <- r$n_sets
  rep_share <- r$replicate_broadly / n; ch_share <- r$ch00_sensitive_prespecified / n
  data.table(theme = r$theme, layer = c("fgsea\n(L2)", "CAMERA\n(L3)", "Replicate\nconsistency", "Leave-one-\nsample-out", "Cook's\nranking", "ch00 /\norganelle", "STV\nliterature", "Phase 7\nclass"),
    status = c("fgsea-significant only",
               "not significant",
               fifelse(rep_share >= 0.75, "consistent / stable", fifelse(rep_share >= 0.5, "partial / minor", "weak / dependent")),
               fifelse(r$loo_max_single_sample_share >= 0.5, "weak / dependent", fifelse(r$loo_max_single_sample_share >= 0.25, "partial / minor", "consistent / stable")),
               fifelse(r$cooks_not_retained == 0, "consistent / stable", "partial / minor"),
               fifelse(ch_share >= 0.5, "weak / dependent", fifelse(r$ch00_sensitive_prespecified > 0 | r$sets_with_org_identical_le_gene > 0, "partial / minor", "consistent / stable")),
               fcase(grepl("^same direction", r$literature_relation_category), "same direction, not significant",
                     grepl("^no comparable", r$literature_relation_category), "no comparable claim", default = "not directly comparable"),
               r$phase7_class),
    label = c(paste0(n, " sets\n+", r$fgsea_positive, " / -", r$fgsea_negative),
              paste0("0/", n, " FDR<0.05\ndir. agree ", r$camera_direction_concordant, "/", n),
              paste0("broadly\n", r$replicate_broadly, "/", n),
              paste0("max ", pct(r$loo_max_single_sample_share), "\n", sub("STV_", "", r$loo_max_sample)),
              paste0(n - r$cooks_not_retained, "/", n, "\nretained"),
              paste0(r$ch00_sensitive_prespecified, "/", n, " sensitive\nLE ", pct(r$le_slot_share_org_identical)),
              fcase(grepl("^same direction", r$literature_relation_category), "same dir.\n(n.s.)", grepl("^no comparable", r$literature_relation_category), "none", default = "acute only"),
              r$phase7_class))
}))
f1[, theme := factor(short[theme], levels = rev(short[themes10]))]
f1[, layer := factor(layer, levels = unique(layer))]
f1_cols <- c("fgsea-significant only" = "#c7d9ec", "not significant" = "#bdbdbd", "consistent / stable" = "#a6d96a",
             "partial / minor" = "#fee08b", "weak / dependent" = "#f46d43", "same direction, not significant" = "#d9ef8b",
             "not directly comparable" = "#e0e0e0", "no comparable claim" = "#f0f0f0", "EXPLORATORY" = "#fdae61", "UNSUPPORTED" = "#969696")
p1 <- ggplot(f1, aes(layer, theme, fill = status)) + geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 2.3, lineheight = 0.85) +
  scale_fill_manual(values = f1_cols, breaks = names(f1_cols), name = NULL) + scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "F1. Integrated evidence by biological theme (Phase 7; descriptive)",
       subtitle = "Gene level (L1): 0 of 18,277 genes at padj < 0.05. fgsea = coordinated rank-level shift; CAMERA (estimated inter-gene correlation): 0 of 215 sets at FDR < 0.05.\nNo theme reaches ROBUST or SUPPORTED WITH CAUTION (requires Tier 1/2 sets).",
       caption = "LOO: largest share of a theme's sets that depend on one sample. ch00: sets losing significance without organelle-identical ch00 genes; LE = share of leading-edge slots that are organelle-identical.") +
  theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), legend.position = "bottom", plot.title.position = "plot",
                                       plot.subtitle = element_text(size = 7.5), plot.caption = element_text(size = 6.5, hjust = 0)) +
  guides(fill = guide_legend(nrow = 2))
ggsave(file.path(fig, "F1_integrated_evidence_overview.png"), p1, width = 11.5, height = 6.8, dpi = 300)

# ---------------------------------------------------------------- F2 fgsea vs CAMERA on identical axes
s[, theme_f := factor(theme, levels = names(theme_cols))]
lim <- ceiling(max(-log10(s$padj_fgsea)))
p2 <- ggplot(s, aes(-log10(padj_fgsea), -log10(camera_FDR), colour = theme_f, shape = concordant)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") + geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  annotate("text", x = lim, y = -log10(0.05), label = "CAMERA FDR = 0.05", hjust = 1, vjust = -0.5, size = 2.8, colour = "grey30") +
  annotate("text", x = -log10(0.05), y = lim, label = "fgsea padj = 0.05", hjust = -0.05, vjust = 1, size = 2.8, colour = "grey30") +
  geom_point(size = 1.8, alpha = 0.85) + coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
  scale_colour_manual(values = theme_cols, name = NULL, drop = TRUE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4), labels = c(`TRUE` = "same direction", `FALSE` = "opposite direction"), name = "CAMERA vs fgsea") +
  labs(x = "fgsea  -log10(padj)", y = "CAMERA  -log10(FDR), estimated inter-gene correlation",
       title = "F2. fgsea and CAMERA agree in direction but not in statistical strength",
       subtitle = sprintf("%d / %d fgsea-significant sets have the same direction in CAMERA; 0 / %d reach CAMERA FDR < 0.05 (all FDR >= %.3f).\nIdentical axes; every point lies on the floor of the plot.",
                          sum(s$concordant), nrow(s), nrow(s), min(s$camera_FDR))) +
  theme_bw(base_size = 9) + theme(legend.position = "right", legend.text = element_text(size = 7), plot.subtitle = element_text(size = 7.5))
ggsave(file.path(fig, "F2_fgsea_vs_camera.png"), p2, width = 10, height = 7.5, dpi = 300)

# ---------------------------------------------------------------- F3 theme summary
cav <- c("Ribosome / translation / RNA processing" = "26/34 sets depend on STV_infected_rep1 (LOO); 0/34 broadly consistent",
         "Photosynthesis / plastid / thylakoid" = "fgsea-only; persists after ch00 caveat (0/25 ch00-sensitive)",
         "ATP synthase / proton transport" = "organelle-identical feature class: 4/5 sets ch00-sensitive",
         "Defence / biotic response" = "4/8 sets depend on STV_free_rep2 (boundary); overlaps hormone/stress",
         "Hormone signalling" = "same direction as STV studies; gene level n.s. (ACO1 padj 0.30)",
         "Endoplasmic reticulum / protein folding" = "4/4 sets depend on STV_free_rep1 (LOO)",
         "Heme / oxidoreductase" = "mixed direction (5 up / 6 down); P450-dominated",
         "Glutathione / detoxification" = "small (3 sets); overlaps water/abiotic and heme",
         "Water deprivation / abiotic stress" = "small (3 sets); overlaps defence and glutathione",
         "Carbohydrate binding / cell wall" = "mixed direction; 2/4 sets depend on STV_free_rep2")
f3 <- m[, .(theme, dir = fcase(fgsea_direction %like% "^positive", paste0("higher rank in STV-infected (", n_sets, "/", n_sets, ")"),
                                fgsea_direction %like% "^negative", paste0("lower rank in STV-infected (", n_sets, "/", n_sets, ")"),
                                default = paste0("mixed (", fgsea_positive, " up / ", fgsea_negative, " down)")),
            dir_class = fcase(fgsea_direction %like% "^positive", "positive", fgsea_direction %like% "^negative", "negative", default = "mixed"),
            phase7_class)]
f3[, caveat := vapply(cav[theme], function(z) paste(strwrap(z, width = 48), collapse = "\n"), character(1))]
f3[, th := factor(short[theme], levels = rev(short[themes10]))]
f3l <- rbind(f3[, .(th, col = "Direction (fgsea NES)", txt = dir, fillv = dir_class)],
             f3[, .(th, col = "Phase 7 class", txt = phase7_class, fillv = phase7_class)],
             f3[, .(th, col = "Major caveat", txt = caveat, fillv = "caveat")])
f3l[, col := factor(col, levels = c("Direction (fgsea NES)", "Phase 7 class", "Major caveat"))]
p3 <- ggplot(f3l, aes(col, th, fill = fillv)) + geom_tile(colour = "white", linewidth = 0.8) + geom_text(aes(label = txt), size = 2.4) +
  scale_fill_manual(values = c(positive = "#f4a582", negative = "#92c5de", mixed = "#e0e0e0", EXPLORATORY = "#fdae61", UNSUPPORTED = "#969696", caveat = "#fafafa"),
                    breaks = c("positive", "negative", "mixed", "EXPLORATORY", "UNSUPPORTED"), name = NULL) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "F3. Biological themes: direction, Phase 7 evidence class and main caveat",
       subtitle = "Directions are coordinated rank-level shifts (fgsea), not activation or suppression. No theme is ROBUST or SUPPORTED WITH CAUTION;\nno fgsea signal survived correlation-aware validation (CAMERA); no gene is significant at padj < 0.05.") +
  theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), legend.position = "bottom", plot.title.position = "plot",
                                       plot.subtitle = element_text(size = 7.5)) +
  theme(axis.text.x.top = element_text(face = "bold"))
p3 <- p3 + theme(plot.margin = margin(5, 15, 5, 5))
ggsave(file.path(fig, "F3_theme_summary.png"), p3, width = 12, height = 5.8, dpi = 300)

# ---------------------------------------------------------------- F4 ethylene case study
eth <- eth[desc_category != "prior-study locus (other annotation)"]
eth[, cat := factor(c("biosynthesis (ACS/ACO)" = "biosynthesis\n(ACS/ACO)", "perception / signalling" = "perception /\nsignalling",
                      "ethylene-responsive (by description; mostly ERF TFs)" = "ethylene-responsive\n(mostly ERF TFs)",
                      "GO-annotated only" = "GO-annotated\nonly")[desc_category],
                    levels = rev(c("biosynthesis\n(ACS/ACO)", "perception /\nsignalling", "ethylene-responsive\n(mostly ERF TFs)", "GO-annotated\nonly")))]
eth[, lab := fcase(gene_id == "Solyc07g049530.3", sprintf("Solyc07g049530.3 ACO1\n(padj %.2f)", padj),
                   gene_id == "Solyc12g009240.1", sprintf("Solyc12g009240.1 ERF (ITAG4.0);\n'PR-protein' in Badeli 2026 (padj %.2f)", padj),
                   default = NA_character_)]
p4a <- ggplot(eth, aes(stat, cat)) + geom_vline(xintercept = 0, colour = "grey50") +
  geom_jitter(data = eth[is.na(lab)], height = 0.15, width = 0, size = 1.4, alpha = 0.6, colour = "grey30") +
  geom_point(data = eth[!is.na(lab)], colour = "#d6604d", size = 2.4) +
  geom_text_repel(data = eth[!is.na(lab)], aes(label = lab), size = 2.4, min.segment.length = 0, nudge_y = 0.35, lineheight = 0.85) +
  labs(x = "DESeq2 Wald statistic (STV-infected vs STV-free)", y = NULL,
       title = "A. Gene level: no ethylene-related gene at padj < 0.05",
       subtitle = sprintf("%d tested ethylene-related genes; %d with negative statistic; minimum padj %.3f", nrow(eth), sum(eth$stat < 0), min(eth$padj, na.rm = TRUE))) +
  theme_bw(base_size = 9)
es2 <- es[!is.na(NES)]
es2[, lab := sprintf("fgsea padj %.3f | CAMERA FDR %.3f | %s | %s", padj_fgsea, camera_FDR, tier, replicate_class)]
es2[, term_name := vapply(term_name, function(z) paste(strwrap(z, width = 22), collapse = "\n"), character(1))]
p4b <- ggplot(es2, aes(NES, reorder(term_name, NES))) + geom_col(fill = "#92c5de", width = 0.6) + geom_vline(xintercept = 0, colour = "grey50") +
  geom_text(aes(x = 0.05, label = lab), hjust = 0, size = 2.4) + xlim(min(es2$NES) * 1.05, 2.6) +
  labs(x = "fgsea NES", y = NULL, title = "B. Pathway level: directional, fgsea-only (not CAMERA-supported)",
       subtitle = "'ethylene biosynthetic process' (GO:0009693) not tested: 4 genes in the ranked universe") +
  theme_bw(base_size = 9)
p4 <- plot_grid(p4a, p4b, ncol = 1, rel_heights = c(1.25, 0.75), align = "v", axis = "l")
p4 <- plot_grid(ggdraw() + draw_label("F4. Ethylene case study: gene-level evidence vs pathway-level exploratory direction", fontface = "bold", size = 11, x = 0.01, hjust = 0),
                p4, ncol = 1, rel_heights = c(0.05, 1))
ggsave(file.path(fig, "F4_ethylene_case_study.png"), p4, width = 10, height = 8, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_22.txt"))
cat("22_phase7_figures.R finished\n")
