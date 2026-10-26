#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 7 (Task 14): descriptive GO redundancy among the 215 fgsea-significant GO sets.
# No test and no new significance threshold. Quantities: Jaccard of full set membership (ranked universe), Jaccard of
# fgsea leading edges, GO ancestor/descendant pairs (is_a + part_of, GO release used in Phase 6), and descriptive modules
# from average-linkage clustering of 1 - leading-edge Jaccard, done separately for positive and negative NES.
# The cut height (0.7) was fixed before inspecting module content; module counts at 0.5 and 0.9 are reported as sensitivity.
# Usage (project root, env stv-r): Rscript scripts/R/21_go_redundancy_modules.R   (after 20_biological_integration.R)
# Outputs: results/biological_integration/go_redundancy_{pairs,summary,set_modules,module_summary}.tsv,
#          figures/biological_integration/F5_go_redundancy_heatmap.png, F5b_go_modules_overview.png, sessionInfo_21.txt
suppressPackageStartupMessages({ library(data.table); library(ontologyIndex); library(pheatmap); library(ggplot2) })
options(warn = 1)
out <- "results/biological_integration"; fig <- "figures/biological_integration"
dir.create(fig, showWarnings = FALSE, recursive = TRUE)
fwt <- function(x, f) fwrite(x, file.path(out, f), sep = "\t", quote = FALSE, na = "NA")
cut_h <- 0.7; cut_sens <- c(0.5, 0.9)
theme_cols <- c("Ribosome / translation / RNA processing" = "#1f78b4", "Photosynthesis / plastid / thylakoid" = "#33a02c",
                "ATP synthase / proton transport" = "#b2df8a", "Defence / biotic response" = "#e31a1c", "Hormone signalling" = "#ff7f00",
                "Endoplasmic reticulum / protein folding" = "#6a3d9a", "Heme / oxidoreductase" = "#b15928",
                "Glutathione / detoxification" = "#fb9a99", "Water deprivation / abiotic stress" = "#a6cee3",
                "Carbohydrate binding / cell wall" = "#cab2d6", "Other (no theme regex matched)" = "#d9d9d9")

s <- fread(file.path(out, "theme_set_level_evidence.tsv"), sep = "\t", quote = "")
stopifnot(nrow(s) == 215)
s[, set_key := paste(collection, go_id, sep = ":")]
gs <- readRDS("results/functional/annotation/go_gene_sets.rds")
le <- fread("results/functional/go_gsea/leading_edge_significant_sets_long.tsv", sep = "\t", quote = "", select = c("collection", "go_id", "gene_id"))
mem <- lapply(seq_len(nrow(s)), function(i) gs$sets[[s$collection[i]]][[s$go_id[i]]])
lem <- split(le$gene_id, paste(le$collection, le$go_id, sep = ":"))[s$set_key]
stopifnot(all(lengths(mem) == s$size), all(lengths(lem) == s$leading_edge_size))
jac <- function(L) {
  n <- length(L); M <- diag(1, n)
  for (i in seq_len(n - 1)) for (j in (i + 1):n) M[i, j] <- M[j, i] <- length(intersect(L[[i]], L[[j]])) / length(union(L[[i]], L[[j]]))
  dimnames(M) <- list(s$set_key, s$set_key); M
}
J_full <- jac(mem); J_le <- jac(lem)

ont <- get_ontology("reference/functional/go/release_2026-08-05/go-basic.obo", propagate_relationships = c("is_a", "part_of"))
anc <- ont$ancestors[s$go_id]
stopifnot(!any(vapply(anc, is.null, logical(1))))

ij <- which(upper.tri(J_full), arr.ind = TRUE)
pairs <- data.table(set_a = s$set_key[ij[, 1]], set_b = s$set_key[ij[, 2]], term_a = s$term_name[ij[, 1]], term_b = s$term_name[ij[, 2]],
                    theme_a = s$theme[ij[, 1]], theme_b = s$theme[ij[, 2]],
                    same_direction = sign(s$NES[ij[, 1]]) == sign(s$NES[ij[, 2]]),
                    jaccard_full = round(J_full[ij], 4), jaccard_le = round(J_le[ij], 4))
pairs[, go_ancestor_relation := mapply(function(a, b) b %in% anc[[a]] || a %in% anc[[b]], s$go_id[ij[, 1]], s$go_id[ij[, 2]])]
fwt(pairs[jaccard_full >= 0.25 | jaccard_le >= 0.25 | go_ancestor_relation][order(-jaccard_le)], "go_redundancy_pairs.tsv")

mods <- list(); hcs <- list()
for (d in c("positive", "negative")) {
  idx <- which(if (d == "positive") s$NES > 0 else s$NES < 0)
  hc <- hclust(as.dist(1 - J_le[idx, idx]), method = "average"); hcs[[d]] <- hc
  cl <- cutree(hc, h = cut_h)
  sz <- sort(table(cl), decreasing = TRUE)
  relab <- setNames(sprintf("%s%02d", if (d == "positive") "P" else "N", seq_along(sz)), names(sz))
  sens <- sapply(cut_sens, function(h) length(unique(cutree(hc, h = h))))
  mods[[d]] <- data.table(set_key = s$set_key[idx], module = unname(relab[as.character(cl)]), n_modules_h0.5 = sens[1], n_modules_h0.9 = sens[2])
}
md <- rbindlist(mods)
s <- merge(s, md, by = "set_key")
fwt(s[, .(collection, go_id, term_name, theme, NES, padj_fgsea, camera_FDR, tier, module)][order(module, padj_fgsea)], "go_redundancy_set_modules.tsv")

msum <- s[, {
  k <- set_key; sub_le <- J_le[k, k, drop = FALSE]; sub_f <- J_full[k, k, drop = FALSE]; ut <- upper.tri(sub_le)
  .(direction = fifelse(NES[1] > 0, "positive", "negative"), n_sets = .N, n_BP = sum(collection == "BP"), n_MF = sum(collection == "MF"),
    n_CC = sum(collection == "CC"), label_term = term_name[which.min(padj_fgsea)], label_go_id = go_id[which.min(padj_fgsea)],
    themes = paste(names(sort(table(theme), decreasing = TRUE)), sort(table(theme), decreasing = TRUE), sep = " x", collapse = "; "),
    union_le_genes = length(unique(unlist(lem[k]))), mean_jaccard_le = if (any(ut)) round(mean(sub_le[ut]), 3) else NA_real_,
    mean_jaccard_full = if (any(ut)) round(mean(sub_f[ut]), 3) else NA_real_, min_padj_fgsea = signif(min(padj_fgsea), 3),
    min_camera_FDR = signif(min(camera_FDR), 3))
}, by = module][order(direction, -n_sets)]
fwt(msum, "go_module_summary.tsv")

ut_all <- upper.tri(J_le)
summ <- data.table(
  metric = c("significant_sets", "set_pairs", "pairs_jaccard_full_ge_0.5", "pairs_jaccard_le_ge_0.5", "pairs_identical_membership",
             "pairs_go_ancestor_relation", "sets_with_any_ancestor_or_descendant_among_215", "modules_h0.7_positive", "modules_h0.7_negative",
             "singleton_modules_h0.7", "modules_h0.5_total", "modules_h0.9_total", "module_cut_height"),
  value = c(nrow(s), nrow(pairs), sum(pairs$jaccard_full >= 0.5), sum(pairs$jaccard_le >= 0.5), sum(pairs$jaccard_full == 1),
            sum(pairs$go_ancestor_relation), length(unique(c(pairs[go_ancestor_relation == TRUE, set_a], pairs[go_ancestor_relation == TRUE, set_b]))),
            sum(msum$direction == "positive"), sum(msum$direction == "negative"), sum(msum$n_sets == 1),
            mods$positive$n_modules_h0.5[1] + mods$negative$n_modules_h0.5[1], mods$positive$n_modules_h0.9[1] + mods$negative$n_modules_h0.9[1], cut_h))
fwt(summ, "go_redundancy_summary.tsv")

# F5: leading-edge Jaccard heatmap, ordered by the within-direction clustering
ord <- c(hcs$positive$labels[hcs$positive$order], hcs$negative$labels[hcs$negative$order])
stopifnot(length(ord) == 215, setequal(ord, s$set_key))
ann <- data.frame(direction = ifelse(s$NES > 0, "positive NES", "negative NES"), theme = s$theme, row.names = s$set_key)
ann_cols <- list(direction = c("positive NES" = "#b2182b", "negative NES" = "#2166ac"), theme = theme_cols[names(theme_cols) %in% s$theme])
png(file.path(fig, "F5_go_redundancy_heatmap.png"), width = 3000, height = 2700, res = 300)
pheatmap(J_le[ord, ord], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
         annotation_row = ann, annotation_col = ann, annotation_colors = ann_cols, annotation_legend = TRUE,
         color = colorRampPalette(c("white", "#fee08b", "#d73027"))(50), breaks = seq(0, 1, length.out = 51), border_color = NA,
         gaps_row = sum(s$NES > 0), gaps_col = sum(s$NES > 0),
         main = "Leading-edge Jaccard similarity, 215 fgsea-significant GO sets (descriptive; no test)", fontsize = 8)
dev.off()

# F5b: modules with >= 3 sets, coloured by theme composition
big <- msum[n_sets >= 3]
comp <- s[module %in% big$module, .N, by = .(module, theme)]
comp <- merge(comp, big[, .(module, label = paste0(module, ": ", substr(label_term, 1, 45), " (", n_sets, ")"), direction, n_sets)], by = "module")
comp[, label := factor(label, levels = unique(comp[order(direction, n_sets)]$label))]
p <- ggplot(comp, aes(x = N, y = label, fill = theme)) + geom_col(width = 0.75) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = theme_cols, name = "Phase 6B theme") +
  labs(x = "Number of fgsea-significant GO sets in module", y = NULL,
       title = "Descriptive GO modules (average linkage on 1 - leading-edge Jaccard, cut 0.7)",
       subtitle = "Modules with >= 3 sets; label = most significant set. Modules are redundancy groups, not independent findings.\nNo set is significant under CAMERA (FDR >= 0.997).") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", legend.text = element_text(size = 7)) + guides(fill = guide_legend(ncol = 2))
ggsave(file.path(fig, "F5b_go_modules_overview.png"), p, width = 10, height = 8, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_21.txt"))
print(summ); print(msum[n_sets >= 3, .(module, n_sets, label_term, themes, mean_jaccard_le)])
cat("\n21_go_redundancy_modules.R finished\n")
