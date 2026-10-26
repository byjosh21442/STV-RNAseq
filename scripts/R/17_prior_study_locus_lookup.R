#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6B, Task 11: descriptive lookup of the loci named by the two prior analyses of GSE137303
# (Fukuhara et al. 2020, via the gene names quoted by Badeli et al. 2026; Badeli et al. 2026 Table 2 / text) in the
# primary Phase 4 DESeq2 results and the Phase 6 leading edges. Not a gene list, no threshold, no new test.
# Loci are matched by the Solyc locus prefix only (ITAG3.0 IDs in the source vs ITAG4.0 IDs here); the ITAG4.0
# description is reported next to the name used by the source, and differences are not resolved here.
# Usage (project root, env stv-r): Rscript scripts/R/17_prior_study_locus_lookup.R
# Inputs : results/deseq2/deseq2_results_all.tsv, results/deseq2/deseq2_results_shrunken.tsv,
#          results/functional/ranking/deseq2_wald_ranking.tsv, reference/ITAG4.0_tx2gene.tsv,
#          results/functional/go_gsea/leading_edge_significant_sets_long.tsv,
#          results/functional_robustness/themes/theme_membership_significant_sets.tsv
# Output : results/functional_robustness/prior_studies/prior_study_locus_lookup.tsv
source("scripts/R/functional_robustness_helpers.R")
out <- file.path(fr_out, "prior_studies")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Loci and names exactly as given by the sources (Badeli et al. 2026, Research Square v1, read 2026-10-26).
src <- data.table(
  locus = c("Solyc12g009240", "Solyc07g049530", "Solyc08g067160", "Solyc04g008330", "Solyc02g086990",
            "Solyc07g007770", "Solyc02g031950", "Solyc01g109010", "Solyc01g095080"),
  name_in_source = c("Pathogenesis-related (PR)-protein", "ACC oxidase 1 (ACO1)", "Uncharacterized protein (DUF247)",
                     "beta-Galactosidase 4-like", "PGIP (Polygalacturonase inhibitor)", "ERF017 (Ethylene-responsive TF)",
                     "Glyoxalase II (GLO2)", "GMP synthase [glutamine-hydrolyzing]", "NAC1"),
  source = c(rep("Badeli 2026 Table 2", 8), "Badeli 2026 text (gene highlighted by Fukuhara 2020)"))

res <- fread("results/deseq2/deseq2_results_all.tsv", sep = "\t", quote = "")
shr <- fread("results/deseq2/deseq2_results_shrunken.tsv", sep = "\t", quote = "",
             select = c("gene_id", "log2FoldChange_apeglm"))
rk  <- fread("results/functional/ranking/deseq2_wald_ranking.tsv", sep = "\t", quote = "", select = c("gene_id", "rank"))
t2g <- fread("reference/ITAG4.0_tx2gene.tsv", sep = "\t", quote = "")
le  <- fread("results/functional/go_gsea/leading_edge_significant_sets_long.tsv", sep = "\t", quote = "",
             select = c("collection", "go_id", "gene_id"))
th  <- fread("results/functional_robustness/themes/theme_membership_significant_sets.tsv", sep = "\t", quote = "",
             select = c("collection", "go_id", "theme", "tier"))
fr_check(nrow(res) == 18277L, "primary DESeq2 table changed")
gene_col <- setdiff(names(t2g), names(t2g)[1])[1]   # second column = gene ID
itag_genes <- unique(t2g[[gene_col]])

prefix <- function(x) sub("\\.[0-9]+$", "", x)
rows <- lapply(seq_len(nrow(src)), function(i) {
  g_itag <- itag_genes[prefix(itag_genes) == src$locus[i]]
  g_res  <- res[prefix(gene_id) == src$locus[i]]
  base <- data.table(locus = src$locus[i], name_in_source = src$name_in_source[i], source = src$source[i],
                     itag40_gene_id = if (length(g_itag)) paste(g_itag, collapse = ";") else NA_character_,
                     tested_18277 = nrow(g_res) > 0)
  if (!nrow(g_res)) return(base)
  g <- g_res$gene_id[1]
  sets <- unique(le[gene_id == g, .(collection, go_id)])
  thm  <- merge(sets, th, by = c("collection", "go_id"))
  cbind(base, data.table(
    itag40_description = g_res$description[1], baseMean = signif(g_res$baseMean[1], 4),
    log2FC_MLE = round(g_res$log2FoldChange_MLE[1], 3), log2FC_apeglm = round(shr[gene_id == g, log2FoldChange_apeglm], 3),
    wald_stat = round(g_res$stat[1], 3), rank_of_18238 = if (g %in% rk$gene_id) rk[gene_id == g, rank] else NA_integer_,
    pvalue = signif(g_res$pvalue[1], 3), padj = signif(g_res$padj[1], 3),
    padj_NA_reason = g_res$padj_NA_reason[1],
    significant_sets_with_gene_in_leading_edge = nrow(sets),
    themes_of_those_sets = if (nrow(thm)) paste(sort(unique(thm$theme)), collapse = "; ") else "",
    tiers_of_those_sets = if (nrow(thm)) paste(sort(unique(thm$tier)), collapse = "; ") else ""))
})
tab <- rbindlist(rows, fill = TRUE)
fr_wt(tab, file.path(out, "prior_study_locus_lookup.tsv"))
print(tab[, .(locus, name_in_source, itag40_gene_id, tested_18277, itag40_description = substr(itag40_description, 1, 45),
              log2FC_MLE, wald_stat, rank_of_18238, padj, significant_sets_with_gene_in_leading_edge)])
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_17.txt"))
cat("17_prior_study_locus_lookup.R finished\n")
