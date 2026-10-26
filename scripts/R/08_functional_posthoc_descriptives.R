#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6, POST-HOC descriptive check (not pre-specified; does not change any test or criterion).
# Share of leading-edge genes located on SL4.0ch00 (unanchored sequence) in each FDR-significant GO set. Motivated by
# the observation (after the run) that chloroplast-type ATP-synthase genes on ch00 recur in many leading edges.
# Usage (project root, env stv-r): Rscript scripts/R/08_functional_posthoc_descriptives.R
# Inputs : results/functional/go_gsea/leading_edge_significant_sets_long.tsv, reference/SGN/ITAG4.0/ITAG4.0_gene_models.gff
# Output : results/functional/go_gsea/posthoc_leading_edge_ch00_share.tsv
suppressPackageStartupMessages(library(data.table))
gff <- fread(cmd = "awk -F'\\t' '$3==\"gene\"' reference/SGN/ITAG4.0/ITAG4.0_gene_models.gff", header = FALSE, sep = "\t", quote = "")
gff[, gene_id := sub("^ID=gene:([^;]+);.*$", "\\1", V9)]
ch00 <- gff[V1 == "SL4.0ch00", gene_id]
le <- fread("results/functional/go_gsea/leading_edge_significant_sets_long.tsv", sep = "\t", quote = "")
le[, on_ch00 := gene_id %in% ch00]
out <- le[, .(leading_edge_size = .N, leading_edge_on_ch00 = sum(on_ch00),
              ch00_genes = paste(gene_id[on_ch00], collapse = ";")), by = .(collection, go_id, term_name, NES, padj)]
out[, fraction_on_ch00 := round(leading_edge_on_ch00 / leading_edge_size, 3)]
setorder(out, -fraction_on_ch00, padj)
fwrite(out, "results/functional/go_gsea/posthoc_leading_edge_ch00_share.tsv", sep = "\t", quote = FALSE)
cat("ITAG4.0 genes on SL4.0ch00:", length(ch00), "\n")
cat("unique leading-edge genes:", uniqueN(le$gene_id), "; on ch00:", uniqueN(le[on_ch00 == TRUE, gene_id]), "\n")
cat("significant sets:", nrow(out), "; with >= 10% leading edge on ch00:", out[fraction_on_ch00 >= 0.10, .N],
    "; with any:", out[leading_edge_on_ch00 > 0, .N], "\n")
