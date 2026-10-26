#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6B, Task 9 (D-038): SL4.0ch00 genes - organelle similarity classes (from
# scripts/bash/ch00_organelle_blast.sh) and the pre-specified sensitivity run without near-identical organelle-type ch00 genes.
# ch00 genes are NOT removed for being unanchored. Rules fixed in docs/functional_robustness_plan.md §7.
# Usage (project root, env stv-r): Rscript scripts/R/15_functional_ch00_investigation.R
# Inputs : results/functional_robustness/ch00/{ch00_genes_gff.tsv,blastn_*.tsv,blastp_*.tsv}, reference/ITAG4.0_tx2gene.tsv,
#          reference/functional/organelle/*.cds_aa.faa, reference/SGN/ITAG4.0/ITAG4.0_gene_models.gff, Phase 6 tables
# Outputs: results/functional_robustness/ch00/*
source("scripts/R/functional_robustness_helpers.R")
run_utc <- fr_check_frozen()
out <- file.path(fr_out, "ch00")

t2g <- fread("reference/ITAG4.0_tx2gene.tsv", sep = "\t")
gff <- fread(cmd = "awk -F'\\t' '$3==\"gene\"' reference/SGN/ITAG4.0/ITAG4.0_gene_models.gff", header = FALSE, sep = "\t", quote = "")
gff <- gff[, .(gene_id = sub("^ID=gene:([^;]+);.*$", "\\1", V9), chromosome = V1, start = V4, end = V5, strand = V7)]
ch00 <- fread(file.path(out, "ch00_genes_gff.tsv"), header = FALSE, col.names = c("gene_id", "start", "end", "strand"))
fr_check(nrow(ch00) == 513L && setequal(ch00$gene_id, gff[chromosome == "SL4.0ch00", gene_id]), "ch00 gene list inconsistent with the GFF")

# organelle protein names from the RefSeq CDS FASTA headers
hdr <- unlist(lapply(Sys.glob("reference/functional/organelle/*.cds_aa.faa"), function(f) grep("^>", readLines(f), value = TRUE)))
prot <- data.table(sseqid = sub("^>(\\S+).*$", "\\1", hdr),
                   organelle_gene = sub("^.*\\[gene=([^]]+)\\].*$", "\\1", hdr),
                   organelle_protein = sub("^.*\\[protein=([^]]+)\\].*$", "\\1", hdr))
prot[organelle_gene == hdr, organelle_gene := NA_character_]

bn <- fread(file.path(out, "blastn_itag40_cdna_vs_organelle_genomes.tsv"), sep = "\t")
bn[, gene_id := t2g$GENEID[match(qseqid, t2g$TXNAME)]]
fr_check(!anyNA(bn$gene_id), "blastn query IDs not in tx2gene")
bn[, organelle := fifelse(grepl("NC_007898", sseqid), "plastid", fifelse(grepl("NC_035963", sseqid), "mitochondrion", "other"))]
setorder(bn, gene_id, -bitscore)
best <- bn[, .SD[1], by = gene_id]
best[, organelle_class := fifelse(pident >= 98 & qcovs >= 80, paste0(organelle, "-type near-identical"),
                          fifelse(pident >= 80 & qcovs >= 50, "organelle-similar, diverged", "no close organelle similarity"))]
bp <- fread(file.path(out, "blastp_ch00_proteins_vs_organelle_proteins.tsv"), sep = "\t")
bp[, gene_id := t2g$GENEID[match(qseqid, t2g$TXNAME)]]
bp <- merge(bp, prot, by = "sseqid", all.x = TRUE)
setorder(bp, gene_id, -bitscore)
bpb <- bp[, .SD[1], by = gene_id][, .(gene_id, blastp_organelle_gene = organelle_gene, blastp_organelle_protein = organelle_protein,
                                     blastp_pident = pident, blastp_query_cov_pct = qcovs)]

res <- fread("results/deseq2/deseq2_results_all.tsv", sep = "\t", quote = "", na.strings = "NA")
rk <- fread("results/functional/ranking/deseq2_wald_ranking.tsv", sep = "\t", quote = "")
le <- fread("results/functional/go_gsea/leading_edge_significant_sets_long.tsv", sep = "\t", quote = "")
le_n <- le[, .(n_significant_sets_in_leading_edge = .N), by = gene_id]

all_genes <- merge(gff, best[, .(gene_id, best_hit = sseqid, organelle, pident, aln_length = length, qlen, qcovs, organelle_class)],
                   by = "gene_id", all.x = TRUE)
all_genes[is.na(organelle_class), organelle_class := "no close organelle similarity"]
all_genes <- merge(all_genes, res[, .(gene_id, baseMean, stat, padj, padj_NA_reason, description)], by = "gene_id", all.x = TRUE)
all_genes[, tested := !is.na(baseMean)][, ranked := gene_id %in% rk$gene_id]
all_genes <- merge(all_genes, le_n, by = "gene_id", all.x = TRUE)
all_genes[is.na(n_significant_sets_in_leading_edge), n_significant_sets_in_leading_edge := 0L]
all_genes[, on_ch00 := chromosome == "SL4.0ch00"]

ch <- merge(all_genes[on_ch00 == TRUE], bpb, by = "gene_id", all.x = TRUE)
setorder(ch, -n_significant_sets_in_leading_edge, -ranked, gene_id)
fr_wt(ch, file.path(out, "ch00_gene_classification.tsv"))
fr_wt(all_genes[grepl("near-identical", organelle_class)][order(chromosome, start)], file.path(out, "organelle_type_genes_all_chromosomes.tsv"))
cls_summary <- all_genes[, .(genes = .N, tested = sum(tested), ranked = sum(ranked),
                             in_any_significant_leading_edge = sum(n_significant_sets_in_leading_edge > 0),
                             median_baseMean_ranked = round(median(baseMean[ranked]), 1)),
                         by = .(location = fifelse(on_ch00, "SL4.0ch00", "SL4.0ch01-12"), organelle_class)]
fr_wt(cls_summary[order(location, organelle_class)], file.path(out, "organelle_class_summary.tsv"))
print(cls_summary[order(location, organelle_class)])

# --- sensitivity: ranking and gene sets without near-identical organelle-type ch00 genes -------------------------------
excl <- ch[ranked == TRUE & grepl("near-identical", organelle_class), gene_id]
fr_wt(data.table(gene_id = excl), file.path(out, "ch00_sensitivity_excluded_genes.tsv"))
stats0 <- setNames(rk$stat, rk$gene_id)
stats1 <- stats0[!names(stats0) %in% excl]
pairs <- fr_go_pairs()
sets1 <- fr_sets(pairs, names(stats1))
gs <- readRDS("results/functional/annotation/go_gene_sets.rds")
cmp <- list(); allr <- list()
for (k in fr_collections) {
  r1 <- fr_fgsea(sets1[[k]], stats1)
  c1 <- fr_compare_primary(k, r1, "without_organelle_type_ch00")
  c1[, n_excluded_in_set := vapply(go_id, function(id) sum(gs$sets[[k]][[id]] %in% excl), 0L)]
  c1[, ch00_sensitive := significant_fdr05 %in% TRUE & status != "retained"]
  c1[, material_warnings := fr_material_warnings(r1)]
  cmp[[k]] <- c1
}
cmp <- rbindlist(cmp)
fr_wt(cmp[order(collection, padj, go_id)], file.path(out, "ch00_sensitivity_vs_primary.tsv"))
summ <- cmp[, .(excluded_genes = length(excl), primary_significant = sum(significant_fdr05 %in% TRUE),
                retained = sum(status == "retained", na.rm = TRUE), supported = sum(status == "supported", na.rm = TRUE),
                lost = sum(status == "lost", na.rm = TRUE), ch00_sensitive = sum(ch00_sensitive, na.rm = TRUE),
                new_significant = sum(new_significant),
                spearman_NES = round(cor(NES, NES_s, method = "spearman", use = "complete.obs"), 4),
                material_warnings = max(material_warnings)), by = collection]
summ[, fraction_ch00_sensitive := round(ch00_sensitive / primary_significant, 4)]
fr_wt(summ, file.path(out, "ch00_sensitivity_summary.tsv"))
print(summ)
fr_check(all(summ$fraction_ch00_sensitive <= 0.5), "ch00 artifact materially affects the findings (> 50% of significant sets ch00-sensitive)")
writeLines(c(paste("run_start_utc:", run_utc), capture.output(sessionInfo())), file.path(out, "sessionInfo_15.txt"))
cat("15_functional_ch00_investigation.R finished\n")
