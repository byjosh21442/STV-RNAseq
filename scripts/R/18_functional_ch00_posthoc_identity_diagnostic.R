#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6B, Task 9 - POST-HOC DIAGNOSTIC (not pre-specified; does not change any classification).
# Motivation (methodological, found when reading the BLAST output, not from pathway names): blastn was run with
# -max_hsps 1, so spliced organelle genes (cDNA split into several HSPs on the organelle genome) get a single-HSP
# query coverage below the pre-specified 80 % even at ~100 % identity and fall into "organelle-similar, diverged".
# This script repeats the pre-specified ch00 sensitivity with a broader exclusion: every ranked SL4.0ch00 gene whose
# best organelle hit has identity >= 98 % and single-HSP coverage >= 50 %. It only checks whether the ch00 STOP
# assessment (> 50 % of the significant sets of a collection ch00-sensitive) would change; the frozen rule stays primary.
# Usage (project root, env stv-r): Rscript scripts/R/18_functional_ch00_posthoc_identity_diagnostic.R
# Inputs : results/functional_robustness/ch00/ch00_gene_classification.tsv, Phase 6 ranking and gene sets
# Outputs: results/functional_robustness/ch00/posthoc_identity98_{excluded_genes,vs_primary,summary}.tsv
source("scripts/R/functional_robustness_helpers.R")
run_utc <- fr_check_frozen()
out <- file.path(fr_out, "ch00")

ch <- fread(file.path(out, "ch00_gene_classification.tsv"), sep = "\t", quote = "")
rk <- fread("results/functional/ranking/deseq2_wald_ranking.tsv", sep = "\t", quote = "")
excl_pre <- fread(file.path(out, "ch00_sensitivity_excluded_genes.tsv"))$gene_id
excl <- ch[ranked == TRUE & !is.na(pident) & pident >= 98 & qcovs >= 50, gene_id]
fr_check(all(excl_pre %in% excl), "post-hoc exclusion must contain the pre-specified exclusion")
fr_wt(ch[gene_id %in% excl, .(gene_id, organelle_class, organelle, pident, qcovs, qlen, aln_length, stat,
                              n_significant_sets_in_leading_edge, blastp_organelle_gene, pre_specified_exclusion = gene_id %in% excl_pre)],
      file.path(out, "posthoc_identity98_excluded_genes.tsv"))

stats0 <- setNames(rk$stat, rk$gene_id)
stats1 <- stats0[!names(stats0) %in% excl]
sets1 <- fr_sets(fr_go_pairs(), names(stats1))
gs <- readRDS("results/functional/annotation/go_gene_sets.rds")
cmp <- rbindlist(lapply(fr_collections, function(k) {
  r1 <- fr_fgsea(sets1[[k]], stats1)
  c1 <- fr_compare_primary(k, r1, "posthoc_without_ch00_identity98")
  c1[, n_excluded_in_set := vapply(go_id, function(id) sum(gs$sets[[k]][[id]] %in% excl), 0L)]
  c1[, ch00_sensitive := significant_fdr05 %in% TRUE & status != "retained"]
  c1[, material_warnings := fr_material_warnings(r1)]
  c1
}))
fr_wt(cmp[order(collection, padj, go_id)], file.path(out, "posthoc_identity98_vs_primary.tsv"))
summ <- cmp[, .(excluded_genes = length(excl), primary_significant = sum(significant_fdr05 %in% TRUE),
                retained = sum(status == "retained", na.rm = TRUE), supported = sum(status == "supported", na.rm = TRUE),
                lost = sum(status == "lost", na.rm = TRUE), ch00_sensitive = sum(ch00_sensitive, na.rm = TRUE),
                new_significant = sum(new_significant),
                spearman_NES = round(cor(NES, NES_s, method = "spearman", use = "complete.obs"), 4),
                material_warnings = max(material_warnings)), by = collection]
summ[, fraction_ch00_sensitive := round(ch00_sensitive / primary_significant, 4)]
summ[, note := "post-hoc diagnostic; not pre-specified; classification unchanged"]
fr_wt(summ, file.path(out, "posthoc_identity98_summary.tsv"))
print(summ)
writeLines(c(paste("run_start_utc:", run_utc), capture.output(sessionInfo())), file.path(out, "sessionInfo_18.txt"))
cat("18_functional_ch00_posthoc_identity_diagnostic.R finished\n")
