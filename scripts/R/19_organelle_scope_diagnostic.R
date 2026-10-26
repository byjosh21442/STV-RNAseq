#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6B closure, organelle assessment (Tasks 7 - 8). POST-HOC, descriptive; no quantification, no test.
# (1) Corrected organelle-identity classification of ITAG4.0 genes from all blastn HSPs: coverage = union of query intervals
#     covered by HSPs with identity >= 98 % (criteria frozen in docs/organelle_requantification_assessment.md section 1.4).
# (2) Organelle coverage of the 177 SL4.0ch00 segments.
# (3) Scope of the ch00/organelle ambiguity from EXISTING outputs only: tested/ranked/leading-edge membership, baseMean,
#     Salmon read share per class (primary quant.sf, read-only), Salmon exact-duplicate removal, theme contribution and the
#     status of each set in the existing Phase 6B ch00 sensitivity runs (pre-specified and post-hoc).
# Usage (project root, env stv-r): Rscript scripts/R/19_organelle_scope_diagnostic.R
# Inputs : results/organelle_assessment/blast/*.tsv (scripts/bash/organelle_multihsp_blast.sh), reference/ITAG4.0_tx2gene.tsv,
#          reference/SGN/ITAG4.0/ITAG4.0_gene_models.gff, results/deseq2/deseq2_results_all.tsv, Phase 6 ranking and leading
#          edges, results/functional_robustness/{themes,ch00}/*.tsv, results/quantification/salmon/*/{quant.sf,aux_info/meta_info.json},
#          results/genome/index/salmon_index/duplicate_clusters.tsv
# Outputs: results/organelle_assessment/*.tsv
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
options(warn = 1)
stopifnot(file.exists("results/organelle_assessment/prespecification/prespecification_checksum.tsv"))
ck <- fread("results/organelle_assessment/prespecification/prespecification_checksum.tsv", colClasses = "character")
stopifnot(identical(sub(" .*$", "", system2("sha256sum", ck$file[1], stdout = TRUE)), ck$sha256[1]))
out <- "results/organelle_assessment"
b <- file.path(out, "blast")
fwt <- function(x, f) fwrite(x, file.path(out, f), sep = "\t", quote = FALSE, na = "NA")
id_min <- 98; cov_near <- 0.80; cov_partial <- 0.50   # frozen section 1.4 (same thresholds as D-038)
cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
          "evalue", "bitscore", "qlen", "slen")
org_of <- function(s) fifelse(grepl("NC_007898", s), "plastid", fifelse(grepl("NC_035963", s), "mitochondrion", "other"))
union_len <- function(a, z) {                         # total length of the union of closed intervals [a, z]
  o <- order(a); a <- a[o]; z <- z[o]; tot <- 0; cs <- a[1]; ce <- z[1]
  for (i in seq_along(a)[-1]) if (a[i] <= ce + 1) ce <- max(ce, z[i]) else { tot <- tot + ce - cs + 1; cs <- a[i]; ce <- z[i] }
  tot + ce - cs + 1
}
cover <- function(h) h[pident >= id_min, .(n_hsp = .N, union_q = union_len(pmin(qstart, qend), pmax(qstart, qend)),
                                           union_s = union_len(pmin(sstart, send), pmax(sstart, send)),
                                           w_ident = round(sum(pident * length) / sum(length), 3),
                                           max_single_hsp = max(abs(qend - qstart) + 1), qlen = qlen[1], slen = slen[1]),
                       by = .(qseqid, sseqid)]

# ---- (1) transcripts -------------------------------------------------------------------------------------------------
h <- fread(file.path(b, "blastn_itag40_cdna_vs_organelle_allhsps.tsv"), col.names = cols)
cv <- cover(h)
cv[, `:=`(organelle = org_of(sseqid), coverage = round(union_q / qlen, 4), single_hsp_coverage = round(max_single_hsp / qlen, 4))]
fwt(cv[order(qseqid, -coverage)], "multihsp_transcript_organelle_coverage.tsv")
best <- cv[order(qseqid, -coverage, -w_ident)][, .SD[1], by = qseqid]
best[, corrected_class := fifelse(coverage >= cov_near, paste0(organelle, "-type near-identical"),
                          fifelse(coverage >= cov_partial, paste0(organelle, "-type partial (50-80 %)"), "below 50 % organelle-identical coverage"))]

t2g <- fread("reference/ITAG4.0_tx2gene.tsv")
gff <- fread(cmd = "awk -F'\\t' '$3==\"gene\"' reference/SGN/ITAG4.0/ITAG4.0_gene_models.gff", header = FALSE, sep = "\t", quote = "")
gff <- gff[, .(gene_id = sub("^ID=gene:([^;]+);.*$", "\\1", V9), chromosome = V1)]
genes <- merge(gff, t2g[, .(gene_id = GENEID, transcript_id = TXNAME)], by = "gene_id", all.x = TRUE)
stopifnot(nrow(genes) == 34075L, !anyDuplicated(genes$gene_id))
genes <- merge(genes, best[, .(transcript_id = qseqid, organelle, coverage, single_hsp_coverage, w_ident, n_hsp, corrected_class)],
               by = "transcript_id", all.x = TRUE)
genes[is.na(corrected_class), corrected_class := "no organelle hit"]
res <- fread("results/deseq2/deseq2_results_all.tsv", sep = "\t", quote = "", select = c("gene_id", "baseMean", "stat", "padj"))
rk  <- fread("results/functional/ranking/deseq2_wald_ranking.tsv", sep = "\t", quote = "", select = "gene_id")
le  <- fread("results/functional/go_gsea/leading_edge_significant_sets_long.tsv", sep = "\t", quote = "", select = c("collection", "go_id", "gene_id"))
genes <- merge(genes, res, by = "gene_id", all.x = TRUE)
genes[, `:=`(tested = !is.na(baseMean), ranked = gene_id %in% rk$gene_id, on_ch00 = chromosome == "SL4.0ch00")]
genes <- merge(genes, le[, .(n_significant_sets_in_leading_edge = .N), by = gene_id], by = "gene_id", all.x = TRUE)
genes[is.na(n_significant_sets_in_leading_edge), n_significant_sets_in_leading_edge := 0L]
c6 <- fread("results/functional_robustness/ch00/ch00_gene_classification.tsv", sep = "\t", quote = "", select = c("gene_id", "organelle_class"))
ph <- fread("results/functional_robustness/ch00/posthoc_identity98_excluded_genes.tsv", sep = "\t", quote = "", select = "gene_id")
pre <- fread("results/functional_robustness/ch00/ch00_sensitivity_excluded_genes.tsv")$gene_id
genes <- merge(genes, c6[, .(gene_id, phase6b_class = organelle_class)], by = "gene_id", all.x = TRUE)
genes[, `:=`(phase6b_prespecified_excluded = gene_id %in% pre, phase6b_posthoc_excluded = gene_id %in% ph$gene_id)]
dup <- fread("results/genome/index/salmon_index/duplicate_clusters.tsv")
genes[, salmon_duplicate_removed := transcript_id %in% dup$DuplicateRef]
genes[, salmon_duplicate_retained_copy := transcript_id %in% dup$RetainedRef]
setorder(genes, -on_ch00, -n_significant_sets_in_leading_edge, gene_id)
# Written for ch00 genes and anchored genes with any >= 98 % organelle HSP only (keeps the table small); all other genes
# are anchored "no organelle hit" and are still counted in class_summary_by_rule.tsv.
fwt(genes[on_ch00 | corrected_class != "no organelle hit"], "gene_organelle_classification_corrected.tsv")

near <- function(x) grepl("near-identical", x)
part <- function(x) grepl("partial", x)
loc <- function(x) fifelse(x, "SL4.0ch00", "SL4.0ch01-12")
cls <- rbind(
  genes[!is.na(phase6b_class) | !on_ch00, .(rule = "Phase 6B pre-specified (single HSP; ch00 only)", class = fifelse(is.na(phase6b_class), "not classified (anchored)", phase6b_class), location = loc(on_ch00), tested, ranked, le = n_significant_sets_in_leading_edge > 0)][location == "SL4.0ch00"],
  genes[, .(rule = "corrected (union of HSPs >= 98 %)", class = corrected_class, location = loc(on_ch00), tested, ranked, le = n_significant_sets_in_leading_edge > 0)])
cls_sum <- cls[, .(genes = .N, tested = sum(tested), ranked = sum(ranked), in_any_significant_leading_edge = sum(le)), by = .(rule, location, class)]
cls_sum <- rbind(cls_sum, data.table(rule = "Phase 6B post-hoc (single HSP identity >= 98 %, coverage >= 50 %; ranked ch00 only)", location = "SL4.0ch00",
  class = "excluded", genes = nrow(ph), tested = nrow(ph), ranked = nrow(ph), in_any_significant_leading_edge = sum(genes[phase6b_posthoc_excluded == TRUE, n_significant_sets_in_leading_edge > 0])))
fwt(cls_sum[order(rule, location, class)], "class_summary_by_rule.tsv")
cross <- genes[on_ch00 == TRUE, .N, by = .(phase6b_class, corrected_class)][order(phase6b_class, corrected_class)]
fwt(cross, "ch00_class_crosstab_phase6b_vs_corrected.tsv")
total_baseMean <- sum(genes[tested == TRUE, baseMean])
bm <- genes[tested == TRUE, .(genes = .N, baseMean_min = round(min(baseMean), 1), baseMean_median = round(median(baseMean), 1),
                              baseMean_max = round(max(baseMean), 1), baseMean_sum = round(sum(baseMean), 0)), by = .(location = loc(on_ch00), corrected_class)]
bm[, share_of_total_baseMean := round(baseMean_sum / total_baseMean, 5)]
fwt(bm[order(location, corrected_class)], "basemean_by_corrected_class.tsv")

# ---- (2) ch00 segments ---------------------------------------------------------------------------------------------------
seg <- fread(file.path(b, "ch00_segments.tsv"))
hs <- fread(file.path(b, "blastn_ch00_segments_vs_organelle_allhsps.tsv"), col.names = cols)
sc <- cover(hs)
sc[, `:=`(organelle = org_of(sseqid), segment = sub("_.*$", "", qseqid))]
sc <- dcast(sc, segment ~ organelle, value.var = c("union_q", "union_s"), fill = 0)
seg <- merge(seg, sc, by = "segment", all.x = TRUE)
for (v in setdiff(names(seg), c("segment", "start", "end", "length"))) set(seg, which(is.na(seg[[v]])), v, 0)
gpos <- fread(cmd = "awk -F'\\t' '$1==\"SL4.0ch00\" && $3==\"gene\"' reference/SGN/ITAG4.0/ITAG4.0_gene_models.gff", header = FALSE, sep = "\t", quote = "")
gpos <- gpos[, .(gene_id = sub("^ID=gene:([^;]+);.*$", "\\1", V9), gstart = V4)]
gpos <- merge(gpos, genes[, .(gene_id, corrected_class, tested)], by = "gene_id")
gpos[, segment := seg$segment[findInterval(gstart, seg$start)]]
gs <- gpos[, .(genes = .N, genes_near_identical = sum(near(corrected_class)), genes_partial = sum(part(corrected_class)), genes_tested = sum(tested)), by = segment]
seg <- merge(seg, gs, by = "segment", all.x = TRUE)
for (v in c("genes", "genes_near_identical", "genes_partial", "genes_tested")) set(seg, which(is.na(seg[[v]])), v, 0L)
pl <- grep("plastid", names(seg), value = TRUE); mi <- grep("mitochondrion", names(seg), value = TRUE)
seg[, `:=`(plastid_identical_fraction_of_segment = round(get(grep("^union_q_plastid", names(seg), value = TRUE)) / length, 4),
           mito_identical_fraction_of_segment = round(get(grep("^union_q_mitochondrion", names(seg), value = TRUE)) / length, 4),
           plastome_fraction_covered = round(get(grep("^union_s_plastid", names(seg), value = TRUE)) / 155461, 4),
           mitogenome_fraction_covered = round(get(grep("^union_s_mitochondrion", names(seg), value = TRUE)) / 446257, 4))]
fwt(seg[order(start)], "ch00_segment_organelle_coverage.tsv")
segsum <- seg[, .(segments = .N, bp = sum(length),
                  bp_plastid_identical = sum(get(grep("^union_q_plastid", names(seg), value = TRUE))),
                  bp_mito_identical = sum(get(grep("^union_q_mitochondrion", names(seg), value = TRUE))),
                  segments_ge50pct_plastid = sum(plastid_identical_fraction_of_segment >= 0.5),
                  segments_ge50pct_mito = sum(mito_identical_fraction_of_segment >= 0.5),
                  segments_ge90pct_organelle = sum(plastid_identical_fraction_of_segment + mito_identical_fraction_of_segment >= 0.9))]
fwt(segsum, "ch00_segment_summary.tsv")

# ---- (3) Salmon read share and theme contribution (existing outputs only) ---------------------------------------------------
smp <- c("STV_free_rep1", "STV_free_rep2", "STV_free_rep3", "STV_infected_rep1", "STV_infected_rep2", "STV_infected_rep3")
grp <- genes[, .(transcript_id, grp = fifelse(on_ch00 & near(corrected_class), "ch00 near-identical",
                                       fifelse(on_ch00 & part(corrected_class), "ch00 partial",
                                       fifelse(on_ch00, "ch00 other", fifelse(near(corrected_class) | part(corrected_class),
                                                                              "ch01-12 near-identical or partial", "ch01-12 other")))))]
rs <- rbindlist(lapply(smp, function(s) {
  q <- fread(file.path("results/quantification/salmon", s, "quant.sf"), select = c("Name", "NumReads"))
  m <- fromJSON(file.path("results/quantification/salmon", s, "aux_info/meta_info.json"))
  q <- merge(q, grp, by.x = "Name", by.y = "transcript_id", all.x = TRUE)
  r <- q[, .(reads = sum(NumReads)), by = grp]
  r[, `:=`(sample = s, fraction_of_mapped = round(reads / m$num_mapped, 6), num_processed = m$num_processed, num_mapped = m$num_mapped,
           decoy_fragments = m$num_decoy_fragments, decoy_fraction_of_processed = round(m$num_decoy_fragments / m$num_processed, 4))]
}))
fwt(rs[order(sample, grp)], "salmon_read_share_by_class.tsv")

th <- fread("results/functional_robustness/themes/theme_membership_significant_sets.tsv", sep = "\t", quote = "",
            select = c("collection", "go_id", "theme", "tier", "ch00_sensitive"))
ph_cmp <- fread("results/functional_robustness/ch00/posthoc_identity98_vs_primary.tsv", sep = "\t", quote = "",
                select = c("collection", "go_id", "ch00_sensitive"))
setnames(ph_cmp, "ch00_sensitive", "posthoc_ch00_sensitive")
lel <- merge(le, genes[, .(gene_id, on_ch00, corrected_class)], by = "gene_id")
le_set <- lel[, .(le_genes = .N, le_org_near = sum(near(corrected_class)), le_org_near_or_partial = sum(near(corrected_class) | part(corrected_class)),
                  le_ch00 = sum(on_ch00)), by = .(collection, go_id)]
ts <- merge(merge(th, le_set, by = c("collection", "go_id"), all.x = TRUE), ph_cmp, by = c("collection", "go_id"), all.x = TRUE)
tc <- ts[, .(sets = .N, sets_with_org_near_in_LE = sum(le_org_near > 0), sets_with_org_near_or_partial_in_LE = sum(le_org_near_or_partial > 0),
             LE_gene_slots = sum(le_genes), LE_slots_org_near_or_partial = sum(le_org_near_or_partial),
             ch00_sensitive_prespecified = sum(ch00_sensitive == TRUE), ch00_sensitive_posthoc = sum(posthoc_ch00_sensitive == TRUE, na.rm = TRUE)), by = theme]
tc[, share_LE_slots_org := round(LE_slots_org_near_or_partial / LE_gene_slots, 4)]
fwt(tc[order(-sets)], "theme_contribution_organelle_identical.tsv")
fwt(ts[order(collection, go_id)], "set_level_organelle_identical_leading_edge.tsv")
allsets <- unique(ts[, .(collection, go_id, le_org_near_or_partial)])
ov <- allsets[, .(significant_sets = .N, sets_with_org_identical_LE_gene = sum(le_org_near_or_partial > 0)), by = collection]
fwt(ov, "collection_overview_organelle_identical.tsv")
print(cls_sum[order(rule, location, class)]); print(cross); print(bm); print(segsum); print(rs[order(grp, sample)]); print(tc[order(-sets)]); print(ov)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_19.txt"))
cat("19_organelle_scope_diagnostic.R finished\n")
