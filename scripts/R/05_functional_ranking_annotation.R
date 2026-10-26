#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6, Tasks 1-4: freeze the DESeq2 Wald-statistic ranking, audit the functional annotation
# resources, build GO gene sets (BP/MF/CC) and audit Plant Reactome tomato pathways.
# NO enrichment test is run here: everything this script reports is available before any GSEA result is seen.
# Usage (project root, env stv-r): Rscript scripts/R/05_functional_ranking_annotation.R
# Inputs : results/deseq2/deseq2_results_all.tsv; results/deseq2/deseq2_objects.rds (direction check only);
#          reference/ITAG4.0_tx2gene.tsv; reference/SGN/ITAG4.0/ITAG4.0_{goterms.txt,proteins.fasta};
#          reference/functional/* (scripts/bash/download_functional_annotation.sh)
# Outputs: results/functional/ranking/*, results/functional/annotation/*, results/functional/plant_reactome_gsea/*,
#          figures/functional/ranking_wald_statistic.png, figures/functional/go_gene_set_sizes.png,
#          results/functional/annotation/go_gene_sets.rds (git-ignored; read by 06_functional_gsea.R)
# Decisions: docs/decisions.md D-023 (ranking), D-024 (GO source), D-025 (GO processing), D-026 (size range),
#            D-028 (Plant Reactome). STOP conditions halt the script with an error.
suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ontologyIndex)
})
options(warn = 1)

# --- pre-specified parameters --------------------------------------------------------------------------------
n_tested_expected <- 18277L
min_size <- 15L    # D-026: GSEA User Guide default; fgsea vignette examples
max_size <- 500L   # D-026
go_roots <- c(BP = "GO:0008150", MF = "GO:0003674", CC = "GO:0005575")
ns_code  <- c(biological_process = "BP", molecular_function = "MF", cellular_component = "CC")
contrast_regex <- "^Wald statistic: condition STV[ _]infected vs STV[ _]free$"

res_tsv   <- "results/deseq2/deseq2_results_all.tsv"
res_rds   <- "results/deseq2/deseq2_objects.rds"
t2g_tsv   <- "reference/ITAG4.0_tx2gene.tsv"
go40_txt  <- "reference/SGN/ITAG4.0/ITAG4.0_goterms.txt"
prot40_fa <- "reference/SGN/ITAG4.0/ITAG4.0_proteins.fasta"
plaza_gz  <- "reference/functional/plaza_dicots_05/go.sly.csv.gz"
go_obo    <- "reference/functional/go/release_2026-08-05/go-basic.obo"
go41_txt  <- "reference/functional/sgn_itag4.1/ITAG4.1_goterms.txt"
prot41_fa <- "reference/functional/sgn_itag4.1/ITAG4.1_proteins.fasta"
pr_tab    <- "reference/functional/plant_reactome/current/gene_ids_by_pathway_and_species.tab"

out_rank <- "results/functional/ranking"
out_ann  <- "results/functional/annotation"
out_pr   <- "results/functional/plant_reactome_gsea"
fig_dir  <- "figures/functional"
for (d in c(out_rank, out_ann, out_pr, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

summ <- data.table(section = character(), metric = character(), value = character())
add <- function(section, metric, value) {
  summ <<- rbind(summ, data.table(section = section, metric = metric, value = paste(value, collapse = ",")))
}
check <- function(ok, msg) if (!isTRUE(ok)) stop("STOP CONDITION: ", msg, call. = FALSE)
wt <- function(x, file) fwrite(x, file, sep = "\t", quote = FALSE, na = "NA")
q5 <- function(x) if (length(x)) paste(quantile(x, c(0, 0.25, 0.5, 0.75, 1), names = FALSE, type = 7), collapse = "/") else "NA"
read_fasta <- function(path) {
  x <- readLines(path)
  h <- startsWith(x, ">")
  id <- sub("^>(\\S+).*$", "\\1", x[h])
  grp <- cumsum(h)[!h]
  seqs <- vapply(split(x[!h], grp), paste, "", collapse = "")
  out <- setNames(rep("", length(id)), id)
  out[as.integer(names(seqs))] <- seqs
  out
}
read_id_go <- function(path) {   # SGN goterms format: ID<TAB>GO:...,GO:...  (GO column may be absent)
  p <- strsplit(readLines(path), "\t", fixed = TRUE)
  data.table(id = vapply(p, `[`, "", 1),
             go = vapply(p, function(v) if (length(v) >= 2) v[2] else "", ""))
}
locus_of <- function(id) sub("[.].*$", "", id)   # diagnostic only (never used to map)

# =================================================================================================================
# TASK 1 - freeze the ranking
# =================================================================================================================
res <- fread(res_tsv, sep = "\t", quote = "", na.strings = "NA")
check(nrow(res) == n_tested_expected, paste("DESeq2 result rows:", nrow(res)))
check(!anyDuplicated(res$gene_id), "duplicated gene IDs in the DESeq2 result")
add("ranking", "input_file", res_tsv)
add("ranking", "genes_tested", nrow(res))
add("ranking", "duplicate_gene_ids", sum(duplicated(res$gene_id)))

# direction: description of the stored Wald statistic, identical values, sign(stat) == sign(MLE log2FC)
obj <- readRDS(res_rds)
stat_desc <- mcols(obj$res)$description[colnames(obj$res) == "stat"]
check(grepl(contrast_regex, stat_desc), paste("unexpected Wald statistic description:", stat_desc))
m <- match(res$gene_id, rownames(obj$res))
check(!anyNA(m), "DESeq2 TSV genes missing from the stored results object")
check(isTRUE(all.equal(obj$res$stat[m], res$stat)), "TSV stat differs from the stored DESeq2 results")
check(isTRUE(all.equal(res$stat, res$log2FoldChange_MLE / res$lfcSE_MLE)), "stat != log2FoldChange_MLE / lfcSE_MLE")
check(all(sign(res$stat) == sign(res$log2FoldChange_MLE)), "sign(stat) differs from sign(log2FC)")
add("ranking", "stat_description_in_DESeq2_object", stat_desc)
add("ranking", "direction_convention", "stat = log2FC_MLE / lfcSE; positive = higher in STV_infected, negative = higher in STV_free")
add("ranking", "sign_stat_equals_sign_log2FC_MLE", TRUE)
rm(obj); invisible(gc())

# universe: a gene is ranked if its Wald statistic is finite AND DESeq2 reported a p-value for it
res[, stat_finite := is.finite(stat)]
res[, exclusion_reason := fifelse(!stat_finite, "stat_NA_or_Inf",
                           fifelse(is.na(pvalue), paste0("pvalue_NA_", fifelse(is.na(padj_NA_reason) | padj_NA_reason == "",
                                                                               "unspecified", padj_NA_reason)), NA_character_))]
excl <- res[!is.na(exclusion_reason)]
rk <- res[is.na(exclusion_reason)]
add("ranking", "stat_NA_or_Inf", sum(!res$stat_finite))
add("ranking", "excluded_genes", nrow(excl))
for (r in unique(excl$exclusion_reason)) add("ranking", paste0("excluded_", r), sum(excl$exclusion_reason == r))
check(all(excl$exclusion_reason == "pvalue_NA_cooks_outlier"), "unexpected exclusion reason (not Cook's outlier)")
setorder(rk, -stat, gene_id)
rk[, rank := seq_len(.N)]
stats_vec <- setNames(rk$stat, rk$gene_id)
ties_fgsea <- sum(duplicated(stats_vec[stats_vec != 0]))   # same definition as fgsea::preparePathwaysAndStats
tie_groups <- rk[, .N, by = stat][N > 1]
add("ranking", "genes_ranked", nrow(rk))
add("ranking", "ties_fgsea_definition", ties_fgsea)
add("ranking", "tie_groups", nrow(tie_groups))
add("ranking", "stat_exactly_zero", sum(rk$stat == 0))
add("ranking", "stat_min", min(rk$stat))
add("ranking", "stat_max", max(rk$stat))
add("ranking", "stat_quantiles_0_25_50_75_100", q5(rk$stat))
add("ranking", "stat_mean_sd", c(mean(rk$stat), sd(rk$stat)))
add("ranking", "positive_stat_genes", sum(rk$stat > 0))
add("ranking", "negative_stat_genes", sum(rk$stat < 0))
for (k in 3:5) {
  add("ranking", paste0("abs_stat_gt_", k, "_observed"), sum(abs(rk$stat) > k))
  add("ranking", paste0("abs_stat_gt_", k, "_expected_N01"), round(nrow(rk) * 2 * pnorm(-k), 2))
}
# position the excluded genes would take (diagnostic only)
all_sorted <- res[stat_finite == TRUE][order(-stat, gene_id)]
excl[, hypothetical_rank_if_included := match(gene_id, all_sorted$gene_id)]

wt(rk[, .(rank, gene_id, stat, log2FoldChange_MLE, lfcSE_MLE, pvalue, padj, baseMean, description)],
   file.path(out_rank, "deseq2_wald_ranking.tsv"))
wt(excl[order(-stat), .(gene_id, stat, log2FoldChange_MLE, lfcSE_MLE, pvalue, padj, padj_NA_reason, maxCooks, baseMean,
                       exclusion_reason, hypothetical_rank_if_included, description)],
   file.path(out_rank, "ranking_excluded_genes.tsv"))
wt(rk[order(-abs(stat))][1:20, .(rank, gene_id, stat, log2FoldChange_MLE, baseMean, padj, description)],
   file.path(out_rank, "ranking_extreme_statistics_top20.tsv"))
if (nrow(tie_groups)) wt(rk[stat %in% tie_groups$stat, .(rank, gene_id, stat, baseMean)], file.path(out_rank, "ranking_ties.tsv"))

png(file.path(fig_dir, "ranking_wald_statistic.png"), width = 2400, height = 1100, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
h <- hist(rk$stat, breaks = 100, freq = FALSE, col = "grey80", border = "grey60",
          main = sprintf("DESeq2 Wald statistic (n = %s)", format(nrow(rk), big.mark = ",")),
          xlab = "Wald statistic (STV_infected vs STV_free)")
curve(dnorm(x), add = TRUE, lwd = 2, col = "#D95F02")
legend("topright", legend = c("observed", "N(0,1) reference"), fill = c("grey80", NA), border = c("grey60", NA),
       lty = c(NA, 1), lwd = c(NA, 2), col = c(NA, "#D95F02"), bty = "n", cex = 0.8)
plot(rk$rank, rk$stat, type = "l", lwd = 1.5, xlab = "Rank (1 = most positive)", ylab = "Wald statistic",
     main = "Ranked list used for preranked GSEA")
abline(h = 0, lty = 2, col = "grey50")
mtext("positive = higher in STV_infected", side = 3, line = 0.2, adj = 0, cex = 0.7)
dev.off()

# =================================================================================================================
# TASK 2 - annotation resources
# =================================================================================================================
t2g <- fread(t2g_tsv, sep = "\t")
check(identical(names(t2g), c("TXNAME", "GENEID")), "tx2gene header")
check(!anyDuplicated(t2g$TXNAME) && !anyDuplicated(t2g$GENEID), "tx2gene is not 1:1")
itag_genes <- t2g$GENEID
check(all(res$gene_id %in% itag_genes), "DESeq2 genes outside ITAG4.0")
add("annotation", "itag40_genes", length(itag_genes))
add("annotation", "itag40_tx2gene_one_to_one", TRUE)

# --- SGN ITAG4.0_goterms.txt (previously documented inconsistency) ---------------------------------------------
g40 <- read_id_go(go40_txt)
g40[, has_go := grepl("GO:", go)]
g40[, exact_itag40_gene := id %in% itag_genes]
g40[, same_locus_other_version := !exact_itag40_gene & locus_of(id) %in% locus_of(itag_genes)]
add("sgn_itag40_goterms", "rows", nrow(g40))
add("sgn_itag40_goterms", "duplicated_ids", sum(duplicated(g40$id)))
add("sgn_itag40_goterms", "ids_with_GO", sum(g40$has_go))
add("sgn_itag40_goterms", "ids_exact_itag40_gene", sum(g40$exact_itag40_gene))
add("sgn_itag40_goterms", "ids_with_GO_exact_itag40_gene", sum(g40$has_go & g40$exact_itag40_gene))
add("sgn_itag40_goterms", "ids_not_exact_but_same_locus_other_version", sum(g40$same_locus_other_version))
add("sgn_itag40_goterms", "ids_with_GO_not_exact_but_same_locus", sum(g40$has_go & g40$same_locus_other_version))
g40_exact_go <- g40[has_go & exact_itag40_gene, id]
g40_locus_go <- unique(locus_of(g40[has_go & same_locus_other_version, id]))

# --- SGN ITAG4.1_goterms.txt (transcript IDs) with protein identity ----------------------------------------------
g41 <- read_id_go(go41_txt)
g41[, has_go := grepl("GO:", go)]
p40 <- read_fasta(prot40_fa)
p41 <- read_fasta(prot41_fa)
check(!anyDuplicated(names(p40)) && setequal(names(p40), t2g$TXNAME), "ITAG4.0 protein IDs != tx2gene transcripts")
g41[, exact_itag40_tx := id %in% t2g$TXNAME]
g41[, identical_protein := exact_itag40_tx & id %in% names(p41) & p40[id] == p41[id]]
g41[is.na(identical_protein), identical_protein := FALSE]
g41[, gene_id := t2g$GENEID[match(id, t2g$TXNAME)]]
add("sgn_itag41_goterms", "rows", nrow(g41))
add("sgn_itag41_goterms", "tx_with_GO", sum(g41$has_go))
add("sgn_itag41_goterms", "tx_with_GO_exact_itag40_tx", sum(g41$has_go & g41$exact_itag40_tx))
add("sgn_itag41_goterms", "tx_with_GO_exact_tx_identical_protein", sum(g41$has_go & g41$identical_protein))
g41_genes <- unique(g41[has_go & identical_protein, gene_id])

# --- PLAZA 5.0 dicots (primary GO source, D-024) ----------------------------------------------------------------
plaza_hdr <- readLines(gzfile(plaza_gz), n = 12)
plaza_meta <- sub("^# ", "", plaza_hdr[startsWith(plaza_hdr, "# ")])
check(any(grepl("annotation source/version : ITAG4.0", plaza_meta, fixed = TRUE)), "PLAZA file is not ITAG4.0-based")
for (l in plaza_meta) add("plaza", "file_header", l)
plaza <- fread(cmd = sprintf("zcat '%s' | grep -v '^# '", plaza_gz), sep = "\t", quote = "", header = TRUE)
setnames(plaza, 1, "gene_id")
check(identical(names(plaza), c("gene_id", "species", "go", "evidence", "go_source", "provider", "comment", "description",
                                "is_propagated_from_child_term")), "unexpected PLAZA columns")
check(all(plaza$species == "sly"), "PLAZA species column not 'sly'")
plaza_genes_any <- unique(plaza$gene_id)
check(all(plaza_genes_any %in% itag_genes), "PLAZA gene IDs outside ITAG4.0 (version mismatch)")
add("plaza", "rows", nrow(plaza))
add("plaza", "genes_any_source", length(plaza_genes_any))
add("plaza", "genes_exact_itag40", sum(plaza_genes_any %in% itag_genes))
add("plaza", "rows_with_NOT_qualifier_text", sum(grepl("NOT", plaza$comment, fixed = TRUE)))
src_tab <- plaza[, .N, by = .(go_source, provider, is_propagated_from_child_term)][order(-N)]
wt(src_tab, file.path(out_ann, "plaza_rows_by_source.tsv"))
prim <- unique(plaza[go_source == "primary" & is_propagated_from_child_term == 0L, .(gene_id, go, evidence, provider)])
ev_tab <- prim[, .N, by = .(provider, evidence)][order(-N)]
wt(ev_tab, file.path(out_ann, "plaza_primary_direct_evidence.tsv"))
prim_pairs <- unique(prim[, .(gene_id, go)])
add("plaza", "primary_direct_rows_dedup_gene_go_evidence_provider", nrow(prim))
add("plaza", "primary_direct_gene_go_pairs", nrow(prim_pairs))
add("plaza", "primary_direct_genes", uniqueN(prim_pairs$gene_id))
rm(plaza); invisible(gc())

# --- Gene Ontology (go-basic, dated release) -------------------------------------------------------------------
obo_head <- readLines(go_obo, n = 3)
add("go_ontology", "file", go_obo)
add("go_ontology", "data_version", sub("^data-version: ", "", grep("^data-version:", obo_head, value = TRUE)))
ont <- get_ontology(go_obo, propagate_relationships = c("is_a", "part_of"), extract_tags = "everything")
first_val <- function(x) if (is.list(x)) vapply(x, function(v) if (length(v)) v[[1]] else NA_character_, "") else as.character(x)
as_lst <- function(x) if (is.list(x)) x else as.list(x)
term_ns <- setNames(unname(ns_code[first_val(ont$namespace)]), ont$id)
obsolete <- setNames(as.logical(ont$obsolete), ont$id)
is_go <- startsWith(ont$id, "GO:")   # get_ontology also returns the [Typedef] stanzas (part_of, regulates, ...)
add("go_ontology", "non_GO_entries_ignored", ont$id[!is_go])
check(!anyNA(term_ns[is_go & !obsolete]), "GO terms with unknown namespace")
alt <- as_lst(ont$alt_id); names(alt) <- ont$id
alt <- alt[lengths(alt) > 0]
alt <- alt[!vapply(alt, function(v) all(is.na(v)), TRUE)]
alt2primary <- setNames(rep(names(alt), lengths(alt)), unlist(alt, use.names = FALSE))
repl <- as_lst(ont$replaced_by); names(repl) <- ont$id
add("go_ontology", "terms", sum(is_go))
add("go_ontology", "terms_obsolete", sum(obsolete[is_go]))
add("go_ontology", "alt_ids", length(alt2primary))
add("go_ontology", "relations_propagated", "is_a,part_of")

map_term <- function(g) {
  if (!startsWith(g, "GO:")) return(c(NA_character_, "not_in_ontology_dropped"))
  if (g %in% ont$id && !obsolete[[g]]) return(c(g, "current"))
  if (g %in% names(alt2primary)) {
    p <- alt2primary[[g]]
    if (!obsolete[[p]]) return(c(p, "alt_id_to_primary"))
  }
  if (g %in% ont$id && obsolete[[g]]) {
    r <- repl[[g]]; r <- r[!is.na(r)]
    if (length(r) == 1 && r %in% ont$id && !obsolete[[r]]) return(c(r, "obsolete_replaced_by"))
    return(c(NA_character_, "obsolete_dropped"))
  }
  c(NA_character_, "not_in_ontology_dropped")
}
used <- unique(prim_pairs$go)
mt <- t(vapply(used, map_term, character(2)))
term_map <- data.table(plaza_go = used, go_mapped = mt[, 1], mapping = mt[, 2])
term_map[, n_gene_pairs := prim_pairs[, .N, by = go][match(term_map$plaza_go, go), N]]
wt(term_map[mapping != "current"][order(mapping, plaza_go)], file.path(out_ann, "go_term_version_mapping.tsv"))
for (s in c("current", "alt_id_to_primary", "obsolete_replaced_by", "obsolete_dropped", "not_in_ontology_dropped")) {
  add("go_processing", paste0("terms_", s), sum(term_map$mapping == s))
  add("go_processing", paste0("gene_term_pairs_", s), sum(term_map[mapping == s, n_gene_pairs]))
}
direct <- unique(prim_pairs[term_map, on = c(go = "plaza_go"), nomatch = 0L][!is.na(go_mapped), .(gene_id, term = go_mapped)])
add("go_processing", "direct_gene_term_pairs_after_mapping", nrow(direct))

# propagation through is_a + part_of (go-basic: no cross-namespace links; verified below)
dterms <- unique(direct$term)
anc <- ont$ancestors[dterms]
check(all(mapply(function(t, a) t %in% a, dterms, anc)), "ontologyIndex ancestors do not include the term itself")
anc_long <- data.table(term = rep(dterms, lengths(anc)), ancestor = unlist(anc, use.names = FALSE))
anc_long[, ns_term := term_ns[term]][, ns_anc := term_ns[ancestor]]
check(all(anc_long$ns_term == anc_long$ns_anc), "GO propagation crosses namespaces")
prop <- unique(direct[anc_long, on = "term", allow.cartesian = TRUE, nomatch = 0L][, .(gene_id, term = ancestor, ns = ns_anc)])
add("go_processing", "propagated_gene_term_pairs_incl_roots", nrow(prop))
prop <- prop[!term %in% go_roots]
add("go_processing", "propagated_gene_term_pairs_excl_roots", nrow(prop))

# --- per-gene mapping table for the tested universe -------------------------------------------------------------
dir_str <- direct[, .(go_terms_direct = paste(sort(term), collapse = ";"), n_go_direct = .N), by = gene_id]
nsc <- dcast(prop[, .N, by = .(gene_id, ns)], gene_id ~ ns, value.var = "N", fill = 0L)
for (k in c("BP", "MF", "CC")) if (!k %in% names(nsc)) nsc[, (k) := 0L]
ranked_set <- rk$gene_id
mp <- res[, .(project_gene_id = gene_id, description)]
mp[, solyc_gene_id := fifelse(project_gene_id %in% plaza_genes_any, project_gene_id, NA_character_)]
mp <- merge(mp, dir_str, by.x = "project_gene_id", by.y = "gene_id", all.x = TRUE, sort = FALSE)
mp <- merge(mp, nsc[, .(gene_id, n_go_BP_propagated = BP, n_go_MF_propagated = MF, n_go_CC_propagated = CC)],
            by.x = "project_gene_id", by.y = "gene_id", all.x = TRUE, sort = FALSE)
for (k in c("n_go_direct", "n_go_BP_propagated", "n_go_MF_propagated", "n_go_CC_propagated")) mp[is.na(get(k)), (k) := 0L]
mp[, mapping_status := fifelse(n_go_direct > 0, "exact_ITAG4.0_id_primary_GO",
                       fifelse(project_gene_id %in% prim_pairs$gene_id, "exact_ITAG4.0_id_primary_GO_all_terms_dropped",
                       fifelse(!is.na(solyc_gene_id), "exact_ITAG4.0_id_non_primary_GO_only_not_used", "no_GO_in_PLAZA")))]
mp[, annotation_source := fifelse(n_go_direct > 0, "PLAZA5.0_dicots_primary_direct", NA_character_)]
mp[, annotation_version := fifelse(n_go_direct > 0, "plaza_dicots_05_2021-09-20|GO_release_2026-08-05", NA_character_)]
mp[, in_ranked_universe := project_gene_id %in% ranked_set]
mp[, sgn_itag40_goterms_status := fifelse(project_gene_id %in% g40_exact_go, "exact_id_with_GO",
                                  fifelse(locus_of(project_gene_id) %in% g40_locus_go, "GO_only_under_other_version_not_used",
                                          "no_GO"))]
mp[, sgn_itag41_goterms_status := fifelse(project_gene_id %in% g41_genes, "exact_tx_identical_protein_with_GO", "no_GO")]
setcolorder(mp, c("solyc_gene_id", "project_gene_id", "description", "go_terms_direct", "n_go_direct",
                  "n_go_BP_propagated", "n_go_MF_propagated", "n_go_CC_propagated", "annotation_source",
                  "annotation_version", "mapping_status", "in_ranked_universe", "sgn_itag40_goterms_status",
                  "sgn_itag41_goterms_status"))
setorder(mp, project_gene_id)
wt(mp, file.path(out_ann, "gene_annotation_mapping_tested_genes.tsv"))
for (s in sort(unique(mp$mapping_status))) {
  add("mapping_status_tested", s, sum(mp$mapping_status == s))
  add("mapping_status_ranked", s, sum(mp$mapping_status == s & mp$in_ranked_universe))
}

# --- coverage table ----------------------------------------------------------------------------------------------
tested_set <- res$gene_id
cov_row <- function(source, genes) {
  data.table(source = source,
             tested_genes_covered = sum(tested_set %in% genes), tested_universe = length(tested_set),
             pct_tested = round(100 * mean(tested_set %in% genes), 2),
             ranked_genes_covered = sum(ranked_set %in% genes), ranked_universe = length(ranked_set),
             pct_ranked = round(100 * mean(ranked_set %in% genes), 2))
}
cov <- rbind(
  cov_row("PLAZA5.0 primary direct GO after GO-version mapping (USED)", direct$gene_id),
  cov_row("PLAZA5.0 BP after propagation (non-root)", prop[ns == "BP", gene_id]),
  cov_row("PLAZA5.0 MF after propagation (non-root)", prop[ns == "MF", gene_id]),
  cov_row("PLAZA5.0 CC after propagation (non-root)", prop[ns == "CC", gene_id]),
  cov_row("PLAZA5.0 any source incl. orthology projections (not used)", plaza_genes_any),
  cov_row("SGN ITAG4.0_goterms exact gene ID (not used)", g40_exact_go),
  cov_row("SGN ITAG4.1_goterms exact transcript + identical protein (not used)", g41_genes))
wt(cov, file.path(out_ann, "go_annotation_coverage.tsv"))

# =================================================================================================================
# TASK 3 - GO gene sets restricted to the ranked universe
# =================================================================================================================
prop_r <- prop[gene_id %in% ranked_set]
sets <- lapply(c(BP = "BP", MF = "MF", CC = "CC"), function(k) {
  x <- prop_r[ns == k]
  split(x$gene_id, x$term)
})
catalog <- rbindlist(lapply(names(sets), function(k) {
  s <- sets[[k]]
  data.table(term = names(s), namespace = k, term_name = first_val(ont$name[names(s)]),
             size_ranked_universe = lengths(s))
}))
catalog[, in_size_range_15_500 := size_ranked_universe >= min_size & size_ranked_universe <= max_size]
memb_key <- unlist(lapply(sets, function(s) vapply(s, function(g) paste(sort(g), collapse = ","), "")), use.names = FALSE)
catalog[, membership_key := memb_key]
catalog[, identical_membership_group_size := .N, by = .(namespace, membership_key)]
catalog[, membership_key := NULL]
setorder(catalog, namespace, -size_ranked_universe, term)
wt(catalog, file.path(out_ann, "go_gene_set_catalog.tsv"))
for (k in names(sets)) {
  ck <- catalog[namespace == k]
  tested_k <- ck[in_size_range_15_500 == TRUE]
  genes_k <- unique(prop_r[ns == k, gene_id])
  genes_tested_k <- unique(unlist(sets[[k]][tested_k$term]))
  add(paste0("go_sets_", k), "terms_with_ranked_genes", nrow(ck))
  add(paste0("go_sets_", k), "terms_in_size_range_15_500", nrow(tested_k))
  add(paste0("go_sets_", k), "terms_below_15", sum(ck$size_ranked_universe < min_size))
  add(paste0("go_sets_", k), "terms_above_500", sum(ck$size_ranked_universe > max_size))
  add(paste0("go_sets_", k), "ranked_genes_annotated_nonroot", length(genes_k))
  add(paste0("go_sets_", k), "pct_ranked_genes_annotated_nonroot", round(100 * length(genes_k) / length(ranked_set), 2))
  add(paste0("go_sets_", k), "ranked_genes_in_tested_sets", length(genes_tested_k))
  add(paste0("go_sets_", k), "pct_ranked_genes_in_tested_sets", round(100 * length(genes_tested_k) / length(ranked_set), 2))
  add(paste0("go_sets_", k), "size_quantiles_all_terms", q5(ck$size_ranked_universe))
  add(paste0("go_sets_", k), "size_quantiles_tested_terms", q5(tested_k$size_ranked_universe))
  add(paste0("go_sets_", k), "tested_terms_sharing_identical_membership", sum(tested_k$identical_membership_group_size > 1))
}
saveRDS(list(sets = sets, ranked_genes = ranked_set, go_data_version = summ[metric == "data_version", value],
             min_size = min_size, max_size = max_size), file.path(out_ann, "go_gene_sets.rds"))

png(file.path(fig_dir, "go_gene_set_sizes.png"), width = 2700, height = 900, res = 200)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
for (k in names(sets)) {
  sz <- catalog[namespace == k, size_ranked_universe]
  hist(log10(sz), breaks = 40, col = "grey80", border = "grey60",
       main = sprintf("GO %s: %d terms (%d in 15-500)", k, length(sz), sum(sz >= min_size & sz <= max_size)),
       xlab = "log10(genes in ranked universe)")
  abline(v = log10(c(min_size, max_size)), lty = 2, col = "#D95F02", lwd = 2)
}
dev.off()

# =================================================================================================================
# TASK 4 - Plant Reactome tomato audit (no GSEA unless validated; D-028)
# =================================================================================================================
pr <- fread(pr_tab, sep = "\t", header = FALSE, quote = "")
check(ncol(pr) == 4, "unexpected Plant Reactome file layout")
setnames(pr, c("pathway_id", "pathway_name", "species", "member_id"))
sly <- unique(pr[species == "Solanum lycopersicum"])
ids <- data.table(member_id = unique(sly$member_id))
ids[, id_class := fifelse(grepl("^Solyc[0-9]{2}g[0-9]{6}[.][0-9]+$", member_id), "Solyc_gene_id", "other_identifier")]
ids[, exact_itag40_gene := member_id %in% itag_genes]
ids[, in_ranked_universe := member_id %in% ranked_set]
ids[, itag40_same_locus_other_version := !exact_itag40_gene & id_class == "Solyc_gene_id" &
                                         locus_of(member_id) %in% locus_of(itag_genes)]
setorder(ids, id_class, member_id)
wt(ids, file.path(out_pr, "plant_reactome_tomato_identifier_audit.tsv"))
pw <- sly[, .(n_members = uniqueN(member_id), n_exact_itag40 = sum(member_id %in% itag_genes),
              n_ranked = sum(member_id %in% ranked_set)), by = .(pathway_id, pathway_name)]
pw[, in_size_range_15_500 := n_ranked >= min_size & n_ranked <= max_size]
setorder(pw, -n_ranked, pathway_id)
wt(pw, file.path(out_pr, "plant_reactome_tomato_pathway_coverage.tsv"))
wt(sly[member_id %in% ranked_set, .(pathway_id, pathway_name, gene_id = member_id)][order(pathway_id, gene_id)],
   file.path(out_pr, "plant_reactome_tomato_exact_mapping_ranked.tsv"))
add("plant_reactome", "file", pr_tab)
add("plant_reactome", "species_in_file", uniqueN(pr$species))
add("plant_reactome", "tomato_rows", nrow(sly))
add("plant_reactome", "tomato_pathways", uniqueN(sly$pathway_id))
add("plant_reactome", "tomato_member_ids", nrow(ids))
add("plant_reactome", "tomato_ids_Solyc_format", sum(ids$id_class == "Solyc_gene_id"))
add("plant_reactome", "tomato_ids_other_format", sum(ids$id_class != "Solyc_gene_id"))
add("plant_reactome", "tomato_ids_exact_itag40", sum(ids$exact_itag40_gene))
add("plant_reactome", "tomato_ids_same_locus_other_version_not_used", sum(ids$itag40_same_locus_other_version))
add("plant_reactome", "ranked_genes_covered", sum(ids$in_ranked_universe))
add("plant_reactome", "pct_ranked_genes_covered", round(100 * sum(ids$in_ranked_universe) / length(ranked_set), 2))
add("plant_reactome", "pathway_size_quantiles_members", q5(pw$n_members))
add("plant_reactome", "pathway_size_quantiles_ranked", q5(pw$n_ranked))
add("plant_reactome", "pathways_with_15_500_ranked_genes", sum(pw$in_size_range_15_500))
add("plant_reactome", "curation_status_tomato", "computationally projected from curated rice pathways (Plant Reactome docs)")

wt(summ, "results/functional/functional_annotation_summary.tsv")
writeLines(capture.output(sessionInfo()), file.path(out_ann, "sessionInfo_05_functional_ranking_annotation.txt"))
cat("05_functional_ranking_annotation.R finished\n")
print(summ, nrows = 500)
