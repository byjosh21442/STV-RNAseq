# STV_GSE137303 - Phase 6B shared helpers (sourced by scripts/R/09 - 16). No analysis is run when this file is sourced.
# fgsea parameters are copied verbatim from scripts/R/06_functional_gsea.R (D-027); GO processing mirrors
# scripts/R/05_functional_ranking_annotation.R (D-025) and is checked against the Phase 6 gene sets before use.
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(BiocParallel)
})
options(warn = 1)

fr_seed        <- 20261025L
fr_seed_diag   <- 20261026L
fr_fdr         <- 0.05
fr_support_fdr <- 0.25
fr_fg_args <- list(minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000, sampleSize = 101,
                   scoreType = "std", gseaParam = 1, BPPARAM = SerialParam())
fr_collections <- c("BP", "MF", "CC")
fr_out <- "results/functional_robustness"
fr_fig <- "figures/functional_robustness"
fr_samples <- c("STV_free_rep1", "STV_free_rep2", "STV_free_rep3",
                "STV_infected_rep1", "STV_infected_rep2", "STV_infected_rep3")
fr_coef <- "condition_STV_infected_vs_STV_free"

fr_check <- function(ok, msg) if (!isTRUE(ok)) stop("STOP CONDITION: ", msg, call. = FALSE)
fr_wt <- function(x, file) fwrite(x, file, sep = "\t", quote = FALSE, na = "NA")
fr_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

# The Phase 6B pre-specification must be frozen and unchanged before any 6B analysis runs.
fr_check_frozen <- function() {
  ck <- fread(file.path(fr_out, "prespecification/prespecification_checksum.tsv"), colClasses = "character")
  for (i in seq_len(nrow(ck))) {
    s <- sub(" .*$", "", system2("sha256sum", ck$file[i], stdout = TRUE))
    fr_check(identical(s, ck$sha256[i]), paste("frozen pre-specification changed:", ck$file[i]))
  }
  now <- fr_now()
  fr_check(now >= max(ck$frozen_utc), "run started before the Phase 6B pre-specification was frozen")
  now
}

# Sample metadata exactly as in scripts/R/02_deseq2.R (labels and order asserted).
fr_coldata <- function() {
  sheet <- read.csv("config/samplesheet.csv", stringsAsFactors = FALSE)
  fr_check(identical(sheet$sample, fr_samples), "samplesheet sample order/labels changed")
  fr_check(identical(sheet$condition, rep(c("STV_free", "STV_infected"), each = 3)), "samplesheet conditions changed")
  data.frame(row.names = sheet$sample, sample = sheet$sample,
             condition = factor(sheet$condition, levels = c("STV_free", "STV_infected")))
}

# Ranking rule of D-023/D-033 applied to any DESeq2 results table: finite statistic (and, unless include_na_p,
# a reported p-value); sorted by decreasing statistic, ties broken by gene ID (as in Phase 6).
fr_ranking <- function(gene_id, stat, pvalue, include_na_p = FALSE) {
  keep <- is.finite(stat) & (include_na_p | !is.na(pvalue))
  x <- data.table(gene_id = gene_id[keep], stat = stat[keep])
  setorder(x, -stat, gene_id)
  setNames(x$stat, x$gene_id)
}

# GO gene-term pairs (direct PLAZA primary terms after version mapping, from the tracked mapping table) propagated over
# is_a + part_of with the same go-basic release; roots excluded. Returns gene_id, term, ns.
fr_go_pairs <- function() {
  mp <- fread("results/functional/annotation/gene_annotation_mapping_tested_genes.tsv", sep = "\t", quote = "",
              colClasses = "character")
  direct <- mp[as.integer(n_go_direct) > 0,
               .(term = unlist(strsplit(go_terms_direct, ";", fixed = TRUE))), by = .(gene_id = project_gene_id)]
  ont <- ontologyIndex::get_ontology("reference/functional/go/release_2026-08-05/go-basic.obo",
                                     propagate_relationships = c("is_a", "part_of"), extract_tags = "everything")
  first_val <- function(x) if (is.list(x)) vapply(x, function(v) if (length(v)) v[[1]] else NA_character_, "") else as.character(x)
  ns_code <- c(biological_process = "BP", molecular_function = "MF", cellular_component = "CC")
  term_ns <- setNames(unname(ns_code[first_val(ont$namespace)]), ont$id)
  dterms <- unique(direct$term)
  fr_check(all(dterms %in% ont$id), "direct GO terms missing from the GO release")
  anc <- ont$ancestors[dterms]
  anc_long <- data.table(term = rep(dterms, lengths(anc)), ancestor = unlist(anc, use.names = FALSE))
  prop <- unique(direct[anc_long, on = "term", allow.cartesian = TRUE, nomatch = 0L][, .(gene_id, term = ancestor)])
  prop[, ns := term_ns[term]]
  prop <- prop[!term %in% c("GO:0008150", "GO:0003674", "GO:0005575")]
  attr(prop, "term_name") <- setNames(first_val(ont$name), ont$id)
  prop
}

fr_sets <- function(pairs, universe) {
  lapply(setNames(fr_collections, fr_collections), function(k) {
    x <- pairs[ns == k & gene_id %in% universe]
    lapply(split(x$gene_id, x$term), sort)
  })
}

fr_same_sets <- function(a, b) {
  identical(sort(names(a)), sort(names(b))) &&
    all(vapply(names(a), function(t) identical(sort(a[[t]]), sort(b[[t]])), TRUE))
}

# fgsea exactly as in Phase 6 (warnings captured, not muffled silently).
fr_fgsea <- function(pathways, stats, seed = fr_seed) {
  warns <- character()
  set.seed(seed)
  res <- withCallingHandlers(do.call(fgseaMultilevel, c(list(pathways = pathways, stats = stats), fr_fg_args)),
                             warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  attr(res, "warnings") <- warns
  res
}
fr_material_warnings <- function(res) {
  w <- attr(res, "warnings")
  sum(!grepl("likely overestimated", w, fixed = TRUE))
}

# Phase 6 primary results (read-only).
fr_phase6 <- function(k) fread(sprintf("results/functional/go_gsea/gsea_go_%s_all.tsv", k), sep = "\t", quote = "")
fr_phase6_seed <- function(k) fread(sprintf("results/functional/go_gsea/gsea_go_%s_seed_diagnostic.tsv", k), sep = "\t", quote = "")

# Compare a sensitivity fgsea run with the Phase 6 primary result for one collection.
fr_compare_primary <- function(k, res, label) {
  p6 <- fr_phase6(k)
  s <- res[, .(go_id = pathway, NES_s = NES, padj_s = padj, size_s = size)]
  m <- merge(p6[, .(go_id, term_name, NES, padj, size, significant_fdr05)], s, by = "go_id", all = TRUE)
  m[, collection := k][, run := label]
  m[, status := fifelse(!significant_fdr05 %in% TRUE, NA_character_,
                fifelse(!is.na(padj_s) & padj_s < fr_fdr & sign(NES_s) == sign(NES), "retained",
                fifelse(!is.na(padj_s) & padj_s < fr_support_fdr & sign(NES_s) == sign(NES), "supported", "lost")))]
  m[, new_significant := !(significant_fdr05 %in% TRUE) & !is.na(padj_s) & padj_s < fr_fdr]
  m[]
}
