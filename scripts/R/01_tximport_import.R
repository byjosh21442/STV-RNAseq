#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 3: validate Salmon outputs, sample metadata and tx2gene, then import gene-level estimates
# with tximport as a smoke test. No differential expression, no filtering, no thresholds, no contrasts.
# Usage (project root, env stv-r): Rscript scripts/R/01_tximport_import.R
# Writes: results/quantification/tximport_validation.tsv, results/quantification/tximport_sessionInfo.txt
suppressPackageStartupMessages(library(tximport))

samples_csv <- "results/quantification/salmon/samples.csv"   # index written by the Nextflow workflow
sheet_csv   <- "config/samplesheet.csv"
tx2gene_tsv <- "reference/ITAG4.0_tx2gene.tsv"
quant_root  <- "results/quantification/salmon"
out_tsv     <- "results/quantification/tximport_validation.tsv"
out_session <- "results/quantification/tximport_sessionInfo.txt"

checks <- data.frame(check = character(), value = character(), status = character())
add <- function(check, value, ok = TRUE) {
  checks[nrow(checks) + 1, ] <<- list(check, paste(value, collapse = ","), if (isTRUE(ok)) "OK" else "FAIL")
}

# --- sample metadata ---------------------------------------------------------------------------------------
wf    <- read.csv(samples_csv, stringsAsFactors = FALSE)
sheet <- read.csv(sheet_csv, stringsAsFactors = FALSE)
add("samples_in_workflow_index", nrow(wf), nrow(wf) == 6)
add("sample_ids_match_samplesheet", setequal(wf$sample, sheet$sample), setequal(wf$sample, sheet$sample))
coldata <- sheet[, c("sample", "condition", "replicate", "geo_gsm", "run_accession")]
wf <- wf[match(coldata$sample, wf$sample), ]
add("condition_matches_workflow_index", identical(coldata$condition, wf$condition), identical(coldata$condition, wf$condition))
add("samples_per_condition", paste(names(table(coldata$condition)), table(coldata$condition), sep = "="))
add("libtype_requested", unique(wf$libtype), length(unique(wf$libtype)) == 1)

# --- quant.sf discovery ------------------------------------------------------------------------------------
files <- file.path(quant_root, coldata$sample, "quant.sf")
names(files) <- coldata$sample
add("quant_sf_found", sum(file.exists(files)), all(file.exists(files)))
same_dir <- normalizePath(dirname(files)) == normalizePath(wf$quant_dir)
add("quant_dir_matches_workflow_index", all(same_dir), all(same_dir))

# --- Salmon run settings recorded in meta_info.json --------------------------------------------------------
meta <- lapply(file.path(dirname(files), "aux_info", "meta_info.json"), jsonlite::fromJSON)
applied <- vapply(meta, function(m) paste(m$library_types, collapse = "+"), character(1))
add("salmon_library_type_applied", unique(applied),
    length(unique(applied)) == 1 && all(applied == wf$libtype | wf$libtype == "A"))
gc_bias <- vapply(meta, function(m) isTRUE(m$gc_bias_correct), logical(1))
add("salmon_gc_bias_correct", unique(gc_bias), length(unique(gc_bias)) == 1)
seq_bias <- vapply(meta, function(m) isTRUE(m$seq_bias_correct), logical(1))
add("salmon_seq_bias_correct", unique(seq_bias), length(unique(seq_bias)) == 1)
add("salmon_version", unique(vapply(meta, `[[`, character(1), "salmon_version")),
    length(unique(vapply(meta, `[[`, character(1), "salmon_version"))) == 1)

# --- tx2gene -----------------------------------------------------------------------------------------------
tx2gene <- read.delim(tx2gene_tsv, stringsAsFactors = FALSE)
add("tx2gene_columns", names(tx2gene), identical(names(tx2gene), c("TXNAME", "GENEID")))
add("tx2gene_rows", nrow(tx2gene))
add("tx2gene_duplicated_TXNAME", sum(duplicated(tx2gene$TXNAME)), !any(duplicated(tx2gene$TXNAME)))
quant_names <- lapply(files, function(f) read.delim(f, stringsAsFactors = FALSE)$Name)
same_names <- all(vapply(quant_names, identical, logical(1), quant_names[[1]]))
add("quant_sf_same_transcripts_all_samples", same_names, same_names)
add("quant_sf_transcripts", length(quant_names[[1]]))
missing_tx <- setdiff(quant_names[[1]], tx2gene$TXNAME)
add("quant_sf_names_missing_from_tx2gene", length(missing_tx), length(missing_tx) == 0)
add("tx2gene_TXNAME_absent_from_quant_sf", length(setdiff(tx2gene$TXNAME, quant_names[[1]])), TRUE)
expected_genes <- unique(tx2gene$GENEID[tx2gene$TXNAME %in% quant_names[[1]]])

# --- tximport (gene level) ---------------------------------------------------------------------------------
# countsFromAbundance = "no" (default) keeps estimated counts plus the average-transcript-length matrix that
# DESeqDataSetFromTximport() turns into a normalization offset (tximport and DESeq2 vignettes).
# ignoreTxVersion = FALSE (default): identifiers already match exactly.
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
add("tximport_countsFromAbundance", txi$countsFromAbundance)
add("gene_matrix_dim", dim(txi$counts))
add("genes_equal_expected", nrow(txi$counts) == length(expected_genes), nrow(txi$counts) == length(expected_genes))
add("gene_ids_equal_expected", setequal(rownames(txi$counts), expected_genes), setequal(rownames(txi$counts), expected_genes))
add("columns_equal_samples", identical(colnames(txi$counts), coldata$sample), identical(colnames(txi$counts), coldata$sample))
add("NA_counts", sum(is.na(txi$counts)), !anyNA(txi$counts))
add("NA_or_nonpositive_length", sum(is.na(txi$length) | txi$length <= 0), !anyNA(txi$length) && all(txi$length > 0))
salmon_sums <- vapply(files, function(f) sum(read.delim(f)$NumReads), numeric(1))
rel_diff <- max(abs(colSums(txi$counts) - salmon_sums) / salmon_sums)
add("max_relative_diff_gene_vs_salmon_total", signif(rel_diff, 3), rel_diff < 1e-9)
add("gene_count_totals", round(colSums(txi$counts)))
add("genes_zero_in_all_samples", sum(rowSums(txi$counts) == 0))

# --- report ------------------------------------------------------------------------------------------------
write.table(checks, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), out_session)
print(checks, right = FALSE, row.names = FALSE)
if (any(checks$status == "FAIL")) {
  stop("tximport validation failed: ", paste(checks$check[checks$status == "FAIL"], collapse = ", "))
}
cat("tximport validation OK ->", out_tsv, "\n")
