#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 6, Task 7: MapMan / Mercator4 preparation (no over-representation analysis; D-029).
# Mercator4 (v8.0 on 2026-10-26) is a web service (https://www.plabipd.de/mercator_main.html); the submission is a
# documented manual step (results/functional/mapman/MERCATOR4_SUBMISSION.md). This script prepares and checks inputs.
# Usage (project root, env stv-r): Rscript scripts/R/07_mapman_mercator_prep.R
# Inputs : reference/SGN/ITAG4.0/ITAG4.0_proteins.fasta (unmodified Mercator4 input), reference/ITAG4.0_tx2gene.tsv,
#          results/deseq2/deseq2_results_shrunken.tsv, results/deseq2/deseq2_results_all.tsv,
#          results/functional/ranking/deseq2_wald_ranking.tsv
# Outputs: results/functional/mapman/{mercator4_input_manifest.tsv, mercator4_input_fasta_record.tsv,
#          mapman_expression_STV_infected_vs_free.txt, mapman_prep_summary.tsv}
suppressPackageStartupMessages(library(data.table))
options(warn = 1)

fasta   <- "reference/SGN/ITAG4.0/ITAG4.0_proteins.fasta"
t2g_tsv <- "reference/ITAG4.0_tx2gene.tsv"
shr_tsv <- "results/deseq2/deseq2_results_shrunken.tsv"
all_tsv <- "results/deseq2/deseq2_results_all.tsv"
rk_tsv  <- "results/functional/ranking/deseq2_wald_ranking.tsv"
out_dir <- "results/functional/mapman"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
check <- function(ok, msg) if (!isTRUE(ok)) stop("STOP CONDITION: ", msg, call. = FALSE)
wt <- function(x, file) fwrite(x, file, sep = "\t", quote = FALSE, na = "NA")
summ <- data.table(metric = character(), value = character())
add <- function(m, v) summ <<- rbind(summ, data.table(metric = m, value = paste(v, collapse = ",")))

# --- FASTA (read only; never modified) ---------------------------------------------------------------------------
x <- readLines(fasta)
h <- startsWith(x, ">")
ids <- sub("^>(\\S+).*$", "\\1", x[h])
seqs <- vapply(split(x[!h], cumsum(h)[!h]), paste, "", collapse = "")
check(length(seqs) == length(ids), "FASTA records without sequence")
check(!anyDuplicated(ids), "duplicated FASTA identifiers")
t2g <- fread(t2g_tsv)
check(setequal(ids, t2g$TXNAME), "FASTA identifiers differ from ITAG4.0 tx2gene transcripts")
term_stop <- endsWith(seqs, "*")
core <- sub("[*]$", "", seqs)
internal <- lengths(regmatches(core, gregexpr("[*]", core)))
bad_chars <- grepl("[^ACDEFGHIKLMNPQRSTVWYXBZUO*]", seqs)
rk <- fread(rk_tsv)
tested <- fread(all_tsv, select = "gene_id")$gene_id
man <- data.table(protein_id = ids, gene_id = t2g$GENEID[match(ids, t2g$TXNAME)], aa_length = nchar(core),
                  terminal_stop = term_stop, internal_stop_count = internal, nonstandard_characters = bad_chars)
man[, tested_in_deseq2 := gene_id %in% tested][, in_ranked_universe := gene_id %in% rk$gene_id]
setorder(man, protein_id)
wt(man, file.path(out_dir, "mercator4_input_manifest.tsv"))
sha <- sub(" .*$", "", system2("sha256sum", fasta, stdout = TRUE))
wt(data.table(file = fasta, sha256 = sha, bytes = file.size(fasta), records = length(ids),
              terminal_stop_records = sum(term_stop), internal_stop_records = sum(internal > 0),
              nonstandard_character_records = sum(bad_chars), sequence_type = "protein",
              source = "SGN ITAG4.0 release (reference/reference_manifest.tsv)"),
   file.path(out_dir, "mercator4_input_fasta_record.tsv"))
add("fasta_file", fasta); add("fasta_sha256", sha); add("fasta_records", length(ids))
add("fasta_records_terminal_stop", sum(term_stop)); add("fasta_records_internal_stop", sum(internal > 0))
add("fasta_records_nonstandard_characters", sum(bad_chars))
add("fasta_ids_equal_tx2gene_transcripts", TRUE)

# --- MapMan expression table (descriptive; D-017 apeglm LFC as display value) ---------------------------------------
shr <- fread(shr_tsv, select = c("gene_id", "log2FoldChange_apeglm"))
al <- fread(all_tsv, select = c("gene_id", "log2FoldChange_MLE", "stat"))
check(nrow(shr) == 18277L && nrow(al) == 18277L && setequal(shr$gene_id, al$gene_id), "DESeq2 tables differ in genes")
ex <- merge(shr, al, by = "gene_id")
check(!anyNA(ex), "NA values in the MapMan expression metrics")
ex[, IDENTIFIER := t2g$TXNAME[match(gene_id, t2g$GENEID)]]
check(!anyNA(ex$IDENTIFIER) && !anyDuplicated(ex$IDENTIFIER), "gene -> protein identifier mapping not 1:1")
ex <- ex[order(IDENTIFIER), .(IDENTIFIER, log2FC_apeglm_infected_vs_free = log2FoldChange_apeglm,
                              log2FC_MLE_infected_vs_free = log2FoldChange_MLE, wald_stat_infected_vs_free = stat)]
fwrite(ex, file.path(out_dir, "mapman_expression_STV_infected_vs_free.txt"), sep = "\t", quote = FALSE)
add("expression_rows", nrow(ex))
add("expression_identifier", "ITAG4.0 protein/transcript ID exactly as in the FASTA (case preserved)")
add("expression_display_column", "log2FC_apeglm_infected_vs_free (positive = higher in STV_infected)")
add("expression_range_apeglm", signif(range(ex$log2FC_apeglm_infected_vs_free), 4))
add("genes_not_tested_omitted", length(ids) - nrow(ex))
wt(summ, file.path(out_dir, "mapman_prep_summary.tsv"))
print(summ)
cat("07_mapman_mercator_prep.R finished\n")
