#!/usr/bin/env Rscript
# STV_GSE137303 - Phase 7 (Tasks 2, 3, 8 - 13): integrated evidence matrix and theme deep dives.
# Descriptive integration of FROZEN outputs only: no new test, no new threshold for significance, no re-run of DESeq2, fgsea
# or CAMERA. Classification follows D-045 as revised by the owner (docs/biological_integration.md section 2): evidence status
# from inferential evidence (DESeq2, CAMERA), UNSUPPORTED only by explicit review, and a separate descriptive robustness
# descriptor. The superseded U1 - U3 rules (versions 1 and 2) are kept as historical columns and are not applied.
# The only new computation on expression values is a descriptive per-gene replicate pattern for the ethylene case study
# (blind VST from deseq2_objects.rds, as in D-036; no p-values).
# Usage (project root, env stv-r): Rscript scripts/R/20_biological_integration.R
# Inputs : results/functional_robustness/{themes,comparison,replicate_consistency,loo,ch00,cooks_sensitivity,prior_studies}/,
#          results/organelle_assessment/{set_level_organelle_identical_leading_edge,gene_organelle_classification_corrected}.tsv,
#          results/functional/go_gsea/leading_edge_significant_sets_long.tsv, results/functional/annotation/go_gene_sets.rds,
#          results/deseq2/{deseq2_results_all,deseq2_results_shrunken}.tsv, results/deseq2/deseq2_objects.rds,
#          results/biological_integration/curated_theme_literature.tsv (manual literature summary; see docs)
# Outputs: results/biological_integration/*.tsv, sessionInfo_20.txt
suppressPackageStartupMessages({ library(data.table); library(DESeq2) })
options(warn = 1)
fr <- "results/functional_robustness"
out <- "results/biological_integration"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
fwt <- function(x, f) fwrite(x, file.path(out, f), sep = "\t", quote = FALSE, na = "NA")
rd <- function(f, ...) fread(f, sep = "\t", quote = "", ...)
short_desc <- function(x) trimws(sub(" \\(AHRD.*$", "", gsub("%2C", ",", x)))
themes10 <- c("Ribosome / translation / RNA processing", "Photosynthesis / plastid / thylakoid", "ATP synthase / proton transport",
              "Defence / biotic response", "Hormone signalling", "Endoplasmic reticulum / protein folding",
              "Heme / oxidoreductase", "Glutathione / detoxification", "Water deprivation / abiotic stress",
              "Carbohydrate binding / cell wall")
samples <- c("STV_free_rep1", "STV_free_rep2", "STV_free_rep3", "STV_infected_rep1", "STV_infected_rep2", "STV_infected_rep3")

# ---------------------------------------------------------------- set-level evidence (215 fgsea-significant sets)
tm  <- rd(file.path(fr, "themes/theme_membership_significant_sets.tsv"))
cmp <- rbindlist(lapply(c("BP", "MF", "CC"), function(k) rd(file.path(fr, "comparison", sprintf("method_comparison_%s.tsv", k)))))
rc  <- rd(file.path(fr, "replicate_consistency/replicate_consistency_all_significant_sets.tsv"))
loo <- rd(file.path(fr, "loo/loo_dependence_by_set.tsv"))
ph  <- rd(file.path(fr, "ch00/posthoc_identity98_vs_primary.tsv"))[significant_fdr05 == TRUE]
chp <- rd(file.path(fr, "ch00/ch00_sensitivity_vs_primary.tsv"))[significant_fdr05 == TRUE]
org <- rd("results/organelle_assessment/set_level_organelle_identical_leading_edge.tsv")
le  <- rd("results/functional/go_gsea/leading_edge_significant_sets_long.tsv")
res <- rd("results/deseq2/deseq2_results_all.tsv")
shr <- rd("results/deseq2/deseq2_results_shrunken.tsv", select = c("gene_id", "log2FoldChange_apeglm"))
rk  <- rd("results/functional/ranking/deseq2_wald_ranking.tsv", select = c("gene_id", "rank"))
gcl <- rd("results/organelle_assessment/gene_organelle_classification_corrected.tsv",
          select = c("gene_id", "corrected_class", "on_ch00"))
stopifnot(nrow(tm) == 215, uniqueN(tm[, .(collection, go_id)]) == 215, all(tm$sig_f))

s <- merge(tm[, .(collection, go_id, theme, term_name, tier, NES, padj_fgsea = padj, camera_direction = Direction,
                  camera_FDR = FDR, replicate_class = classification, loo_free_rep1_dependent = free1_dependent,
                  ch00_sensitive_prespecified = ch00_sensitive, cooks_status, direction)],
           cmp[, .(collection, go_id, size, leading_edge_size, camera_correlation = Correlation, camera_VIF = VIF,
                   concordant, FDR_camera_cor001_descriptive = FDR_cor001, fgsea_seed_unstable)],
           by = c("collection", "go_id"), all.x = TRUE)
s <- merge(s, rc[, .(collection, go_id, pairs_in_expected_direction, r_drop_STV_free_rep1)], by = c("collection", "go_id"), all.x = TRUE)
dep_cols <- paste0("dependent_on_", samples)
s <- merge(s, loo[, c("collection", "go_id", samples, dep_cols), with = FALSE], by = c("collection", "go_id"), all.x = TRUE)
setnames(s, samples, paste0("loo_status_drop_", samples))
s <- merge(s, ph[, .(collection, go_id = go_id, ch00_sensitive_posthoc = ch00_sensitive)], by = c("collection", "go_id"), all.x = TRUE)
s <- merge(s, org[, .(collection, go_id, le_genes_org = le_genes, le_org_identical = le_org_near_or_partial)],
           by = c("collection", "go_id"), all.x = TRUE)
stopifnot(nrow(s) == 215, !anyNA(s$tier), !anyNA(s$replicate_class), !anyNA(s[[dep_cols[1]]]), !anyNA(s$ch00_sensitive_posthoc),
          !any(s$cooks_status == "lost"))
s[, n_samples_dependent := rowSums(.SD), .SDcols = dep_cols]
setorder(s, theme, collection, padj_fgsea)
fwt(s, "theme_set_level_evidence.tsv")

# ---------------------------------------------------------------- representative leading-edge genes per theme
lt <- merge(le[, .(collection, go_id, gene_id, stat, padj_gene, description)], tm[, .(collection, go_id, theme)], by = c("collection", "go_id"))
lt <- merge(lt, gcl, by = "gene_id", all.x = TRUE)
lt[is.na(corrected_class), `:=`(corrected_class = "no organelle hit", on_ch00 = FALSE)]
gt <- lt[, .(n_theme_sets_in_LE = uniqueN(paste(collection, go_id)), stat = stat[1], padj_gene = padj_gene[1],
             description = short_desc(description[1]), organelle_identical = grepl("near-identical|partial", corrected_class[1]),
             on_ch00 = on_ch00[1]), by = .(theme, gene_id)]
gt[, abs_stat := abs(stat)]
setorder(gt, theme, -n_theme_sets_in_LE, -abs_stat)
gt[, abs_stat := NULL]
gt[, rank_in_theme := seq_len(.N), by = theme]
fwt(gt[rank_in_theme <= 10], "theme_representative_genes.tsv")

# ---------------------------------------------------------------- curated literature (manual input, documented)
lit <- rd(file.path(out, "curated_theme_literature.tsv"))
stopifnot(setequal(lit$theme, themes10))

# ---------------------------------------------------------------- theme matrix + Phase 7 classification (rules: docs section 2)
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
mat <- rbindlist(lapply(themes10, function(th) {
  x <- s[theme == th]; n <- nrow(x); g <- gt[theme == th]
  dep_share <- sapply(dep_cols, function(cc) mean(x[[cc]]))
  rep_top <- g[rank_in_theme <= 5, paste0(gene_id, " (", substr(description, 1, 40), "; ", n_theme_sets_in_LE, " sets; stat ",
                                          fmt(stat), "; padj ", fmt(padj_gene, 3), ")", collapse = " | ")]
  data.table(
    theme = th, n_sets = n, sets_BP = sum(x$collection == "BP"), sets_MF = sum(x$collection == "MF"), sets_CC = sum(x$collection == "CC"),
    top_go_terms = x[order(padj_fgsea)][1:min(3, n), paste0(go_id, " ", term_name, collapse = "; ")],
    all_go_ids = paste(x$go_id, collapse = ";"),
    fgsea_positive = sum(x$NES > 0), fgsea_negative = sum(x$NES < 0),
    fgsea_direction = fifelse(all(x$NES > 0), "positive (higher rank in STV-infected)", fifelse(all(x$NES < 0), "negative (lower rank in STV-infected)", "mixed")),
    NES_min = round(min(x$NES), 3), NES_max = round(max(x$NES), 3),
    fgsea_padj_min = signif(min(x$padj_fgsea), 3), fgsea_padj_median = signif(median(x$padj_fgsea), 3), fgsea_sets_padj_lt_0.01 = sum(x$padj_fgsea < 0.01),
    camera_up = sum(x$camera_direction == "Up"), camera_down = sum(x$camera_direction == "Down"),
    camera_direction_concordant = sum(x$concordant), camera_sets_FDR_lt_0.05 = sum(x$camera_FDR < 0.05),
    camera_FDR_min = signif(min(x$camera_FDR), 3),
    camera_secondary_cor001_sets_FDR_lt_0.05_descriptive = sum(x$FDR_camera_cor001_descriptive < 0.05),
    tier3 = sum(x$tier == "Tier 3"), tier4 = sum(x$tier == "Tier 4"), tier1_or_2 = sum(x$tier %in% c("Tier 1", "Tier 2")),
    le_genes_unique = nrow(g), le_genes_padj_lt_0.05 = sum(g$padj_gene < 0.05, na.rm = TRUE), le_gene_padj_min = signif(min(g$padj_gene, na.rm = TRUE), 3),
    representative_le_genes = rep_top,
    replicate_broadly = sum(x$replicate_class == "broadly condition-consistent"), replicate_partially = sum(x$replicate_class == "partially consistent"),
    replicate_ambiguous = sum(x$replicate_class == "ambiguous"), replicate_strongly_free_rep1 = sum(x$replicate_class == "strongly influenced by STV_free_rep1"),
    loo_sets_dependent_any_sample = sum(x$n_samples_dependent > 0),
    loo_dependent_STV_free_rep1 = sum(x$dependent_on_STV_free_rep1),
    loo_max_single_sample_share = round(max(dep_share), 3),
    loo_max_sample = sub("dependent_on_", "", names(dep_share)[which.max(dep_share)]),
    loo_dependent_sets_by_sample = {d <- sort(colSums(x[, dep_cols, with = FALSE]), decreasing = TRUE); d <- d[d > 0]
                                    if (length(d)) paste0(sub("dependent_on_", "", names(d)), " ", d, collapse = "; ") else "none"},
    cooks_not_retained = sum(x$cooks_status != "retained"),
    ch00_sensitive_prespecified = sum(x$ch00_sensitive_prespecified), ch00_sensitive_posthoc = sum(x$ch00_sensitive_posthoc),
    sets_with_org_identical_le_gene = sum(x$le_org_identical > 0),
    le_slot_share_org_identical = round(sum(x$le_org_identical) / sum(x$le_genes_org), 4),
    free_rep1_material = mean(x$dependent_on_STV_free_rep1) >= 0.5 | any(x$replicate_class == "strongly influenced by STV_free_rep1"),
    phase6b_class = "EXPLORATORY")
}))
mat <- merge(mat, lit, by = "theme", sort = FALSE)
# Superseded D-045 rules, kept as history (section 2.4): v1 = first coded (U1 pre-specified run only, U2 LOO only);
# v2 = Phase 7 commit d29a849 (U1 both runs, U2 LOO and weak replicate consistency). Neither is applied any more.
mat[, `:=`(rule_U1_v1_prespecified_only = ch00_sensitive_prespecified / n_sets >= 0.5,
           rule_U1_ch00 = ch00_sensitive_prespecified / n_sets >= 0.5 & ch00_sensitive_posthoc / n_sets >= 0.5,
           rule_U2_v1_loo_only = loo_max_single_sample_share >= 0.5,
           rule_U2_single_sample = loo_max_single_sample_share >= 0.5 & replicate_broadly / n_sets < 0.5,
           rule_U3_tier4 = tier4 / n_sets >= 0.5)]
mat[, superseded_phase7_class_v1 := fifelse(rule_U1_v1_prespecified_only | rule_U2_v1_loo_only | rule_U3_tier4, "UNSUPPORTED", phase6b_class)]
mat[, superseded_phase7_class_v2 := fifelse(rule_U1_ch00 | rule_U2_single_sample | rule_U3_tier4, "UNSUPPORTED", phase6b_class)]
mat[, superseded_phase7_class_v2_reason := fcase(
  rule_U1_ch00, paste0("U1: ", ch00_sensitive_prespecified, "/", n_sets, " sets ch00-sensitive (pre-specified) and ", ch00_sensitive_posthoc, "/", n_sets, " (post-hoc): reference/feature composition"),
  rule_U2_single_sample, paste0("U2: ", round(100 * loo_max_single_sample_share), " % of sets dependent on ", loo_max_sample, " in LOO and only ", replicate_broadly, "/", n_sets, " broadly replicate-consistent"),
  rule_U3_tier4, paste0("U3: ", tier4, "/", n_sets, " sets Tier 4"),
  rule_U2_v1_loo_only, paste0("EXPLORATORY with caveat: ", round(100 * loo_max_single_sample_share), " % of sets dependent on ", loo_max_sample, " in LOO (boundary), but ", replicate_broadly, "/", n_sets, " broadly replicate-consistent"),
  default = "fgsea-only (Tier 3); direction concordant with CAMERA; no CAMERA-significant set; no ROBUST/SUPPORTED route without Tier 1/2")]
stopifnot(all(mat$tier1_or_2 == 0), all(mat$camera_sets_FDR_lt_0.05 == 0), all(mat$le_genes_padj_lt_0.05 == 0))

# D-045 as revised by the owner. Evidence status: without Tier 1/2 sets (CAMERA support) no theme can be ROBUST or SUPPORTED
# WITH CAUTION; UNSUPPORTED only after an explicit review showing that the pattern itself lacks support or arises primarily
# from a technical feature class (section 2.3), never from a fixed share of LOO, ch00 or other robustness diagnostics.
d045_unsupported <- "ATP synthase / proton transport"
rng <- function(x, d = 2) paste0(fmt(min(x), d), " - ", fmt(max(x), d))
cx <- merge(chp[, .(collection, go_id, NES, NES_pre = NES_s, padj_pre = padj_s, ch_pre = ch00_sensitive)],
            ph[, .(collection, go_id, NES_post = NES_s, padj_post = padj_s, ch_post = ch00_sensitive)], by = c("collection", "go_id"))
cx <- merge(cx, tm[, .(collection, go_id, theme)], by = c("collection", "go_id"))
cx_txt <- cx[ch_pre & ch_post, .(txt = paste0(.N, " sets lose fgsea support without the organelle-identical features in both exclusion runs: NES ",
                                              rng(NES), " -> ", rng(NES_pre), " (pre-specified; padj ", rng(padj_pre), ") and ", rng(NES_post),
                                              " (post-hoc; padj ", rng(padj_post), ")")), by = theme]
keep_txt <- s[!(ch00_sensitive_prespecified & ch00_sensitive_posthoc), .(txt = paste0(go_id, " ", term_name, collapse = "; ")), by = theme]
mat[, evidence_status := fifelse(theme %in% d045_unsupported, "UNSUPPORTED", "EXPLORATORY")]
mat[, evidence_status_basis := paste0("Inferential: CAMERA ", camera_sets_FDR_lt_0.05, "/", n_sets, " sets at FDR < 0.05 (minimum FDR ", camera_FDR_min,
                                      "); ", le_genes_padj_lt_0.05, " leading-edge genes at padj < 0.05. Exploratory: ", n_sets,
                                      " fgsea sets at padj < 0.05 (Tier 3: ", tier3, "; Tier 4: ", tier4, "); CAMERA direction concordant ",
                                      camera_direction_concordant, "/", n_sets, ".")]
mat[theme %in% d045_unsupported, evidence_status_basis := paste0(
  evidence_status_basis, " Explicit review: the pattern is carried primarily by the organelle-identical feature class; ",
  cx_txt$txt[match(theme, cx_txt$theme)], ". Not ch00-sensitive: ", keep_txt$txt[match(theme, keep_txt$theme)], " (broader term).")]
# Robustness descriptor: a descriptive label only (not a significance category; it never changes evidence_status). The cut
# points are labelling conventions matching the F1 colour bands; the counts behind each label are in robustness_basis.
mat[, robustness_descriptor := fcase(
  rule_U1_ch00, "technically sensitive",
  loo_max_single_sample_share > 0.5 & replicate_broadly / n_sets < 0.5, "strongly sample-sensitive",
  replicate_broadly / n_sets < 0.75, "mixed robustness",
  loo_max_single_sample_share >= 0.5 | loo_sets_dependent_any_sample / n_sets > 0.5, "moderately sample-sensitive",
  default = "relatively stable")]
mat[, robustness_basis := paste0("Replicate consistency (D-036): ", replicate_broadly, "/", n_sets, " broadly, ", replicate_partially, " partially, ",
                                 replicate_ambiguous, " ambiguous. LOO (D-037): ", loo_sets_dependent_any_sample, "/", n_sets,
                                 " sets depend on one sample (by sample: ", loo_dependent_sets_by_sample, "). ch00 (D-038, D-043): ",
                                 ch00_sensitive_prespecified, "/", n_sets, " sets sensitive (pre-specified), ", ch00_sensitive_posthoc, "/", n_sets,
                                 " (post-hoc); organelle-identical leading-edge slots ", round(100 * le_slot_share_org_identical, 1),
                                 " %. Cook's ranking (D-033): ", n_sets - cooks_not_retained, "/", n_sets, " retained, 0 lost.")]
mat[, theme := factor(theme, levels = themes10)]; setorder(mat, theme)
fwt(mat, "theme_evidence_matrix.tsv")

# ---------------------------------------------------------------- leading-edge overlap between themes (Task 12, 14)
th_all <- c(themes10, "Other (no theme regex matched)")
le_th <- lapply(th_all, function(th) unique(lt[theme == th, gene_id])); names(le_th) <- th_all
ov <- CJ(theme_a = th_all, theme_b = th_all)[theme_a < theme_b]
ov[, `:=`(n_a = lengths(le_th[theme_a]), n_b = lengths(le_th[theme_b]))]
ov[, shared_le_genes := mapply(function(a, b) length(intersect(le_th[[a]], le_th[[b]])), theme_a, theme_b)]
ov[, jaccard_le := round(shared_le_genes / mapply(function(a, b) length(union(le_th[[a]], le_th[[b]])), theme_a, theme_b), 4)]
ov[, overlap_coefficient := round(shared_le_genes / pmin(n_a, n_b), 4)]
setorder(ov, -jaccard_le)
fwt(ov, "theme_le_overlap.tsv")

# ---------------------------------------------------------------- ethylene case study (Task 8)
gs <- readRDS("results/functional/annotation/go_gene_sets.rds")
eth_go <- c("GO:0009693" = "ethylene biosynthetic process", "GO:0009873" = "ethylene-activated signaling pathway",
            "GO:0071369" = "cellular response to ethylene stimulus", "GO:0009723" = "response to ethylene")
eth_members <- lapply(names(eth_go), function(g) gs$sets$BP[[g]]); names(eth_members) <- names(eth_go)
tested <- res[, .(gene_id, description, baseMean, log2FoldChange_MLE, stat, pvalue, padj)]
tested <- merge(tested, shr, by = "gene_id", all.x = TRUE)
tested <- merge(tested, rk, by = "gene_id", all.x = TRUE)
rx <- c(biosynthesis = "1-aminocyclopropane-1-carboxylate (synthase|oxidase)|\\bACC (synthase|oxidase)",
        perception_signalling = "[Ee]thylene receptor|[Ee]thylene response sensor|[Ee]thylene[- ]insensitive|EIN3|EIL[0-9]|EIL-|CTR1|[Cc]onstitutive triple response|EIN3-binding F-box",
        ERF_ethylene_responsive = "[Ee]thylene[- ]responsive|[Ee]thylene responsive")
tested[, sd := short_desc(description)]
tested[, desc_category := fcase(grepl(rx["biosynthesis"], sd), "biosynthesis (ACS/ACO)",
                                grepl(rx["perception_signalling"], sd), "perception / signalling",
                                grepl(rx["ERF_ethylene_responsive"], sd), "ethylene-responsive (by description; mostly ERF TFs)",
                                default = NA_character_)]
tested[, sd := NULL]
for (g in names(eth_go)) tested[, (g) := gene_id %in% eth_members[[g]]]
prior <- rd(file.path(fr, "prior_studies/prior_study_locus_lookup.tsv"))
prior_eth <- prior[grepl("ACO|ethylene|ERF|EREB|NAC|PR", paste(name_in_source, itag40_description))]
tested[, prior_study_name := prior_eth$name_in_source[match(gene_id, prior_eth$itag40_gene_id)]]
eth <- tested[!is.na(desc_category) | `GO:0009693` | `GO:0009873` | `GO:0071369` | `GO:0009723` | !is.na(prior_study_name)]
eth[is.na(desc_category), desc_category := fifelse(!is.na(prior_study_name) & !(`GO:0009693` | `GO:0009873` | `GO:0071369` | `GO:0009723`),
                                                     "prior-study locus (other annotation)", "GO-annotated only")]
eth[, go_ethylene_sets := apply(.SD, 1, function(v) paste(names(eth_go)[v], collapse = ";")), .SDcols = names(eth_go)]
eth[, (names(eth_go)) := NULL]
le_eth <- le[go_id %in% names(eth_go), .(eth_le_sets = paste(go_id, collapse = ";")), by = gene_id]
le_any <- lt[, .(n_significant_le_sets = uniqueN(paste(collection, go_id)), le_themes = paste(sort(unique(theme)), collapse = "; ")), by = gene_id]
eth <- merge(eth, le_eth, by = "gene_id", all.x = TRUE)
eth <- merge(eth, le_any, by = "gene_id", all.x = TRUE)
eth[is.na(n_significant_le_sets), n_significant_le_sets := 0L]
# descriptive per-gene replicate pattern (blind VST; no test)
obj <- readRDS("results/deseq2/deseq2_objects.rds")
vs <- assay(obj$vsd_blind)
stopifnot(identical(colnames(vs), samples))
v <- vs[eth$gene_id, , drop = FALSE]
fr_i <- 1:3; in_i <- 4:6
sgn <- sign(eth$log2FoldChange_MLE)
pairs <- sapply(seq_len(nrow(v)), function(i) mean(outer(v[i, in_i], v[i, fr_i], "-") * sgn[i] > 0))
d_all <- rowMeans(v[, in_i, drop = FALSE]) - rowMeans(v[, fr_i, drop = FALSE])
d_wo1 <- rowMeans(v[, in_i, drop = FALSE]) - rowMeans(v[, 2:3, drop = FALSE])
eth[, `:=`(pairs_in_direction_of_9 = round(pairs * 9), r_drop_STV_free_rep1 = round(d_wo1 / d_all, 3),
           free_rep1_is_extreme_free = apply(v[, fr_i, drop = FALSE], 1, function(z) which.max(abs(z - median(z))) == 1))]
vt <- as.data.table(round(v, 3)); setnames(vt, paste0("vst_", samples)); eth <- cbind(eth, vt)
eth[, `:=`(description = short_desc(description), gene_level_significant = !is.na(padj) & padj < 0.05)]
eth[, desc_category := factor(desc_category, levels = c("biosynthesis (ACS/ACO)", "perception / signalling",
                                                        "ethylene-responsive (by description; mostly ERF TFs)", "GO-annotated only",
                                                        "prior-study locus (other annotation)"))]
setorder(eth, desc_category, stat)
stopifnot(!any(eth$gene_level_significant))
fwt(eth, "ethylene_case_study.tsv")
# prior-study ethylene loci that were not tested (filtered or absent from ITAG4.0)
fwt(prior_eth[, .(locus, name_in_source, source, itag40_gene_id, tested_18277, itag40_description, padj_NA_reason)],
    "ethylene_prior_study_loci.tsv")
es <- rbind(s[go_id %in% names(eth_go)][, record := "tested set (15 - 500 genes)"],
            data.table(collection = "BP", go_id = "GO:0009693", term_name = eth_go[["GO:0009693"]],
                       size = length(eth_members[["GO:0009693"]]), record = "not tested (< 15 genes in the ranked universe)"), fill = TRUE)
fwt(es, "ethylene_go_sets.tsv")

# ---------------------------------------------------------------- ribosome / translation composition (Task 9)
rib <- s[theme == "Ribosome / translation / RNA processing"]
rib[, subtheme := fcase(grepl("ribosome biogenesis|rRNA|ribonucleoprotein complex biogenesis|preribosome|nucleol", term_name), "ribosome biogenesis / rRNA processing",
                        grepl("tRNA|RNA modification|methylation", term_name), "tRNA / RNA modification",
                        grepl("ribosom|structural constituent", term_name), "ribosome structure (cytosolic or organellar)",
                        grepl("translation|peptide|amide biosynthetic", term_name), "translation / peptide biosynthesis",
                        default = "other RNA processing")]
fwt(rib[, .(collection, go_id, term_name, subtheme, NES, padj_fgsea, camera_FDR, tier, replicate_class,
            loo_dependent_STV_free_rep1 = dependent_on_STV_free_rep1, n_samples_dependent, le_org_identical)], "ribosome_sets_by_subtheme.tsv")
rg <- unique(lt[theme == "Ribosome / translation / RNA processing", .(gene_id, description = short_desc(description), stat, padj_gene, corrected_class)])
rg[, gene_category := fcase(grepl("60S ribosomal|40S ribosomal|[Rr]ibosomal protein L[0-9]|[Rr]ibosomal protein S[0-9]", description) & !grepl("50S|30S|chloroplast|mitochond", description), "cytosolic ribosomal protein",
                            grepl("50S ribosomal|30S ribosomal|chloroplast|mitochond", description) & grepl("ribosom", description), "organellar ribosomal protein (plastid/mitochondrial)",
                            grepl("[Nn]ucleol|rRNA|snoRNA|[Ff]ibrillarin|pre-rRNA|ribosome biogenesis|WD40|[Hh]elicase|[Pp]seudouridine|[Nn]op[0-9]", description), "biogenesis / RNA processing factor",
                            grepl("[Tt]ranslation|[Ee]longation factor|tRNA", description), "translation factor / tRNA-related",
                            default = "other")]
fwt(rg[, .(n_genes = .N, median_stat = round(median(stat), 3), organelle_identical = sum(grepl("near-identical|partial", corrected_class))), by = gene_category][order(-n_genes)],
    "ribosome_le_gene_composition.tsv")

# ---------------------------------------------------------------- photosynthesis / plastid vs ch00 (Task 10)
ps <- s[theme == "Photosynthesis / plastid / thylakoid", .(collection, go_id, term_name, NES, padj_fgsea, camera_FDR, tier, replicate_class,
        loo_dependent_STV_free_rep1 = dependent_on_STV_free_rep1, n_samples_dependent, ch00_sensitive_prespecified, ch00_sensitive_posthoc,
        le_genes = le_genes_org, le_org_identical, le_org_share = round(le_org_identical / le_genes_org, 4))]
setorder(ps, padj_fgsea)
fwt(ps, "photosynthesis_sets_organelle_check.tsv")

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_20.txt"))
cat("\n== matrix\n"); print(mat[, .(theme, n_sets, fgsea_direction, NES_min, NES_max, camera_direction_concordant, tier4, replicate_broadly,
                                   loo_dependent_STV_free_rep1, loo_max_single_sample_share, loo_max_sample, ch00_sensitive_prespecified,
                                   le_slot_share_org_identical, superseded_phase7_class_v1, superseded_phase7_class_v2, evidence_status,
                                   robustness_descriptor)])
cat("\n== ethylene genes by category\n"); print(eth[, .N, by = desc_category])
cat("\n== ribosome LE composition\n"); print(fread(file.path(out, "ribosome_le_gene_composition.tsv")))
cat("\n== top theme LE overlaps\n"); print(head(ov, 12))
cat("\n20_biological_integration.R finished\n")
