# Results by process — STV_GSE137303

One section per process, with the key verified numbers and links to the detailed documents. Phases 1 – 3 are summarized retrospectively from their committed outputs.

## Primary result (Phase 4; D-013 – D-019)

- **18,277 genes tested** after the prespecified filter (count ≥ 10 in ≥ 3 samples).
- **0 genes with padj < 0.05** (Benjamini–Hochberg, `alpha = 0.05`).
- **0 genes with padj < 0.05 and |log2FC| ≥ 1** (MLE log2FC).
- No statistically significant differential expression was detected under the prespecified DESeq2 model and FDR threshold in this six-sample dataset. This is **not** proof that STV has no transcriptomic effect: with n = 3 biological replicates per condition the power to detect moderate changes is limited.
- Sample QC shows STV_free_rep1 separated from the other samples (D-018 report); it is examined in the next phase without changing this analysis.

## Results by process

| # | Process | Key result | Detail |
|---|---|---|---|
| 1 | Dataset and metadata (GEO GSE137303, PRJNA565137) | 6 single-end 75-nt leaf libraries, cv. M82, 1-month-old, 3 STV-free + 3 STV-infected; HiSeq 2500 | [dataset_audit.md](dataset_audit.md) |
| 2 | Raw-data acquisition (ENA) | 6/6 files MD5 = ENA; reads 14,318,362 – 19,862,460 (SeqKit = ENA = SRA) | [raw_data_validation.md](raw_data_validation.md) |
| 3 | Read QC (SeqKit, FastQC, MultiQC) | Q30 97 % in all; mean per-base Q 33.5 – 38.8; adapters ≤ 0.21 %; per-base content and duplication FAIL in all six (typical of RNA-seq). No trimming (D-009) | [raw_qc_assessment.md](raw_qc_assessment.md) |
| 4 | Reference (SGN SL4.0 + ITAG4.0) | 34,075 genes/mRNAs; decoy-aware Salmon index (k = 31) | [reference_validation.md](reference_validation.md), [reference_decision.md](reference_decision.md) |
| 5 | Quantification (Salmon 2.7.0, `SR`, `--gcBias`) | Mapping 87.95 – 90.26 %; decoy 8.45 – 10.57 %; SR-compatible ratio 0.9924 – 0.9943; no diagnostics | [salmon_quantification.md](salmon_quantification.md) |
| 6 | Gene-level import (tximport 1.38.2) | 33,976 genes × 6; gene totals = Salmon totals; 9,808 genes with zero counts in all samples | [tximport_preparation.md](tximport_preparation.md) |
| 7 | Pre-filter (D-014) | count ≥ 10 in ≥ 3 samples → 18,277 genes retained (53.79 %) | `results/deseq2/deseq2_summary.tsv` |
| 8 | DESeq2 model (D-013, D-015) | `~ condition`, reference STV_free; size factors (counts/avgTxLength) 0.848 – 1.212; parametric dispersion trend (asymptDisp 0.0835, extraPois 5.16); median final dispersion 0.088; 0 non-converged genes | [decisions.md](decisions.md) D-013 – D-015 |
| 9 | Testing (D-016) | 18,238 of the 18,277 tested genes have p-values (39 Cook's-filtered, cutoff 18.0); independent-filtering threshold at quantile 0 %; **padj < 0.05: 0 genes; padj < 0.05 and \|log2FC\| ≥ 1: 0 genes**; raw p < 0.05: 192 (descriptive only, never used as a gene list) | `results/deseq2/deseq2_results_all.tsv` |
| 10 | Shrinkage (apeglm, D-017) | Completed; all genes converged | `results/deseq2/deseq2_results_shrunken.tsv` |
| 11 | Sample QC (D-018) | VST PCA top 500: PC1 49.28 %, PC2 27.17 %; no separation by condition; hclust k = 2 = {STV_free_rep1} vs rest | `figures/deseq2/` |
| 12 | Figures and plotting dependencies (D-019, D-020) | 8 PNG figures; pheatmap 1.0.13 and ggrepel 0.9.8 approved as visualization-only | `figures/deseq2/`, [software_versions.md](software_versions.md) |
