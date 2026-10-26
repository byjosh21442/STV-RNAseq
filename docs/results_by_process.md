# Results by process — STV_GSE137303

One section per process, with the key verified numbers and links to the detailed documents. Phases 1 – 4 are summarized retrospectively from their committed outputs. Phase 5 is a **diagnostic checkpoint** that did not change the primary analysis.

## Primary result (Phase 4; D-013 – D-019)

- **18,277 genes tested** after the prespecified filter (count ≥ 10 in ≥ 3 samples).
- **0 genes with padj < 0.05** (Benjamini–Hochberg, `alpha = 0.05`).
- **0 genes with padj < 0.05 and |log2FC| ≥ 1** (MLE log2FC).
- No statistically significant differential expression was detected under the prespecified DESeq2 model and FDR threshold in this six-sample dataset. This is **not** proof that STV has no transcriptomic effect.

## QC caveat (Phase 5; owner decision D-022)

STV_free_rep1 is consistently separated in expression-space diagnostics, but no independent technical failure or documented experimental covariate justifies its removal. All six samples are kept; the six-sample analysis remains the primary analysis.

## Limitations kept for later discussion

- **n = 3 biological replicates per condition.** Power to detect moderate changes is limited: with 3 replicates, Schurch et al. (2016) report that only 20 – 40 % of significantly DE genes are detected. Only the 0-gene outcome at the prespecified thresholds is reported; no relaxed or nominal-p gene list is used.
- **Observed sample structure.** STV_free_rep1 inflates within-group variance in the STV-free group. Each STV-free sample is closest to the STV-infected sample with the same replicate number, but this is undocumented and is not modelled.
- **Undocumented design factors.** Plant identity, growth batch, extraction date and library batch are unknown. The STV-free and STV-infected plants may come from different M82 sub-lines (a secondary, unverified statement in the 2026 reanalysis).
- **No leave-one-out run.** A sensitivity analysis without STV_free_rep1 was not run (D-022). If ever run, it would have n = 2 STV-free samples and could not replace the primary analysis.

## Results by process

| # | Process | Key result | Detail |
|---|---|---|---|
| 1 | Dataset and metadata (GEO GSE137303, PRJNA565137) | 6 single-end 75-nt leaf libraries, cv. M82, 1-month-old, 3 STV-free + 3 STV-infected; HiSeq 2500. **No documented pairing/blocking factor** | [dataset_audit.md](dataset_audit.md), [sample_structure_investigation.md](sample_structure_investigation.md) §2 |
| 2 | Raw-data acquisition (ENA) | 6/6 files MD5 = ENA; reads 14,318,362 – 19,862,460 (SeqKit = ENA = SRA) | [raw_data_validation.md](raw_data_validation.md) |
| 3 | Read QC (SeqKit, FastQC, MultiQC) | Q30 97 % in all; mean per-base Q 33.5 – 38.8; adapters ≤ 0.21 %; per-base content and duplication FAIL in all six (typical of RNA-seq). No trimming (D-009) | [raw_qc_assessment.md](raw_qc_assessment.md) |
| 4 | Sequencing layout (read names, Phase 5) | All six libraries on instrument SN679, run 271, flowcell HN3JNBCXY, lanes 1 + 2 (lane-1 fraction 0.496 – 0.498): **one sequencing batch** | `results/qc/sample_structure/fastq_header_lanes.tsv` |
| 5 | Reference (SGN SL4.0 + ITAG4.0) | 34,075 genes/mRNAs; decoy-aware Salmon index (k = 31) | [reference_validation.md](reference_validation.md), [reference_decision.md](reference_decision.md) |
| 6 | Quantification (Salmon 2.7.0, `SR`, `--gcBias`) | Mapping 87.95 – 90.26 %; decoy 8.45 – 10.57 %; SR-compatible ratio 0.9924 – 0.9943; no diagnostics | [salmon_quantification.md](salmon_quantification.md) |
| 7 | Gene-level import (tximport 1.38.2) | 33,976 genes × 6; gene totals = Salmon totals; 9,808 genes with zero counts in all samples | [tximport_preparation.md](tximport_preparation.md) |
| 8 | Pre-filter (D-014) | count ≥ 10 in ≥ 3 samples → 18,277 genes retained (53.79 %) | `results/deseq2/deseq2_summary.tsv` |
| 9 | DESeq2 model (D-013, D-015) | `~ condition`, reference STV_free; size factors (counts/avgTxLength) 0.848 – 1.212; parametric dispersion trend (asymptDisp 0.0835, extraPois 5.16); median final dispersion 0.088; 0 non-converged genes | [decisions.md](decisions.md) D-013 – D-015 |
| 10 | Testing (D-016) | 18,238 of the 18,277 tested genes have p-values (39 Cook's-filtered, cutoff 18.0); independent-filtering threshold at quantile 0 %; **padj < 0.05: 0 genes; padj < 0.05 and \|log2FC\| ≥ 1: 0 genes**; raw p < 0.05: 192 (descriptive only, never used as a gene list) | `results/deseq2/deseq2_results_all.tsv` |
| 11 | Shrinkage (apeglm, D-017) | Completed; all genes converged | `results/deseq2/deseq2_results_shrunken.tsv` |
| 12 | Sample QC (D-018) | VST PCA top 500: PC1 49.28 %, PC2 27.17 %; no separation by condition; hclust k = 2 = {STV_free_rep1} vs rest. Figures drawn with pheatmap/ggrepel/ggplot2 (visualization only, D-020) | `figures/deseq2/` |
| 13 | STV label support (Phase 5) | 0 STV reads in free libraries; 304 / 386 / 580 in infected (15.3 – 34.4 per million; 97 – 98 % of LC429302.1 covered). Labels supported | [sample_structure_investigation.md](sample_structure_investigation.md) §3 |
| 14 | Correlation / distance / PCA sensitivity (Phase 5) | STV_free_rep1 most remote (median distance 93.4 vs 68.9 – 76.0; median within-condition Pearson 0.928 vs 0.941 – 0.961); the same configuration in top 100/500/1000/all-gene PCA; every sample's nearest neighbour is from the other condition; same-number proximity observed but not documented | §4 |
| 15 | Cook's diagnostics (Phase 5) | STV_free_rep1 has the largest Cook's distance in 29.9 % of genes (expected 16.7 %); above the cutoff in 33 of the 39 filtered genes (top sample in 19). Filtering unchanged | §5 |
| 16 | Robust PCA (Phase 5) | `rrcov` defaults are degenerate at n = 6 (`PcaGrid` k = n; `PcaHubert` h = n); subtask stopped; no flag used | §6 |
| 17 | Upstream QC synthesis (Phase 5) | **No unique technical failure** of STV_free_rep1: small edge values (mapping −1.08 points, decoy +0.93 points vs next sample); same run/flowcell/lanes; all QC passed | §7 |
| 18 | GEO processed table — historical structure check (Phase 5) | Not used for DE inference. STV_free_rep1 most remote in both FPKM and Cuffnorm "count" columns; same nearest-neighbour pattern | §8 |
| 19 | External reanalysis (Phase 5) | Badeli et al. 2026 (preprint) independently saw STV−_rep1 separated on PC1; no outlier handling | [external_reanalysis_comparison.md](external_reanalysis_comparison.md) |
| 20 | Decision matrix and owner decision (Phase 5) | **Option A approved (D-022): all six samples kept**; the six-sample DESeq2 analysis remains primary. B not applicable (no documented covariate). C (leave-one-out) not run at this checkpoint: n_free would be 2 | §10 |
