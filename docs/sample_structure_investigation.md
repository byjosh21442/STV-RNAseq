# Sample-structure investigation — STV_GSE137303 (Phase 5 diagnostic checkpoint)

**Status:** COMPLETE. **Owner decision 2026-10-25 (D-022): Option A approved.** All six samples are kept, and the six-sample DESeq2 analysis remains the primary differential-expression analysis. STV_free_rep1 is not removed, and no leave-one-out analysis was run at this checkpoint. **The primary DESeq2 analysis (D-013 – D-019) is unchanged:**

- all six samples are retained;
- the design is `~ condition`, with no pairing or blocking factor;
- significance thresholds are unchanged;
- no functional enrichment was run;
- no Phase 4 output was overwritten.

The Phase 5 scripts write only to:

- `results/qc/sample_structure/`
- `results/qc/stv_signal/`
- `figures/sample_structure/`

They read the stored Phase 4 objects. Before anything else runs, `03_sample_structure_diagnostics.R` checks that it reproduces the Phase 4 VST distances and top-500 PCA coordinates (tolerance 1e-8); both checks passed.

## Key messages (for presentation)

1. **Metadata.** No documented pairing/blocking factor was identified in the public metadata examined. All six libraries were sequenced together, on the same instrument run (SN679, run 271), flowcell (HN3JNBCXY) and lanes (1 + 2, ≈ 50/50).
2. **Labels.** The STV-free / STV-infected labels are supported by the reads:
   - 0 STV reads in each free library;
   - 304 – 580 STV reads in each infected library, covering 97 – 98 % of the authors' own STV genome (LC429302.1).
3. **Structure.** `STV_free_rep1` is the most remote sample in every expression-based diagnostic:
   - VST correlation, distance and PCA, for every gene set tested;
   - Cook's-distance burden;
   - the authors' own GEO processed table;
   - the independent 2026 reanalysis.
4. **Same-number proximity.** Each STV-free sample's nearest neighbour is the STV-infected sample with the same replicate number. This holds in our VST data and in the authors' GEO table. It is **observed, not documented**, so it is not used as a pairing factor.
5. **Technical failure.** Upstream QC shows **no unique technical failure** of `STV_free_rep1`. It is at the edge of several metrics, but every difference is small:
   - mapping 87.95 % vs 89.03 – 90.26 %;
   - decoy 10.57 % vs 8.45 – 9.64 %;
   - GC 42.8 % vs 43.8 – 44.2 %.

   All six libraries pass the same QC checks, and `STV_infected_rep1` is usually the next-closest sample.
6. **Robust PCA.** With n = 6 and default parameters, robust PCA (`rrcov`) is mathematically degenerate, so the subtask was stopped (§6).
7. **Owner decision (D-022): A — keep all six samples; the primary analysis stands.** Option B (documented covariate) is not possible because no such covariate exists. Option C (leave-one-out) was not run at this checkpoint; it would leave n = 2 STV-free samples and could only ever be a labelled sensitivity analysis.

**Primary result:** 18,277 genes tested after the prespecified filter; 0 genes with padj < 0.05; 0 genes with padj < 0.05 and |log2FC| ≥ 1. No statistically significant differential expression was detected under the prespecified DESeq2 model and FDR threshold in this six-sample dataset. This is not proof that STV has no transcriptomic effect.

**QC caveat:** STV_free_rep1 is consistently separated in expression-space diagnostics, but no independent technical failure or documented experimental covariate justifies its removal.

---

## 1. Method → source map

Each diagnostic below is used only as far as its source supports it. Full citations are in [methods_references.md](methods_references.md) §17 – §24.

| Diagnostic (this checkpoint) | Source and what it actually says/does | How it is used here | Limits |
|---|---|---|---|
| VST (blind) sample distances, PCA (top 500), hierarchical clustering | DESeq2 vignette 3.22, "Data quality assessment by sample clustering and visualization" and "Approach to count outliers"; Love et al. 2014; rnaseqGene workflow (Love et al. 2015) | Phase 4 objects re-read; distances and PCA reproduced exactly | Exploratory; no test or threshold |
| PCA with other gene sets (top 100, 1000, all 18,277 retained genes) | `plotPCA` default `ntop = 500` (DESeq2); varying the gene set was requested by the owner. No publication prescribes these set sizes | **Sensitivity diagnostic only**; the primary QC PCA remains top 500 (D-018) | No set is "selected" |
| Pearson correlation; median correlation of each replicate with the others in its condition | Gierliński et al. 2015 (yeast, 48 replicates): median Pearson correlation is one of three bad-replicate criteria. Upton et al. 2023 (plant RNA-seq review): replicates with poor correlation can be identified, and *higher replication* supports removing them | Statistic reported per sample | Gierliński's cut-offs were derived from 48 replicates, so no cut-off is applied at n = 3 |
| Spearman correlation | Requested by the owner as a rank-based complement to Pearson. No external threshold is applied | Reported alongside Pearson | Descriptive only |
| Per-sample Cook's distance distribution; `boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)`; count of genes where each sample has the largest Cook's distance; count above the cutoff | DESeq2 vignette, "Approach to count outliers". Suggests a PCA plot and this boxplot "to see if one sample is consistently higher than others" when there are many outliers. Flagging needs ≥ 3 replicates; replacement needs ≥ 7. Cutoff = 99 % quantile of F(p, m − p) | Read from the existing fit. **Cook's filtering is not changed**; the 39 filtered genes are tabulated, not reinterpreted | Replacement is impossible at n = 3 per group |
| Robust PCA (`PcaGrid`, `PcaHubert`), DESeq2 `rlog` input, default parameters | Chen et al. 2020 (BMC Bioinformatics). `PcaGrid` on `rlog` had 100 % sensitivity and specificity in their tests. Smallest tested set: 3 + 3 baseline plus 1 added sample (n = 7); `vst` input and `PcaHubert` gave false positives. Algorithms: Croux et al. 2007 (`PcaGrid`), Hubert et al. 2005 (ROBPCA = `PcaHubert`), Todorov & Filzmoser 2009 (`rrcov`) | Feasibility assessed; output classified; **no removal** | Degenerate at n = 6 (§6) |
| Plant RNA-seq precedent: a replicate flagged by PCA/distance | Khorramdelazad et al. 2018 (lentil, 3 replicates per group): VST, distance matrix and PCA; sample 4C "clustered with the 6 hpi samples … identified as an outlier and was excluded" | Shows the practice exists in plant RNA-seq | Their outlier clustered with a *different experimental group*, which is not the case here |
| Caution against removing biological outliers | Chen et al. 2020: "removing biological outliers will result in underestimation of the natural biological variance and will increase the risk of spurious conclusions". Their KO-3 was removed after records explained it (freezer delay) | Weighed in the decision matrix (§10) | — |
| Replicate number and power | Schurch et al. 2016: with 3 replicates, 20 – 40 % of significantly DE genes are detected; ≥ 6 replicates recommended; "a significant minority chance that one or more replicates within each condition should be rejected" | Weighed for option C (n_free = 2) | Yeast data |
| Batch effects | Upton et al. 2023: randomize library preparation and sequencing runs | Read names checked for run/flowcell/lane | — |
| Paired/blocked designs | DESeq2 vignette: paired or blocked designs (e.g. `~ subject + condition`) use a *known* subject/batch column in `colData`. The vignette offers no procedure to infer pairing from expression | No pairing factor introduced | — |
| STV read detection | Salmon (Patro et al. 2017) competitive selective alignment against host cDNA + STV genomes. LC429302.1 (Fukuhara/Tabara/Takahashi, strain M82); KX949574.1 (Alcalá-Briseño et al. 2017). Upton et al. 2023: non-polyadenylated viral RNA is under-detected by mRNA-focused libraries | Label support only | Abundance not interpreted |

## 2. Task 1 — metadata, batch and pairing

**Sources examined:**

- the GEO SOFT record (`data/metadata/sources/GSE137303_brief.soft.txt`);
- BioSample XML (`PRJNA565137_biosamples.xml`, retrieved 2026-10-26T04:39:41Z);
- BioProject XML;
- SRA RunInfo / experiment packages and ENA read-run records;
- the original Illumina read names inside the FASTQ headers (`scripts/bash/fastq_header_audit.sh` → `results/qc/sample_structure/fastq_header_lanes.tsv`);
- the GEO series design text.

The Fukuhara et al. 2020 full text is subscription-only and could not be accessed, so its methods are **UNKNOWN** here.

### 2a. Per-sample identifiers (all DOCUMENTED)

| Field | STV_free_rep1 | STV_free_rep2 | STV_free_rep3 | STV_infected_rep1 | STV_infected_rep2 | STV_infected_rep3 |
|---|---|---|---|---|---|---|
| Condition (GEO `virus infection`) | STV-free | STV-free | STV-free | STV-infected | STV-infected | STV-infected |
| Replicate (GEO title) | rep1 | rep2 | rep3 | rep1 | rep2 | rep3 |
| GEO sample | GSM4074979 | GSM4074980 | GSM4074981 | GSM4074982 | GSM4074983 | GSM4074984 |
| GEO title | Leaf, STV-free, rep1 | Leaf, STV-free, rep2 | Leaf, STV-free, rep3 | Leaf, STV-infected, rep1 | Leaf, STV-infected, rep2 | Leaf, STV-infected, rep3 |
| GEO description = original ID (also the original SRA file name `N.fastq`) | Sample 12 | Sample 432 | Sample 433 | Sample 22 | Sample 412 | Sample 422 |
| BioSample | SAMN12735328 | SAMN12735327 | SAMN12735326 | SAMN12735325 | SAMN12735324 | SAMN12735323 |
| SRA sample | SRS5376422 | SRS5376423 | SRS5376424 | SRS5376425 | SRS5376426 | SRS5376427 |
| SRA experiment | SRX6834734 | SRX6834735 | SRX6834736 | SRX6834737 | SRX6834738 | SRX6834739 |
| SRA run | SRR10102791 | SRR10102792 | SRR10102793 | SRR10102794 | SRR10102795 | SRR10102796 |
| Reads (SRA = ENA = SeqKit) | 18,080,482 | 16,288,249 | 14,318,362 | 19,230,359 | 16,858,900 | 19,862,460 |
| Reads lane 1 / lane 2 (read names) | 8,989,348 / 9,091,134 | 8,104,250 / 8,183,999 | 7,110,159 / 7,208,203 | 9,553,168 / 9,677,191 | 8,368,492 / 8,490,408 | 9,859,069 / 10,003,391 |

### 2b. Shared attributes and candidate blocking factors

| Field | Value (same for all six unless stated) | Status | Source |
|---|---|---|---|
| Project | GEO GSE137303; BioProject PRJNA565137; SRA study SRP221436 | DOCUMENTED | GEO, BioProject, SRA |
| Cultivar | M82 | DOCUMENTED | GEO characteristics, BioSample |
| Age | 1-month-old | DOCUMENTED | GEO characteristics, BioSample |
| Tissue | leaf | DOCUMENTED | GEO characteristics, BioSample |
| Instrument | Illumina HiSeq 2500 (instrument ID SN679) | DOCUMENTED | GEO/SRA model; instrument ID from read names |
| Run / flowcell / lanes | Run 271, flowcell HN3JNBCXY, lanes 1 and 2 (lane-1 fraction 0.4964 – 0.4976), 64 tiles per lane per sample | DOCUMENTED (read names) | `fastq_header_lanes.tsv` |
| Submission | SRA961017 (center GEO); BioSamples submitted 2019-09-12; SRA LoadDate 2019-09-12; `LibraryName` empty | DOCUMENTED | SRA RunInfo, BioSample XML |
| Extraction protocol | Trizol (identical text for all six) | DOCUMENTED (protocol only) | GEO `extract_protocol` |
| Library kit | TruSeq stranded mRNA (identical text) | DOCUMENTED (protocol only) | GEO `extract_protocol` |
| Extraction date | not reported | **UNKNOWN** | — |
| Library-prep batch | not reported | **UNKNOWN** | — |
| Sequencing batch | one run and flowcell, both lanes for every library (multiplexed together) | DOCUMENTED: **no sequencing-batch difference** | read names |
| Plant identity | GEO overall design: "three 1-month-old STV-infected and three 1-month-old STV-free tomato plants of cv. M82" → one plant per library is INFERRED; no plant IDs | INFERRED / UNKNOWN | GEO series design |
| Plant line / seed lot | The 2026 reanalysis states that the STV-free line was a cured M82 sub-line and the infected line an M82 sub-line naturally harbouring STV, citing Fukuhara et al. 2020 | **Secondary, unverified** (primary text not accessible) | Badeli et al. 2026 |
| Growth batch / greenhouse position / sampling time | not reported | **UNKNOWN** | — |
| Pairing / blocking | none reported in GEO, BioSample, BioProject or SRA | **No documented pairing/blocking factor was identified in the public metadata examined.** | all of the above |

**Original IDs.** The GEO descriptions are two-digit codes for both rep1 samples ("12", "22") and three-digit codes beginning with 4 for rep2/rep3 ("432", "433", "412", "422"). The strings are DOCUMENTED, but their **meaning is UNKNOWN**: they may encode plant, tray, date or nothing. They are **not** used as a covariate. The only way to resolve them is to ask the submitters (GEO contributors: Fukuhara, Tabara, Koiwa, Takahashi).

**Replicate numbers are not treated as pairs.** Replicate numbers appear only in GEO titles. Nothing documents that, for example, `STV_free_rep2` and `STV_infected_rep2` share a plant, tray, date or library batch.

## 3. Task 2 — STV signal (label support only)

**Script:** `scripts/bash/stv_signal_check.sh` (env `stv-rnaseq`). It uses a separate index and separate outputs; the host quantification under `results/quantification/salmon/` is not read or modified.

| Item | Value |
|---|---|
| Primary STV reference | **LC429302.1** "Southern tomato virus M82 genomic RNA, complete genome", 3,376 nt, strain M82, Japan; submitted 2018-10-16 by Fukuhara, Tabara & Takahashi (the GSE137303 authors). SHA-256 `5134f44b…61cf`. <https://www.ncbi.nlm.nih.gov/nuccore/LC429302.1> |
| Independent cross-check | **KX949574.1** "Southern tomato virus isolate Florida, complete genome", 3,433 nt (Alcalá-Briseño et al. 2017, *Genome Announc* 5(7):e01374-16; used by the 2026 reanalysis). SHA-256 `83a78d49…0b50`. <https://www.ncbi.nlm.nih.gov/nuccore/KX949574.1> |
| Tool | Salmon 2.7.0 (selective alignment) |
| Index | `salmon index -t reference/STV_diagnostic/itag4_cdna_plus_stv.fa -i reference/STV_diagnostic/salmon_itag4_cdna_plus_stv_k31 -k 31 -p 16`, where the FASTA = ITAG4.0 cDNA + LC429302.1 + KX949574.1 (competitive; no decoys) |
| Quantification (per sample) | `salmon quant -i reference/STV_diagnostic/salmon_itag4_cdna_plus_stv_k31 -l U -r data/raw/<run>.fastq.gz -p 16 -o results/qc/stv_signal/salmon/<sample> --writeMappings >(awk … keep STV records)`. `-l U` is deliberate, so reads from either strand of the dsRNA genome count. The expected `-l U` vs observed `SR` warning was logged. |
| Outputs | `results/qc/stv_signal/stv_signal_summary.tsv`, `stv_reference_provenance.tsv`, `stv_signal_method.txt` (tracked); `stv_mappings/*.stv.sam` (420 KB, git-ignored by the existing `*.sam` rule; regenerated by the script) |

| Sample | Reads processed | STV reads (Salmon NumReads) | of which LC429302.1 / KX949574.1 | Fraction of reads | Per million | LC429302.1 covered | Primary alignments rev / fwd |
|---|---:|---:|---|---:|---:|---:|---|
| STV_free_rep1 | 18,080,482 | 0 | 0 / 0 | 0 | 0 | 0 % | 0 / 0 |
| STV_free_rep2 | 16,288,249 | 0 | 0 / 0 | 0 | 0 | 0 % | 0 / 0 |
| STV_free_rep3 | 14,318,362 | 0 | 0 / 0 | 0 | 0 | 0 % | 0 / 0 |
| STV_infected_rep1 | 19,230,359 | 386 | 386 / 0 | 0.0020 % | 20.07 | 97.93 % | 386 / 0 |
| STV_infected_rep2 | 16,858,900 | 580 | 548.6 / 31.4 | 0.0034 % | 34.40 | 98.16 % | 579 / 1 |
| STV_infected_rep3 | 19,862,460 | 304 | 304 / 0 | 0.0015 % | 15.31 | 97.01 % | 302 / 2 |

**Interpretation (label support only).** The signal is all-or-nothing and matches the published labels, so the STOP condition was not triggered. Reads are assigned almost entirely to the authors' own strain sequence. The low absolute counts are consistent with Upton et al. (2023)'s general caution about mRNA-focused libraries; they are not interpreted. STV abundance is **not** used biologically or as a covariate.

## 4. Task 3 — correlation, distance and PCA sensitivity (existing VST)

**Script:** `scripts/R/03_sample_structure_diagnostics.R` (env `stv-r`).

**Input:** the stored `vst(blind = TRUE)` of the 18,277 retained genes.

**Outputs:**

- `results/qc/sample_structure/correlation_{pearson,spearman}_vst.tsv`
- `pairwise_sample_distances_vst.tsv`
- `per_sample_distance_correlation_summary.tsv`
- `pca_sensitivity_{variance,coordinates}.tsv`
- figures: `figures/sample_structure/correlation_heatmap_{pearson,spearman}_vst.png` and `pca_sensitivity_vst.png`

### 4a. Correlation matrices

Pearson is shown in the upper triangle and Spearman in the lower triangle.

| | free_rep1 | free_rep2 | free_rep3 | inf_rep1 | inf_rep2 | inf_rep3 |
|---|---:|---:|---:|---:|---:|---:|
| free_rep1 | 1 | 0.9290 | 0.9263 | 0.9532 | 0.9255 | 0.9425 |
| free_rep2 | 0.9203 | 1 | 0.9548 | 0.9411 | **0.9840** | 0.9546 |
| free_rep3 | 0.9176 | 0.9506 | 1 | 0.9685 | 0.9623 | **0.9838** |
| inf_rep1 | 0.9461 | 0.9354 | 0.9642 | 1 | 0.9506 | 0.9603 |
| inf_rep2 | 0.9160 | **0.9824** | 0.9579 | 0.9453 | 1 | 0.9607 |
| inf_rep3 | 0.9353 | 0.9502 | **0.9813** | 0.9547 | 0.9557 | 1 |

### 4b. Per-sample summary

Euclidean distance is on VST over all retained genes.

| Sample | Median distance to all 5 | Median to same condition (2) | Median to other condition (3) | Nearest sample (distance) | Median Pearson, same condition | Median Spearman, same condition |
|---|---:|---:|---:|---|---:|---:|
| STV_free_rep1 | **93.36** (rank 1) | **93.54** | 82.17 | STV_infected_rep1 (76.02) | **0.9277** | **0.9189** |
| STV_free_rep2 | 74.53 | 83.94 | 74.53 | STV_infected_rep2 (45.16) | 0.9419 | 0.9354 |
| STV_free_rep3 | 68.94 | 84.12 | 62.47 | STV_infected_rep3 (43.74) | 0.9405 | 0.9341 |
| STV_infected_rep1 | 76.02 | 74.66 | 76.02 | STV_free_rep3 (62.47) | 0.9554 | 0.9500 |
| STV_infected_rep2 | 70.55 | 74.93 | 68.94 | STV_free_rep2 (45.16) | 0.9557 | 0.9505 |
| STV_infected_rep3 | 70.55 | 70.28 | 74.53 | STV_free_rep3 (43.74) | 0.9605 | 0.9552 |

**Pairwise ranking.** Of the 15 pairs, the two closest are free_rep3–inf_rep3 (43.74) and free_rep2–inf_rep2 (45.16). Next comes free_rep3–inf_rep1 (62.47). Within-condition pairs rank 5, 6, 7, 10, 13 and 14. The three most distant pairs all involve free_rep1:

- free_rep1–free_rep2: 93.36
- free_rep1–free_rep3: 93.72
- free_rep1–inf_rep2: 96.54

**Observations.**

- Every sample's nearest neighbour comes from the *other* condition.
- Each free sample's nearest neighbour has the same replicate number (rep1→rep1, rep2→rep2, rep3→rep3). Infected rep1 and rep3 are both nearest to free_rep3.
- The same-number proximity is a descriptive observation only. It is **not documented**, so it is not a pairing factor.
- Without free_rep1, within-condition similarity is comparable between groups: free_rep2–free_rep3 Pearson 0.9548 vs infected pairs 0.9506 – 0.9607.

### 4c. PCA gene-set sensitivity

This is a **sensitivity diagnostic** and was not used to choose a PCA.

| Gene set | PC1 % | PC2 % | PC3 % | Most remote sample (distance to centroid, all PCs) | Next most remote |
|---|---:|---:|---:|---|---|
| top 100 | 52.35 | 23.94 | 17.99 | STV_free_rep1 (20.70) | STV_infected_rep1 (15.49) |
| **top 500 (primary QC, D-018)** | 49.28 | 27.17 | 16.66 | STV_free_rep1 (33.17) | STV_free_rep2 (24.34) |
| top 1000 | 47.96 | 27.95 | 16.18 | STV_free_rep1 (39.28) | STV_free_rep2 (29.36) |
| all 18,277 retained | 42.96 | 29.55 | 16.38 | STV_free_rep1 (64.44) | STV_free_rep2 (49.68) |

The configuration is the same in all four gene sets:

- free_rep1 is the extreme of PC1.
- free_rep2 and inf_rep2 lie together at the opposite end.
- free_rep3, inf_rep3 and inf_rep1 group on PC2.
- PC3 separates inf_rep1 from inf_rep3.
- No PC1/PC2 view separates the two conditions.

The result does not depend on the gene set.

## 5. Task 4 — Cook's distances (existing fit, filtering untouched)

**Outputs:**

- `results/qc/sample_structure/cooks_per_sample_summary.tsv`
- `cooks_filtered_genes_influential_sample.tsv` (39 rows; Cook's distance, raw counts and normalized counts for each sample)
- figure: `figures/sample_structure/cooks_log10_boxplot.png`

The cutoff is `qf(0.99, 2, 4)` = 18.00. There are no NA or zero Cook's distances.

| Sample | Median Cook's | 99th pct | Max | Genes above cutoff | Genes where this sample has the max Cook's (expected ≈ 3,046 = 1/6) | Cook's-filtered genes (39) where this sample is most influential |
|---|---:|---:|---:|---:|---:|---:|
| STV_free_rep1 | **0.197** | **4.88** | 34.97 | **33** | **5,463 (29.9 %)** | **19** |
| STV_free_rep2 | 0.150 | 2.64 | 18.47 | 2 | 2,371 (13.0 %) | 2 |
| STV_free_rep3 | 0.136 | 2.57 | 24.25 | 1 | 2,844 (15.6 %) | 1 |
| STV_infected_rep1 | 0.115 | 2.65 | 32.30 | 12 | 2,868 (15.7 %) | 9 |
| STV_infected_rep2 | 0.132 | 2.16 | 27.90 | 1 | 2,560 (14.0 %) | 1 |
| STV_infected_rep3 | 0.089 | 2.01 | 30.54 | 10 | 2,171 (11.9 %) | 7 |

**Reading (vignette criterion "one sample consistently higher than others").** STV_free_rep1 has the highest Cook's distances across the distribution (median, 99th percentile) and the highest burden:

- It has the largest Cook's distance in 29.9 % of genes, about 1.8× the 1/6 expected if samples were exchangeable.
- It exceeds the cutoff in 33 of the 39 filtered genes: as the top sample in 19, and alongside inf_rep1 (8) or inf_rep3 (6) in 14.

The boxplot shows a modest upward shift of the whole free_rep1 distribution, not a separate population. The 39 genes are tabulated only. They remain Cook's-filtered (`padj = NA`) as in Phase 4 and are not reinterpreted.

## 6. Task 5 — robust PCA: feasibility assessed, subtask STOPPED

**Package:**

- `rrcov` 1.7-7 (with `robustbase` 0.99-7 and `pcaPP` 2.0-5), installed **only** into the isolated `stv-r` environment;
- 4 new packages, no updates or removals;
- `envs/stv-r.yml` and the lock file were updated.

**Method:** following Chen et al. 2020:

- `PcaGrid(X)` and `PcaHubert(X)` with defaults (`crit.pca.distances = 0.975`; `PcaHubert` alpha = 0.75);
- X = samples × 18,277 genes;
- primary input `rlog(dds, blind = TRUE)` (the input Chen et al. recommend), with the stored VST as a secondary input;
- seed 20261025.

**Outputs:**

- `results/qc/sample_structure/robust_pca_rrcov_classification.tsv`
- `figures/sample_structure/robust_pca_outlier_maps.png`

| Fit | k chosen | h (ROBPCA subset) | Degeneracy | rrcov flag |
|---|---:|---:|---|---|
| `PcaGrid` on rlog | 6 | — | **k = n = 6**: no orthogonal space, so all orthogonal distances ≈ 1e-13 | STV_free_rep1 (score distance 4.23 > 3.80) |
| `PcaGrid` on VST | 6 | — | k = n | STV_free_rep1 (5.81 > 3.80) |
| `PcaHubert` on rlog | 3 | **6** | **h = n**: the "robust" covariance step uses every sample, so it has no resistance to an outlier | none |
| `PcaHubert` on VST | 3 | 6 | h = n | none |

**Why the subtask was stopped.** Neither method behaves as a robust outlier detector at n = 6 with default parameters:

- With `PcaGrid`, the automatic choice gives k = 6, so the six samples span the whole PCA subspace. The classification then rests only on a χ²(6) score-distance cutoff estimated from six points.
- With `PcaHubert` (ROBPCA), the subset size h equals n, so nothing is down-weighted.

Chen et al. validated these methods with at least 7 samples (3 + 3 + 1 added). Tuning k, alpha or the cutoff after seeing the data would be parameter selection driven by the result, so it was not done. **The `PcaGrid` flag on STV_free_rep1 is therefore not used as evidence** and does not authorize removal.

## 7. Task 6 — six-sample upstream QC

**Output:** `results/qc/sample_structure/upstream_qc_six_samples.tsv`, which also lists the minimum and maximum sample for each metric.

| Metric (source) | free_rep1 | free_rep2 | free_rep3 | inf_rep1 | inf_rep2 | inf_rep3 |
|---|---:|---:|---:|---:|---:|---:|
| Raw reads (SeqKit = ENA) | 18,080,482 | 16,288,249 | 14,318,362 | 19,230,359 | 16,858,900 | 19,862,460 |
| Q30 bases % (SeqKit) | 97 | 97 | 97 | 97 | 97 | 97 |
| Mean GC % (SeqKit) | **42.81** | 44.03 | 44.14 | 43.80 | 44.16 | 43.78 |
| FastQC per-sequence GC | PASS | PASS | WARN | WARN | WARN | WARN |
| Max adapter % (FastQC, poly-A) | **0.209** | 0.110 | 0.144 | 0.191 | 0.110 | 0.156 |
| Duplicate reads % (FastQC/MultiQC) | 68.80 | 65.41 | 64.80 | 69.48 | 66.64 | 67.46 |
| Overrepresented sequences (n) | 3 | 0 | 3 | 2 | 2 | 2 |
| Run:flowcell / lanes | SN679:271:HN3JNBCXY / 1+2 | same | same | same | same | same |
| Lane-1 fraction | 0.4972 | 0.4976 | 0.4966 | 0.4968 | 0.4964 | 0.4964 |
| Salmon mapping % | **87.95** | 90.24 | 89.89 | 89.03 | 90.26 | 89.97 |
| Salmon decoy % of processed | **10.57** | 8.45 | 8.81 | 9.64 | 8.47 | 8.69 |
| Salmon SR-compatible fragment ratio | 0.99237 | 0.99394 | 0.99396 | 0.99331 | 0.99427 | 0.99362 |
| Salmon strand-mapping bias | 0.0213 | 0.0168 | 0.0172 | 0.0209 | 0.0165 | 0.0177 |
| Library size (tximport estimated counts) | 15,902,183 | 14,698,030 | 12,871,350 | 17,121,521 | 15,217,015 | 17,871,112 |
| Genes detected (count > 0) | 21,995 | 21,655 | 21,723 | 21,826 | 21,561 | 22,262 |
| Size factor (counts / avgTxLength) | 1.077 | 0.999 | 0.848 | 1.031 | 0.924 | 1.212 |
| PC1 / PC2 (top 500) | −29.70 / −13.33 | 17.73 / −15.50 | 7.76 / 16.96 | −12.79 / 8.39 | 17.83 / −7.65 | −0.83 / 11.12 |
| Median distance, all / same condition | **93.36 / 93.54** | 74.53 / 83.94 | 68.94 / 84.12 | 76.02 / 74.66 | 70.55 / 74.93 | 70.55 / 70.28 |
| Cook's: above cutoff / is max / of 39 | **33 / 5,463 / 19** | 2 / 2,371 / 2 | 1 / 2,844 / 1 | 12 / 2,868 / 9 | 1 / 2,560 / 1 | 10 / 2,171 / 7 |
| STV reads (per million) | 0 | 0 | 0 | 386 (20.1) | 580 (34.4) | 304 (15.3) |

**Technical-failure assessment.**

- **Not extreme:** STV_free_rep1 is in the middle for depth, library size, size factor and detected genes. Its Q30, lanes, run and flowcell are identical to the other libraries.
- **At the edge:** it is at the edge for GC, adapter content, mapping, decoy fraction, SR-compatible ratio and strand bias. The margins are small (largest: mapping −1.08 percentage points and decoy +0.93 points vs the next sample), and STV_infected_rep1 is the next sample on most of these metrics.
- **No failed check:** every library passed the same raw-QC and Salmon checks. No FastQC module fails for free_rep1 alone; the two FAIL modules, per-base content and duplication, fail in all six (typical for RNA-seq; see `docs/raw_qc_assessment.md`).

**Conclusion:** no unique technical failure of STV_free_rep1 was identified. Its distinctness is visible only at the expression level, and the public metadata cannot tell biological from undocumented technical causes. Ordinary biological variability is **not** labelled a failure here.

## 8. Task 7 — historical processed-data structure check (not used for differential-expression inference)

**Script:** `scripts/R/04_geo_historical_structure_check.R`.

**Input:** `data/processed/original_geo/GSE137303_count.txt.gz`, MD5 verified against its sidecar. This is Cuffnorm 2.2.1 output on SL2.50 with NCBI RefSeq identifiers. The column → sample mapping was verified in [geo_counts_audit.md](geo_counts_audit.md).

**Why it is usable for structure only.** The FPKM columns are documented by GEO as Cuffnorm FPKM. The "count" columns are INFERRED to be library-size-normalized, not raw, so neither table is valid DESeq2 input (D-005). Both are used only to describe sample-to-sample structure:

- `log2(x + 1)` on the 23,838 rows that are non-zero in at least one sample;
- Pearson/Spearman correlation, Euclidean distance, and PCA on the top 500 variable rows.

Every output file carries the header "Historical processed-data structure check - not used for differential-expression inference". Outputs are in `results/qc/sample_structure/geo_historical/`; figures are `figures/sample_structure/geo_historical_{pca,spearman_fpkm}.png`.

| Question | GEO FPKM | GEO "count" |
|---|---|---|
| Is free_rep1 unusual? | **Yes**: largest median distance to all others (101.9 vs 79.9 – 88.6) and to its condition (102.5); extreme on PC1 (−32.9; PC1 40.0 %) | **Yes**: largest median distance (149.6 vs 129.9 – 140.3) |
| Same-number proximity? | **Yes**: nearest neighbours free_rep1→inf_rep1, free_rep2→inf_rep2, free_rep3→inf_rep3. Closest pairs: free_rep2–inf_rep2 (Spearman 0.9827), free_rep3–inf_rep3 (0.9815) | **Yes**: same nearest neighbours |

The authors' own pipeline (TopHat2/Bowtie2, SL2.50, Cuffnorm) shows the same structure as ours (Salmon, SL4.0/ITAG4.0, tximport, VST). **The structure is therefore not an artefact of our quantification pipeline.** This is a statement about data structure only.

## 9. Task 8 — external reanalysis (structure only)

See [external_reanalysis_comparison.md](external_reanalysis_comparison.md). Badeli et al. 2026 (preprint) report:

- PC1 (~58 %) "separated one of the STV-free samples (labeled 'STV−_rep1') away from all other samples";
- PC2 (~20 %) set apart STV+_rep2;
- no clustering by infection status.

They did not apply an outlier test, remove a sample or run a sensitivity analysis; their model was single-factor. Their gene-level findings are **not** used here.

## 10. Task 9 — decision matrix

| Option | Evidence for | Evidence against / constraints | Assessment |
|---|---|---|---|
| **A. Keep all six (default; primary analysis)** | No documented technical failure; no documented covariate; DESeq2 already limits single-sample influence through Cook's filtering (39 genes); Chen et al. 2020 warn that removing biological outliers underestimates biological variance | free_rep1 adds variance to the STV-free group, which lowers power. At n = 3, outlier replacement is impossible (needs ≥ 7) | **Recommended; approved by the owner (D-022). The primary result stands: 0 genes at padj < 0.05.** |
| **B. Documented covariate only** | — | **None exists.** All six share run/flowcell/lanes; extraction date, library batch, plant identity, growth batch and pairing are UNKNOWN. Replicate numbers, PCA clusters and the same-number proximity are not documented and are **prohibited** as covariates. The original-ID pattern ("12"/"22" vs "4xx") has undocumented meaning | **Not applicable.** The only possible action is to ask the submitters what the original IDs and replicate numbers mean |
| **C. Leave-one-out excluding STV_free_rep1 (NOT executed)** | free_rep1 is the most remote sample in: VST correlation/distance; all four PCA gene sets; Cook's burden; GEO FPKM; GEO counts; the 2026 reanalysis | No upstream QC failure or record explaining it (unlike Chen's KO-3). The `PcaGrid` flag is degenerate (§6). **n_free would become 2**: the DESeq2 vignette says at least 3 replicates are needed for Cook's flagging, and power falls further (Schurch et al. 2016). Upton et al. 2023 tie removal to *higher* replication. Excluding a sample after seeing the results risks selective analysis | **Not recommended as a replacement.** If the owner authorizes it, it must be labelled "Sensitivity analysis (leave-one-out, STV_free_rep1 excluded)", reported next to the primary analysis, and must never replace it |

**Recommendation: preserve the primary analysis only (A).** No documented covariate exists to investigate (B). A leave-one-out run (C) is not needed for scientific validity.

**Owner decision (2026-10-25, D-022): Option A approved.** Rationale recorded by the owner:

- STV_free_rep1 is consistently the most distant expression sample;
- however, no unique technical failure was identified;
- the published infection label is supported by the independent STV-read diagnostic;
- the same structure is visible in the historical GEO processed data;
- an independent GSE137303 reanalysis also observed separation of STV−_rep1;
- no documented batch, block or pairing factor exists;
- removing the sample would reduce the STV-free group to n = 2.

The leave-one-out sensitivity analysis was not run at this checkpoint.

**Limitations kept for later discussion:**

- n = 3 biological replicates per condition (limited power; Schurch et al. 2016);
- the observed sample structure (STV_free_rep1 separation; undocumented same-replicate-number proximity);
- unknown plant, growth, extraction and library batches;
- an unverified secondary report of separate M82 sub-lines.

## 11. Files

| Kind | Path |
|---|---|
| Scripts | `scripts/bash/fastq_header_audit.sh`, `scripts/bash/stv_signal_check.sh`, `scripts/R/03_sample_structure_diagnostics.R`, `scripts/R/04_geo_historical_structure_check.R` |
| Tables (to be tracked on commit) | `results/qc/sample_structure/*.tsv`, `results/qc/sample_structure/geo_historical/*.tsv`, `results/qc/stv_signal/*.tsv`, `results/qc/stv_signal/stv_signal_method.txt` |
| Run logs and sessionInfo | `results/qc/sample_structure/sample_structure_run.log`, `sample_structure_sessionInfo.txt`, `geo_historical/geo_historical_run.log`, `geo_historical_sessionInfo.txt` |
| Figures | `figures/sample_structure/` (7 PNG files) |
| Provenance | `data/metadata/sources/PRJNA565137_{biosamples,bioproject}.xml`, `STV_{LC429302.1,KX949574.1}.gb`, `phase5_retrieved_utc.txt` |
| Git-ignored | `reference/STV/`, `reference/STV_diagnostic/` (1.1 GB index), `results/qc/stv_signal/salmon/`, `results/qc/stv_signal/stv_mappings/*.stv.sam`, `results/deseq2/deseq2_objects.rds` |
