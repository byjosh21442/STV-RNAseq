# External reanalysis comparison — methods and sample structure only

**Scope (owner instruction):** compare methods and the sample-structure / QC observations of the 2026 reanalysis with this project. **Its genes, pathways, enrichment results and biological conclusions are not used**, and nothing here changes the primary analysis.

## Source

Badeli, G., Kaboosi, K., & Nasrollanejad, S. (2026). *Reanalysis of the tomato transcriptome revealed host responses to Southern tomato virus infection* [Preprint, Research Square, "In Review"; posted 2026-02-17]. <https://doi.org/10.21203/rs.3.rs-8823240/v1> (CC BY 4.0).

- **Access:** HTML full text <https://www.researchsquare.com/article/rs-8823240/v1.html>, read 2026-10-25; metadata checked on Crossref.
- **Status:** **not peer reviewed**; cited as a methods comparison only.

## Methods side by side

| Step | Badeli et al. 2026 (as stated) | This project |
|---|---|---|
| Input | GSE137303 raw FASTQ, 6 SE libraries | Same (ENA FASTQ, MD5-verified) |
| Trimming | Trimmomatic 0.39, default settings | None (D-003/D-009: QC passed; Salmon soft-clips) |
| Reference | Ensembl SL3.0 + STV KX949574 as an extra "chromosome" | SGN SL4.0 decoys + ITAG4.0 cDNA (D-001). STV was checked separately: LC429302.1 primary, KX949574.1 cross-check |
| Alignment / quantification | Rsubread 2.13.3 `align` (genome) + `featureCounts` (gene level, unique, multi-mappers excluded); 34,075 gene features | Salmon 2.7.0 selective alignment, `-l SR`, `--gcBias` → tximport 1.38.2 (`avgTxLength` offsets); 33,976 genes in tximport |
| Pre-filter | genes with < 10 total counts across all 6 samples removed | count ≥ 10 in ≥ 3 samples (D-014): 18,277 retained |
| Normalization | DESeq2 median-ratio size factors | DESeq2 with tximport normalization factors |
| Model | DESeq2 1.42.0, single-factor `condition` (STV+ vs STV−), Wald test | DESeq2 1.50.2, `~ condition`, Wald test (D-013/D-015) |
| Significance | BH FDR < 0.05 **and** \|log2FC\| ≥ 1; exploratory FDR < 0.10 / < 0.20 | padj < 0.05 (primary); stringent padj < 0.05 and \|MLE log2FC\| ≥ 1 (D-016); no relaxed threshold |
| Shrinkage | apeglm | apeglm (D-017) |
| Sample QC visualization | VST; PCA on top 500 variable genes; Euclidean distance heatmap on all genes | VST (blind); PCA top 500; Euclidean distances on all retained genes; complete-linkage clustering (D-018). This checkpoint adds correlation, gene-set sensitivity, Cook's and robust-PCA diagnostics |
| Outlier handling | None reported (no outlier test, no removal, no sensitivity analysis) | None applied; diagnostics only (this checkpoint) |
| Design wording | Calls the experiment a "paired design", but the model is single-factor; no pairing variable is described | No documented pairing; none used |
| Software environment | R 4.5.1, Windows 10 | R 4.5.3, conda `stv-r`, WSL Ubuntu |

**Internal inconsistency in the source.** The preprint lists R 4.5.1 together with DESeq2 1.42.0, a Bioconductor 3.18 / R 4.3 release. This is recorded, not resolved.

## Sample-structure observations

| Observation | Badeli et al. 2026 | This project |
|---|---|---|
| Separation of STV_free_rep1 | PC1 (~58 %) "separated one of the STV-free samples (labeled 'STV−_rep1') away from all other samples … this particular replicate had an atypical expression profile" | Same sample: PC1 49.3 % (top 500), and the most remote sample in every gene set, correlation and distance metric and in Cook's burden ([sample_structure_investigation.md](sample_structure_investigation.md) §4 – §5) |
| Second axis | PC2 (~20 %) separated STV+_rep2 "somewhat from the others" | PC2 27.2 %: free_rep2 and inf_rep2 together vs free_rep3, inf_rep3 and inf_rep1 |
| Grouping by condition | "No clear clustering of samples by STV status". Distance clustering: "One STV-infected sample clustered with two STV-free samples, while the other two STV-infected samples formed another cluster with the remaining STV-free sample" | No grouping by condition. Every sample's nearest neighbour is from the other condition; hclust k = 2 gives {free_rep1} vs the rest |
| Contradiction inside the source | The Fig. 4 caption says "Replicates tend to cluster by condition", which contradicts the Results text above | — |
| Sequencing / mapping QC | Text says "no obvious outlier in library preparation or sequencing quality". Table 1 lists STV−_rep1 with the lowest value (80.6 %) in a percentage column (range 80.6 – 89.9 %). The column heading was not captured by text extraction, so the metric is not asserted | No unique technical failure. free_rep1 has the lowest Salmon mapping (87.95 %) and highest decoy fraction (10.57 %), with small margins (§7 of the investigation) |
| STV reads | < 0.01 % of reads aligned to STV (KX949574) | 0.0015 – 0.0034 % in infected libraries and 0 in free libraries (LC429302.1 + KX949574.1) |
| Explanation offered | "The variability could stem from slight differences in physiological state" (speculative) | Not explained; the public metadata cannot separate biological from undocumented technical causes |

## What the comparison establishes

1. An independent pipeline (different genome build, aligner, counting and normalization) found the **same** STV_free_rep1 separation. So did the authors' own Cuffnorm table (historical check, investigation §8). The structure is a property of the data, not of one pipeline.
2. The reanalysis also kept all six samples with a single-factor model and found no gene at FDR < 0.05 with \|log2FC\| ≥ 1. This matches our primary result (0 genes at padj < 0.05), but agreement in gene counts is not used as evidence for any decision here.
3. Claims in the reanalysis that come from Fukuhara et al. 2020, a source we could not access, are treated as **secondary and unverified**. These include: STV status confirmed by dsRNA extraction and RT-PCR; separate cured and naturally infected M82 sub-lines; plants grown simultaneously.
