# agriGO v2 feasibility — STV_GSE137303 (Phase 6, Task 8; D-030)

**Conclusion: C — agriGO is not suitable for this dataset/version and should not be used.** No agriGO job was run. agriGO was not forced into the analysis because it appeared in the original proposal.

Evidence table: [`results/functional/agrigo/agrigo_feasibility_evidence.tsv`](../results/functional/agrigo/agrigo_feasibility_evidence.tsv). Sources:
- agriGO v2 home, manual and download pages, accessed 2026-10-26;
- Tian et al. 2017 (doi:10.1093/nar/gkx382);
- Kim & Volsky 2005 (doi:10.1186/1471-2105-6-144).

## 1. Tomato identifier support

agriGO v2 (version 2.0; site "Last maintain: 2021-04-04"; species news up to 2026-04-10) lists these *Solanum lycopersicum* backgrounds:

| Background | GO version | Annotated IDs | Source |
|---|---|---|---|
| Transcript ID / gene ID (ITAG3.2) | 2016 | 20,036 | Sol Genomics Network |
| **Transcript ID (ITAG4.0 version)** | 2016 | 25,407 | **"User defined"** |
| ITAG2.4 loci / transcripts | 2010 | 19,662 | SGN |
| ITAG2.3 loci (Phytozome) | 2010 | 16,941 | Phytozome |
| Gramene R50, PLAZA 3.0, UniProt 2016 | 2016 | 20,933 – 26,757 | various |

## 2. Compatibility with ITAG4.0

- An ITAG4.0 transcript-ID background exists. It is labelled "User defined", and its annotation source, evidence codes and filtering are not documented.
- Its GO is the agriGO 2016 dataset. agriGO's ontology tables are GO Consortium 2009-10 and **2016-8** (45,110 terms). Our gene sets use GO release 2026-08-05 (`data-version 2026-07-26`; 48,340 terms).
- **Identifier mapping is not the obstacle.** Every ITAG4.0 gene in our ranking has exactly one transcript (`reference/ITAG4.0_tx2gene.tsv`, 1:1), so a gene → transcript ID conversion is exact and defensible. Option **B** ("suitable only after a defensible ID mapping") therefore does not describe the problem.

## 3. SEA requirements

- SEA needs a query **gene list** (plus a background). With **0 FDR-significant genes**, the only possible lists would be nominal-p or relaxed-threshold lists, which the analysis plan forbids.
- **Not authorized; not run.**

## 4. PAGE requirements and assumptions

- **Input**: identifiers plus ≥ 1 numeric column, "fold change (FC) or log2-transformed FC value (latter preferred)". This is compatible in format: we could supply 18,238 genes with apeglm or MLE log2FC.
- **Method**: parametric Z-score of the mean value of a GO set against all genes, using the Central Limit Theorem (Kim & Volsky 2005).
  - This is a competitive test that treats genes as independent — the same limitation as preranked GSEA with gene-set permutation (Tamayo et al. 2016; Wu & Smyth 2012). PAGE adds no protection against inter-gene correlation.
  - It tests the same question (GO-set shift in the infected-vs-free comparison) as the primary fgsea analysis.
- **Its input metric** is a fold change. The primary ranking uses the Wald statistic (sign plus precision; D-023). Running PAGE on MLE log2FC would weight noisy low-count estimates, and running it on the Wald statistic departs from the tool's documented input.

## 5. Can a custom annotation or background fix the version mismatch reproducibly?

- **Annotation**: PAGE accepts an own GO annotation file "limited to 4MB". Our PLAZA-based direct annotations (103,133 gene–term pairs) would fit. However, agriGO places uploaded terms on its own 2016-8 ontology:
  - terms created or redefined after 2016 cannot be resolved;
  - agriGO's handling of `alt_id`/obsolete terms and its propagation rules are undocumented.
  
  The resulting sets would differ from the 2026 GO sets used for the primary analysis, without a documented way to reconcile them.
- **Reproducibility**:
  - agriGO is web-only;
  - there is no versioned or programmatic interface;
  - job results are deleted after a few days (the manual states "THREE DAYS", then "After 7 days all information of the job will be deleted");
  - the server reports intermittent maintenance.
  
  A run could not be re-executed from the project environment.

## 6. Decision

- **A** (PAGE methodologically and technically suitable): **No**:
  - outdated internal ontology;
  - undocumented ITAG4.0 background;
  - non-reproducible web execution;
  - no added statistical question.
- **B** (suitable only after a defensible ID mapping): **No**. The ID mapping is already exact; the blocking problems are ontology version, provenance and reproducibility.
- **C** (not suitable for this dataset/version; should not be used): **Yes**.

If a parametric gene-set test is wanted later as a complement to GSEA, it should be run locally and reproducibly on the same documented gene sets. That would be a new owner decision; nothing of this kind was done in Phase 6.
