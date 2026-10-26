# Phase 7 — Biological interpretation (STV infection, tomato M82 leaves)

Status: Phase 7, 2026-10-26. This is observational transcriptomic evidence from 3 STV-free and 3 STV-infected plants.
Supporting documents:
- [biological_integration.md](biological_integration.md) (evidence matrix, theme detail, ethylene case study, GO redundancy);
- [literature_biological_context.md](literature_biological_context.md) (STV biology, sources, comparisons with earlier
  analyses).

**Question.** What transcriptomic differences, if any, can be defensibly associated with STV infection in M82 tomato
leaves? Which apparent patterns remain exploratory because of statistical, replicate-level or methodological uncertainty?

## 1. Primary finding (gene level)

- DESeq2 tested 18,277 genes (prespecified model `~ condition`; count ≥ 10 in ≥ 3 samples).
- **No statistically significant differential expression was detected under the prespecified DESeq2 model and FDR
  threshold in this six-sample dataset: 0 genes at padj < 0.05, and 0 genes at padj < 0.05 with |log2FC| ≥ 1.**
- This is the central result. It is not proof that STV has no transcriptomic effect: with three replicates per group,
  power to detect moderate changes is limited.

## 2. Functional patterns (fgsea; coordinated rank-level shifts)

Preranked GSEA on the full DESeq2 Wald ranking found 215 GO sets at padj < 0.05 (BP 120, MF 51, CC 44). These are small,
coordinated shifts of many genes along the ranking. No single gene drives them, and no leading-edge gene is significant.
The main directions:
- **Higher rank in STV-infected leaves:** ribosome / translation; plastid / photosynthesis; ATP synthase; purine
  ribonucleotide metabolism.
- **Lower rank in STV-infected leaves:** defence response; ethylene / hormone signalling; glutathione; water deprivation /
  abiotic and oxidative stress; ER membrane components.
- **Mixed:** heme / oxidoreductase; cell wall / carbohydrate.

These patterns do not override the gene-level null result.

## 3. Robustness (correlation-aware validation and replicate diagnostics)

- **CAMERA** (estimated inter-gene correlation) found **0 of 215 sets at FDR < 0.05**. All 215 have the same direction as
  in fgsea. The fgsea patterns are therefore coordinated tendencies whose statistical significance is **method-dependent**
  and did not survive correlation-aware competitive testing. Correlation is highest in the ribosome sets (median VIF 37).
- **Replicate consistency.** 144 of 215 sets are broadly condition-consistent and none is driven mainly by STV_free_rep1.
  The ribosome theme has no broadly consistent set.
- **Leave-one-sample-out.**
  - No collection exceeds 50 % single-sample dependence (maximum 47.7 %).
  - At theme level, the ribosome sets depend mostly on STV_infected_rep1 (26/34) and the ER sets entirely on
    STV_free_rep1 (4/4).
- **Cook's ranking choice.** Immaterial: 210 of 215 sets retained, 5 supported, 0 lost.
- **ch00 / organelle-identical features.** ATP synthase (4/5 sets) and the purine-ribonucleotide module depend on them.
  Photosynthesis / plastid sets do not.
- **GO redundancy.** The 215 sets collapse to about 40 – 100 descriptive modules (73 at the chosen cut). They are not 215
  findings.

## 4. Biological themes

No theme is ROBUST or SUPPORTED WITH CAUTION: that would require sets validated by CAMERA (Tier 1 or 2), and there are
none.

| Theme | Direction in STV-infected | Class | Defensible statement |
|---|---|---|---|
| Photosynthesis / plastid | higher rank | EXPLORATORY | A plastid / photosynthesis-related directional tendency. It persists conceptually after the ch00 caveat (no ch00-sensitive set; about 1 % organelle-identical leading-edge slots) but is fgsea-only |
| Defence response | lower rank | EXPLORATORY (boundary caveat) | A defence-related lower-rank tendency, broadly replicate-consistent (7/8). Half the sets depend on STV_free_rep2 in LOO |
| Hormone / ethylene | lower rank | EXPLORATORY | Lower-ranked ethylene-related sets, the same direction as both earlier analyses of these data. Not significant at gene level (ACO1 padj 0.297) |
| Glutathione; water / abiotic stress | lower rank | EXPLORATORY | Small themes that overlap the defence signal; not independent |
| Heme / oxidoreductase; cell wall | mixed | EXPLORATORY | No coherent direction |
| Ribosome / translation | higher rank | **UNSUPPORTED** | Strongest fgsea signal, but driven mostly by one infected replicate and never broadly consistent |
| Endoplasmic reticulum | lower rank | **UNSUPPORTED** | All sets depend on the outlier sample STV_free_rep1 |
| ATP synthase / proton transport | higher rank | **UNSUPPORTED** | Sensitive to reference / feature composition (organelle-identical feature class) |

Notes:
- The defence, ethylene, glutathione and water-stress tendencies form one partly overlapping, lower-ranked stress- and
  defence-related pattern. They should be read as one exploratory pattern, not as four findings.
- A higher rank of ribosomal transcripts says nothing about ribosome amount or protein synthesis.

## 5. Comparison with Fukuhara et al. (2020)

- Same six libraries. The original used SL2.50 / ITAG2.4 (GEO), TopHat2 / Bowtie2 and Cuffnorm FPKM; its statistical
  method is not available in the accessible record.
- The original reported down-regulation of ethylene biosynthesis and signalling genes.
- This reanalysis (SL4.0 / ITAG4.0, Salmon, tximport, DESeq2, fgsea, CAMERA) agrees on the **direction** for ethylene-related
  genes and GO sets. It finds **no gene at padj < 0.05**.
- Reference, annotation, quantification, normalization, statistical framework, n = 3 and sample heterogeneity could each
  contribute to the difference; none is shown to be the cause. The original is not "wrong"; the conclusions differ in
  statistical strength, not in direction.

## 6. Comparison with Badeli et al. (2026, preprint)

- An independent reanalysis of the same data (SL3.0, Rsubread / featureCounts, DESeq2).
- It also finds **0 DEGs**, with ACO1 lower but not significant (padj 0.225; here 0.297). It also finds STV_free_rep1 as
  the PC1 outlier.
- Gene-level conclusions agree.
- This project adds a pathway-level layer: fgsea direction, which is not CAMERA-supported. It makes no mechanistic claim,
  whereas Badeli et al. call lower ethylene a "plausible explanation" of the growth phenotype.

## 7. STV biological context

- STV (*Amalgavirus lycopersici*) is a dsRNA amalgavirus with no known virions, transmitted only vertically through seed.
  It is systemic, highly prevalent, genetically very uniform, and asymptomatic in single infection.
- Reported host effects:
  - taller, higher-yielding M82 plants (Fukuhara et al.);
  - altered miRNAs with few viral siRNAs (Elvira-González et al. 2020);
  - synergism with CMV and PepMV in mixed infections (Elvira-González et al. 2021).
- An absence of gene-level differential expression is **consistent with** this persistent, asymptomatic lifestyle. It
  contrasts with the large reprogramming seen in tomato leaves infected by acute viruses (ToBRFV, TYLCV / TYLCTHV, TSWV,
  PepMV), which are not directly comparable.
- No primary study addresses STV and photosynthesis, translation or defence directly. The pathway tendencies here are
  hypotheses from this dataset, not established STV biology.

## 8. Technical limitations

- n = 3 vs 3; limited power for moderate effects.
- Sample heterogeneity:
  - STV_free_rep1 is separated in expression space (kept; D-022);
  - several themes depend on single samples in LOO (STV_free_rep1, STV_free_rep2, STV_infected_rep1).
- No gene-level FDR significance; nothing is "differentially expressed".
- Pathway significance is method-dependent: fgsea (gene-permutation null) vs CAMERA (correlation-aware). All themes are
  fgsea-only.
- GO redundancy: overlapping terms are not independent evidence.
- ch00 caveat: organelle-identical features drive the ATP-synthase and purine-ribonucleotide signals. Sequence identity is
  not genomic origin.
- Reference and annotation differences against the original study (SL2.50 / ITAG2.4) and Badeli et al. (SL3.0).
- Observational design. No documented pairing or batch factors. The STV-free and STV-infected plants may come from
  different M82 sub-lines (unverified statement in Badeli et al.), so line effects cannot be separated from infection.
- Single tissue (leaf) and one developmental stage (1 month).
- No hormone, protein or physiological measurements.
- MapMan / Mercator4 not yet available (site under maintenance).

## 9. Integrated interpretation

In M82 tomato leaves, STV infection is **not associated with any statistically significant gene-level expression
difference** in this six-sample dataset. This is consistent with the persistent, asymptomatic lifestyle described for STV,
but does not exclude subtle effects that three replicates cannot detect.

At the pathway level, the ranking shows coordinated **directional tendencies**:
- plastid / photosynthesis-related sets rank higher in STV-infected leaves;
- a partly shared set of defence-, ethylene-, glutathione- and water/abiotic-stress-related terms ranks lower.

The ethylene direction agrees with both earlier analyses of these data. These tendencies are mostly replicate-consistent
but **did not survive correlation-aware validation**. They are biologically plausible but statistically uncertain,
exploratory hypotheses for an adequately replicated experiment.

The ribosome / translation, ER and ATP-synthase signals are **not supported**: they depend on a single sample or on
organelle-identical features.

## 10. What cannot be concluded

The data do **not** support any of the following:
- STV *causes*, *induces*, *activates*, *suppresses* or *inhibits* any gene, pathway or process.
- Any individual gene is differentially expressed, including ACO1, the ERF / PR locus Solyc12g009240 and WRKY genes.
- Ethylene biosynthesis or signalling is significantly altered, or ethylene levels are lower in STV-infected plants.
- STV suppresses plant defence or immunity.
- STV increases photosynthesis or chloroplast activity.
- STV increases ribosome production or protein synthesis, or hijacks host translation.
- STV affects ER function or protein folding.
- The ATP-synthase signal reflects a nuclear transcriptional response. Organelle-identical features are not shown to be
  organellar transcripts either.
- The transcriptome explains the reported growth, yield or germination phenotypes, or demonstrates mutualism.
- STV has no effect at all (absence of evidence with n = 3).
- Any result generalizes to other cultivars, tissues, developmental stages or field conditions.
- Any conclusion depends on MapMan / Mercator4, which has not been run.
