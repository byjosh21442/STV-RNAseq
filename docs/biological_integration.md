# Biological integration: evidence matrix by theme (Phase 7)

Status: Phase 7, 2026-10-26; classification revised the same day under the owner's D-045 decision (§2) and approved as
final by the owner; **Phase 7 is closed**. This phase
integrates **frozen** results. It runs no new test of significance and no new DESeq2, fgsea or CAMERA run, and it does not
change any threshold or sample. Machine-readable version:
`results/biological_integration/theme_evidence_matrix.tsv`. Interpretation: [phase7_biological_interpretation.md](phase7_biological_interpretation.md).
Literature: [literature_biological_context.md](literature_biological_context.md).

## 1. Evidence hierarchy (kept separate throughout)

| Level | Evidence | Frozen result |
|---|---|---|
| L1 | Primary gene-level inference: DESeq2, 3 STV-free vs 3 STV-infected, 18,277 tested genes | **0 genes at padj < 0.05; 0 stringent DEGs** |
| L2 | Preranked GSEA (fgsea) on the DESeq2 Wald ranking | 215 significant GO sets (BP 120, MF 51, CC 44): coordinated rank-level shifts |
| L3 | Correlation-aware competitive test (CAMERA, estimated inter-gene correlation) | **0 of 215 sets at FDR < 0.05**; direction agrees with fgsea in 215 of 215 |
| L4 | Replicate consistency, leave-one-sample-out (LOO), Cook's ranking choice, STV_free_rep1, SL4.0ch00 / organelle-identical features | Phase 6B diagnostics plus the corrected organelle classification (D-044) |
| L5 | External literature | Used only after L1 – L4 and only for plausibility. It never chooses thresholds, sets, genes, samples or classes |

L2 never overrides L1. Every theme below is built from L2 sets and then weighed against L3 and L4.

## 2. Classification framework (D-045, revised by the owner)

Every statement keeps three levels apart:
- **Inferential.** DESeq2 at gene level (0 of 18,277 genes at padj < 0.05) and CAMERA at set level (0 of 215 sets at
  FDR < 0.05). None of the functional themes received correlation-aware FDR support.
- **Exploratory.** fgsea on the frozen DESeq2 Wald ranking (215 GO sets at padj < 0.05). These are coordinated
  rank-level tendencies.
- **Robustness.** Replicate consistency (D-036), leave-one-sample-out (LOO, D-037), Cook's ranking choice (D-033) and the
  ch00 / organelle-identical diagnostics (D-038, D-043, D-044). They describe how strongly each tendency depends on single
  samples or on feature composition. They are not significance tests.

Each theme therefore has two separate labels: an **evidence status** and a descriptive **robustness descriptor**.

### 2.1 Evidence status (dimension 1)

| Status | Meaning in this project | Themes |
|---|---|---|
| ROBUST | Convincing support from several internal analyses, including inferential validation (CAMERA-supported Tier 1 sets; D-035, D-039) | none |
| SUPPORTED WITH CAUTION | Meaningful inferential support (Tier 1 or 2), with substantial limitations | none |
| EXPLORATORY | A biologically interpretable fgsea pattern in the frozen results that lacks independent inferential validation | 9 |
| UNSUPPORTED | Only after an explicit review shows that the pattern itself lacks adequate support, is contradicted by stronger evidence, or arises primarily from a technical feature class. **Never** because a fixed share of LOO, ch00 or other diagnostics fails | 1 (ATP synthase, §2.3) |

CAMERA supports no set (Tier 1 = Tier 2 = 0), so no theme can be ROBUST or SUPPORTED WITH CAUTION.

### 2.2 Robustness descriptor (dimension 2; descriptive only)

The descriptor summarizes the frozen diagnostics of each theme:
- B = share of sets "broadly condition-consistent" (D-036);
- M = largest share of sets that depend on the same sample in LOO;
- A = share of sets that depend on at least one sample in LOO;
- the number of sets that are ch00-sensitive in both exclusion runs.

Labels are assigned in this order:
1. **technically sensitive**: at least half of the sets lose fgsea support without the organelle-identical features in
   both exclusion runs;
2. **strongly sample-sensitive**: M > 50 % and B < 50 %;
3. **mixed robustness**: 50 % ≤ B < 75 % (components of the theme differ in robustness);
4. **moderately sample-sensitive**: M ≥ 50 % or A > 50 %;
5. **relatively stable**: otherwise (B ≥ 75 %, M < 50 %, A ≤ 50 %).

Disclosure and limits:
- The labels were defined after all diagnostics were known. The cut points are labelling conventions; they are the colour
  bands already used in figure F1. They are **not** statistical boundaries, and a descriptor never changes an evidence
  status.
- Themes have 3 – 34 sets, so one set can move a theme across a cut. Boundary cases: defence and cell wall (M exactly
  50 %), glutathione (A = 2/3 with 3 sets) and hormone (A = 2/5).
- With three replicates per group, each LOO refit keeps only two replicates in one group. Gene-level sensitivity is
  already low with three replicates (Schurch et al. 2016), so some loss of fgsea support in five-sample refits is expected
  and is not by itself evidence of an artefact. D-037 counts a set as dependent on a sample only when dropping either other
  same-condition sample keeps it supported.
- The counts behind every label are in the matrix (`robustness_basis`).

**Owner approval (final, 2026-10-26; D-049).** The descriptors and their cut points are approved as implemented and are not
to be tuned further. The labels are descriptive. The thresholds are project conventions defined post hoc. They are not
inferential thresholds, they do not determine statistical significance, and they do not override DESeq2 or CAMERA.
Borderline classifications such as defence and cell wall are accepted because the underlying quantitative diagnostics stay
visible in the matrix.

### 2.3 Explicit reviews: ATP synthase and endoplasmic reticulum

**ATP synthase / proton transport: UNSUPPORTED — technically sensitive (final, owner-approved 2026-10-26).**
- *Question.* Does the existing evidence show only technical sensitivity, or that the apparent pattern is primarily
  explained by organelle-identical / ch00 features?
- *Evidence (frozen D-038 / D-043 tables).*
  - Four of the five sets are ATP-synthase-specific: GO:0006754, GO:0045259, GO:0016469 and GO:0046933 (22 – 59 genes).
  - Organelle-identical features are 4 – 8 of their 7 – 24 leading-edge genes (33 – 57 %).
  - Removing the 5 pre-specified features per set lowers NES from 1.61 – 1.88 to 1.18 – 1.38 (padj 0.29 – 0.60). The
    post-hoc rule (5 – 8 features) gives 1.08 – 1.31 (padj 0.38 – 0.72). NES falls by 21 – 37 %, and none of the four sets
    remains fgsea-significant under either rule.
  - The fifth set, "ATP metabolic process" (GO:0046034; 99 genes), is not ch00-sensitive (NES 1.81 → 1.65 / 1.58;
    padj 0.016 / 0.042). It is a broad metabolic term and does not carry the ATP-synthase-specific pattern.
- *Decision.* Removing about five features from sets of 22 – 59 genes removes the fgsea signal of every
  ATP-synthase-specific set under both rules. This is more than sensitivity: the apparent pattern is carried primarily by
  the organelle-identical feature class. UNSUPPORTED is kept, now by explicit review rather than by the superseded rule U1.
- *Scope.* This does not show that organellar ATP-synthase transcripts are unchanged. Sequence identity does not establish
  genomic origin, the library is poly(A)-selected, and no organelle-aware re-quantification was run (D-044, owner decision).

**Endoplasmic reticulum: EXPLORATORY — strongly sample-sensitive (final, owner-approved 2026-10-26).**
- *Question.* Is the pattern merely strongly sample-sensitive, or do the replicate diagnostics leave insufficient evidence
  for a meaningful group-level pattern?
- *Evidence for a group-level pattern.*
  - All four sets have all 9 infected–free sample pairs in the expected direction (P = 9/9).
  - Without STV_free_rep1, about half of the infected–free difference remains (r₁ 0.48 – 0.50).
  - No set is in the "strongly influenced by STV_free_rep1" class (r₁ < 0.25, D-036).
  - CAMERA direction agrees in 4/4 sets.
- *Sensitivity.*
  - All four sets lose fgsea support when STV_free_rep1 is left out (they stay supported when STV_free_rep2 or STV_free_rep3
    is left out).
  - Three of the four also lose it without STV_infected_rep2.
  - Only 1/4 is broadly replicate-consistent; 3/4 are partially consistent.
  - Three of the four sets are near-duplicates (272 – 274 genes; leading edges of 101 – 102 genes), so the theme is
    effectively two signals (module N07).
- *Decision.* The diagnostics show a consistent group-level direction whose size depends strongly on the outlier sample.
  They do not show that the pattern is absent, contradicted or technical. The theme is therefore EXPLORATORY and strongly
  sample-sensitive; it is the weakest of the exploratory themes. Its leading edge consists of sterol and phospholipid
  enzymes, so "protein folding" is not supported by the content.

### 2.4 History of D-045 (preserved)

| Version | Defined | Rules (UNSUPPORTED if any holds) | Rationale | Results already visible | Classes |
|---|---|---|---|---|---|
| v1, first coded | 2026-10-26, about 16:30 UTC, in `20_biological_integration.R` before its first successful output | U1: ≥ 50 % of sets ch00-sensitive in the pre-specified run. U2 v1: ≥ 50 % of sets depend on the same sample in LOO. U3: ≥ 50 % Tier 4 | Conservative, project-defined screen: a theme whose fgsea support rests mostly on one sample, on the ch00 feature class or on contradictory methods should not be read as a finding | All Phase 6B outputs and summaries (ch00 counts, LOO dependence by collection and on STV_free_rep1); the theme-level LOO shares were first seen in the v1 output | UNSUPPORTED: ribosome, ATP synthase, defence, ER, cell wall. EXPLORATORY: the other 5 |
| v2, Phase 7 (post hoc) | About 16:33 UTC the same day, after the v1 matrix had been viewed and before any interpretation was written | U1: ≥ 50 % in **both** exclusion runs. U2 v2: U2 v1 **and** < 50 % of sets broadly replicate-consistent. U3 unchanged | The Phase 7 instruction not to treat any single diagnostic as definitive | The v1 theme matrix, including the classes above | UNSUPPORTED: ribosome, ER, ATP synthase. EXPLORATORY: the other 7. Committed in `d29a849`; owner review requested |
| Owner revision | 2026-10-26, after external methodological review | No automatic rule. Evidence status from inferential evidence; UNSUPPORTED only by explicit review (§2.3); diagnostics give a descriptive robustness descriptor (§2.2) | Neither U2 version has an external basis as an inferential threshold, and robustness diagnostics are not significance tests | Everything above | §2.5 |

Notes:
- U2 v1 was a conservative project-defined rule. U2 v2 was introduced post hoc during Phase 7.
- Neither has a validated external basis as a formal inferential threshold. No pathway-level leave-one-out method paper
  was found in Phase 6B (`methods_references.md` §38), and 50 % is a project convention. The same applies to the 50 %
  cut points of U1 and U3.
- The methodological review that motivated separating inferential evidence from robustness diagnostics is documented in
  `methods_references.md` §47.
- The flags `rule_U1_v1_prespecified_only`, `rule_U1_ch00`, `rule_U2_v1_loo_only`, `rule_U2_single_sample` and
  `rule_U3_tier4`, and the classes `superseded_phase7_class_v1` and `superseded_phase7_class_v2`, stay in the matrix as
  history. They are not applied.

### 2.5 Classification under the revised D-045

| Theme | v1 class | v2 class | Evidence status | Robustness descriptor |
|---|---|---|---|---|
| Ribosome / translation | UNSUPPORTED | UNSUPPORTED | EXPLORATORY | strongly sample-sensitive |
| Photosynthesis / plastid | EXPLORATORY | EXPLORATORY | EXPLORATORY | mixed robustness |
| ATP synthase / proton transport | UNSUPPORTED | UNSUPPORTED | **UNSUPPORTED** (explicit review) | technically sensitive |
| Defence response | UNSUPPORTED | EXPLORATORY | EXPLORATORY | moderately sample-sensitive |
| Hormone / ethylene | EXPLORATORY | EXPLORATORY | EXPLORATORY | relatively stable |
| Endoplasmic reticulum | UNSUPPORTED | UNSUPPORTED | EXPLORATORY | strongly sample-sensitive |
| Heme / oxidoreductase | EXPLORATORY | EXPLORATORY | EXPLORATORY | moderately sample-sensitive |
| Glutathione | EXPLORATORY | EXPLORATORY | EXPLORATORY | moderately sample-sensitive |
| Water / abiotic stress | EXPLORATORY | EXPLORATORY | EXPLORATORY | relatively stable |
| Cell wall / carbohydrate | UNSUPPORTED | EXPLORATORY | EXPLORATORY | moderately sample-sensitive |

The revision changes the role of LOO, not the data, and it does not rescue any pathway. No theme gained inferential support:
ribosome / translation and ER move from UNSUPPORTED to EXPLORATORY only because sample dependence is now described rather
than used as a rejection threshold.

## 3. Integrated evidence matrix

Abbreviations:
- **+** / **−**: sets ranked higher / lower in STV-infected leaves.
- **LOO max**: largest share of the theme's sets that depend on one sample.
- **ch00**: sets ch00-sensitive in the pre-specified / post-hoc runs.
- **org. LE**: share of leading-edge gene slots that belong to the organelle-identical feature class (corrected rule,
  D-044).
- **Rep. b/p/a**: broadly / partially / ambiguous replicate classes.

Common to all ten themes:
- No representative leading-edge gene is individually significant.
- No set is significant under CAMERA.
- CAMERA direction agrees with fgsea for every set.

| Theme | Sets (BP/MF/CC) | Direction | NES range | Min fgsea padj | CAMERA min FDR | LE genes (min gene padj) | Rep. b/p/a | LOO max (sample); free_rep1 | Cook's not retained | ch00 pre/post; org. LE | Evidence status | Robustness descriptor |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Ribosome / translation | 34 (15/8/11) | + 34 | 1.41 – 2.92 | 5.0 × 10⁻³⁰ | 0.997 | 569 (1.00) | 0/32/2 | **76 % (STV_infected_rep1)**; 6 | 0 | 0/1; 0.2 % | EXPLORATORY | strongly sample-sensitive |
| Photosynthesis / plastid | 25 (8/0/17) | + 25 | 1.59 – 2.44 | 4.2 × 10⁻¹¹ | 0.997 | 290 (1.00) | 15/10/0 | 32 % (STV_free_rep2; also 32 % STV_infected_rep1); 1 | 1 | 0/0; 1.1 % | EXPLORATORY | mixed robustness |
| ATP synthase / proton transport | 5 (2/1/2) | + 5 | 1.61 – 1.88 | 3.3 × 10⁻³ | 0.997 | 60 (1.00) | 5/0/0 | 60 % (STV_free_rep2); 0 | 0 | **4/4; 28.6 %** | **UNSUPPORTED** (explicit review, §2.3) | technically sensitive |
| Defence response | 8 (8/0/0) | − 8 | −1.99 – −1.75 | 1.2 × 10⁻⁸ | 0.999 | 150 (0.193) | 7/1/0 | 50 % (STV_free_rep2; also 50 % STV_infected_rep1); 1 | 0 | 0/0; 0 | EXPLORATORY | moderately sample-sensitive |
| Hormone / ethylene | 5 (5/0/0) | − 5 | −2.02 – −1.69 | 4.0 × 10⁻³ | 0.999 | 45 (1.00) | 4/1/0 | 20 %; 1 | 0 | 1/0; 0 | EXPLORATORY (1 Tier 4 set) | relatively stable |
| Endoplasmic reticulum | 4 (0/0/4) | − 4 | −1.57 – −1.49 | 1.9 × 10⁻⁴ | 0.997 | 188 (1.00) | 1/3/0 | **100 % (STV_free_rep1)**; 4 | 0 | 0/0; 0 | EXPLORATORY (explicit review, §2.3) | strongly sample-sensitive |
| Heme / oxidoreductase | 11 (2/8/1) | + 5 / − 6 | −1.83 – 1.95 | 3.0 × 10⁻⁴ | 0.997 | 196 (0.297) | 9/2/0 | 36 % (STV_free_rep1); 4 | 0 | 1/0; 0.7 % | EXPLORATORY | moderately sample-sensitive |
| Glutathione | 3 (2/1/0) | − 3 | −1.95 – −1.60 | 1.4 × 10⁻³ | 0.998 | 56 (1.00) | 3/0/0 | 33 %; 1 | 1 | 0/0; 0 | EXPLORATORY | moderately sample-sensitive |
| Water / abiotic stress | 3 (3/0/0) | − 3 | −2.17 – −1.62 | 1.6 × 10⁻⁶ | 0.999 | 124 (1.00) | 3/0/0 | 33 %; 1 | 0 | 0/0; 0 | EXPLORATORY | relatively stable |
| Cell wall / carbohydrate | 4 (2/2/0) | + 2 / − 2 | −1.87 – 1.63 | 1.7 × 10⁻⁴ | 0.998 | 136 (0.297) | 4/0/0 | 50 % (STV_free_rep2); 0 | 1 | 1/1; 0 | EXPLORATORY | moderately sample-sensitive |

The full matrix (`theme_evidence_matrix.tsv`) also contains:
- all GO IDs, the top three terms, fgsea medians and CAMERA up/down counts;
- the descriptive count of sets at FDR < 0.05 in the secondary CAMERA run with `inter.gene.cor = 0.01`. This run is not
  used for any class (D-034);
- the LOO dependence by sample, `evidence_status_basis`, `robustness_basis` and `primary_caveat`;
- the superseded D-045 flags and classes (§2.4);
- the five representative leading-edge genes, the STV and tomato–virus literature, the interpretation and the caveats.

Representative genes are chosen by a fixed rule: most theme sets in whose leading edge the gene appears, then |Wald
statistic|.

Set-level detail for all 215 sets: `theme_set_level_evidence.tsv`.

The 113 fgsea-significant sets that match no theme regex ("Other") are not turned into new themes. They are summarized as
descriptive GO modules in §6.

## 4. Theme by theme

### 4.1 Ribosome / translation (Task 9): EXPLORATORY, strongly sample-sensitive

- **What the sets show.** 34 sets rank higher in STV-infected leaves. Composition:
  - 12 ribosome-structure sets (cytosolic and organellar);
  - 8 ribosome-biogenesis / rRNA-processing sets;
  - 9 translation / peptide-biosynthesis sets;
  - 4 tRNA / RNA-modification sets;
  - 1 other.

  This is the strongest fgsea signal: the "ribosome" set has padj 5 × 10⁻³⁰.
- **Leading-edge genes (569).**
  - 166 cytosolic ribosomal proteins (60S/40S);
  - 55 organellar ribosomal proteins (50S/30S; plastid or mitochondrial). One is in the organelle-identical feature class.
  - 121 translation-factor / tRNA-related genes and 51 biogenesis factors.
  - 176 genes outside these description classes (`ribosome_le_gene_composition.tsv`).

  The top genes are ribosomal proteins with Wald statistics of 0.4 – 1.1: many small shifts. No gene is close to
  significance (all gene padj = 1.00).
- **Why the evidence is weak.**
  - Inferential: CAMERA FDR ≥ 0.997 in all 34 sets.
  - These sets carry the strongest inter-gene correlation of all themes (median estimated correlation 0.43; median VIF
    36.6). This is exactly the situation in which a gene-permutation null overstates significance.
  - Robustness: **no** set is broadly replicate-consistent (0/34; 32 partially, 2 ambiguous).
  - 26 of 34 sets depend on **STV_infected_rep1** in LOO (29/34 lost when it is dropped; 33/34 retained when
    STV_infected_rep2 is dropped). Phase 6B documented the dependence at collection level (CC 47.7 %). Phase 7 shows that
    it is concentrated in this theme.
  - Under the revised D-045 this concentration is described (strongly sample-sensitive) and no longer used as a rejection
    threshold. The tendency stays exploratory, and it is the least robust of the larger themes.
- **Biological plausibility.**
  - A higher transcript rank of ribosomal genes is not evidence of more ribosomes or more protein synthesis. STV is not
    shown to affect host translation.
  - In the acute-virus literature, translation-related transcripts change in both directions. Examples: TSWV in tomato
    lowered protein-biosynthesis genes (Prigigallo et al. 2025). The NIK1–LIMYB immune pathway represses ribosomal-protein
    genes in Arabidopsis (Zorzatto et al. 2015).
  - Neither study concerns a persistent virus, and neither is comparable with this dataset.

### 4.2 Photosynthesis / plastid (Task 10): EXPLORATORY, mixed robustness

- **Pattern.** 25 sets rank higher in STV-infected leaves:
  - plastid envelope, thylakoid and photosynthetic membrane;
  - light-harvesting genes (chlorophyll a-b binding proteins dominate the leading edges).

  15 of 25 are broadly replicate-consistent.
- **ch00 caveat.**
  - No set is ch00-sensitive under either exclusion rule.
  - Organelle-identical features occupy 1.1 % of the theme's leading-edge slots (maximum 6.1 % in any set).
  - 19 of 25 sets contain none (`photosynthesis_sets_organelle_check.tsv`).
  - The broader photosynthesis / plastid tendency therefore **persists conceptually after the ch00 caveat**. It does not
    depend on organelle-identical gene models, unlike ATP synthase.
- **Limits.**
  - CAMERA FDR ≥ 0.997.
  - The core "photosynthesis" module (P07, 7 sets) depends on STV_infected_rep1 in LOO (7/7).
  - The plastid-envelope module (P03, 12 sets) is broadly consistent (12/12), with dependence mostly on STV_free_rep2
    (4/12).
- **Literature.** Acute tomato viruses lower photosynthesis-gene expression (PepMV, TYLCTHV, TSWV) and damage chloroplasts
  (ToCV). This dataset shows the opposite rank direction, which is not comparable. STV causes no leaf ultrastructural
  change (Elvira-González et al. 2020).
- **Verdict.** An exploratory, directionally consistent tendency with mixed robustness: the plastid-envelope component is
  broadly replicate-consistent, while the core photosynthesis component depends on STV_infected_rep1. It is not a
  demonstrated change in photosynthesis.

### 4.3 ATP synthase / proton transport (Task 11): UNSUPPORTED, technically sensitive

- 5 sets rank higher in STV-infected leaves.
- Organelle-identical features fill **28.6 %** of their leading-edge slots, and all 5 sets contain them.
- 4 of 5 sets lose significance without them, under both exclusion rules.
- The largest positive GO module, purine-ribonucleotide biosynthesis (P01, 19 sets, mostly "Other"), is an extension of
  this caveat: 17/19 sets are ch00-sensitive and all 19 contain organelle-identical leading-edge genes.
- **This theme is sensitive to reference / feature composition** (D-038, D-044). It is not used as biological evidence for
  STV. "Organelle-identical feature class" describes sequence identity only. It does not establish organellar origin.
- The explicit D-045 review (§2.3) shows that the apparent pattern is carried primarily by this feature class: without it,
  none of the four ATP-synthase-specific sets remains fgsea-significant (NES falls by 21 – 37 %).

### 4.4 Defence and stress cluster (Task 12)

| Theme | Direction | Rep. broadly | LOO max | Representative leading-edge genes | Status; robustness |
|---|---|---|---|---|---|
| Defence response | − 8/8 | 7/8 | 50 % STV_free_rep2 (and 50 % STV_infected_rep1) | allene oxide cyclase, NB-ARC, EIX receptor 1, acidic chitinase, NIM1-like 2 | EXPLORATORY; moderately sample-sensitive |
| Water / abiotic stress | − 3/3 | 3/3 | 33 % | catalase, WRKY6, WRKY39, WRKY42, WRKY33 | EXPLORATORY; relatively stable |
| Glutathione | − 3/3 | 3/3 | 33 % | glutathione S-transferases | EXPLORATORY; moderately sample-sensitive |
| Heme / oxidoreductase | 5 up / 6 down | 9/11 | 36 % STV_free_rep1 | cytochrome P450s, allene oxide synthase | EXPLORATORY; moderately sample-sensitive |
| Cell wall / carbohydrate | 2 up / 2 down | 4/4 | 50 % STV_free_rep2 | wall-associated receptor kinases, LRK10-like | EXPLORATORY; moderately sample-sensitive |

- **Not independent signals.** Leading-edge Jaccard between themes (`theme_le_overlap.tsv`):
  - glutathione–water 0.17 (26 shared genes; overlap coefficient 0.46);
  - glutathione–heme 0.13;
  - heme–water 0.12;
  - defence–hormone 0.11 (20 shared genes; 44 % of the hormone leading edge).

  The negative GO modules mix these themes. For example, module N01 "peroxidase activity" (8 sets) spans the water,
  heme, glutathione and "Other" themes.
- **Reading.** These are one partly shared, lower-ranked stress/defence-related tendency spread over overlapping GO terms,
  not five findings.
- **Direction relative to other studies.**
  - It matches the non-significant PR-locus trend in Badeli et al. (preprint).
  - It is opposite to the defence induction reported for acute tomato viruses (ToBRFV, TYLCTHV, TSWV, TYLCV, PepMV),
    which are not comparable.
- **Limits.**
  - No leading-edge gene is significant (minimum padj 0.193, Solyc12g009240.1).
  - The cell-wall and heme themes have no coherent direction.

### 4.5 Hormone / ethylene (Task 8): EXPLORATORY, relatively stable

See §5.

### 4.6 Endoplasmic reticulum (Task 13): EXPLORATORY, strongly sample-sensitive

- 4 CC sets rank lower in STV-infected leaves: ER, ER subcompartment, ER membrane, and the nuclear outer membrane–ER
  network.
- **Single-sample dependence.**
  - All 4 are *lost* when STV_free_rep1 is dropped and *retained* when STV_free_rep2 or STV_free_rep3 is dropped.
  - In the replicate diagnostic, all 9 infected–free pairs are ordered as expected, but removing STV_free_rep1 halves the
    infected–free difference (r₁ = 0.48 – 0.50). Only 1 of 4 sets is broadly consistent.
- **Leading edge.** The top genes are sterol and phospholipid biosynthesis enzymes (Δ7-sterol-C5-desaturase,
  emopamil-binding protein, Δ14-sterol reductase, phosphatidylserine synthase), not folding chaperones. "ER / protein
  folding" overstates the content.
- **Verdict.** STV_free_rep1 is the documented outlier sample (Phase 5, D-022). It is kept in every primary analysis; this
  classification removes nothing. The first Phase 7 rules (v1 and v2) classed the theme UNSUPPORTED. Under the revised
  D-045 it is **EXPLORATORY and strongly sample-sensitive** (explicit review in §2.3): the group-level direction is
  consistent across all sample pairs, but its size and fgsea support depend strongly on STV_free_rep1. It is the weakest of
  the exploratory themes.

### 4.7 Remaining themes

Heme / oxidoreductase, glutathione, water / abiotic stress and cell wall / carbohydrate: see §4.4.

## 5. Ethylene case study (Task 8)

Outputs:
- `results/biological_integration/ethylene_case_study.tsv` (120 tested genes);
- `ethylene_go_sets.tsv`;
- `ethylene_prior_study_loci.tsv`;
- figure F4.

Gene selection (fixed rule; no threshold):
- tested genes annotated (propagated) to GO:0009693, GO:0009873, GO:0071369 or GO:0009723;
- tested genes whose ITAG4.0 description names ACC synthase / oxidase, an ethylene receptor, CTR1, EIN2/EIN3/EIL/EBF, or
  "ethylene-responsive";
- loci named by Fukuhara et al. / Badeli et al.

**Gene level (statistical evidence).**
- **No ethylene-related gene has padj < 0.05.**
- 11 ACS/ACO genes, 15 perception / signalling genes, 54 "ethylene-responsive" (mostly ERF) genes and 39 GO-annotated-only
  genes were tested. 74 of the 120 have a negative statistic.
- **ACO1 (Solyc07g049530.3).**
  - MLE log2FC −1.81, apeglm −1.56, Wald −4.06 (rank 18,236 of 18,238), padj **0.297**.
  - Lower in all 9 infected–free sample pairs (blind VST). Free: 9.77 / 9.66 / 8.53; infected: 7.52 / 7.48 / 8.24.
  - Removing STV_free_rep1 keeps 86 % of the difference.
  - STV_free_rep3 is the lowest free sample, close to STV_infected_rep3.
  - In PLAZA 5.0 it is annotated only to the untested "ethylene biosynthetic process" set, so it cannot be in the leading
    edge of the three tested ethylene sets. It is in 8 other significant leading edges.
- **Solyc12g009240.1.**
  - ITAG4.0: "Ethylene-responsive transcription factor"; Badeli et al.: "PR-protein". It is not annotated to any ethylene
    GO set.
  - Most negative statistic of the ranking (−4.41), padj **0.193**; lower in all 9 pairs.
- **Receptors / CTR1 / EIN3-EIL / EBF.** Small statistics (|Wald| ≤ 0.9), all padj 1.00.
- **Not tested.** ERF017 (Solyc07g007770) and NAC1 (Solyc01g095080), named by Badeli et al., were removed by the count
  filter.

**Pathway level (exploratory directional evidence).**
- "ethylene-activated signaling pathway" NES −2.01, padj 0.004.
- "cellular response to ethylene stimulus" NES −2.02, padj 0.004.
- "response to ethylene" NES −1.69, padj 0.041.
- All three: Tier 3, broadly replicate-consistent, CAMERA FDR 0.999.
- "ethylene biosynthetic process" (GO:0009693) has 4 genes in the universe and was not tested.

**Comparison.**
- Fukuhara et al. (abstract) report down-regulation of ethylene biosynthesis and signalling genes.
- Badeli et al. (preprint) report ACO1 lower but padj 0.225, and no DEG.
- This reanalysis agrees on **direction** and on the **absence of gene-level significance** (as in Badeli).

**Conclusion.** Ethylene-related transcripts and GO sets show a consistent **lower-rank tendency** in STV-infected leaves.
This is pathway-level, exploratory, directional evidence, and it did not survive correlation-aware validation. Ethylene
biosynthesis or signalling is **not** shown to be significantly altered, and no inference on ethylene levels is possible
(no hormone measurement).

## 6. GO redundancy (Task 14)

Method (descriptive; `scripts/R/21_go_redundancy_modules.R`):
- Jaccard of full set membership (ranked universe) and of fgsea leading edges for all 23,005 pairs of the 215 sets.
- GO ancestor/descendant pairs (is_a + part_of; the GO release used in Phase 6).
- Modules from average-linkage clustering of 1 − leading-edge Jaccard, separately for positive and negative NES, cut at
  0.7. The cut was fixed before inspecting module content.

| Measure | Value |
|---|---|
| Pairs with full-membership Jaccard ≥ 0.5 | 329 (16 identical) |
| Pairs with leading-edge Jaccard ≥ 0.5 | 353 |
| Ancestor/descendant pairs | 369, involving 179 of 215 sets |
| Modules at cut 0.7 | 73 (44 positive, 29 negative; 34 singletons) |
| Modules at cut 0.5 / 0.9 (sensitivity) | 103 / 44 |

Largest modules (`go_module_summary.tsv`, figure F5b):
- **Positive.**
  - P01 purine-ribonucleotide biosynthesis (19 sets; organelle-identical / ch00-sensitive, §4.3);
  - P02 purine-ribonucleotide metabolism (14);
  - P03 plastid envelope (12);
  - P04 pigment biosynthesis (7);
  - P05 translation-factor activity (7);
  - P06 ribosome (7);
  - P07 photosynthesis (7);
  - P08 hexose metabolism (6).
- **Negative.**
  - N01 peroxidase activity (8);
  - N02 carbohydrate binding (6);
  - N03 response to external biotic stimulus (6);
  - N04 lipid modification (6);
  - N05 heme binding (5);
  - N06 response to acid chemical (5);
  - N07 endoplasmic reticulum (4);
  - N08 ethylene-activated signalling (4).

**Reading.** The 215 significant sets collapse to about 40 – 100 redundancy groups, depending on the cut. The ten themes
are not ten independent findings, and the number of significant GO terms says nothing about the size of the biological
effect. No module is significant under CAMERA.

## 7. Figures (Task 15; `figures/biological_integration/`)

| Figure | Content | Note |
|---|---|---|
| F1 `F1_integrated_evidence_overview.png` | Ten themes × evidence layers (fgsea, CAMERA, replicate, LOO, Cook's, ch00, literature), then evidence status and robustness descriptor in separate columns | CAMERA cells are grey ("not significant"); fgsea cells use a neutral colour; the robustness palette is separate from the status palette. Regenerated for the D-045 revision |
| F2 `F2_fgsea_vs_camera.png` | −log10 fgsea padj vs −log10 CAMERA FDR on identical axes | All points lie on the floor (CAMERA FDR ≥ 0.997); direction agreement 215/215 |
| F3 `F3_theme_summary.png` | Direction, evidence status, robustness descriptor and main caveat (`primary_caveat`) per theme | Directions are rank-level shifts. Regenerated for the D-045 revision |
| F4 `F4_ethylene_case_study.png` | Gene-level Wald statistics vs pathway-level NES | Shows the gene-level null beside the pathway tendency. Point jitter and `ggrepel` label placement use the project seed `20261025` (D-050); two regenerations from identical inputs gave the same SHA-256 (`run_log.md` entry 139) |
| F5 `F5_go_redundancy_heatmap.png`, F5b `F5b_go_modules_overview.png` | Leading-edge Jaccard heatmap; module composition | Descriptive only |

## 8. Outputs and reproducibility

| Script | Outputs (`results/biological_integration/`) |
|---|---|
| `scripts/R/20_biological_integration.R` | `theme_evidence_matrix.tsv`, `theme_set_level_evidence.tsv`, `theme_representative_genes.tsv`, `theme_le_overlap.tsv`, `ethylene_case_study.tsv`, `ethylene_go_sets.tsv`, `ethylene_prior_study_loci.tsv`, `ribosome_sets_by_subtheme.tsv`, `ribosome_le_gene_composition.tsv`, `photosynthesis_sets_organelle_check.tsv`, `sessionInfo_20.txt` |
| `scripts/R/21_go_redundancy_modules.R` | `go_redundancy_pairs.tsv`, `go_redundancy_summary.tsv`, `go_redundancy_set_modules.tsv`, `go_module_summary.tsv`, F5, F5b, `sessionInfo_21.txt` |
| `scripts/R/22_phase7_figures.R` | F1 – F4, `sessionInfo_22.txt` |
| `scripts/bash/literature_search_europepmc.sh` | `literature/europepmc_*.tsv` (metadata only) |
| manual (documented) | `curated_theme_literature.tsv`: literature summaries, interpretation and caveats per theme |

Run order: 20 → 21 → 22 (env `stv-r`). Inputs are frozen Phase 4 – 6B and D-044 outputs. Git-ignored inputs:
- `go_gene_sets.rds` and `deseq2_objects.rds` (regenerable; SHA-256 in `reproducibility_audit.md`);
- the GO OBO file (SHA-256 in `annotation_sources_manifest.tsv`).

D-045 revision (2026-10-26): scripts 20 and 22 were re-run on the unchanged inputs. Only `theme_evidence_matrix.tsv`, F1
and F3 changed. All other outputs of script 20, F2 and every quantitative column of the matrix are byte-identical (`run_log.md`
entry 136).
