# Functional robustness plan — Phase 6B (pre-specified, frozen before any CAMERA result)

Status: **pre-specification**. Sections 1 – 10 were written and frozen (SHA-256 snapshot in
`results/functional_robustness/prespecification/`) before CAMERA, the Cook's-inclusion run, the leave-one-out runs or the
ch00 sensitivity run were executed. Section 11 is appended after the runs and does not change sections 1 – 10.

Scope: Phase 6B **validates** the Phase 6 preranked GSEA (fgsea) results. It does not tell a biological story. The primary
DESeq2 analysis (D-013 – D-019), the six-sample design (D-022), the FDR threshold and the Phase 6 fgsea results stay
unchanged. No DEG list is built, no nominal-p over-representation analysis is run. Decisions D-033 – D-041 in
[decisions.md](decisions.md).

## 1. Two different hypotheses (not statistically identical tests)

| | fgsea (Phase 6, primary) | CAMERA (Phase 6B, validation) |
|---|---|---|
| Input | DESeq2 Wald statistic per gene (negative-binomial GLM, tximport offsets) | voom log-CPM + precision weights, limma linear model |
| Null hypothesis | Genes of the set are randomly distributed along the ranked list (preranked, gene-permutation; weighted KS-type running sum) | The mean moderated statistic of the set equals that of the other genes, after inflating its variance for inter-gene correlation (competitive; Goeman & Bühlmann 2007 terminology) |
| Inter-gene correlation | Ignored (gene-permutation null) | Estimated per set (VIF = 1 + (m − 1) ρ̄) |
| Output | ES, NES, p, BH padj, leading edge | Direction (Up/Down), p, BH FDR, estimated correlation |

Agreement between the two is evidence of robustness to the test and to the expression model. Disagreement is reported as
**method dependence**, not as a falsification of either result.

## 2. Ranking decision D-023 (Cook's distance) — final rule

Evidence and reasoning: [cooks_gsea_decision.md](cooks_gsea_decision.md). Evidence is ambiguous, therefore:

- **Primary (unchanged):** the Phase 6 ranking of 18,238 genes (finite Wald statistic **and** DESeq2 p-value reported;
  the 39 Cook's-flagged genes excluded).
- **Sensitivity (labelled, never primary):** fgsea on all 18,277 genes with a finite Wald statistic (Cook's-flagged genes
  included). Gene sets are rebuilt from the tracked mapping table + the same GO release with the same propagation and
  15 – 500 size rule. Identity check: restricting the rebuilt sets to the 18,238 primary genes must reproduce
  `go_gene_sets.rds` exactly; otherwise STOP (GO mapping change).
- Same fgsea call as Phase 6 (`fgseaMultilevel`, eps 0, nPermSimple 10000, sampleSize 101, scoreType std, gseaParam 1,
  seed 20261025, SerialParam, BH within collection).
- Interpretation (fixed now): a primary-significant set is **Cook's-stable** if padj < 0.05 with the same NES sign in the
  inclusion run. Sets significant only in the inclusion run are listed as sensitivity-only and are never added to the
  primary result. Which sample carries the maximum Cook's distance for each of the 39 genes is reported descriptively.

## 3. CAMERA implementation (Tasks 4 – 5; all choices fixed before any CAMERA output exists)

| Item | Choice | Source |
|---|---|---|
| Counts | `tximport(type = "salmon", tx2gene = reference/ITAG4.0_tx2gene.tsv, countsFromAbundance = "lengthScaledTPM")` from the same six `quant.sf` files | tximport vignette, "limma-voom" section: limma-voom does not use the offset matrix, so scaled counts from abundances are recommended |
| Samples | all six; order and labels from `config/samplesheet.csv` (checked against the DESeq2 colData) | D-022 |
| Genes / universe | the 18,238 primary ranked genes (same as fgsea); no `filterByExpr` re-filtering (the D-014 prefilter is already applied; `filterByExpr` is reported as a diagnostic only) | owner instruction: keep the existing prefilter |
| Normalization | `DGEList` on the 18,238-gene matrix (library size = column sums of retained genes), `calcNormFactors(method = "TMM")` | limma/edgeR docs; RNAseq123 workflow |
| Transformation | `voom(y, design)` (log-CPM + observation-level precision weights) | Law et al. 2014; `?voom` |
| Design | `model.matrix(~ condition)`, factor levels STV_free (reference), STV_infected; **no** batch, pair or surrogate terms | D-013 |
| Contrast | coefficient 2 (`conditionSTV_infected`): positive = higher in STV_infected (same sign convention as the DESeq2 Wald statistic) | D-013 |
| Gene sets | identical lists to Phase 6 (`go_gene_sets.rds`, 15 – 500 genes in the 18,238 universe), indices by `ids2indices` | D-024 – D-026 |
| Test (primary) | `camera(v, index, design, contrast = 2, inter.gene.cor = NA, allow.neg.cor = FALSE, use.ranks = FALSE, trend.var = FALSE, sort = FALSE)` — correlation **estimated per set** | `?camera`: with `inter.gene.cor = NA` "camera gives rigorous error rate control for all sample sizes and all gene sets"; strict mode chosen because this is a validation |
| Test (secondary, descriptive) | identical call with the package default `inter.gene.cor = 0.01` | `?camera`: default 0.01 is "a useful compromise between strict error rate control and interpretable gene set rankings"; RNAseq123 workflow |
| Multiplicity | BH within each collection (one `camera` call per collection), FDR < 0.05 | same as Phase 6 |
| Equivalence diagnostics | same samples/labels/reference/contrast sign; Spearman ρ between the limma moderated t (`eBayes`) and the DESeq2 Wald statistic over the 18,238 genes is reported; ρ ≤ 0 → STOP (model not equivalent) | — |
| Reported per set | NGenes, Direction, PValue, FDR, estimated Correlation, VIF = 1 + (NGenes − 1)·max(Correlation, 0) | `?camera` |

A result of zero CAMERA-significant sets is valid. No CAMERA parameter is changed after results are seen.

## 4. fgsea vs CAMERA comparison and tiers (Task 6; defined before any set identity is examined)

Per collection (BP, MF, CC), over the identical tested sets: significant in fgsea (padj < 0.05), in CAMERA (primary,
FDR < 0.05), in both, fgsea-only, CAMERA-only; direction agreement (fgsea NES > 0 ↔ CAMERA "Up"); overlap fraction
(|both| / |union|) and |both| / |fgsea-significant|; Spearman correlation of the signed significance
sign × (−log10 p) between methods over all tested sets.

Definitions: *concordant* = same direction. *Supportive* = the other method has FDR/padj < 0.25 in the concordant
direction (GSEA User Guide exploratory cut-off). *Opposing* = the other method has FDR/padj < 0.25 in the opposite
direction. *fgsea seed-unstable* = significance call differs between seed 20261025 and the Phase 6 diagnostic seed
20261026 (`gsea_go_*_seed_diagnostic.tsv`).

| Tier | Rule (applied in this order) |
|---|---|
| Tier 4 — discordant/unstable | significant in ≥ 1 method AND (the other method is *opposing* OR the set is significant in both with opposite directions OR [fgsea padj < 0.05, CAMERA neither significant nor supportive, and fgsea seed-unstable]) |
| Tier 1 — robust across methods | FDR < 0.05 in both, concordant, and not fgsea seed-unstable |
| Tier 2 — supported | significant in one method and *supportive* in the other; also Tier-1 sets that are fgsea seed-unstable |
| Tier 3 — fgsea only | fgsea padj < 0.05, CAMERA not supportive and not opposing (method-dependent; **not** called false) |
| Tier 3C — CAMERA only | CAMERA FDR < 0.05, fgsea not supportive and not opposing (extension of Tier 3, method-dependent) |

Sets significant in neither method get no tier. Tier definitions are not changed after pathway names are seen.

STOP ("major CAMERA contradiction"): in any collection ≥ 10 % of the sets significant in both methods have opposite
directions, or the between-method Spearman correlation of signed significance is < 0. Zero CAMERA hits is not a STOP.

## 5. Replicate-consistency diagnostic (Task 7; descriptive, no inferential p-values)

- Data: `vsd_blind` (blind VST, D-018) from `results/deseq2/deseq2_objects.rds`; all six samples; no sample hidden.
- Genes per set: the Phase 6 fgsea leading edge for fgsea-significant sets; the full gene set (labelled) for CAMERA-only
  sets (CAMERA has no leading edge).
- Score: each gene's VST values are z-scored across the six samples; sample score = mean z over the set's genes.
- Quantities: expected sign s (sign of NES; CAMERA direction for CAMERA-only sets); Δ_all = mean(infected) − mean(free);
  P = fraction of the 9 infected–free sample pairs ordered in the expected direction; r₁ = Δ without STV_free_rep1 / Δ_all
  (z-scores kept from all six samples); r for dropping each other sample is reported as a control.
- Classification (fixed): (1) s·Δ_all ≤ 0 → **ambiguous**; (2) else r₁ < 0.25 → **strongly influenced by STV_free_rep1**;
  (3) else P ≥ 8/9 and r₁ ≥ 0.5 → **broadly condition-consistent**; (4) else P ≥ 6/9 → **partially consistent**;
  (5) else **ambiguous**.
- Computed for every set significant in either method. Figures only for the selection in §9.
- STOP ("pathways predominantly driven by one sample"): > 50 % of Tier 1 + Tier 2 sets in a collection classified
  "strongly influenced by STV_free_rep1".

## 6. Leave-one-sample-out influence at the pathway level (Task 8; sensitivity only)

Support: case-deletion influence (Cook 1977), which is exactly what DESeq2's per-gene Cook's distance measures (Love et al.
2014); replicate subsampling to assess result stability in RNA-seq (Schurch et al. 2016, gene level). No published
pathway-level leave-one-out precedent was found, so this is a **labelled sensitivity analysis**, not six new primary
analyses, and it is **never** used to justify removing a sample.

- For each of the six samples: drop it; `DESeqDataSetFromTximport` on `txi_final.rds` (countsFromAbundance "no") with
  the remaining five samples; genes fixed to the 18,277 primary tested genes (no re-filtering); `DESeq(test = "Wald",
  fitType = "parametric", minReplicatesForReplace = 7)`; `results(name = "condition_STV_infected_vs_STV_free")`; ranking =
  finite statistic and non-NA p-value (same rule as primary); sets restricted to that run's ranked genes, 15 – 500; same
  fgsea call and seed.
- With two replicates left in one group, DESeq2 does not flag Cook's outliers for that group (≥ 3 replicates required);
  this is a known limitation of the reduced runs.
- Per primary-significant set and run: *retained* (padj < 0.05, same sign), *supported* (padj < 0.25, same sign), *lost*
  (otherwise). Per run: counts and Spearman ρ of NES vs primary over commonly tested sets.
- *Sample-s-dependent* set: supported in both runs that drop the other two samples of the same condition, but lost when
  s is dropped. STOP if, for any sample, s-dependent sets are > 50 % of the primary-significant sets of a collection.
- CAMERA is not re-run per dropped sample.

## 7. SL4.0ch00 investigation (Task 9)

Facts available before the run: SGN SL4.0 release notes describe chromosome 0 as "152 unplaced contigs"; earlier builds
separated organellar sequences out of the nuclear assembly and removed completely mitochondrial contigs. The SL4.0 FASTA
contains ch00 – ch12 only (no plastome or mitogenome sequence); ch00 = 9,643,250 bp with 513 ITAG4.0 genes (45 ranked).

> **Correction (2026-10-26, after Phase 6B; D-044; the pre-specified text above is kept unchanged).** "No plastome or
> mitogenome sequence" is accurate only for **separate, labelled** records. SL4.0ch00 contains the complete plastome and
> mitogenome sequence at ≥ 98 % identity
> ([organelle_requantification_assessment.md](organelle_requantification_assessment.md) §4.3). The pre-specified rules and
> the results of §11 are not affected.

- Organelle references: tomato plastome **NC_007898.3** (155,461 bp) and mitogenome **NC_035963.1** (446,257 bp), NCBI
  RefSeq, downloaded with SHA-256 recorded (git-ignored under `reference/functional/organelle/`).
- Justification for BLAST: the recurring ch00 leading-edge genes carry plastid-encoded protein names; whether their
  sequence is (near-)identical to the organelle genome determines whether reads from organellar transcripts can be
  assigned to them. BLAST+ from the existing `rachis-qiime2-2026.4` conda env (version recorded).
- blastn of each ranked ch00 gene's ITAG4.0 cDNA against both organelle genomes; blastp of its protein against the
  organelle-encoded proteins (support only).
- Classes (fixed): **plastid-type near-identical** (best blastn hit to NC_007898.3 with identity ≥ 98 % and query
  coverage ≥ 80 %); **mitochondrion-type near-identical** (same thresholds, NC_035963.1); **organelle-similar, diverged**
  (identity ≥ 80 % and coverage ≥ 50 % but below near-identical); **nuclear unanchored** (otherwise). ch00 genes are not
  removed merely because they are unanchored.
- Sensitivity: near-identical organelle-type ch00 genes are removed from the ranking and the gene sets; fgsea re-run
  with the Phase 6 call. A primary-significant set is *ch00-sensitive* if padj ≥ 0.05 or the NES sign flips.
- STOP ("ch00 artifact that materially affects the findings"): > 50 % of the primary-significant sets of any collection
  are ch00-sensitive.

## 8. Themes and final classification (Tasks 10, 14)

Themes are keyword groups over GO term names, taken from the Phase 6 top themes (post-hoc grouping of already-seen
Phase 6 results, fixed here before CAMERA). A set may belong to several themes; unmatched sets are "other".

| Theme | Case-insensitive regular expression on term name |
|---|---|
| Ribosome / translation / RNA processing | `ribosom\|translation\|rRNA\|ribonucleoprotein\|preribosome\|nucleolus\|tRNA` |
| Photosynthesis / plastid / thylakoid | `photosynth\|plastid\|chloroplast\|thylakoid\|photosystem\|chlorophyll\|light.harvesting\|stroma` |
| ATP synthase / proton transport | `ATP synthase\|proton-transporting\|proton transmembrane\|ATP biosynthetic\|ATP metabolic` |
| Defence / biotic response | `defense\|immune\|response to (other organism\|external biotic\|biotic\|bacterium\|fungus\|virus\|oomycete\|insect\|wounding)\|chitin` |
| Hormone signalling | `ethylene\|jasmonic\|salicylic\|abscisic\|hormone\|auxin\|gibberellin\|cytokinin` |
| Water deprivation / abiotic stress | `water deprivation\|osmotic\|salt\|desiccation\|drought\|cold\|heat\|hypoxia\|oxidative stress\|abiotic` |
| Endoplasmic reticulum / protein folding | `endoplasmic reticulum\|Golgi\|unfolded protein\|protein folding` |
| Glutathione / detoxification | `glutathione\|detoxification\|toxin\|xenobiotic` |
| Heme / oxidoreductase | `heme\|tetrapyrrole\|monooxygenase\|oxidoreductase` |
| Carbohydrate binding / cell wall | `carbohydrate binding\|cell wall\|polysaccharide\|glucan\|pectin\|xyloglucan` |

Classification per theme (over its sets significant in ≥ 1 method):

- **ROBUST**: ≥ 1 Tier 1 set; ≥ 50 % of its Tier 1 + 2 sets "broadly condition-consistent"; no Tier 1 + 2 set "strongly
  influenced by STV_free_rep1"; ≤ 25 % of its fgsea-significant sets STV_free_rep1-dependent in the leave-one-out run;
  ≥ 1 Tier 1 set not ch00-sensitive.
- **SUPPORTED WITH CAUTION**: ≥ 1 Tier 1 or Tier 2 set, but a ROBUST criterion fails.
- **EXPLORATORY**: only Tier 3, 3C or 4 sets.

Wording: "higher/lower ranked in STV-infected leaves", "positive/negative enrichment"; never "activation",
"suppression" or causal language. Literature (priority: STV studies > tomato–virus transcriptomics > persistent/
endornavirus-like viruses > tomato defence/photosynthesis/ribosome responses to viruses > broad plant–virus) is used for
plausibility only, with the directness of each source stated.

## 9. Figures (Task 13; selection fixed now)

Only figures that answer: (Q1) robust across methods? (Q2) replicate-consistent? (Q3) sensitive to STV_free_rep1?
(Q4) evidence for/against a theme?

| Figure | Content | Question |
|---|---|---|
| F1 `method_agreement_signed_significance.png` | fgsea vs CAMERA signed −log10 p, all tested sets, coloured by tier, no labels | Q1 |
| F2 `replicate_consistency_<collection>.png` | six-sample scores (free_rep1 marked) for: up to 6 Tier 1 sets per collection ordered by CAMERA FDR, then fgsea padj, then GO ID; plus the 2 strongest fgsea-only sets (Tier 3/4) per collection by fgsea padj, labelled "diagnostic example" | Q2, Q3 |
| F3 `loo_retention.png` | primary-significant sets retained/supported/lost per dropped sample and collection | Q3 |
| F4 `ch00_sensitivity_nes.png` | NES primary vs after excluding organelle-type ch00 genes, coloured by theme "photosynthesis/plastid" and "ATP synthase" | Q4 |
| F5 `cooks_sensitivity_nes.png` | NES primary vs Cook's-inclusion run | Q1 (ranking robustness) |

No other pathway-specific figure.

## 10. STOP conditions (Phase 6B)

The CAMERA model cannot be made equivalent (§3); GO mapping changes (§2 identity check or CAMERA sets ≠ Phase 6 sets);
any primary Phase 6 output changes (`sha256sum -c results/functional_robustness/phase6_outputs_sha256.txt`); the Cook's
decision cannot be sensitivity-tested; a major CAMERA contradiction (§4); pathways predominantly driven by one sample
(§5, §6); a ch00 artifact that materially affects the findings (§7); metadata/sample labels become questionable (sample
order/labels checked in every script); any choice would have to be made after seeing pathway identities. On a STOP:
nothing is committed; the condition is reported.

---

## 11. Outcome (appended 2026-10-26 after all runs; §1 – §10 unchanged, frozen copy sha256 `e7b9b261…`)

All scripts ran after the freeze (2026-10-26T07:36:28Z; every R script re-checks the frozen snapshot with
`fr_check_frozen()`). Logs are in the git-ignored `logs/phase6b_*.log`; every number below is in a tracked table.

### 11.1 STOP assessment

| STOP condition (§10) | Check | Result |
|---|---|---|
| CAMERA model not equivalent | Spearman ρ(limma moderated t, DESeq2 Wald) over 18,238 genes; same samples, labels, reference, sign | ρ = **0.9626**, sign agreement 0.9327 → not triggered |
| GO mapping changed | Rebuilt sets restricted to the 18,238 genes = `go_gene_sets.rds`; CAMERA indices built from the same sets | identical → not triggered |
| Primary Phase 6 outputs changed | `sha256sum -c results/functional_robustness/phase6_outputs_sha256.txt` (76 files), before commit | all OK → not triggered |
| Cook's decision not testable | Inclusion run executed (all 18,277 finite-statistic genes) | run → not triggered |
| Major CAMERA contradiction | ≥ 10 % of jointly significant sets discordant, or signed-significance ρ < 0 | 0 jointly significant sets; ρ = 0.923 / 0.914 / 0.893 (BP/MF/CC) → not triggered |
| Driven by one sample (§5) | > 50 % of Tier 1 + 2 sets "strongly influenced by STV_free_rep1" | no Tier 1 + 2 sets; 0 of 215 sets in that class → not triggered |
| Driven by one sample (§6) | sample-s-dependent sets > 50 % of a collection | maximum **47.7 %** (CC, STV_infected_rep1) → not triggered (close to the threshold; reported) |
| ch00 artifact (§7) | ch00-sensitive sets > 50 % of a collection | 12.5 % / 9.8 % / 9.1 % → not triggered |
| Labels questionable | sample order/labels asserted in scripts 10 – 16 | all checks passed |
| Choice after seeing pathway identities | — | none: tiers, classes, themes and figure rules were applied as frozen; one post-hoc diagnostic (§11.3) is labelled and does not change any class |

### 11.2 Results against the pre-specified rules

- **Helper validation.** The 6B helper reproduces the Phase 6 fgsea results (max |Δp| 5.6 × 10⁻¹⁶, max |ΔNES| 5.3 × 10⁻¹⁵,
  identical leading edges; `helper_validation_phase6.tsv`). Script 09 reproduces the Phase 6 leading-edge summary exactly.
- **Cook's (§2).** Inclusion run: of the 215 primary-significant sets, 210 are Cook's-stable and 5 stay at padj < 0.25 with the
  same sign (BP 4, CC 1); none is lost; 10 sets are significant only in the inclusion run (sensitivity-only, not added).
  Spearman ρ(NES) ≥ 0.994. STV_free_rep1 carries the largest Cook's distance for 19 of the 39 genes.
- **CAMERA (§3).** Primary (`inter.gene.cor = NA`): **0 sets at FDR < 0.05 in BP, MF and CC**; smallest FDR 0.997 – 0.999;
  0 sets at FDR < 0.25. Estimated inter-gene correlation: median 0.046 / 0.039 / 0.077 over all tested sets, 0.056 / 0.065 /
  0.091 over fgsea-significant sets (median VIF 5.6 / 5.4 / 14.5). Secondary (`inter.gene.cor = 0.01`, descriptive): 40 / 7 / 26
  significant sets, 72 of them also fgsea-significant. limma itself: 0 genes at adj. P < 0.05 (descriptive).
- **Tiers (§4).** Tier 1 = 0, Tier 2 = 0, Tier 3C = 0. **Tier 3 = 210** (BP 115, MF 51, CC 44); **Tier 4 = 5** (BP; fgsea-only,
  CAMERA not supportive, seed-unstable). CAMERA direction agrees with the NES sign for 100 % of the 215 fgsea-significant sets
  and for 89 % of all tested sets.
- **Replicate consistency (§5).** 144 sets broadly condition-consistent, 69 partially consistent, 2 ambiguous, 0 strongly
  influenced by STV_free_rep1.
- **Leave-one-out (§6).** Six reduced fits ran (0 material fgsea warnings). Retained + supported per dropped sample ranges from
  34 % (CC, STV_free_rep2) to 99 % (BP, STV_free_rep3). Sample-dependent sets: highest for
  STV_infected_rep1 (CC 47.7 %, MF 37.3 %, BP 26.7 %) and STV_free_rep2 (CC 43.2 %, BP 35.0 %); STV_free_rep1: BP 13.3 %,
  MF 37.3 %, CC 11.4 %.
- **ch00 (§7).** 186 of 513 ch00 genes are near-identical to an organelle genome (71 plastid-type, 115 mitochondrion-type);
  22 of them are ranked and were excluded. ch00-sensitive sets: BP 15, MF 5, CC 4. The 9 lost sets (BP 3, MF 2, CC 4) each
  contain 5 of the 22 excluded genes. Details: [ch00_functional_investigation.md](ch00_functional_investigation.md).
- **Themes (§8).** Every theme is **EXPLORATORY**: no theme has a Tier 1 or Tier 2 set, which the ROBUST and SUPPORTED WITH
  CAUTION classes require. `themes/theme_evidence_classification.tsv`.
- **Figures (§9).** F1 – F5 produced as specified; F2 contains only the "diagnostic example" panels because Tier 1 is empty.

### 11.3 Corrections, deviations and post-hoc additions (none changes a pre-specified rule or class)

1. **Fact in §7:** "45 ranked" ch00 genes is wrong; the tracked ranking contains **47** ch00 genes (checked by GFF join and by
   the `Solyc00g` prefix). The rule does not depend on this number.
2. **Class label:** scripts write "no close organelle similarity" for the §7 class "nuclear unanchored"; the same label is
   used for anchored chromosomes, where "unanchored" would be wrong. Same rule.
3. **Scope of BLAST:** blastn was run on all ITAG4.0 cDNAs (one transcript per gene), not only on ranked ch00 genes, so that
   organelle-type genes on ch01 – ch12 could be counted descriptively. The sensitivity used only ch00 genes, as specified.
4. **BLAST coverage limitation (found in the BLAST output, not from pathway names):** blastn used `-max_hsps 1`, so a spliced
   organelle gene whose cDNA aligns as several HSPs gets a single-HSP coverage < 80 % even at ~100 % identity. 129 of the 136
   ch00 genes classed "organelle-similar, diverged" have identity ≥ 98 %. The frozen rule was kept. A **post-hoc diagnostic**
   (`scripts/R/18_functional_ch00_posthoc_identity_diagnostic.R`, labelled as such) excluded every ranked ch00 gene with identity
   ≥ 98 % and coverage ≥ 50 % (42 genes): ch00-sensitive sets BP 19 (15.8 %), MF 3 (5.9 %), CC 4 (9.1 %). The ch00 STOP
   assessment does not change. The BLAST+ warning "Examining 5 or more matches is recommended" is harmless here (two subject
   sequences).
5. **Post-hoc lookup for Task 11:** `scripts/R/17_prior_study_locus_lookup.R` reports the loci named by the prior analyses in the
   primary results (descriptive; no threshold; not a gene list).
6. **Git:** `camera/voom_limma_gene_statistics.tsv` (1.9 MB) is git-ignored and regenerated by script 11 (D-042).
7. **Reference year:** the frozen Cook's document cites Korotkevich et al. as 2021 (latest bioRxiv version); `methods_references.md`
   §28 cites the first posted version (2016). Same DOI (10.1101/060012). Cook (1977) resolves on Crossref to the canonical DOI
   10.1080/00401706.1977.10489493; 10.2307/1268249 is an alias of the same record.

### 11.4 Reading of the outcome (technical, no biological interpretation)

Under the pre-specified rules, the Phase 6 fgsea signals are **method-dependent**: they are robust to the Cook's ranking choice,
mostly replicate-consistent and partly sensitive to single-sample removal, but none reaches significance under CAMERA with
estimated inter-gene correlation. CAMERA agrees with fgsea on direction and broadly on ranking, which points to a difference in
significance calibration (gene-permutation null vs correlation-adjusted variance) rather than to opposite signals. All themes
remain **EXPLORATORY**: hypotheses for future, adequately replicated experiments, not findings.
