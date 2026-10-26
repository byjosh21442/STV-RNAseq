# Phase 6B decision audit (D-033 – D-044)

Written 2026-10-26 (Phase 6B closure). This is an audit only. No decision was re-run or changed. The recommendation labels
are the agent's technical assessment. Final status is an **owner decision**.

## Labels

| Label | Meaning |
|---|---|
| APPROVE | Sound as implemented. Keep it as the project rule |
| APPROVE WITH CAVEAT | Keep it. Reports must carry the stated caveat |
| REVISE | The rule or its documentation should change before Phase 7 |
| OWNER DECISION REQUIRED | Not a technical question, or a choice between equally defensible options |

## Audit criteria

- **Timing.**
  - *Pre-specified*: the rule is in `functional_robustness_plan.md` §1 – §10, frozen 2026-10-26T07:36:28Z (snapshot in
    `results/functional_robustness/prespecification/`) before any Phase 6B result.
  - *Post-hoc*: written after results were seen.
- **Affects primary inference.** Whether the decision changes the Phase 4 gene-level result or the Phase 6 primary fgsea
  outputs. The Phase 6 outputs are verified byte-identical: 76 / 76 SHA-256.
- **Reproducible.** Whether the outputs can be regenerated from tracked code. Inputs are either tracked or git-ignored with
  a recorded SHA-256.
  - The git-ignored R intermediates are `txi_final.rds` (script 01), `deseq2_objects.rds` (script 02) and
    `go_gene_sets.rds` (script 05). They are read by scripts 10, 11, 13, 14, 15 and 18. Their SHA-256 at closure is recorded in
    [reproducibility_audit.md](reproducibility_audit.md).

## Summary table

| ID | Question | Decision (short) | Timing | Result (before → after / primary vs sensitivity) | Affects primary inference? | Reproducible? | Post-hoc bias risk | External support | Recommendation |
|---|---|---|---|---|---|---|---|---|---|
| D-033 | Keep the 39 Cook's-filtered genes out of the GSEA ranking? | Primary unchanged (18,238); sensitivity with all 18,277 | Pre-specified | 210 / 215 sets stable, 5 supported, 0 lost, ρ(NES) ≥ 0.994 | No | Yes (`10_functional_cooks_sensitivity.R`) | Low: the rule was fixed before results and not chosen by set count | Mixed: Love et al. 2014 favour exclusion in small designs; fgsea leaves the choice open | **APPROVE** |
| D-034 | Which correlation-adjusted test validates fgsea? | CAMERA on voom/TMM; primary `inter.gene.cor = NA`; 0.01 descriptive | Pre-specified | Primary: 0 sets at FDR < 0.05 (all collections). Descriptive 0.01: BP 40, MF 7, CC 26 | No (Phase 6 outputs unchanged), but it determines the tiers and themes | Yes (`11_functional_camera.R`; inputs: `quant.sf`, SHA-256 checked) | Low: the choice was frozen before any CAMERA output | Strong: limma 3.66.0 `?camera` says that with `NA` "camera gives rigorous error rate control for all sample sizes and all gene sets", and that with 0.01 it "will rank biologically interpretable sets more highly"; Wu & Smyth 2012 | **APPROVE WITH CAVEAT** |
| D-035 | How to grade agreement between fgsea and CAMERA? | Tiers 1 – 4, 3C | Pre-specified | Tier 1 0, Tier 2 0, Tier 3 210, Tier 4 5 | No | Yes (`12_functional_method_comparison.R`) | Low | GSEA User Guide (FDR 0.25 exploratory); otherwise project-defined | **APPROVE WITH CAVEAT** |
| D-036 | Are pathway signals carried by all replicates? | P (pairs in direction) and r₁ (dependence on STV_free_rep1); descriptive classes | Pre-specified | 144 broadly, 69 partially, 2 ambiguous, 0 rep1-driven | No | Yes (`13_functional_replicate_consistency.R`; blind VST from regenerable `deseq2_objects.rds`) | Low | None specific; thresholds are project-defined | **APPROVE WITH CAVEAT** |
| D-037 | Do calls depend on any single sample? | 6 five-sample refits, sensitivity only | Pre-specified | Max sample-dependent share 47.7 % (CC, STV_infected_rep1); STV_free_rep1 13.3 / 37.3 / 11.4 % | No; explicitly never used to drop a sample | Yes (`14_functional_loo_sensitivity.R`; `txi_final.rds` regenerable by script 01) | Low | Principle-level only (Cook 1977; DESeq2 Cook's; Schurch et al. 2016). No pathway-level precedent found | **APPROVE WITH CAVEAT** |
| D-038 | Are ch00 leading-edge genes organelle copies, and do results depend on them? | Single-HSP BLAST classes; exclusion sensitivity | Pre-specified | 186 / 513 near-identical; ch00-sensitive 12.5 / 9.8 / 9.1 %; ATP-synthase sets depend | No | Yes (`ch00_organelle_blast.sh`, `15_*.R`; BLAST+ from an external env, version recorded) | Low for the rule. **Known limitation:** `-max_hsps 1` under-calls coverage. The corrected rule gives 334 near-identical (D-044). The wording "no organelle genome in the index" was inaccurate and is now corrected | SGN release notes; BLAST+ manual | **APPROVE WITH CAVEAT** |
| D-039 | How to group sets into themes and label confidence? | Fixed regexes over GO names; ROBUST / SUPPORTED WITH CAUTION / EXPLORATORY | Pre-specified (before CAMERA) | All 10 themes EXPLORATORY; 113 sets in "Other" | No | Yes (`16_functional_robustness_figures_themes.R`) | Low. Disclosed: the "Hormone" regex contains *ethylene*, a word seen in Phase 6 results and in the Fukuhara et al. 2020 abstract | None; transparent rule-based grouping | **APPROVE WITH CAVEAT** |
| D-040 | Which figures? | F1 – F5 only, selection rules fixed | Pre-specified | 5 figure groups | No | Yes (`16_*.R`) | Low | — | **APPROVE** |
| D-041 | Leading-edge summary computed by an untracked helper | Reproduce with tracked `09_functional_leading_edge_summary.R` | Procedural (incident) | Exact reproduction (2,858 genes, all summary values equal) | No | Yes | None: the result was reproduced, not changed | — | **APPROVE** |
| D-042 | Track the 1.95 MB voom gene table? | Git-ignored; regenerated by script 11 | **Post-hoc** (after results) | — | No | Yes (script 11), given the SHA-verified `quant.sf` | Low (storage policy, not science) | — | **OWNER DECISION REQUIRED** (Task 2 below) |
| D-043 | (1) Prior-study loci in our results? (2) Does the `-max_hsps 1` limitation change the ch00 STOP? | Scripts 17 and 18, labelled post-hoc | **Post-hoc** | (1) Descriptive lookup; (2) max ch00-sensitive 15.8 % (BP), no STOP | No; no class changed | Yes (`17_*.R`, `18_*.R`) | Medium in principle (written after results); mitigated by labelling and by no class change | — | **APPROVE WITH CAVEAT** |
| D-044 (new) | Is an organelle-aware re-quantification needed? | Assessment only; recommendation A | Criteria pre-specified (frozen §1, 15:18:09Z); findings post-hoc | Organelle genomes already inside ch00 (decoy + targets); max dependence 15.8 % | No | Yes (`organelle_multihsp_blast.sh`, `19_*.R`, `ch00_ncbi_provenance.sh`) | Low for the decision rule (frozen); findings descriptive | Salmon docs; Srivastava et al. 2020; nf-core/rnaseq 3.27.0 usage | **OWNER DECISION REQUIRED** |

## Rationale and caveats per decision

**D-033 — APPROVE.**
- The primary ranking is unchanged and the alternative is tested openly.
- The result shows that the choice is immaterial: 0 sets lost.

**D-034 — APPROVE WITH CAVEAT.**
- `inter.gene.cor = NA` is the mode that limma documents as giving rigorous error control. It was frozen before any CAMERA
  output, so switching now would be post-hoc.
- Caveat 1: the "0 CAMERA-significant sets" conclusion depends on this pre-specified choice. The descriptive 0.01 run
  finds 40 / 7 / 26 sets. Reports must state both and must not promote the 0.01 run.
- Caveat 2 (minor): the call uses `allow.neg.cor = FALSE`. limma's pre-3.29.6 estimated mode used `TRUE`. With FALSE,
  negative correlations do not reduce the variance inflation factor, which is slightly more conservative. The median
  estimated correlation in fgsea-significant sets is positive (0.056 – 0.091), so the practical effect is expected to be
  small. This was not re-run.

**D-035 — APPROVE WITH CAVEAT.** The tiers inherit the D-034 choice. "Tier 3" means method-dependent, not false. It must not
be reported as "not enriched".

**D-036 — APPROVE WITH CAVEAT.** The class thresholds (8/9, 6/9, r₁ 0.25 / 0.5) are project-defined and descriptive, with no
external calibration. Do not report them as tests.

**D-037 — APPROVE WITH CAVEAT.**
- Support is principle-level only.
- Five-sample fits leave n = 2 in one group.
- The result (max 47.7 % sample-dependent) is close to the 50 % STOP line and must be reported as such.

**D-038 — APPROVE WITH CAVEAT.**
- The pre-specified class under-calls organelle identity because of `-max_hsps 1`. Under the corrected union-of-HSP rule,
  334 ch00 genes are near-identical, against 186 (`organelle_requantification_assessment.md` §4).
- The statement that the index had "no organelle genome" was inaccurate. ch00 contains the complete plastome and
  mitogenome. Correction notes were added; the original text is kept.
- The STOP assessment and the ATP-synthase caveat are unchanged.
- 7 ranked organelle-identical genes are outside every existing sensitivity run. Three of them are in 2 significant leading
  edges, at low abundance.

**D-039 — APPROVE WITH CAVEAT.**
- Regex themes are coarse: 113 of 215 sets fall in "Other".
- The theme labels inherit D-034 / D-035.
- The ethylene disclosure must stay attached to any statement on hormone signalling.

**D-040, D-041 — APPROVE.** No issue found.

**D-042 — OWNER DECISION REQUIRED.** See Task 2 below.

**D-043 — APPROVE WITH CAVEAT.**
- Both scripts were written after the results and are labelled as post-hoc.
- Neither changed a class.
- Script 18 is now complemented by the corrected classification of D-044, which supersedes it as the most complete
  description of organelle identity. Script 18 still stands as the record of what was checked in Phase 6B.

**D-044 — OWNER DECISION REQUIRED.** Accept recommendation A (no organelle-aware re-quantification)? Also two optional
follow-ups, both owner decisions:
- adopt the corrected "organelle-identical feature class" as a reporting label;
- run a gene-exclusion fgsea sensitivity with the corrected class.

## Task 2 — D-042: the voom gene table

| Question | Answer |
|---|---|
| File | `results/functional_robustness/camera/voom_limma_gene_statistics.tsv`: 18,238 rows × `gene_id, logFC, AveExpr, t, P.Value, adj.P.Val`; 1,952,291 bytes; SHA-256 `a10edfcb0e7080f3c7e20f8627a5491f746a4f344807be5030f6e56692d09893`; gzip -9 = 718,243 bytes |
| Primary or secondary? | **Secondary.** It is the gene-level layer of the CAMERA validation model (D-034). It is not the project's gene-level result, which is DESeq2 |
| Regenerable? | Yes: `scripts/R/11_functional_camera.R` from the six final `quant.sf` (SHA-256 in `results/quantification/quant_sf_sha256.tsv`) and the tracked tx2gene / universe. The `quant.sf` are git-ignored; losing them would require the Salmon run again |
| Needed to audit CAMERA results? | **No.** The set-level results and diagnostics (ρ = 0.963, sign agreement 0.933, TMM, prior df) are in tracked `camera_summary.tsv` and `camera_go_*_all.tsv` |
| Unique information? | Yes, but secondary: per-gene limma logFC / t / P under voom. These are not in any tracked table. They would be needed only for a gene-level limma audit |
| Repository impact | Public repository; > 1 MB. Under D-032, the owner approved only five Phase 6 tables above 1 MB |
| Options | (a) **REGENERATE ON DEMAND** (current state; keep ignored; record the SHA-256 above as the reference); (b) **TRACK** as is (1.95 MB; needs owner approval); (c) track a gzip copy (0.72 MB; format change, needs owner approval) |
| Agent recommendation | **(a) REGENERATE ON DEMAND.** The table is secondary and fully regenerable, and nothing reported depends on it alone. Recording its SHA-256 lets a regenerated copy be verified. **Owner decision.** |

D-042 in `decisions.md` keeps its original text. The recommendation is appended as an update.

## Owner decisions (2026-10-26, after this audit)

- **D-044 approved: Option A**, no new organelle-aware Salmon re-quantification.
- **D-042 resolved: REGENERATE ON DEMAND.** The voom table stays untracked. Script 11 and the reference SHA-256 above are
  kept.
- **No gene-exclusion fgsea** with the corrected organelle-identical class. The existing sensitivity analyses suffice
  for the ATP-synthase caveat.
- The label **"organelle-identical feature class"** may be used when technically appropriate. It describes sequence
  identity only and must not imply a proven biological or genomic origin.
- The APPROVE and APPROVE WITH CAVEAT assessments above stand; their caveats carry into Phase 7.
