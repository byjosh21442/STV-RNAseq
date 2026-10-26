# Organelle-aware re-quantification — assessment (Phase 6B closure)

Status: **assessment only.** No new Salmon index was built, no FASTQ was quantified, and tximport, DESeq2, fgsea and CAMERA
were not re-run. The primary analysis (Phases 3 – 6B) is frozen.

## 1. Pre-analysis question and decision criteria (frozen 2026-10-26T15:17:12Z, before the detailed assessment)

This section was written before the SGN/ITAG4.0 provenance research, the Salmon reference reconstruction, the method
literature search and the corrected BLAST diagnostic. Its SHA-256 snapshot is in
`results/organelle_assessment/prespecification/`. Later sections do not change it.

**Disclosure.** Some quantities were already known when this section was written: the Phase 6B ch00 sensitivity results
(`functional_robustness_plan.md` §11). The criteria below concern reference correctness and the dependence of conclusions.
None of them depends on obtaining more DEGs or more significant pathways.

### 1.1 Question

> Does the current SL4.0/ITAG4.0 reference and Salmon decoy-aware quantification create ambiguity between nuclear/unplaced
> ch00 loci and chloroplast/mitochondrial transcripts large enough to justify a secondary, organelle-aware quantification
> sensitivity analysis?

The biological target of the project remains the **nuclear leaf transcriptome response** (STV-infected vs STV-free). An
organelle-aware branch, if ever run, would be a **sensitivity branch** and would never silently replace the primary analysis.

### 1.2 Evidence dimensions (assessed only from documentation and existing outputs)

| ID | Dimension | What is assessed | Source type |
|---|---|---|---|
| E1 | Reference completeness | Whether the index used for the primary quantification contains the tomato plastome and mitogenome, as quantification targets or as decoys | Index build commands, logs, reference manifests |
| E2 | Sequence ambiguity | How many ITAG4.0 gene models are (near-)identical to organelle sequence (criteria §1.4), and how many of these are expressed, i.e. in the 18,277 tested genes | BLAST diagnostics (existing + corrected) |
| E3 | Annotation provenance | What SGN/ITAG4.0 documentation says about ch00, organelle sequence and organellar genes; whether origin (organellar contig, NUPT/NUMT, contamination) can be established | Official documentation, primary assembly papers |
| E4 | Effect on conclusions | (a) any gene-level result; (b) the share of Phase 6 significant sets that depend on organelle-identical gene models; (c) whether any **theme classification** or any planned interpretation depends on them | Existing Phase 6/6B tables |
| E5 | Method guidance | Whether official Salmon / tximport / nf-core documentation or primary methods papers say how organelle genomes should be represented (targets, decoys or neither) when the biological target is nuclear | Official documentation, primary literature |
| E6 | Reproducibility and cost | Whether a branch could be built from versioned public accessions with tracked code and kept separate from the primary results, at a feasible compute and storage cost | Accessions, existing scripts |

### 1.3 Decision rules (applied in this order)

- **STOP (report to owner, no recommendation)** if any of the following holds:
  - official documentation shows the current architecture is incorrect, not merely incomplete (STOP 6 of the instruction);
  - organelle-identical gene models carry a gene-level result: a gene with padj < 0.05, of which there are none;
  - the decision could only be made by running the new quantification (STOP 11).
- **C — strongly recommended before final interpretation** if organelle sequence is absent from the index (E1) **and** at
  least one of:
  1. organelle-identical gene models make ≥ 50 % of the significant sets of any collection sensitive. This is the Phase 6B
     ch00 STOP threshold, applied to the corrected classification;
  2. a theme that the project intends to interpret is shown to depend on them, for themes not already explained by the
     ch00 caveat;
  3. E5 sources state that organelle sequence **should** be included for nuclear RNA-seq of plant leaves with
     poly(A)-selected libraries.
- **B — optional sensitivity quantification scientifically useful** if organelle sequence is absent from the index (E1)
  **and** expressed organelle-identical gene models exist (E2 ≥ 1 in the tested set) **and** their influence is limited
  (E4 below the C thresholds) **and** a reproducible branch is feasible (E6).
- **A — no new quantification needed** otherwise. That means organelle sequence is already represented, or no expressed
  gene model is organelle-identical, or the residual ambiguity cannot affect any conclusion that the project will draw.

The expected number of significant genes or sets in a re-quantified analysis is **not** a criterion. Neither is its
direction.

### 1.4 Corrected sequence-identity criteria (fixed now, for Task 8)

These criteria are used only if a corrected BLAST diagnostic is run. The thresholds are the ones pre-specified in Phase 6B
(D-038). Only the way coverage is computed changes:
- **Coverage across HSPs.** For each ITAG4.0 transcript and each organelle genome, coverage is the **union of query
  intervals** covered by all HSPs with identity ≥ 98 % (e-value ≤ 1e-20). It is expressed as a fraction of the transcript
  length. Identity is the alignment-length-weighted identity of those HSPs.
- **Near-identical (organelle-type).** Union coverage ≥ 80 % with weighted identity ≥ 98 %, for the organelle genome with
  the larger coverage.
- **Partial organelle-identical.** Union coverage ≥ 50 % and < 80 % at ≥ 98 % identity. It is reported separately and is
  never merged with near-identical.
- **Short-segment rule.** A gene is never called organelle-type because of a single short HSP. The union coverage criterion
  requires ≥ 50 % or ≥ 80 % of the transcript length, so a short HSP alone cannot qualify.
- **Interpretation.** Identity and coverage describe **sequence similarity only**. They do not establish genomic placement,
  annotation correctness or biological origin (organellar transcript vs nuclear insertion).

### 1.5 Scope metrics to report (Task 7; from existing outputs only)

- ch00 genes and transcripts.
- How many are in the tested set (18,277), the ranking (18,238) and any significant leading edge.
- Counts by organelle-similarity class, under the Phase 6B rule, the Phase 6B post-hoc rule and the corrected rule.
- baseMean distribution by class.
- For the photosynthesis/plastid, ATP-synthase/proton-transport and other themes: leading-edge membership of
  organelle-identical genes, and set status under the existing ch00 sensitivity runs.
- Whether any conclusion other than ATP synthase depends on them.

No new threshold is introduced for these descriptive metrics.
