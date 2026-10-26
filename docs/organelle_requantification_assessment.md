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

---

Sections 2 – 11 were written after §1 was frozen (2026-10-26T15:18:09Z). Everything in them is **post-hoc and descriptive**
(D-044). No quantification, tximport, DESeq2, fgsea or CAMERA run was performed.

## 2. Current Salmon reference (reconstruction from code and logs; Task 5)

```text
ITAG4.0_cDNA.fasta (SGN; 34,075 transcripts, one per gene)          S_lycopersicum_chromosomes.4.00.fa.gz (SGN SL4.0;
        |                                                            13 records: SL4.0ch00 + SL4.0ch01-ch12, 782,520,033 bp)
        |                                                                    |                     |
        +---------------- gentrome.fa = cDNA + genome ---------------------+             decoys.txt = the 13 headers
                                         |
             salmon index -t gentrome.fa --decoys decoys.txt -k 31 --ramLimit 6     (Salmon 2.7.0)
                                         |
      33,976 transcript targets (99 exact duplicates removed, keep_duplicates = false) + 13 decoy sequences
                                         |
      salmon quant --libType SR --fldMean 250 --fldSD 25 --gcBias  ->  quant.sf (6 samples)  ->  tximport  ->  DESeq2
```

| Item | Value | Source |
|---|---|---|
| Index build | `workflow/modules/salmon_index.nf`: gentrome = gunzip(ITAG4.0 cDNA) + gunzip(SL4.0 genome); decoys = genome FASTA headers; `ext.args` in `workflow/nextflow.config` | code |
| Index metadata | k = 31; num_refs 33,989; first_decoy_index 33,976; num_decoys 13 (SL4.0ch00 – ch12); keep_duplicates false; Salmon 2.7.0 | `results/genome/index/salmon_index/info.json` (git-ignored) |
| Index log | 99 exact-duplicate transcripts removed (91 Solyc00g); 19 poly-A tails clipped; 11 non-ACGT bases replaced | index log; `duplicate_clusters.tsv` |
| Quantification | `--libType SR --fldMean 250 --fldSD 25 --gcBias`, 6 threads; mapping 87.95 – 90.26 %; decoy fragments 1.26 – 1.91 M per sample (8.45 – 10.57 % of processed) | `meta_info.json`, `docs/salmon_quantification.md` |
| Separate organelle records | **None.** No plastome or mitogenome FASTA record was added, either as target or as decoy | code |
| Organelle sequence inside the index | **Yes, complete** (this assessment, §4.3): SL4.0ch00 (a decoy) contains the full plastome and the full mitogenome at ≥ 98 % identity, and 334 ch00 transcripts (targets) are near-identical to them | §4 |

## 3. Provenance of SL4.0ch00 (Task 4)

| Source (accessed 2026-10-26) | Finding |
|---|---|
| SGN `build_4.00/ReleaseNotes.txt` | "De novo Canu assembly (N50 5.5 Mb) followed by Hi-C scaffolding (12 chromosomes and unplaced contigs)"; "Only 152 unplaced contigs in Chr 00 compared to 4,374 in SL3.0". No statement on organelle screening for SL4.0 |
| Same file, older builds | SL2.40: "4 contigs that were identified as completely mitochondrial have been removed." SL1.50: "Organellar sequences were separated (and are thus not included in this data set)." |
| SGN SL4.0 FASTA (sha256 in `reference/reference_manifest.tsv`) | 13 records; ch00 = 9,643,250 bp with 176 gaps of exactly 100 N → 177 gap-free segments (9,625,650 bp) |
| ITAG4.0 GFF / README | All 513 ch00 genes have source `maker_ITAG` (ordinary gene models; 409 are new Solyc00g5xxxxx IDs, 104 carry older Solyc00g0xx – 3xx IDs). README: 34,075 gene models, 29,281 preserved, 4,794 novel; no organelle statement |
| Hosmani et al. 2019, bioRxiv 10.1101/767764 (abstract, Crossref) | SL4.0 assembled de novo from PacBio reads, Hi-C scaffolded; "almost all scaffolds have been anchored". The full text could not be retrieved; no organelle statement was available |
| NCBI GenBank SL4.0 = **GCA_000188115.5** (Datasets v2 sequence report; `provenance/ncbi_GCA_000188115.5_sequence_report.tsv`) | 12 chromosomes (CM001064.4 – CM001075.4) + **140 unplaced scaffolds** (8,795,813 bp). No organelle molecule is listed |
| Exact-sequence comparison SGN ch00 vs GenBank unplaced scaffolds (`scripts/bash/ch00_ncbi_provenance.sh`) | 164 of 177 ch00 segments are present base-for-base in GenBank (8,354,466 bp). This includes all 4 segments that are ≥ 90 % plastid-identical. For example, ch00 segment 15 (228,710 bp; contains the complete plastome) is identical to **AEKE04000081.1 (Contig195)**. **Absent** from GenBank: 12 of the 17 segments that are ≥ 90 % mitochondrion-identical (832,437 bp; 143 genes, 5 tested), plus segment 1 (438,747 bp; 4 genes, 0 tested; not organelle-like; probably split differently, not investigated) |
| Ensembl Plants (REST `info/assembly/solanum_lycopersicum_gca000188115v5cm`) | SL4.0 (GCA_000188115.5), golden path 781,672,596 bp (= the GenBank total); top-level regions use the GenBank accessions (CM001064.4 …). Data release 80 |
| Precedent (Ensembl REST `info/assembly/arabidopsis_thaliana`) | TAIR10 includes **Mt** (366,924 bp) and **Pt** (154,478 bp) as top-level sequences. A plant reference that represents organelles as explicit, labelled molecules |

**What this establishes.**
- (i) ch00 is the SL4.0 set of unplaced contigs.
- (ii) The official SL4.0 assembly carries organelle-identical unplaced contigs as ordinary nuclear-labelled sequence. This
  holds for the SGN release and, for the plastid-identical contigs, for the GenBank release.
- (iii) ITAG4.0 annotated ordinary gene models on them.
- (iv) The GenBank release omits most of the mitochondrion-identical contigs that SGN keeps. The sources accessed do not
  state why.

**What it does not establish.** Whether a given ch00 contig is an organellar molecule assembled into the nuclear set
(assembly contamination) or a nuclear insertion (NUPT/NUMT). A full-length plastome copy in one 228,710 bp contig is more
consistent with an assembled organellar molecule than with an insertion. That is an inference, not a documented fact.

## 4. Corrected sequence-identity diagnostic (Task 8; criteria frozen in §1.4)

### 4.1 Method

`scripts/bash/organelle_multihsp_blast.sh`:
- BLAST+ 2.17.0+ megablast, e-value ≤ 1e-20, **all HSPs** kept (no `-max_hsps 1`).
- Queries: the 34,075 ITAG4.0 cDNAs and the 177 ch00 segments.
- Subjects: NC_007898.3 and NC_035963.1. Their SHA-256 is checked against the Phase 6B manifest.

`scripts/R/19_organelle_scope_diagnostic.R` computes the union-of-HSP coverage at ≥ 98 % identity and classifies each gene.

### 4.2 Gene classes (corrected rule vs Phase 6B single-HSP rule)

| Location | Corrected class | Genes | Tested | Ranked | In ≥ 1 significant leading edge |
|---|---|---|---|---|---|
| SL4.0ch00 | plastid-type near-identical | 155 | 38 | 38 | 13 |
| SL4.0ch00 | mitochondrion-type near-identical | 179 | 5 | 5 | 1 |
| SL4.0ch00 | partial (50 – 80 %; 3 plastid, 2 mito) | 5 | 0 | 0 | 0 |
| SL4.0ch00 | below 50 % organelle-identical coverage | 6 | 0 | 0 | 0 |
| SL4.0ch00 | no organelle hit | 168 | 4 | 4 | 1 |
| SL4.0ch01 – 12 | plastid-type near-identical | 42 | 5 | 5 | 2 |
| SL4.0ch01 – 12 | plastid-type partial | 19 | 1 | 1 | 1 |
| SL4.0ch01 – 12 | mitochondrion-type near-identical / partial | 46 / 20 | 0 | 0 | 0 |
| SL4.0ch01 – 12 | below 50 % | 69 | 16 | 16 | 2 |

- The corrected rule calls **334** ch00 genes near-identical, against 186 under the Phase 6B single-HSP rule. Every Phase 6B
  near-identical gene stays near-identical.
- 148 genes that Phase 6B called "diverged" (126) or "no close similarity" (22) are near-identical when all HSPs are
  combined: their single best HSP covered only part of the transcript (`ch00_class_crosstab_phase6b_vs_corrected.tsv`).
- Every gene excluded by the Phase 6B post-hoc rule (42) is corrected near-identical.
- Ranked organelle-identical genes **not** covered by any existing sensitivity run: 7. These are Solyc00g500052.1 on ch00,
  5 anchored near-identical genes (Solyc09g061460.3, Solyc12g038080.1, Solyc09g055840.1, Solyc09g064400.2,
  Solyc10g061887.1) and 1 anchored partial (Solyc02g021290.1). Three of them sit in significant leading edges, in two sets:
  - CC GO:0009579 (thylakoid): Solyc09g061460.3;
  - MF GO:0046906 (tetrapyrrole binding): Solyc12g038080.1 and Solyc02g021290.1.

  Their baseMean is 8.1 – 151.8 and none has padj < 0.05.

### 4.3 Organelle content of ch00 itself

- **Segments.** 177 segments total:
  - 21 are ≥ 90 % organelle-identical: 4 plastid, 17 mito; 1,545,876 bp; 300 genes, 37 tested.
  - 3 are intermediate: 232,117 bp; 45 genes.
  - 153 are < 10 % organelle-identical: 7,847,657 bp; 168 genes, 4 tested.
- **Identical bp.** 570,377 bp of ch00 are plastid-identical and 1,131,400 bp mito-identical (≥ 98 %).
- **Organelle genome coverage.** The union over all segments covers **155,461 / 155,461 bp of the plastome** and
  **446,257 / 446,257 bp of the mitogenome** at ≥ 98 % identity. Segment 15 alone (228,710 bp) is 100 % plastid-identical
  and covers the whole plastome. Segment 20 (164,455 bp) is 100 % mito-identical.
- **Consequence.** The primary Salmon index already contains both organelle genomes: as decoy sequence (ch00) and, for the
  annotated parts, as transcript targets (ch00 gene models).

## 5. Scope from existing outputs (Task 7)

- **Expression.**
  - ch00 near-identical tested genes: 43, baseMean 9 – 1,276.
  - Plastid-type genes: median baseMean 34.7; they carry 0.046 % of the summed baseMean of all tested genes.
  - Anchored organelle-like tested genes: 6 near-identical or partial, baseMean 8 – 152.
- **Reads (primary `quant.sf`, read-only).**
  - ch00 near-identical transcripts: 6,495 – 8,975 reads per sample, 0.036 – 0.061 % of mapped reads.
  - Anchored organelle-like transcripts: 241 – 413 reads, 0.0015 – 0.0026 %.
  - Decoy-assigned fragments are 8.45 – 10.57 % of processed reads. Salmon does not report which decoy they fell on, so the
    organelle-derived part of that fraction is unknown. Measuring it would need a new run, and the decision does not
    depend on it (§7).
- **Exact duplicates.** The 99 transcripts removed at indexing include 88 ch00 organelle-identical ones (62 plastid-type,
  26 mito-type), 2 anchored plastid-type and 9 with no organelle hit. Reads matching them are counted on the retained
  identical copy.
- **Gene level.** No gene has padj < 0.05, so no gene-level result exists that could depend on these models.
- **Set level** (`collection_overview_organelle_identical.tsv`).
  - Significant sets with ≥ 1 organelle-identical leading-edge gene: BP 28 / 120, CC 12 / 44, MF 8 / 51. Containing such a
    gene is not the same as depending on it.
  - Measured dependence (Phase 6B): at most **15.8 %** of a collection's significant sets are ch00-sensitive (BP, post-hoc
    rule). All other values are ≤ 12.5 %.
- **Themes** (`theme_contribution_organelle_identical.tsv`):

| Theme | Sets | Sets with organelle-identical LE gene | Share of LE gene slots | ch00-sensitive (pre-spec. / post-hoc) |
|---|---|---|---|---|
| ATP synthase / proton transport | 5 | 5 | 28.6 % | 4 / 4 |
| Photosynthesis / plastid / thylakoid | 25 | 6 | 1.1 % | 0 / 0 |
| Ribosome / translation / RNA processing | 34 | 5 | 0.2 % | 0 / 1 (no excluded gene in that set) |
| Heme / oxidoreductase | 11 | 1 | 0.7 % | 1 / 0 |
| Other (no theme) | 113 | 31 | 3.0 % | 17 / 20 |
| Defence, hormone, ER, glutathione, water, carbohydrate | 27 | 0 | 0 % | hormone 1 / 0, carbohydrate 1 / 1 (no excluded gene) |

Only the **ATP-synthase / proton-transport** theme depends on organelle-identical gene models. This dependence was already
documented and caveated in Phase 6B (`ch00_functional_investigation.md` §5).

## 6. Method guidance (Task 6)

| Source | Claim used here |
|---|---|
| Salmon documentation, "Preparing transcriptome indices (mapping-based mode)", <https://salmon.readthedocs.io/en/latest/salmon.html> (accessed 2026-10-26) | "We generally recommend that you build a decoy-aware transcriptome file". Option 2: use "the entire genome of the organism as the decoy sequence". Decoys "mitigate potential spurious mapping of reads that actually arise from some unannotated genomic locus that is sequence-similar to an annotated transcriptome". **No statement on organelles** |
| Srivastava et al. 2020, *Genome Biol* 21:239, doi:10.1186/s13059-020-02151-8 (Crossref-verified; full text PMC7487471) | Decoy rule: if the best decoy score is strictly greater than the best transcript score, the fragment is discarded; ties go to transcripts. The full-genome decoy (SAF) gave the highest concordance with the oracle. Duplicate transcripts are generally discarded at indexing. **No statement on organelles** |
| Patro et al. 2017, *Nat Methods* 14(4):417, doi:10.1038/nmeth.4197 (citation verified through the Crossref reference list of Srivastava 2020; direct Crossref call failed) | Salmon model (already cited for D-002) |
| nf-core/rnaseq 3.27.0 usage, <https://nf-co.re/rnaseq/latest/docs/usage/> (accessed 2026-10-26) | "the entire genome of the organism is used by default for the decoy-aware transcriptome … as is recommended"; "It is recommended to provide the most complete reference genome for your species"; Salmon collapses exact-duplicate transcripts by default. **No organelle-specific guidance** |
| Ensembl TAIR10 (precedent, not guidance) | Arabidopsis reference includes Mt and Pt as labelled top-level molecules |

No official source found states that organelle genomes **should** be added as targets or decoys for nuclear RNA-seq of
plant leaves with poly(A) selection. The guidance is to use the complete genome as decoy, which the project did. The tomato
SL4.0 genome happens to carry its organelle sequence inside ch00, unlabelled.

**Expected behaviour.** The following is inferred from the documented decoy rule and was **not tested**. For reads from
organellar transcripts:
1. Reads inside an annotated ch00 organelle-identical gene model align equally well to the transcript and to the ch00
   decoy. The tie goes to the transcript, so they are counted on that ch00 gene.
2. Reads outside annotated models (intergenic organelle regions, UTR extensions) align only to the decoy and are
   discarded.
3. Reads shared with an anchored near-identical gene are split by the EM.

Adding a separate organelle FASTA as **decoy** would not change 1 – 3, because identical decoy sequence is already present.
Adding organelle transcripts as **targets** would mainly relabel case 1: reads would be shared between new organelle
targets and identical ch00 models, or the identical copies would be collapsed as duplicates. Only case 3 could move reads
away from the 6 anchored organelle-like tested genes.

## 7. Assessment against the frozen criteria (§1.3)

| Dimension | Finding |
|---|---|
| E1 reference completeness | Organelle sequence **is represented**. The plastome and mitogenome are complete at ≥ 98 % inside the ch00 decoy, and 334 ch00 target transcripts are near-identical to them. They are not represented as **labelled** organelle records |
| E2 sequence ambiguity | Expressed organelle-identical models: 43 on ch00 + 6 anchored (near-identical or partial) among the 18,277 tested genes |
| E3 annotation provenance | ch00 = unplaced contigs; organelle-identical contigs are part of the official SL4.0 assembly (SGN; plastid contigs also in GenBank). Contamination vs NUPT/NUMT cannot be established from the sources |
| E4 effect on conclusions | (a) no gene-level result; (b) max 15.8 % of a collection ch00-sensitive (< 50 %); (c) only ATP synthase depends, already caveated. 7 ranked genes are not covered by a sensitivity run; 3 are in 2 leading edges, at low abundance |
| E5 method guidance | Whole-genome decoy recommended; no source requires organelle records |
| E6 reproducibility / cost | A branch would be feasible from versioned accessions (NC_007898.3, NC_035963.1) with tracked code, but per §6 it would mainly relabel reads already assigned to ch00 models |

Decision rules, in order:
- **STOP — not triggered.**
  - The architecture follows the documented recommendation (full-genome decoy), so it is not incorrect.
  - There is no padj < 0.05 gene.
  - The decision did not require a new quantification.
- **C — not met.** Its precondition, organelle sequence absent from the index (E1), is false. The secondary conditions are
  also not met: 15.8 % < 50 %, the only dependent theme is already caveated, and no E5 source requires inclusion.
- **B — not met.** Its precondition, organelle sequence absent (E1), is false.
- **A — applies.** "Organelle sequence is already represented."

## 8. Recommendation

**Option A — no new quantification needed** (recommendation only; owner decision pending, D-044). No
`organelle_sensitivity_plan.md` was written, because §1.3 requires one only for B or C.

Caveats that remain valid without a new quantification:
1. **Feature class.** Counts on the 334 ch00 near-identical gene models most likely reflect organellar transcripts
   captured despite poly(A) selection, not nuclear gene expression. For the 88 anchored near-identical and 39 anchored
   partial models, organellar origin of the reads is possible but less certain (a nuclear insertion could also be
   transcribed). For Phase 7 interpretation, the agent suggests a documented label: "organelle-identical feature class,
   corrected rule" (`gene_organelle_classification_corrected.tsv`). This is a reporting choice and needs no
   re-quantification. **Owner decision.**
2. **ATP synthase / proton transport** must not be interpreted as a nuclear transcriptional response. This is unchanged
   from Phase 6B.
3. **Seven ranked organelle-identical genes are outside the existing sensitivity runs.** Three of them sit in the leading
   edges of 2 significant sets. A gene-exclusion fgsea re-run with the corrected class would close this gap. It is a new
   analysis, so it was **not run**. **Owner decision** whether it is wanted before Phase 7.
4. **Reference source.** If the index were ever rebuilt from GenBank / Ensembl SL4.0 instead of SGN, 12 mito-identical
   contigs (832 kb) would be missing from the decoys. The SGN FASTA (SHA-256 in the manifest) must remain the reference.

**Owner decision (2026-10-26):**
- **Option A approved** (D-044): no new organelle-aware Salmon re-quantification.
- The gene-exclusion fgsea of caveat 3 is **not** run.
- The label "organelle-identical feature class" (caveat 1) may be used when technically appropriate. It states sequence
  identity only: sequence identity ≠ genomic origin. Wording in later phases must not imply that these features are
  proven organellar transcripts.

## 9. Correction of Phase 6B wording

Phase 6B stated that "the FASTA contains no plastome or mitogenome sequence" and that the Salmon index "had no organelle
genome, either as transcripts or as decoys". The relevant places are `ch00_functional_investigation.md` §1 and §5,
`functional_robustness_plan.md` §7 and `methods_references.md` §40.

The accurate statement: the reference has **no separate, labelled** organelle record, but SL4.0ch00 contains complete
plastome and mitogenome sequence (§4.3). This **refines** the Phase 6B text and does not contradict its conclusions.
Phase 6B §5 already allowed for "organellar contigs that remained in the nuclear assembly".

Correction notes were added at each place. The original sentences are kept.

## 10. Limitations

- Identity and coverage describe sequence similarity only. They do not establish genomic origin.
- The behaviour described in §6 is inferred from Salmon's documented rules. It was not tested by re-quantification.
- The organelle-derived share of decoy-assigned reads is unknown.
- The Hosmani et al. full text was not available; only the abstract was used.
- The reason why GenBank omits 12 mito-identical contigs is undocumented in the sources accessed.
- ch00 segment 1 has no exact GenBank counterpart. This was not investigated: it carries no organelle-identical sequence
  and no tested gene.

## 11. Outputs and reproducibility

| Script (post-hoc) | Outputs (`results/organelle_assessment/`) |
|---|---|
| Freeze of §1 (manual helper; command in `run_log.md`) | `prespecification/organelle_assessment_question_frozen.md`, `prespecification_checksum.tsv` |
| `scripts/bash/organelle_multihsp_blast.sh` (env with BLAST+ 2.17.0+) | `blast/blast_version.txt`, `blast/ch00_segments.tsv`, `blast/blastn_itag40_cdna_vs_organelle_allhsps.tsv`, `blast/blastn_ch00_segments_vs_organelle_allhsps.tsv`, `blast/blast_run_summary.tsv` |
| `scripts/R/19_organelle_scope_diagnostic.R` (env `stv-r`) | `multihsp_transcript_organelle_coverage.tsv`, `gene_organelle_classification_corrected.tsv` (ch00 genes + anchored genes with any ≥ 98 % organelle HSP; 709 rows), `class_summary_by_rule.tsv`, `ch00_class_crosstab_phase6b_vs_corrected.tsv`, `basemean_by_corrected_class.tsv`, `ch00_segment_organelle_coverage.tsv`, `ch00_segment_summary.tsv`, `salmon_read_share_by_class.tsv`, `theme_contribution_organelle_identical.tsv`, `set_level_organelle_identical_leading_edge.tsv`, `collection_overview_organelle_identical.tsv`, `sessionInfo_19.txt` |
| `scripts/bash/ch00_ncbi_provenance.sh` (network: NCBI Datasets v2 + E-utilities; FASTA kept in a temp dir) | `provenance/ncbi_GCA_000188115.5_sequence_report.tsv`, `provenance/ch00_segment_vs_ncbi_genbank.tsv`, `provenance/ch00_vs_ncbi_summary.tsv` |

Run order: `organelle_multihsp_blast.sh` → `19_organelle_scope_diagnostic.R` → `ch00_ncbi_provenance.sh`. Inputs that are
git-ignored (`reference/`, `results/genome/`, `results/quantification/salmon/`) are verified by SHA-256 where recorded.
