# SL4.0ch00 functional investigation (Phase 6B, Task 9; D-038)

**Question.** Phase 6 found that several recurring leading-edge genes are chloroplast-type ATP-synthase subunits annotated on
the unanchored pseudo-chromosome SL4.0ch00 (post-hoc check, `functional_analysis_plan.md` §8.3). Are these genes copies of
organelle sequences? Do the Phase 6 pathway results depend on them?

The rules were fixed before the run ([functional_robustness_plan.md](functional_robustness_plan.md) §7). The outcome is
summarized in §11 of that plan. This document gives the technical detail and makes no biological interpretation.

## 1. What SL4.0ch00 is (SGN documentation)

Source: SGN `build_4.00/ReleaseNotes.txt` (<https://solgenomics.net/ftp/tomato_genome/assembly/build_4.00/ReleaseNotes.txt>).
It was retrieved 2026-10-26; Last-Modified 2019-09-13; SHA-256 `857189ab…` in `reference/reference_manifest.tsv`.

| Release | Statement (verbatim) | Relevance |
|---|---|---|
| SL4.0 | "De novo Canu assembly (N50 5.5 Mb) followed by Hi-C scaffolding (12 chromosomes and unplaced contigs)." "Only 152 unplaced contigs in Chr 00 compared to 4,374 in SL3.0" | ch00 = unplaced contigs |
| SL2.40 | "4 contigs that were identified as completely mitochondrial have been removed." | Earlier builds screened organellar contigs |
| SL1.50 | "Organellar sequences were separated (and are thus not included in this data set)." | Same |

The SL4.0 notes do not say whether organelle-derived contigs remain in ch00. In the SL4.0 FASTA used here, ch00 is
9,643,250 bp and carries 513 ITAG4.0 genes. The FASTA contains no plastome or mitogenome sequence. Our Salmon index
(D-001, D-010) therefore had no organelle genome, either as transcripts or as decoys. Reads from organellar transcripts could be
assigned to any ITAG4.0 transcript that is (near-)identical to organelle sequence.

> **Correction (2026-10-26, Phase 6B closure, D-044).** The FASTA has no **separate, labelled** plastome or mitogenome
> record. However, a multi-HSP BLAST of the 177 gap-free ch00 segments shows that SL4.0ch00 itself contains the **complete**
> plastome and mitogenome sequence at ≥ 98 % identity. One 228,710 bp segment covers the whole plastome. The Salmon index
> therefore contained organelle sequence: as decoy (ch00) and, for the annotated parts, as ch00 transcript targets. The
> conclusions of this document are unchanged. Details are in
> [organelle_requantification_assessment.md](organelle_requantification_assessment.md) §4.3.

## 2. Methods

| Step | Detail |
|---|---|
| Organelle references | Tomato plastome **NC_007898.3** (*S. lycopersicum* chloroplast, complete genome) and mitogenome **NC_035963.1** (*S. lycopersicum* TGRC:LA1421 mitochondrion, complete genome), NCBI RefSeq via E-utilities `efetch` (genome FASTA + `fasta_cds_aa`), 2026-10-26. Files are git-ignored under `reference/functional/organelle/`; SHA-256 and first header are in `results/functional_robustness/ch00/organelle_reference_manifest.tsv` (87 plastid and 39 mitochondrial proteins). `scripts/bash/download_organelle_reference.sh` |
| BLAST+ | 2.17.0+ from the pre-existing conda env `rachis-qiime2-2026.4` (read-only use; nothing installed). `scripts/bash/ch00_organelle_blast.sh` |
| blastn | `-task megablast`, query = all ITAG4.0 cDNAs (one per gene), subject = both organelle genomes, `-evalue 1e-20 -max_target_seqs 2 -max_hsps 1`, tabular output with `pident`, `qcovs`. 1,061 of 34,075 queries have a hit |
| blastp | Proteins of the 513 ch00 genes vs the 126 organelle-encoded proteins, `-evalue 1e-10 -max_target_seqs 1 -max_hsps 1` (support only) |
| Classes (pre-specified) | Best blastn hit by bit score. **plastid-/mitochondrion-type near-identical**: identity ≥ 98 % and query coverage ≥ 80 %. **organelle-similar, diverged**: identity ≥ 80 % and coverage ≥ 50 %. Otherwise **no close organelle similarity** (the plan's "nuclear unanchored" for ch00 genes). Genes are not removed for being unanchored |
| Sensitivity (pre-specified) | Remove ranked ch00 genes in a near-identical class from the ranking and the gene sets. Re-run fgsea with the Phase 6 call and seed. A primary-significant set is *ch00-sensitive* if padj ≥ 0.05 or its NES sign flips. STOP if > 50 % of a collection is ch00-sensitive. `scripts/R/15_functional_ch00_investigation.R` |
| Post-hoc diagnostic (not pre-specified) | Broader exclusion: every ranked ch00 gene with identity ≥ 98 % and coverage ≥ 50 % (§4). `scripts/R/18_functional_ch00_posthoc_identity_diagnostic.R` |

The BLAST+ warning "Examining 5 or more matches is recommended" refers to `-max_target_seqs` < 5. With only two subject
sequences (blastn) every possible subject is still reported. For blastp, only the best protein is reported, which is the
intended use.

## 3. Results (pre-specified rule)

### 3.1 Organelle similarity of ITAG4.0 genes

| Location | Class | Genes | Tested (18,277) | Ranked (18,238) | In ≥ 1 significant leading edge | Median baseMean (ranked) |
|---|---|---|---|---|---|---|
| SL4.0ch00 | plastid-type near-identical | 71 | 19 | 19 | 5 | 146.7 |
| SL4.0ch00 | mitochondrion-type near-identical | 115 | 3 | 3 | 1 | 14.0 |
| SL4.0ch00 | organelle-similar, diverged | 136 | 20 | 20 | 8 | 29.9 |
| SL4.0ch00 | no close organelle similarity | 191 | 5 | 5 | 1 | 17.2 |
| SL4.0ch01 – 12 | plastid-type near-identical | 32 | 3 | 3 | 1 | 9.1 |
| SL4.0ch01 – 12 | mitochondrion-type near-identical | 32 | 0 | 0 | 0 | — |
| SL4.0ch01 – 12 | organelle-similar, diverged | 448 | 4 | 4 | 2 | 90.9 |
| SL4.0ch01 – 12 | no close organelle similarity | 33,050 | 18,223 | 18,184 | 2,840 | 201.1 |

- 186 of the 513 ch00 genes (36 %) are near-identical to an organelle genome. Only 22 of them pass the expression prefilter.
- 47 ch00 genes are ranked. The frozen plan stated 45, which was a counting error (plan §11.3).
- The 15 ch00 genes in significant leading edges are 5 plastid-type near-identical, 1 mitochondrion-type near-identical, 8
  diverged and 1 with no close similarity.
- On the anchored chromosomes, 3 tested genes are plastid-type near-identical: Solyc09g055840.1, Solyc09g064400.2 and
  Solyc12g038080.1 (1 leading-edge membership). These are reported only; the pre-specified sensitivity covers ch00 genes only.

Table: `results/functional_robustness/ch00/organelle_class_summary.tsv`. Per gene:
- `ch00_gene_classification.tsv`: all 513 ch00 genes, with the blastn and blastp best hits, the DESeq2 statistic and leading-edge counts;
- `organelle_type_genes_all_chromosomes.tsv`.

### 3.2 Recurring ch00 leading-edge genes

| Gene | Class | blastn identity / coverage (%) | blastp best organelle protein | Wald stat | Significant sets with gene in leading edge |
|---|---|---|---|---|---|
| Solyc00g500064.1 | plastid-type near-identical | 100 / 100 | atpB | 0.96 | 34 |
| Solyc00g500140.1 | plastid-type near-identical | 99.3 / 100 | atpE | 1.21 | 34 |
| Solyc00g500323.1 | plastid-type near-identical | 100 / 100 | atpA | 1.86 | 34 |
| Solyc00g500332.1 | plastid-type near-identical | 100 / 100 | atpA | 1.06 | 34 |
| Solyc00g500322.1 | organelle-similar, diverged | 100 / 78 | atpF | 0.92 | 31 |
| Solyc00g500321.1 | organelle-similar, diverged | 100 / 67 | atpH | 0.57 | 9 |
| Solyc00g500335.1 | organelle-similar, diverged | 100 / 59 | atpI | 0.52 | 8 |
| Solyc00g500212.1 | plastid-type near-identical | 100 / 100 | rps4 | 0.53 | 6 |
| Solyc00g500021.1 | organelle-similar, diverged | 98.8 / 77 | psbE | 0.34 | 5 |
| Solyc00g500049.1 | organelle-similar, diverged | 100 / 60 | psbK | 0.45 | 5 |
| Solyc00g500020.1 | organelle-similar, diverged | 100 / 64 | petA | 0.54 | 4 |
| Solyc00g500358.1 | mitochondrion-type near-identical | 100 / 100 | atp1 | 0.35 | 3 |

The other 3 of the 15 have one membership each: Solyc00g500023.1 (diverged, clpP), Solyc00g500054.1 (diverged, psbC) and
Solyc00g500353.1 (no close similarity).

**Annotation note.** For some ch00 genes, the ITAG4.0 functional description (AHRD) and the best blastp organelle protein differ.
For example:
- Solyc00g500064.1 is described as "Ycf15", but its protein is 100 % identical to plastid atpB over 99 % of its length, and its
  cDNA (1,497 nt) matches the plastome over its full length.
- Solyc00g500056.1 is described as "DNA-directed RNA polymerase subunit alpha" but has psaB as its best hit.

Both fields are reported side by side. Nothing is re-annotated, and GO sets are not changed.

### 3.3 Pre-specified sensitivity (22 near-identical ch00 genes removed)

| Collection | Primary-significant | Retained | Supported | Lost | ch00-sensitive (%) | New significant | Spearman ρ (NES) |
|---|---|---|---|---|---|---|---|
| BP | 120 | 105 | 12 | 3 | 15 (12.5 %) | 3 | 0.998 |
| MF | 51 | 46 | 3 | 2 | 5 (9.8 %) | 0 | 0.994 |
| CC | 44 | 40 | 0 | 4 | 4 (9.1 %) | 0 | 0.995 |

- **STOP not triggered** (maximum 12.5 % < 50 %).
- The 9 lost sets each contain 5 excluded genes:
  - BP: ATP biosynthetic process, proton motive force-driven ATP synthesis, purine-containing compound biosynthetic process.
  - MF: proton channel activity, proton-transporting ATP synthase activity (rotational mechanism).
  - CC: proton-transporting ATP synthase complex, proton-transporting two-sector ATPase complex, monoatomic ion channel complex,
    cation channel complex.
- Most "supported" BP sets are nucleotide/purine biosynthesis or metabolism terms that contain the same genes.
- Four sets with no excluded gene moved from padj < 0.05 to 0.05 – 0.07 because the ranking changed:
  - regulation of jasmonic acid mediated signaling pathway;
  - pectin metabolic process;
  - lipid carrier activity;
  - oxidoreductase activity acting on the aldehyde or oxo group of donors.

  These are borderline results, not ch00 effects.
- Photosynthesis/plastid/thylakoid and ribosome/translation sets are **not** ch00-sensitive (0 of 25 and 0 of 34;
  `themes/theme_membership_significant_sets.tsv`). Figure F4: `figures/functional_robustness/ch00_sensitivity_nes.png`.

## 4. Limitation found after the run and post-hoc diagnostic

blastn used `-max_hsps 1`. For a spliced organelle gene (for example atpF, clpP, rpl2, ndhA), the cDNA aligns to the organelle
genome as several HSPs. Only the best HSP is kept, so coverage stays below 80 % even at ~100 % identity. As a result, **129 of
the 136 "diverged" ch00 genes have identity ≥ 98 %**; 20 of them are ranked, including atpF, atpH and atpI (31, 9 and 8
leading-edge memberships).

The pre-specified class was not changed. A labelled post-hoc diagnostic excluded every ranked ch00 gene with identity ≥ 98 %
and coverage ≥ 50 % (42 genes, of which 20 are beyond the pre-specified exclusion):

| Collection | ch00-sensitive (post-hoc) | Lost | Supported | New significant |
|---|---|---|---|---|
| BP | 19 (15.8 %) | 3 | 16 | 12 |
| MF | 3 (5.9 %) | 2 | 1 | 0 |
| CC | 4 (9.1 %) | 4 | 0 | 0 |

The same 9 sets are lost. The ch00 STOP assessment does not change. Outputs:
`results/functional_robustness/ch00/posthoc_identity98_{excluded_genes,vs_primary,summary}.tsv`.

## 5. Technical reading

- ch00 contains many gene models whose sequences are (near-)identical to the tomato plastome or mitogenome. This analysis
  cannot tell organellar contigs that remained in the nuclear assembly apart from nuclear insertions of organellar DNA. Both
  would attract reads from organellar transcripts, because the reference has no organelle genome.
  *Correction (D-044): the reference has no labelled organelle record, but ch00 contains complete organelle sequence; see the
  note in §1.*
- The read counts of these genes may therefore reflect organellar transcripts captured despite poly(A) selection. They should
  not be read as expression of nuclear genes.
- In the Phase 6 results, the sets that contain these genes depend on them: ATP-synthase / proton-transport terms are lost,
  and nucleotide / purine biosynthesis terms fall to 0.05 ≤ padj < 0.25. Under either exclusion rule, the other sets that
  change status contain none of the excluded genes. They are borderline sets whose padj moved just above 0.05 when the
  ranking changed; under the post-hoc rule this includes one ribosome-theme set (tRNA methylation).
- This is a **technical caveat**. It is not evidence for or against any biological effect of STV.
- What could be done later (owner decision, not done here):
  - add the organelle genomes to the Salmon index as transcripts or decoys;
  - or treat organelle-identical gene models as a separate feature class.

  Either would be a change to the quantification and therefore a new analysis.

  *Update (D-044): the first option was assessed without re-quantification in
  [organelle_requantification_assessment.md](organelle_requantification_assessment.md). The recommendation is A, no new
  quantification, pending owner decision. The second option is a reporting choice that needs no re-quantification.*

## 6. Outputs

`results/functional_robustness/ch00/`:
- `organelle_reference_manifest.tsv`, `blast_version.txt`, `blast_run_summary.tsv`;
- `blastn_itag40_cdna_vs_organelle_genomes.tsv`, `blastp_ch00_proteins_vs_organelle_proteins.tsv`, `ch00_genes_gff.tsv`;
- `ch00_gene_classification.tsv`, `organelle_type_genes_all_chromosomes.tsv`, `organelle_class_summary.tsv`;
- `ch00_sensitivity_{excluded_genes,vs_primary,summary}.tsv`;
- `posthoc_identity98_*.tsv`, `sessionInfo_15.txt`, `sessionInfo_18.txt`.

Figure: `figures/functional_robustness/ch00_sensitivity_nes.png`.

References:
- Camacho et al. (2009), BLAST+;
- Kahlau et al. (2006), tomato plastome sequence;
- NCBI RefSeq records NC_007898.3 and NC_035963.1;
- SGN release notes.

See [methods_references.md](methods_references.md) §39 – §40.
