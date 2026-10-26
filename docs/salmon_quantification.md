# Salmon index and quantification (Phase 3)

All six single-end samples were quantified against a decoy-aware SL4.0 + ITAG4.0 index with **Salmon 2.7.0** through the Nextflow workflow. The library type came from `-l A` auto-detection. Run: `exotic_ardinghelli` (Nextflow 26.10.3, session `7785bce3…`), 2026-10-25 20:43–20:49 local time, 5 min 59 s wall time (6 FastQC re-runs + index + 6 quant + MultiQC). Technical QC: **no anomaly detected**.

## 1. Documentation basis (checked before running)

| Source | What it established |
|---|---|
| `salmon index --help`, `salmon quant --help` of the installed 2.7.0 | Option names and defaults (`-k` 31, `-d/--decoys`, `--keepDuplicates` off, `--fldMean` 250, `--fldSD` 25, `-l` default `A`, deterministic two-phase quantification default since 2.6.0) |
| Salmon 2.x docs <https://combine-lab.github.io/salmon> (quick start, library types, mapping modes, CLI, output formats) and `MIGRATION.md` | Salmon 2.x is a Rust rewrite with the same workflow and `quant.sf`/`lib_format_counts.json`/`meta_info.json` formats; the index must be built by 2.x. For decoy-aware indexing, append the genome to the FASTA and pass `-d decoys.txt` (decoys last). Selective alignment is the default mode. Single-end library types are `U`/`SF`/`SR`. |
| Release notes v2.7.0 | Index format v2 (rebuild required); no CLI option changes; quantification byte-identical to 2.6.0 |
| Salmon 1.x docs (readthedocs, labelled 1.12.0) | Decoy-aware index "generally recommend[ed]", full-genome decoys option; k = 31 "works well for reads of 75bp or longer"; for single-end reads the fragment-length distribution "cannot be estimated" and `--fldMean/--fldSD` set it, affecting effective lengths and TPM |
| nf-core/rnaseq 3.27.0 + nf-core/modules `salmon/index`, `salmon/quant` | Genome decoys built from all genome sequence names; Salmon 2.7.0 pinned; no extra fragment-length or bias flags by default |

Accessed 2026-10-25/26; URLs in `docs/methods_references.md` and `docs/pipeline_sources.md`.

## 2. Index

| Item | Value |
|---|---|
| Method | Decoy-aware (full-genome decoys): ITAG4.0 cDNA followed by the complete SL4.0 genome ("gentrome"); `decoys.txt` = the 13 genome sequence names |
| Transcript FASTA | `reference/SGN/ITAG4.0/ITAG4.0_cDNA.fasta` (34,075 records) |
| Decoy FASTA | `reference/SGN/SL4.0/S_lycopersicum_chromosomes.4.00.fa.gz` (13 sequences, 782.5 Mb) |
| Decoy construction | `gunzip -cf genome \| grep '^>' \| cut -d ' ' -f 1 \| sed 's/^>//' > decoys.txt`; `cat <(transcripts) <(genome) > gentrome.fa` (module `workflow/modules/salmon_index.nf`, same as nf-core) |
| Command (rendered) | `salmon index --threads 8 -t gentrome.fa --decoys decoys.txt -i salmon_index -k 31 --ramLimit 6` |
| k-mer / minimizer | k = 31 (Salmon default; docs: adequate for ≥ 75 bp reads); m = 19 (auto) |
| Other parameters | duplicates collapsed (default; `keep_duplicates: false`); poly-A clipping on (default). `--ramLimit 6` (GiB) only caps the external minimizer sort memory on this 15 GiB machine and does not change the index. |
| Build statistics | 638,498,617 distinct canonical k-mers; 14,607,175 unitigs; poly-A tails clipped from 19 transcripts; **99 exact-duplicate transcripts removed** (listed in `duplicate_clusters.tsv`); 0 short (< k) transcripts |
| Targets | 33,976 transcripts + 13 decoys (`num_refs` 33,989, `first_decoy_index` 33,976) |
| Index hashes (`info.json`) | `seq_hash` 4952695c…29985, `name_hash` 86755cea…f4adcb, `decoy_seq_hash` abfafdc3…3d46 |
| Runtime / resources | 3 min 29 s wall time, 359 % CPU, **peak RSS 4.7 GB**, peak vmem 9.9 GB (Nextflow trace) |
| Size on disk | 2.1 GB |
| Warnings | **none** (0 `WARN`/`ERROR` lines in the 102-line build log) |
| Location (git-ignored) | `results/genome/index/salmon_index/`, plus `decoys.txt` and `salmon_index.log` in `results/genome/index/` |

The 99 collapsed duplicates fall into 57 retained sequences: 89 on SL4.0ch00 (unanchored contigs) and 10 on chromosomes 2, 3, 9 and 11. Their ITAG descriptions are mostly organelle-type genes (cytochrome c oxidase subunit 2, NADH dehydrogenase subunit 4, photosystem II D1, maturase K, chloroplast ribosomal proteins …). This is a technical observation, not a biological claim. Salmon assigns all reads for such a sequence to the retained ID, so the removed IDs do not appear in `quant.sf`.

## 3. Quantification parameters

Rendered command (per sample, `workflow/modules/salmon_quant.nf`):

```bash
salmon quant --threads 6 --index salmon_index --libType A -r <sample>.fastq.gz --fldMean 250 --fldSD 25 -o <sample>
```

| Parameter | Value | Justification |
|---|---|---|
| Mapping mode | selective alignment (default) | Salmon docs recommend selective alignment with a decoy-aware index; `--validateMappings` is a no-op in 2.x |
| `--libType` | `A` | Owner instruction for Phase 3 and D-004: strandedness is not documented, so it is inferred per sample |
| `-r` | single-end reads | All six runs are SINGLE 75 nt |
| `--fldMean 250 --fldSD 25` | Salmon defaults, passed explicitly | See §3.1 |
| Bias flags | none (`--seqBias`/`--gcBias`/`--posBias` off) | nf-core default. Fragment-GC correction models whole fragments, whose length distribution cannot be observed in single-end data (only the prior is available). The DESeq2 vignette recommends `--gcBias` "unless you are certain that your data do not contain such bias", so this remains an **owner decision (D-012, Proposed)**. |
| Quantification mode | deterministic two-phase (Salmon 2.x default) | Byte-identical results across runs and thread counts (verified, §6) |
| `--incompatPrior` | 0 (default) | Strand-incompatible mappings are discarded once the library type is resolved |
| Resources | 6 CPUs, 2 tasks in parallel | Peak RSS 2.2–2.3 GB per task (index loaded); 20–26 s per sample |

### 3.1 Single-end fragment-length parameters

- Salmon cannot estimate a fragment-length distribution from single-end reads. `--fldMean/--fldSD` define it, and it enters the effective transcript lengths (and hence TPM and the tximport length offset).
- No source records the library insert size: GEO and ENA give none, and SRA RunInfo reports `InsertSize = 0` and `InsertDev = 0` for all six runs. No Bioanalyzer trace is available.
- The Salmon defaults (mean 250, SD 25) were therefore kept and passed **explicitly**, so the choice is visible in the command line and in `meta_info.json` (`frag_length_source: prior`, mean 250.0, SD 25.0 in all six samples).
- The same assumed distribution applies to every sample, so it does not create between-sample differences. It can affect absolute effective lengths, mostly for transcripts shorter than about 250 nt.
- The values are workflow parameters (`--salmon_fld_mean`, `--salmon_fld_sd`). Changing them is a documented sensitivity option, not a silent edit.

## 4. Quantification QC

Script: [`scripts/bash/summarize_salmon.sh`](../scripts/bash/summarize_salmon.sh). It reads `aux_info/meta_info.json`, `lib_format_counts.json` and `quant.sf` and writes [`results/quantification/salmon_qc_summary.tsv`](../results/quantification/salmon_qc_summary.tsv) and [`results/quantification/salmon_libformat_summary.tsv`](../results/quantification/salmon_libformat_summary.tsv).

| Sample | Condition | SRR | Input reads | Mapped | Mapping rate (%) | Decoy fragments (% input) | Inferred libtype | Quantified transcripts | Transcripts with reads > 0 | Salmon |
|---|---|---|---:|---:|---:|---:|---|---:|---:|---|
| STV_free_rep1 | STV_free | SRR10102791 | 18,080,482 | 15,902,183 | 87.95 | 1,911,775 (10.57) | SR | 33,976 | 22,019 | 2.7.0 |
| STV_free_rep2 | STV_free | SRR10102792 | 16,288,249 | 14,698,030 | 90.24 | 1,375,926 (8.45) | SR | 33,976 | 21,665 | 2.7.0 |
| STV_free_rep3 | STV_free | SRR10102793 | 14,318,362 | 12,871,350 | 89.89 | 1,261,143 (8.81) | SR | 33,976 | 21,744 | 2.7.0 |
| STV_infected_rep1 | STV_infected | SRR10102794 | 19,230,359 | 17,121,521 | 89.03 | 1,853,818 (9.64) | SR | 33,976 | 21,839 | 2.7.0 |
| STV_infected_rep2 | STV_infected | SRR10102795 | 16,858,900 | 15,217,015 | 90.26 | 1,427,474 (8.47) | SR | 33,976 | 21,591 | 2.7.0 |
| STV_infected_rep3 | STV_infected | SRR10102796 | 19,862,460 | 17,871,112 | 89.97 | 1,726,957 (8.69) | SR | 33,976 | 22,274 | 2.7.0 |

Checks (all passed):
- `num_processed` equals the raw read count (SeqKit/ENA) in every sample.
- The same Salmon version, index (33,976 targets + 13 decoys) and transcript set appear in all six `quant.sf`.
- The sum of `NumReads` equals `num_mapped`.
- `meta_info.json` `diagnostics` is empty in every sample, i.e. Salmon raised none of its low-mapping, strand-bias, library-type-mismatch or incompatible-fraction warnings.
- 0 `WARN`/`ERROR` lines in any task log.

### 4.1 Library format (`lib_format_counts.json`, all six inspected)

| Sample | expected_format | compatible_fragment_ratio | Wrong-strand-only fragments (discarded) | SF placements | SR placements | U | strand_mapping_bias |
|---|---|---:|---:|---:|---:|---:|---:|
| STV_free_rep1 | SR | 0.99237 | 122,256 | 346,296 | 15,902,183 | 0 | 0.0213 |
| STV_free_rep2 | SR | 0.99394 | 89,553 | 250,587 | 14,698,030 | 0 | 0.0168 |
| STV_free_rep3 | SR | 0.99396 | 78,196 | 225,687 | 12,871,350 | 0 | 0.0172 |
| STV_infected_rep1 | SR | 0.99331 | 115,256 | 366,177 | 17,121,521 | 0 | 0.0209 |
| STV_infected_rep2 | SR | 0.99427 | 87,749 | 255,624 | 15,217,015 | 0 | 0.0165 |
| STV_infected_rep3 | SR | 0.99362 | 114,779 | 322,172 | 17,871,112 | 0 | 0.0177 |

- **Orientation is consistent:** all six libraries are detected as **SR** (single-end, reads from the reverse strand of the transcript).
- 99.24–99.43 % of the fragments judged by the strand filter are strand-compatible.
- Only 1.7–2.1 % of placements are observed in the sense orientation (`strand_mapping_bias`).
- This matches the expectation inferred in Phase 1 from the reported kit (TruSeq stranded mRNA → `SR`); the expectation is now **measured**, not assumed (D-004).
- Wrong-strand-only fragments (0.55–0.76 % of the judged fragments) were discarded (`--incompatPrior 0`). Their fraction is similar across samples, so no library stands out. Following Salmon's docs, this fraction also bounds possible double-stranded (e.g. genomic DNA) carry-over.

### 4.2 Mapping rate by condition (descriptive only)

| Condition | n | Mean (%) | Min (%) | Max (%) |
|---|---:|---:|---:|---:|
| STV_free | 3 | 89.36 | 87.95 | 90.24 |
| STV_infected | 3 | 89.75 | 89.03 | 90.26 |

- The ranges overlap. The lowest rate (STV_free_rep1, 87.95 %) coincides with the highest decoy share (10.57 %).
- No formal test was made and no biological meaning is attached. Mapping rate is a technical property of the libraries and the reference.

**Verdict:** no abnormal sample, orientation or mapping pattern, so nothing blocked continuing to the tximport import check.

## 5. MultiQC

`results/qc/multiqc_salmon/multiqc_report.html` (git-ignored) aggregates FastQC + Salmon; MultiQC 1.35 parsed all six Salmon directories (`multiqc_salmon.txt`). The Phase 2 raw-QC report in `results/qc/multiqc_raw/` was not overwritten (`--multiqc_dir qc/multiqc_salmon`).

## 6. Library-type setting: `-l A` vs explicit `-l SR` (D-010, owner decision)

D-004 foresees fixing `-l` to the observed type once all six agree, and nf-core/rnaseq also quantifies with the inferred type. As a check, quantification was repeated with `--salmon_libtype SR` using the same index (run `sick_archimedes`, output to `/tmp` only; `results/` untouched). Comparison script: [`scripts/bash/compare_salmon_runs.sh`](../scripts/bash/compare_salmon_runs.sh) → [`results/quantification/libtype_A_vs_SR_comparison.tsv`](../results/quantification/libtype_A_vs_SR_comparison.tsv).

| Sample | `num_mapped` A = SR | Eq. classes A / SR | Transcripts with \|ΔNumReads\| > 0.001 | Max \|ΔNumReads\| | Effective-length differences | Pearson r (log1p NumReads) |
|---|---|---|---:|---:|---:|---:|
| STV_free_rep1 | yes | 26,308 / 26,166 | 621 | 635.0 | 0 | 0.99969 |
| STV_free_rep2 | yes | 26,026 / 25,880 | 720 | 985.4 | 0 | 0.99972 |
| STV_free_rep3 | yes | 26,006 / 25,864 | 680 | 1,377.0 | 0 | 0.99965 |
| STV_infected_rep1 | yes | 26,310 / 26,167 | 690 | 1,017.6 | 0 | 0.99967 |
| STV_infected_rep2 | yes | 25,852 / 25,709 | 661 | 1,044.9 | 0 | 0.99968 |
| STV_infected_rep3 | yes | 26,960 / 26,796 | 774 | 2,218.9 | 0 | 0.99959 |

Facts established:
- `-l A` itself is **reproducible**: rerunning one sample outside Nextflow with 6 and 3 threads gave byte-identical `quant.sf`.
- Explicit `SR` and auto-detected `SR` map exactly the same fragments and produce identical `lib_format_counts.json`.
- They differ only in the equivalence-class structure (≈ 0.55 % fewer classes with explicit SR). As a result, about 2 % of transcripts (621–774) receive different read assignments, in both directions (e.g. STV_free_rep1: 314 higher with SR, 307 lower).
- The changed transcripts are **not** enriched for opposite-strand overlapping genes (2.6 % vs 4.3 % in an unchanged control set), so the mechanism is not established.

**Decision needed (D-010):** keep the `-l A` results as the primary quantification (current state, as instructed for Phase 3), or re-run the final quantification with explicit `SR` (`strandedness=reverse` in the samplesheet), as D-004 and nf-core practice suggest. The tximport/DESeq2 input should be fixed before any differential-expression run.

## 7. Limitations

- **Fragment length is assumed** (§3.1), not measured.
- **99 duplicate transcripts are merged** into 57 retained IDs (Salmon default). Their counts are attributed to the retained gene IDs at gene level. Mostly SL4.0ch00 organelle-type models.
- **Only the 34,075 ITAG4.0 protein-coding models are quantified.** The separately released ITAG4.0 lncRNA annotation is not part of the cDNA FASTA. Reads from unannotated loci are absorbed by the genome decoys (8.5–10.6 % of reads) instead of being misassigned.
- **ITAG4.0 has one transcript per gene,** so the transcript-level step adds no isoform resolution. tximport still provides the average-length offset.
- **No sequence/GC bias correction** (D-012 open).
- **Salmon 2.x is a 2026 Rust rewrite.** The method papers describe the C++ implementation, and the developers report equivalent outputs (`MIGRATION.md`, release notes).

## 8. Reproduce

```bash
conda activate stv-rnaseq
bash scripts/bash/download_reference.sh
bash scripts/bash/validate_reference.sh
nextflow run workflow/main.nf -resume --run_quant true \
  --transcript_fasta reference/SGN/ITAG4.0/ITAG4.0_cDNA.fasta \
  --genome_fasta reference/SGN/SL4.0/S_lycopersicum_chromosomes.4.00.fa.gz \
  --multiqc_title 'STV_GSE137303 QC + Salmon (SL4.0 / ITAG4.0)' --multiqc_dir qc/multiqc_salmon
bash scripts/bash/summarize_salmon.sh
# re-use the published index instead of rebuilding: --salmon_index results/genome/index/salmon_index
```
