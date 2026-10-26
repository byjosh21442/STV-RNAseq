# Nextflow pipeline — `workflow/`

Small DSL2 workflow for this six-sample, **single-end** dataset. Architecture modelled on nf-core/rnaseq 3.27.0 and the Nextflow RNAseq training (sources: `docs/pipeline_sources.md`); it is not a copy of nf-core/rnaseq.

Status (Phase 3): samplesheet validation, FastQC, MultiQC, the decoy-aware **Salmon index (SL4.0 + ITAG4.0) and Salmon quantification of all six samples have run** on real data (run `exotic_ardinghelli`, 5 min 59 s; `docs/salmon_quantification.md`). Phase 2 raw QC: `docs/raw_qc_assessment.md`.

## Diagram

```mermaid
flowchart LR
    S[config/samplesheet.csv] --> V{validate rows<br/>SINGLE-end only}
    V --> R[(data/raw/*.fastq.gz)]
    R --> F[FASTQC]
    F --> M[MULTIQC]
    R -. run_quant=true .-> Q[SALMON_QUANT<br/>-l from strandedness]
    T[ITAG4.0 cDNA] -. no prebuilt index .-> I[SALMON_INDEX<br/>genome decoys]
    G[SL4.0 genome] -.-> I
    I -.-> Q
    I -.-> O4[results/genome/index]
    Q -.-> M
    Q -.-> X[quantification/salmon/samples.csv<br/>for tximport]
    F --> O1[results/qc/fastqc_raw]
    M --> O2[results/qc/multiqc_raw or multiqc_salmon]
```

Dashed = optional branch (`--run_quant true`), executed on real data in Phase 3. Trimming is intentionally absent (D-003/D-009): it would be inserted between the raw FASTQ and Salmon only after an approved QC justification.

## Files

| Path | Content |
|---|---|
| `workflow/main.nf` | Typed `params` block, samplesheet parsing/validation functions, entry workflow, `output {}` block |
| `workflow/modules/fastqc.nf` | `FASTQC` (label `process_low`); reports renamed to the sample ID |
| `workflow/modules/multiqc.nf` | `MULTIQC` (label `process_single`) |
| `workflow/modules/salmon_index.nf` | `SALMON_INDEX` (label `process_high`); full-genome decoys if `--genome_fasta` given; also emits `decoys.txt` and `salmon_index.log` |
| `workflow/modules/salmon_quant.nf` | `SALMON_QUANT` (label `process_medium`); single-end `-r`, `--libType` from samplesheet |
| `workflow/nextflow.config` | Output dir, resource caps/labels, Salmon `ext.args` (k-mer, `--ramLimit`, SE fragment length), reports, `test_stub` profile |
| `workflow/tests/stub/` | Stub samplesheet + `make_fixtures.sh` that generates synthetic fixtures (4 reads/sample, toy FASTA); generated files are git-ignored so no FASTQ is committed |
| `workflow/bootstrap/` | Toy validation workflow (see `docs/nextflow_validation.md`) |

## Running

From the project root with `conda activate stv-rnaseq`:

```bash
# Raw QC (what Phase 2 ran)
nextflow run workflow/main.nf -resume

# Wiring test incl. Salmon, no tools executed
bash workflow/tests/stub/make_fixtures.sh
nextflow run workflow/main.nf -profile test_stub -stub -w work/stub -output-dir /tmp/stv_stub

# Phase 3: decoy-aware index + quantification (approved reference D-001)
nextflow run workflow/main.nf -resume --run_quant true \
  --transcript_fasta reference/SGN/ITAG4.0/ITAG4.0_cDNA.fasta \
  --genome_fasta reference/SGN/SL4.0/S_lycopersicum_chromosomes.4.00.fa.gz \
  --multiqc_title 'STV_GSE137303 QC + Salmon (SL4.0 / ITAG4.0)' --multiqc_dir qc/multiqc_salmon

# Re-quantify with the published index (e.g. explicit library type, D-010)
# nextflow run workflow/main.nf -resume --run_quant true \
#   --salmon_index results/genome/index/salmon_index --salmon_libtype SR
```

## Parameters (`params` block in `main.nf`)

| Param | Default | Meaning |
|---|---|---|
| `input` | `config/samplesheet.csv` | Samplesheet (relative to `project_dir`) |
| `project_dir` | parent of `workflow/` | Base for relative paths |
| `run_quant` | `false` | Enable Salmon |
| `salmon_index` | `''` | Prebuilt index (skips `SALMON_INDEX`) |
| `transcript_fasta` / `genome_fasta` | `''` | Reference files; **no defaults on purpose** (passed on the command line, D-001) |
| `salmon_libtype` | `''` | Global override (e.g. `A`, `SR`); empty = per-sample mapping `auto→A, unstranded→U, forward→SF, reverse→SR` |
| `salmon_kmer` | `31` | `salmon index -k` |
| `salmon_fld_mean` / `salmon_fld_sd` | `250` / `25` | Single-end fragment-length distribution for `salmon quant` (D-011) |
| `multiqc_title` | `STV_GSE137303 raw FASTQ QC` | Report title |
| `fastqc_dir`, `multiqc_dir`, `salmon_dir`, `salmon_index_dir`, `pipeline_info_dir` | `qc/fastqc_raw`, `qc/multiqc_raw`, `quantification/salmon`, `genome/index`, `pipeline_info` | Publish paths under `outputDir` (`results/`) |

Biological assumptions (strandedness, reference) live in the samplesheet or parameters, never in process code.

## Samplesheet validation

Before any task is submitted, `readSamplesheet()` checks required columns, unique IDs, `library_layout == SINGLE`, no `fastq_2`, allowed `strandedness`, non-empty `condition`, `.fastq.gz`/`.fq.gz` suffix and that each FASTQ exists. Any violation stops the run with an explicit error.

## Outputs (`outputDir = results`, publish mode `copy`)

| Output | Path | Git |
|---|---|---|
| FastQC HTML + ZIP | `results/qc/fastqc_raw/<sample>_fastqc.{html,zip}` | ignored |
| MultiQC | `results/qc/multiqc_raw/` (Phase 2 raw QC) and `results/qc/multiqc_salmon/` (Phase 3, FastQC + Salmon) | ignored (regenerable); small TSV summaries in `results/qc/*.tsv` are tracked |
| Salmon | `results/quantification/salmon/<sample>/` + `samples.csv` index | ignored; QC summaries `results/quantification/*.tsv` tracked |
| Salmon index | `results/genome/index/salmon_index/` (2.1 GB), `decoys.txt`, `salmon_index.log` | ignored |
| Run reports | `results/pipeline_info/` (report, timeline, trace, `software_versions.tsv`) | ignored |

## Resources (personal WSL2 machine: 15 GiB RAM, 20 threads)

`executor.cpus = 12`, `executor.memory = 12 GB`, `process.resourceLimits = [cpus: 12, memory: 12 GB]`; labels: `process_single` 1 CPU / 2 GB, `process_low` 1 CPU / 1 GB (FastQC; `--memory` = 80 % of task memory), `process_medium` 6 CPU / 6 GB (Salmon quant), `process_high` 8 CPU / 12 GB (Salmon index). `errorStrategy = 'terminate'`.

Measured in Phase 3 (Nextflow trace): `SALMON_INDEX` 3 min 29 s, peak RSS 4.7 GB (with `--ramLimit 6`); `SALMON_QUANT` 20–26 s per sample, peak RSS 2.2–2.3 GB, two tasks in parallel; FastQC 54–73 s, 0.32–0.34 GB.

## Differences from nf-core/rnaseq (and why)

| nf-core/rnaseq 3.27.0 | This workflow | Reason |
|---|---|---|
| Many steps: fq lint, trimming (Trim Galore/fastp default), rRNA removal, STAR/HISAT2, RSeQC, Qualimap, dupRadar, BBSplit, contamination screening… | FASTQC → (SALMON_INDEX → SALMON_QUANT) → MULTIQC | Six samples, two-week project; approved architecture D-002; no trimming without justification (D-003) |
| Paired- and single-end | Single-end only; PAIRED rows rejected | Verified dataset layout; explicit failure is safer than silent support |
| nf-schema plugin + JSON schema for samplesheet/params | Plain Nextflow functions + typed `params` block | No plugin download; enough for 6 rows |
| `strandedness=auto` → Salmon on 1 M subsampled reads, then thresholds 0.8 / 0.1, then quantification with the inferred type | `auto` → Salmon `-l A` on all reads; `lib_format_counts.json` of all six inspected (D-004: SR in 6/6); whether to re-run with explicit `SR` is D-010 | Simpler; inspection is a human decision step |
| Containers / per-module conda envs | One pinned conda env on PATH (`stv-rnaseq`); same tool versions as nf-core modules | Docker not available in WSL; single lockable env |
| `publishDir` via `conf/modules.config` | Workflow `output {}` block | Nextflow docs: workflow outputs are intended to replace `publishDir` |
| tximport/DESeq2 QC inside the pipeline | Done in R scripts (`scripts/R/`, env `stv-r`) | Keeps statistics outside the workflow; easier to review |
| Salmon index not published unless requested (`--save_reference`) | Index always published to `results/genome/index/` when built | Re-use with `--salmon_index` without rebuilding |
| Versions via topic channel → MultiQC | Same topic-channel pattern → `pipeline_info/software_versions.tsv` | — |
