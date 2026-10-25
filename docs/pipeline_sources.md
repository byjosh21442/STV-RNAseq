# Pipeline reference sources

All accessed **2026-10-25**. These are the authoritative references used to design and implement `workflow/`. Where our implementation differs, the reason is recorded in `docs/nextflow_pipeline.md` ("Differences from nf-core/rnaseq").

## 1. nf-core/rnaseq (principal architecture reference — read, not run)

| What | URL | Used for |
|---|---|---|
| Release 3.27.0 (latest, published 2026-09-23) | <https://github.com/nf-core/rnaseq/releases> (GitHub API `releases/latest`) | Version context |
| Usage docs 3.27.0 | <https://nf-co.re/rnaseq/3.27.0/docs/usage> | Samplesheet columns, `strandedness=auto` inference (Salmon on 1 M reads, thresholds 0.8 / 0.1), Salmon decoy-aware index, `--pseudo_aligner salmon --skip_alignment`, featureCounts not used for quantification since v3.0 |
| `nextflow.config` (master) | <https://raw.githubusercontent.com/nf-core/rnaseq/master/nextflow.config> | `resourceLimits`, `outputDir`, `manifest.nextflowVersion '!>=25.10.4'` |
| nf-core/modules `fastqc/main.nf`, `environment.yml` | <https://github.com/nf-core/modules/tree/master/modules/nf-core/fastqc> | Sample-ID renaming via symlink, `--memory` = task memory / cpus, eval+topic versions; pins fastqc 0.13.0 |
| nf-core/modules `multiqc/main.nf`, `environment.yml` | <https://github.com/nf-core/modules/tree/master/modules/nf-core/multiqc> | `*multiqc_report.html` / `*_data` outputs; pins multiqc 1.35 |
| nf-core/modules `salmon/index`, `salmon/quant` | <https://github.com/nf-core/modules/tree/master/modules/nf-core/salmon> | Genome-decoy index, `-r` for single-end, `lib_format_counts.json`, stub blocks; pins salmon 2.7.0 |

## 2. Nextflow (official documentation and training)

| What | URL | Used for |
|---|---|---|
| Installation | <https://docs.seqera.io/nextflow/install> | Java 17–26 requirement; self-install `curl -s https://get.nextflow.io \| bash`; executable placed in a directory already on PATH |
| Updating Nextflow | <https://docs.seqera.io/nextflow/updating-nextflow> | Stable releases in April/October; `NXF_VER` pinning |
| Releases | <https://github.com/nextflow-io/nextflow/releases> | 26.10.3 = latest stable (2026-10-20) |
| Migrating to 26.10 | <https://docs.seqera.io/nextflow/migrations/26-10> | `executor.cpus/memory` semantics, console output on stderr, `report/timeline/trace.directory`, `when:` deprecated |
| Strict syntax | <https://docs.seqera.io/nextflow/strict-syntax> | Strict parser is the default since 26.04: no `for` loops / top-level statements, `channel` (lower case), explicit closure params, CLI params not type-converted for legacy params |
| Typed parameters | <https://docs.seqera.io/nextflow/typed-parameters> | `params { name: Type = default }` block (no feature flag needed) |
| Language syntax reference | <https://docs.seqera.io/nextflow/reference/syntax> | Params block, output block, functions |
| Workflows / workflow outputs | <https://docs.seqera.io/nextflow/workflow> | `publish:` + `output {}` (replaces `publishDir`), index files, `outputDir` |
| Caching and resuming | <https://docs.seqera.io/nextflow/cache-and-resume> | Task hash, `.nextflow/cache`, `-resume` semantics |
| Training — Nextflow for RNAseq (training 3.8.2) | <https://training.nextflow.io/latest/nf4_science/rnaseq/> , <https://training.nextflow.io/latest/nf4_science/rnaseq/02_single-sample/> | Module layout (`modules/*.nf`), `params` block, `publish:`/`output {}`, `test` profile, FastQC → MultiQC aggregation |
| Example pipeline rnaseq-nf | <https://github.com/nextflow-io/rnaseq-nf> | Minimal `index / quant / fastqc / multiqc` module structure |

## 3. Salmon

| What | URL | Used for |
|---|---|---|
| Salmon docs ("latest", page labelled 1.12.0) | <https://salmon.readthedocs.io/en/latest/salmon.html> | Decoy-aware index, `-l` before `-r`, single-end `--fldMean/--fldSD`, `-l A` |
| Library types | <https://salmon.readthedocs.io/en/latest/library_type.html> | Single-end types `U`, `SF`, `SR` |
| Releases | <https://github.com/COMBINE-lab/salmon/releases> | Latest release v2.8.0 (2026-09-25); nf-core pins 2.7.0 |

Note: the readthedocs "latest" page still reports 1.12.0 while bioconda/GitHub ship 2.x. Before any real quantification, options must be re-checked against `salmon quant --help` of the installed 2.7.0.

## 4. tximport and DESeq2 (later phases)

| What | URL | Used for |
|---|---|---|
| tximport vignette (1.40.0, Bioconductor 3.23) | <https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html> | `tximport(files, type = "salmon", tx2gene)`; workflow publishes `quantification/salmon/samples.csv` index for it |
| DESeq2 vignette (1.52.0) | <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html> | Salmon → tximport → `DESeqDataSetFromTximport` is the recommended path |

## 5. QC tools

| What | URL |
|---|---|
| FastQC | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/> (bioconda 0.13.0) |
| MultiQC | <https://docs.seqera.io/multiqc> (bioconda 1.35) |
| ENA file access | FASTQ URLs = `https://` + ENA `fastq_ftp` field from the ENA portal API record saved in Phase 1 |
