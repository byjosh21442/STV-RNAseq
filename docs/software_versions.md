# Software versions

Recorded 2026-10-25 on WSL2 Ubuntu 26.04 (`~/projects/STV_GSE137303`); updated in Phase 3 (R environment, Salmon used on real data). No global conda configuration, no other conda environment (`base`, `rachis-qiime2-2026.4`) and no shell start-up file was modified.

## Project environment `stv-rnaseq`

- Specification: [`envs/stv-rnaseq.yml`](../envs/stv-rnaseq.yml) (pinned versions, channels `conda-forge` + `bioconda`, `nodefaults`).
- Exact resolved lock (152 packages, linux-64, with MD5): [`envs/stv-rnaseq.linux-64.lock.txt`](../envs/stv-rnaseq.linux-64.lock.txt).
- Created with `CONDA_CHANNEL_PRIORITY=strict conda env create -f envs/stv-rnaseq.yml` — strict priority applied to that command only via environment variable; global `channel_priority` stays `flexible`.
- Rebuild exactly: `conda create --name stv-rnaseq --file envs/stv-rnaseq.linux-64.lock.txt`.
- Size on disk: 1.8 GB.

| Software | Version (build) | Source / channel | Purpose | Status in Phase 2 |
|---|---|---|---|---|
| OpenJDK | 21.0.10 (`h1602c4f_22`) | conda-forge | Java runtime for Nextflow (needs 17–26) and FastQC | used |
| FastQC | 0.13.0 (`hdfd78af_0`) | bioconda | Raw read QC (same version as nf-core/modules) | used |
| MultiQC | 1.35 (`pyhdfd78af_2`) | bioconda | QC aggregation (same version as nf-core/modules) | used |
| Salmon | 2.7.0 (`hfa8f182_0`) | bioconda | Decoy-aware index + quantification (same version as nf-core/modules) | used (Phase 3) |
| SeqKit | 2.14.0 (`hb192632_0`) | bioconda | Read counts / length stats for sanity checks | used |
| pigz | 2.8 (`h421ea60_2`) | conda-forge | gzip integrity tests | used |
| Python | 3.15.1 | conda-forge | MultiQC runtime; FastQC/Salmon summary parsers, reference validation, tx2gene | used |

## Nextflow

| Item | Value |
|---|---|
| Version | **26.10.3 build 14229** (latest stable on 2026-10-25) |
| Method | Official self-install (`curl -s https://get.nextflow.io \| bash`) via [`scripts/bash/install_nextflow.sh`](../scripts/bash/install_nextflow.sh) |
| Launcher location | `~/miniconda3/envs/stv-rnaseq/bin/nextflow` (a directory already on PATH when the env is active — option allowed by the install docs) |
| Framework cache | `~/.nextflow/framework/26.10.3/` |
| Pinning | `NXF_VER=26.10.3` and `NXF_OPTS=-Xms1g -Xmx4g` stored as env-scoped variables (`conda env config vars`), active only inside `stv-rnaseq` |
| Runtime | Groovy 4.0.33 on OpenJDK 21.0.10 |
| Syntax | Strict parser (default since 26.04); `nextflow lint` passes on all workflow files |

Differences from the Phase 1 plan: Java comes from conda-forge inside the project env instead of SDKMAN (keeps Java isolated with the rest of the toolchain and avoids editing `~/.bashrc`); Nextflow is still installed with the official self-install script.

## System tools used (not installed by the project)

| Tool | Version | Purpose |
|---|---|---|
| curl | 8.18.0 | ENA HTTPS download with resume (`-C -`) |
| md5sum | uutils coreutils 0.8.0 | MD5 validation |
| gzip / zcat | 1.14 | Reading the first FASTQ records |
| git | 2.53.0 | Version control (repository-local identity, see `docs/run_log.md`) |

## R environment `stv-r` (Phase 3)

- Specification: [`envs/stv-r.yml`](../envs/stv-r.yml); exact lock (139 packages): [`envs/stv-r.linux-64.lock.txt`](../envs/stv-r.linux-64.lock.txt).
- Created with `CONDA_CHANNEL_PRIORITY=strict conda env create -f envs/stv-r.yml` (180 s; 1.7 GB). Rebuild exactly: `conda create --name stv-r --file envs/stv-r.linux-64.lock.txt`.
- Bioconductor 3.22 builds (newest on bioconda on 2026-10-26; Bioconductor 3.23 is current upstream).

| Software | Version | Source / channel | Purpose | Status |
|---|---|---|---|---|
| R | 4.5.3 (2026-03-11) | conda-forge | Runtime | used |
| tximport | 1.38.2 | bioconda | Gene-level import of Salmon estimates | used (validation/import only) |
| DESeq2 | 1.50.2 | bioconda | Differential expression | installed, **not run** |
| jsonlite | 2.0.0 | conda-forge | Salmon JSON metadata | used |

## Not installed (by design)

STAR, HISAT2, sra-tools, trimming tools, Docker integration, Apptainer.

## Activate

```bash
conda activate stv-rnaseq
nextflow -version   # 26.10.3
conda activate stv-r  # R scripts in scripts/R/
```
