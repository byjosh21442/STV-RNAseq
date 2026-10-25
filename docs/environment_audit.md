# Environment audit

- Date (UTC): 2026-10-25
- Method: read-only commands; nothing installed or modified. Re-run with `bash scripts/bash/env_audit.sh` from the project root inside WSL.
- Raw output is summarized below (not committed as a log).

## Phase 2 update (2026-10-25)

- Active project moved (D-006): `git clone` of the Phase 1 repository to **`/home/cutza/projects/STV_GSE137303`** (ext4). Branch `main`, history intact (`a12dd22`); the Windows copy `/mnt/c/Users/Cutza/Projects/Proyecto RNASeq Genómica` is unchanged and now behind (its `origin` role: local remote of the WSL clone, nothing pushed).
- Re-audit from the WSL path: same host values as below (15 GiB RAM, 20 threads, 942 GB free on `/` before downloads).
- Installed: conda env `stv-rnaseq` (OpenJDK 21.0.10, FastQC 0.13.0, MultiQC 1.35, Salmon 2.7.0, SeqKit 2.14.0, pigz 2.8) and Nextflow 26.10.3 — details in `docs/software_versions.md`. Global conda config and other envs unchanged.
- WSL git: repository-local `user.name`/`user.email` set in the WSL clone only (same identity as the Windows git used for the Phase 1 commit); global git config untouched.
- VS Code is still attached to the Windows folder; files for the WSL repo are written through a Windows staging folder and copied in with CRLF→LF conversion (see `docs/run_log.md`). Opening `~/projects/STV_GSE137303` with VS Code Remote-WSL would remove this indirection.
- Sections below are the original Phase 1 audit (kept for provenance; "MISSING" entries for Java/Nextflow/FastQC/MultiQC/Salmon are now resolved).

## Host / WSL

| Item | Value | Status |
|---|---|---|
| Windows | 10.0.26200.9711 (from `wsl --version`) | OK |
| WSL | 2.7.11.0, default version 2 | OK |
| WSL kernel | 6.18.33.2-microsoft-standard-WSL2 | OK |
| Distributions | `Ubuntu` (default, running, WSL2); `docker-desktop` (stopped) | OK |
| Linux distribution | Ubuntu 26.04 LTS (resolute) | OK |
| CPU | Intel Core 9 270H, x86_64, 10 cores / 20 threads (`nproc` = 20) | OK |
| RAM visible to WSL | 15 GiB total, ~13 GiB available; swap 4 GiB | OK for Salmon; tight for STAR genome indexing (see `pipeline_design.md`) |
| Disk, WSL native `/` | 1007 G, 942 G free | OK |
| Disk, `/mnt/c` | 953 G, 755 G free | OK |

## Tools (inside WSL, interactive shell so `~/.bashrc` conda init applies)

| Tool | Status | Version / path |
|---|---|---|
| git | FOUND | 2.53.0 (`/usr/bin/git`). No global `user.name`/`user.email` in WSL |
| git (Windows host) | FOUND | 2.55.0.windows.3, identity configured |
| conda | FOUND | 26.5.3 (`~/miniconda3`), solver `libmamba`, channels: `conda-forge` only, `channel_priority: flexible`, `base` auto-activated |
| conda envs | — | `base`, `rachis-qiime2-2026.4` (unrelated; not touched) |
| mamba / micromamba | MISSING | conda's libmamba solver is available |
| Java | MISSING | required by Nextflow and FastQC |
| Nextflow | MISSING | — |
| R / Rscript | MISSING | needed later for tximport/DESeq2 |
| Python | FOUND | 3.14.6 (conda `base`) — do not install project tools into `base` |
| Docker | NOT USABLE in WSL | Docker Desktop CLI is on the Windows PATH but reports *"could not be found in this WSL 2 distro … activate the WSL integration in Docker Desktop settings"* |
| Apptainer / Singularity | MISSING | — |
| curl / wget | FOUND | curl 8.18.0, GNU Wget 1.25.0 |
| pigz, sra-tools (`prefetch`, `fasterq-dump`), FastQC, MultiQC, Salmon, STAR, HISAT2, featureCounts | MISSING | expected at this phase |

## Flags

### 1. Project directory is on the Windows filesystem (`/mnt/c`) — RECOMMEND MOVING

Current path (WSL view): `/mnt/c/Users/Cutza/Projects/Proyecto RNASeq Genómica`

Problems for RNA-seq/Nextflow work:

- `/mnt/c` is a drvfs/9P mount: much slower I/O than ext4 for many small files (Nextflow `work/`, Salmon indexes, conda envs).
- All files appear as `-rwxrwxrwx`; permission and file-locking semantics differ from Linux.
- The path contains **spaces** and a **non-ASCII character** (`ó`). Nextflow task scripts, some bioinformatics tools and R packages are known to be fragile with such paths.
- Windows tools write CRLF line endings (observed here: files created from the Windows side had CRLF and had to be normalized before `bash` could run them). `.gitattributes` now enforces LF in the repository.

**Recommendation:** keep this repository as the source of truth, but make the *active* working copy a clone on the native WSL filesystem, e.g.

```bash
mkdir -p ~/projects
git clone "/mnt/c/Users/Cutza/Projects/Proyecto RNASeq Genómica" ~/projects/STV_GSE137303
```

FASTQ files, references, indexes and Nextflow `work/` must live on the WSL filesystem. Nothing was moved or deleted automatically.

### 2. VS Code is not attached to WSL

The workspace is opened through a Windows path and the integrated terminal is Windows PowerShell 5.1; WSL commands were executed through `wsl.exe`. After cloning to `~/projects/STV_GSE137303`, open it from WSL (`code ~/projects/STV_GSE137303`) so that terminals and tools run natively in Ubuntu.

### 3. WSL git has no identity

Commits from WSL need `git config user.name/user.email` (repository-local is sufficient). The initial commit of this phase was made with the Windows git identity.

## Nextflow: not installed — recommended installation (NOT executed)

Source: Nextflow installation docs, <https://docs.seqera.io/nextflow/install> (redirected from nextflow.io/docs/latest/install.html), retrieved 2026-10-25.

- Requirements: Bash ≥ 3.2 and **Java 17 or later (up to 26)**. Java < 17 is unsupported since Nextflow 25.04.
- Officially recommended Java route: SDKMAN!, latest LTS Temurin (user-level, no sudo):
  ```bash
  curl -s https://get.sdkman.io | bash   # then open a new terminal
  sdk install java
  java -version
  ```
- Officially recommended Nextflow route: self-installing package into a user directory:
  ```bash
  curl -s https://get.nextflow.io | bash
  chmod +x nextflow
  mkdir -p "$HOME/.local/bin" && mv nextflow "$HOME/.local/bin/"
  nextflow info
  ```
  Pin the version used by the project (e.g. `export NXF_VER=<version>`) and record it in `docs/run_log.md`.
- Alternative documented by Nextflow: `conda create --name nf-env bioconda::nextflow` — the docs warn this may give outdated versions and Java conflicts.
- Per-tool software (FastQC, MultiQC, Salmon) should come from isolated, versioned environments: per-process conda environments (requires `bioconda` + `conda-forge` channels specified in the env files, without changing the global conda config) **or** containers (requires enabling Docker Desktop WSL integration — a user action in Docker Desktop settings).

Status: **next action, pending approval** (substantial dependency installation is a stop condition).
