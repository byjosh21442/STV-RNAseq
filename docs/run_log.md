# Run log

Everything actually executed in this phase, in order. Host: Windows 10.0.26200 + WSL2 Ubuntu 26.04; commands were launched from Windows PowerShell 5.1 through `wsl.exe` (VS Code was not attached to WSL). Date: 2026-10-25.

## Phase 1 — environment audit (read-only)

1. `wsl.exe --status; wsl.exe -l -v; wsl.exe --version` → WSL 2.7.11.0, Ubuntu (default, WSL2), docker-desktop (stopped).
2. `bash scripts/bash/env_audit.sh` (via `wsl.exe -d Ubuntu --cd <project> -e bash -ic …`) → results in `docs/environment_audit.md`.
3. `conda env list`, `conda config --show channels channel_priority solver`, `git config --global --get user.name|user.email` (WSL and Windows), `ls ~/projects` → read-only.

## Phase 2 — dataset and GEO file audit

4. `bash scripts/bash/fetch_metadata.sh` → GEO SOFT, NCBI SRA runinfo + XML (6 UIDs), ENA read_run TSV, GEO suppl listing into `data/metadata/sources/` (retrieved 2026-10-25T20:12:59Z).
5. `bash scripts/bash/fetch_geo_counts.sh` → downloaded `GSE137303_count.txt.gz` (1,344,528 bytes, MD5 `3800ddd08f44eba1f062f034cce4d329`).
   - First run exited 1 after printing the header: `set -o pipefail` + `zcat | head` SIGPIPE. Fixed by using `set -eu`; re-run succeeded (file not re-downloaded).
6. PowerShell `[xml]` parsing of the SRA XML → per-run layout, reads-per-spot and original submitted file names (`12.fastq` … `433.fastq`).
7. `bash scripts/bash/build_samples_metadata.sh` → `data/metadata/samples_metadata.tsv`, `config/samplesheet.csv`.
   - First run produced empty GEO fields: GEO SOFT is served with CRLF, so GSM keys did not join. Fixed by stripping `\r` in awk; consistency checks changed to hard failures. Re-run: both cross-source checks OK.
8. `bash scripts/bash/verify_geo_counts_annotation.sh` → streamed NCBI `GCF_000188115.3_SL2.50_genomic.gff.gz` (≈8 MB, not saved): 30,389/30,389 IDs match, 30,292 exact (ID, GeneID) pairs.
9. One ad-hoc `wsl.exe … bash -c "zcat … | grep -F '\"GeneID'"` hung on PowerShell quoting; terminal killed; repeated with a PowerShell `GZipStream` read (70 rows with `"GeneID:…,miRBase:…"`).

## Phase 5/6 — reference and pipeline research (no large downloads)

10. `curl` directory listings / READMEs: SGN `tomato_genome/annotation/` (ITAG2.4, 3.2, 4.0, 4.1), `tomato_genome/assembly/` (build_4.00 = current_build), `genomes/Solanum_lycopersicum/` (Heinz1706, M82, id_conversion); NCBI `genomes/all/GCF/000/188/115/`; Ensembl Plants `current/fasta|gtf/solanum_lycopersicum_gca000188115v5cm/`.
11. NCBI Datasets API `genome/accession/…/dataset_report` and `genome/taxon/4081/dataset_report?filters.assembly_source=refseq&filters.assembly_version=all_assemblies` (JSON).
12. Streamed first lines only (nothing saved) of the Ensembl Plants SL4.0 GTF and of the NCBI SL2.50 GFF to read their headers.
13. Web documentation read: Nextflow install, nf-core/rnaseq 3.27.0 usage, FastQC, MultiQC, Salmon (library types, index/quant), tximport, DESeq2, Cuffnorm, agriGO v2.
    - Failed/unavailable (recorded, not worked around): Ensembl REST `info/genomes|assembly/solanum_lycopersicum` (HTTP 400 / "Can not find internal name"); Ensembl Plants species pages (redirect/HTTP 500); Mercator4 site (Service Unavailable).

## Phase 3/4/7/8 — repository

14. Created directories (`workflow/`, `scripts/R/`, `data/raw/`, `reference/`, `results/{qc,quantification,deseq2,functional}/`, `figures/`, `envs/`, `docs/`, `config/`) with `.gitkeep` in empty ones; existing content untouched.
15. Files written from VS Code on Windows get CRLF; each new `*.sh` was normalized with `sed -i 's/\r$//'` (only files created in this session). `.gitattributes` enforces LF in the repository; `data/metadata/sources/**` is stored byte-exact.
16. Nextflow is **not installed** → no `workflow/main.nf` / `nextflow.config` created; installation documented in `docs/environment_audit.md` (Phase 7 fallback path).
17. Git: `git init` + local commit with the Windows git identity (WSL git has no identity configured). No remote created, nothing pushed.

## Not executed (stop conditions)

FASTQ download, reference download, any software installation, trimming, Salmon/STAR/HISAT2, tximport, DESeq2, biological interpretation, GitHub push.
