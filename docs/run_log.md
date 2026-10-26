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

## Phase 1 — not executed (stop conditions)

FASTQ download, reference download, any software installation, trimming, Salmon/STAR/HISAT2, tximport, DESeq2, biological interpretation, GitHub push.

---

# Phase 2 (2026-10-25) — approvals D-002…D-007, D-008 partial; D-001 not approved

All commands below ran in WSL Ubuntu in `~/projects/STV_GSE137303` (launched from PowerShell via `wsl.exe`). Because VS Code is attached to the Windows folder and cannot open `\\wsl.localhost` paths, new/edited files were authored in a Windows staging folder (`C:\Users\Cutza\.stv_staging`, outside both repositories) and copied into the WSL repo with CRLF→LF conversion; helper/ad-hoc scripts lived in `C:\Users\Cutza\.stv_tools` (not part of the project).

## 2A — native WSL working copy

18. `git clone /mnt/c/Users/Cutza/Projects/Proyecto\ RNASeq\ Gen*mica ~/projects/STV_GSE137303` → `a12dd22` on `main`, working tree clean, 36 tracked files, `pwd` = `/home/cutza/projects/STV_GSE137303`. (First attempt failed: PowerShell stripped the inner quotes → `fatal: repository '/mnt/c/Users/Cutza/Projects/Proyecto' does not exist`; retried with escaped spaces.) Windows copy untouched.

## Sources

19. Read current docs (Nextflow install/updating/strict syntax/typed params/workflow outputs/cache-resume/26.10 migration, Nextflow RNAseq training, nf-core/rnaseq 3.27.0, nf-core modules fastqc/multiqc/salmon) and GitHub release lists (Nextflow, Salmon, nf-core/rnaseq) and anaconda.org versions → `docs/pipeline_sources.md`.

## 2B — software

20. `bash scripts/bash/env_audit.sh` from the WSL path (same values as Phase 1).
21. `CONDA_CHANNEL_PRIORITY=strict conda env create -f envs/stv-rnaseq.yml` → env `stv-rnaseq`; `conda list --explicit --md5 > envs/stv-rnaseq.linux-64.lock.txt` (152 packages). Global conda config verified unchanged (`channel_priority: flexible`, `channels: conda-forge`); `rachis-qiime2-2026.4` untouched.

## 2C — Nextflow

22. `bash scripts/bash/install_nextflow.sh` (inside `stv-rnaseq`): official self-install with `NXF_VER=26.10.3`, launcher moved to the env's `bin/`, env-scoped `NXF_VER`/`NXF_OPTS`. `nextflow -version` → 26.10.3 build 14229. A cosmetic non-zero exit was reported once by PowerShell because conda printed a stderr warning ("overwriting environment variables"); `nextflow info` exits 0.
23. Bootstrap workflow: `nextflow lint` clean; run 1 → 4 tasks submitted; run 2 `-resume` → 4 tasks cached (`docs/nextflow_validation.md`).

## 2D/2E — FASTQ acquisition

24. Throughput probe (20 s each, data discarded) on a second ENA connection while the main download ran: FTP 11.6 kB/s, HTTPS 16.8 kB/s → bandwidth already saturated by one transfer; downloads kept sequential.
25. `bash scripts/bash/download_fastq.sh` started 16:56 local time (free disk before: 939 GB on `/`). Observed rate ≈ 0.5–0.8 MB/s.

## 2J — workflow skeleton (built while the download ran)

26. `workflow/main.nf`, `workflow/modules/*.nf`, `workflow/nextflow.config`; `nextflow lint` → 8 files, no errors.
27. Stub wiring test (`-profile test_stub -stub`) → exit 0, 6 tasks incl. `SALMON_INDEX`/`SALMON_QUANT` stubs; libtype mapping `auto→A`, `reverse→SR` confirmed.
28. Real FastQC + MultiQC on synthetic fixtures: first run failed in `MULTIQC` (MultiQC 1.35 writes `multiqc_report_data/`, not `multiqc_data/`); fixed output glob `*_data`; re-run exit 0; `-resume` → all 3 tasks cached.
29. Synthetic fixture FASTQ/FASTA are regenerated by `workflow/tests/stub/make_fixtures.sh` and git-ignored (no FASTQ committed).
30. Checkpoint commit `b9d03cc` in the WSL repo with repository-local identity (same as the Windows identity used for `a12dd22`); pre-commit checks: no staged file > 500 KB, no FASTQ/FASTA/BAM/secrets, `work/`, `.nextflow/`, `data/raw/` ignored.

## 2E — download result

31. `download_fastq.sh` finished at 19:34 local time (9,489 s, exit 0): **6/6 PASS** (size + ENA MD5). Repeated `curl: (56) OpenSSL SSL_read: unexpected eof while reading` from ENA were recovered by `-C -` resume. For SRR10102795/SRR10102796 the first cycle's 5 resumes were exhausted and the second cycle continued the same `.part`. No `.corrupt` files. Disk: `/` 939 GB → 934 GB free; `data/raw` 5.4 GiB (`docs/raw_data_validation.md`).

## 2F — sanity checks

32. `bash scripts/bash/check_raw_fastq.sh` (116 s): gzip OK, structure OK (10,000 records), SeqKit read counts = ENA for all six, all reads 75 nt → `results/qc/raw_fastq_sanity.tsv`.
33. Read-header tabulation (first 200,000 reads/file): all from `SN679:HN3JNBCXY` lane 1.

## 2G/2H — FastQC + MultiQC (through the workflow)

34. `nextflow run workflow/main.nf` (run `kickass_poitras`, 86 s, exit 0): `FASTQC ×6` (1 CPU / 1 GB each, peak RSS ≈ 0.31–0.34 GB, 56–75 s), `MULTIQC` (4.8 s). All six samples present in MultiQC. `-resume` (`wise_hodgkin`): 7/7 cached.
35. Order change vs the instruction list: the workflow skeleton (2J) was built before 2G so that raw QC on real data ran once, through the reproducible workflow, instead of ad-hoc commands followed by a second run.

## 2I — QC interpretation

36. `bash scripts/bash/summarize_fastqc.sh` → `results/qc/raw_qc_summary.tsv`, `results/qc/raw_overrepresented.tsv`; GC-distribution shape and base-composition checks computed read-only from the FastQC ZIPs.
37. `docs/raw_qc_assessment.md`: conclusion **A — no preprocessing appears necessary**; logged as D-009 (Proposed). No FASTQ was modified.

## Phase 2 — not executed (stop conditions)

Reference download, Salmon index, Salmon on real reads, tximport, DESeq2, trimming/filtering, DEG thresholds, biological interpretation, GitHub push.
