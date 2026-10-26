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

# Phase 3 — reference, Salmon quantification, tximport preparation, publication

Owner approvals at the start of Phase 3: D-001 = SL4.0 + ITAG4.0 (ITAG4.1 considered, not selected), D-009 (no trimming). Times are local (UTC−6); the reference manifest records UTC dates (2026-10-26).

## 3A — reference acquisition

38. SGN listings inspected (`assembly/build_4.00/`, `annotation/ITAG4.0_release/`): no checksum files are published.
39. `bash scripts/bash/download_reference.sh` (111 s, exit 0): 9 files, 335 MB, sizes = server `Content-Length`, `gzip -t` OK, SHA-256 → `reference/reference_manifest.tsv`. CDS, lncRNA and repeat GFFs were deliberately not downloaded.

## 3B — reference validation

40. Format inspection: cDNA headers = GFF mRNA `Name`; GFF `ID` = `gene:`/`mRNA:` + Name; genome `SL4.0ch00`–`12`.
41. `bash scripts/bash/validate_reference.sh`, first run: exit 1. My strict `Solyc##g######.v(.t)` name pattern flagged 15 genes with an extra version level. Inspection showed these are genuine upstream names (e.g. `Solyc01g008710.5.1`), unique and consistent with the cDNA. The check was refined (unparseable names and multi-model loci stay material; the extra version level became a NOTE). Re-run: exit 0, no material inconsistency (`docs/reference_validation.md`). `ITAG4.0_goterms.txt` version mismatches were recorded as a NOTE (GO not used now).

## 3C — Salmon documentation and workflow update

42. `salmon index --help`, `salmon quant --help` (2.7.0); Salmon 2.x docs, MIGRATION.md and release notes v2.7.0 read; SRA RunInfo checked for insert size (`InsertSize = 0` for all runs).
43. Workflow changes: `SALMON_INDEX` emits and publishes `decoys.txt` + `salmon_index.log` (tee with `pipefail`) and the index under `results/genome/index/`; typed params `salmon_kmer`, `salmon_fld_mean`, `salmon_fld_sd`, `salmon_index_dir`; `ext.args` for index (`-k 31 --ramLimit 6`) and quant (`--fldMean 250 --fldSD 25`). `nextflow lint workflow/` → 8 files, no errors; stub test (`irreverent_monod`) → exit 0, index/decoys/log published.

## 3D — index + quantification

44. `nextflow run workflow/main.nf -resume --run_quant true --transcript_fasta … --genome_fasta … --multiqc_title 'STV_GSE137303 QC + Salmon (SL4.0 / ITAG4.0)' --multiqc_dir qc/multiqc_salmon` (run `exotic_ardinghelli`, 5 min 59 s, exit 0). FastQC re-ran because `-resume` picks the latest session, which was the stub test; outputs are identical in content to Phase 2 and the Phase 2 MultiQC report was not overwritten. `SALMON_INDEX` 3 min 29 s, peak RSS 4.7 GB, 0 warnings, 99 duplicates collapsed, 19 poly-A tails clipped. `SALMON_QUANT` ×6: 20–26 s each, peak RSS 2.2–2.3 GB, 0 warnings. Log: `logs/phase3_salmon_run.log` (git-ignored).
45. `bash scripts/bash/summarize_salmon.sh` → `results/quantification/salmon_qc_summary.tsv`, `salmon_libformat_summary.tsv`: all six `SR`, mapping 87.95–90.26 %, compatible ratio ≥ 0.9924, no diagnostics, `ISSUES: none`.
46. Library-type check: `--salmon_index results/genome/index/salmon_index --salmon_libtype SR -output-dir /tmp/stv_sr_check` (run `sick_archimedes`, 1 min 7 s; `results/` untouched). Outputs differ from `-l A` for 621–774 transcripts per sample with identical mapping counts. A direct `salmon quant -l A` rerun of STV_free_rep1 with 6 and 3 threads was byte-identical to the pipeline output (deterministic). Opposite-strand-overlap test: no enrichment. `bash scripts/bash/compare_salmon_runs.sh … results/quantification/libtype_A_vs_SR_comparison.tsv`. Logged as D-010 (owner decision).

## 3E — tx2gene and R environment

47. `bash scripts/bash/build_tx2gene.sh` → `reference/ITAG4.0_tx2gene.tsv` (34,075 rows); validation vs cDNA, index duplicates and all six `quant.sf`: exit 0 (0 unexplained mismatches; 99 absent IDs = collapsed duplicates).
48. `conda search` for DESeq2/tximport/r-base (bioconda newest: DESeq2 1.50.2, tximport 1.38.2, R 4.5). `CONDA_CHANNEL_PRIORITY=strict conda env create -f envs/stv-r.yml` (180 s) → R 4.5.3, tximport 1.38.2, DESeq2 1.50.2; lock `envs/stv-r.linux-64.lock.txt` (139 packages). Global `channel_priority` stayed `flexible`; `rachis-qiime2-2026.4` untouched.
49. `Rscript scripts/R/01_tximport_import.R` (env `stv-r`): all checks OK, gene matrix 33,976 × 6, totals = Salmon totals → `results/quantification/tximport_validation.tsv`, `tximport_sessionInfo.txt`. No DESeq2 call.

## 3F — references and documentation

50. DOIs verified via the Crossref API (14 + 2 articles/preprints) and DataCite (nf-core/rnaseq concept DOI); original study metadata via NCBI E-utilities (PMID 31620899 → Fukuhara et al. 2020, *Arch Virol*). Versioned Bioconductor 3.22 vignettes, ENA download docs and SeqKit docs checked.
51. New docs: `docs/reference_validation.md`, `docs/salmon_quantification.md`, `docs/tximport_preparation.md`, `docs/methods_references.md`, `reference/README.md`. Updated: `decisions.md` (D-001/D-009 approved, D-004 outcome, D-010–D-012), `reference_decision.md`, `pipeline_design.md`, `pipeline_sources.md`, `software_versions.md`, `nextflow_pipeline.md`, `README.md` (rewritten).

## 3G — publication audit

52. Pre-commit audit:
    - 34 paths to add/modify; largest new file `reference/ITAG4.0_tx2gene.tsv` (1.2 MB); only other tracked file > 1 MB is the Phase 1 GEO table (1.3 MB).
    - No FASTQ/FASTA/GFF/BAM/index/`work/` paths are tracked.
    - Secret-pattern scan (tokens, keys, passwords): no hits.
    - Ignore rules confirmed for `reference/SGN/`, `results/genome/`, `results/quantification/salmon/`, `results/qc/multiqc_salmon/`, `logs/`, `data/raw/`, `work/`.
    - Non-secret personal data already in history: the GEO contact e-mail inside the archived SOFT record, local paths with the user name in audit docs and in `sessionInfo`, and the repo-local commit e-mail. Left as is.
53. GitHub: `gh` 2.102.0 is installed and authenticated on Windows (account `byjosh21442`, keyring) and absent in WSL. Remote creation and push follow in the commit/push step (hash and verification in the Phase 3 report).
54. `gh repo create STV-RNAseq --public --description "…"` (no README/.gitignore/license). `origin` changed from the Windows copy to `https://github.com/byjosh21442/STV-RNAseq.git`. `git push -u origin main` used Windows `gh.exe` as a one-off credential helper; nothing was stored in the repo config. Commit `a1ec6ae`; remote `main` = local HEAD.

## 3H — final quantification (owner approvals D-010: explicit `SR`; D-012: `--gcBias`, no `--seqBias`)

55. Configuration:
    - `build_samples_metadata.sh` now writes the measured strandedness (`reverse`, documented as the only non-record value). Re-run: `config/samplesheet.csv` and `data/metadata/samples_metadata.tsv` changed only in the strandedness column; the cross-source checks were still OK.
    - New typed param `salmon_gc_bias = true` adds `--gcBias` through `ext.args`.
    - `summarize_salmon.sh` takes `[salmon_dir] [suffix]` and reports the bias flags.
    - `01_tximport_import.R` checks the library type, bias flags and Salmon version from `meta_info.json`.
    - `nextflow lint` → no errors.
56. Archive (no deletion):
    - `results/quantification/salmon/` → `results/quantification/archive/salmon_libA_exploratory/`. Its `samples.csv` paths were rewritten to the archive location; the original is kept as `samples.original_paths.csv`. `quant.sf` SHA-256 values are unchanged from `a1ec6ae`.
    - `results/qc/multiqc_salmon/` → `results/qc/multiqc_salmon_libA_exploratory/`.
    - `tximport_validation.tsv` → `tximport_validation_libA_exploratory.tsv` (`git mv`).
    - The `/tmp` SR-without-GC check run was copied to `results/quantification/archive/salmon_SR_noGC_check/`.
57. `nextflow run workflow/main.nf -resume --run_quant true --salmon_index results/genome/index/salmon_index --multiqc_title '… final …' --multiqc_dir qc/multiqc_salmon` (run `sleepy_mcclintock`, 2 min 25 s, exit 0). Log: `logs/phase3_final_salmon_run.log`.
    - Rendered command per sample: `--libType SR … --fldMean 250 --fldSD 25 --gcBias`.
    - `SALMON_QUANT` ×6: 27–29 s each, 2.4 GB peak RSS, 0 warnings/errors.
    - FastQC re-ran because the sample metadata map (strandedness) is part of its input hash.
58. `bash scripts/bash/summarize_salmon.sh` (final): requested = detected = applied `SR` and `gc_bias_correct=true` in 6/6; mapping 87.95–90.26 %; no diagnostics; `ISSUES: none`.
59. `bash scripts/bash/summarize_salmon.sh results/quantification/archive/salmon_libA_exploratory _libA_exploratory`: the regenerated exploratory tables equal the tables committed in `a1ec6ae`, apart from the 2 new bias columns.
60. Comparisons (`compare_salmon_runs.sh`):
    - Exploratory vs final: identical mapping counts; 1,468–1,593 transcripts change `NumReads`; ~20k effective lengths change; log-scale Pearson 0.99915–0.99929.
    - SR-without-GC vs final: 1,424–1,562 transcripts change.
    - Effective-length checks: consistent shift across samples (median final/exploratory ratio 0.277–0.285). The 3,149 transcripts adjusted in only some samples carry 0.047 % of reads.
    - DESeq2 1.50.2 source checked: `avgTxLength` is divided by the per-gene geometric mean before use as a normalization factor.
    - `results/quantification/quant_sf_sha256.tsv` written for all three runs.
61. `bash scripts/bash/build_tx2gene.sh`: OK against the final `quant.sf`; tx2gene unchanged (only the sample order in the validation table changed).
62. `Rscript scripts/R/01_tximport_import.R` (stv-r) on the final run: all checks OK. 33,976 genes × 6; gene totals = Salmon totals (15,902,183; 14,698,030; 12,871,350; 17,121,521; 15,217,015; 17,871,112); 9,808 genes with zero counts in all samples. No DESeq2 call.
63. Love, Hogenesch & Irizarry (2016), DOI 10.1038/nbt.3682, verified via Crossref. Docs updated: `salmon_quantification.md`, `tximport_preparation.md`, `decisions.md` (D-004 completed, D-010/D-012 approved), `methods_references.md`, `nextflow_pipeline.md`, `pipeline_design.md`, `samplesheet_schema.md`, `README.md`.

## Phase 3 — not executed (stop conditions)

DESeq2 object construction or testing, filtering, DEG thresholds, contrasts, biological interpretation, functional enrichment, `--seqBias`.

## Phase 4 — DESeq2 primary analysis (owner approvals D-013 – D-019)

64. `stv-r` environment: `apeglm=1.32.0`, `r-pheatmap=1.0.13`, `r-ggrepel=0.9.8`, `r-ggplot2` and `r-rcolorbrewer` were added with `CONDA_CHANNEL_PRIORITY=strict … --override-channels -c conda-forge -c bioconda`. `envs/stv-r.yml` was updated and the lock regenerated (139 → 152 packages).
65. `scripts/R/01_tximport_import.R` now saves `results/quantification/txi_final.rds` (git-ignored, 2.2 MB) after all checks pass. It was re-run with the same validation results as #62.
66. `Rscript scripts/R/02_deseq2.R`:
    - The first run stopped on its own assertion: `results(name = …)` writes the coefficient description with spaces, not underscores.
    - The assertion now checks the description returned by `results(contrast = c("condition", "STV_infected", "STV_free"))`.
    - The second run finished with exit 0. Log: `results/deseq2/deseq2_run.log`.
67. Results:
    - 18,277 / 33,976 genes retained.
    - Size factors 0.848 – 1.212; parametric dispersion trend (asymptDisp 0.0835, extraPois 5.16, dispPriorVar 0.776); median final dispersion 0.088; 105 dispersion outliers; 0 non-converged genes.
    - 18,238 genes tested; independent filtering at quantile 0 %.
    - **padj < 0.05: 0; stringent: 0**; raw p < 0.05: 192 (descriptive only).
    - 39 genes Cook's-filtered (cutoff 18.0); apeglm converged for all genes.
    - VST PCA top 500: PC1 49.28 %, PC2 27.17 %. Hierarchical clustering k = 2 = {STV_free_rep1} vs the rest.
    - Outputs: `results/deseq2/`, `results/deseq2/plot_data/`, `figures/deseq2/` (8 PNG files); `results/deseq2/deseq2_objects.rds` (git-ignored, 12.4 MB).
68. `docs/decisions.md`: D-013 – D-019 recorded as implemented.
69. Owner approvals given at the Phase 5 review (2026-10-25):
    - D-020: pheatmap 1.0.13 and ggrepel 0.9.8 approved as visualization-only dependencies.
    - D-021: Git result policy. The four Phase 4 TSVs > 1 MB are tracked; FASTQ, indexes, Salmon outputs, `work/`, RDS and SAM files stay ignored.
    - Docs updated: `README.md` (status, primary result), `docs/software_versions.md`, `docs/tximport_preparation.md`, `docs/methods_references.md` (§8, §15, §17, §18), `docs/results_by_process.md` (new).
    - Phase 4 was committed separately from, and before, Phase 5; hashes are in the Phase 5 handoff report.
