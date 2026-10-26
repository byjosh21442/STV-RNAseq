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

## Phase 5 — sample-structure diagnostic checkpoint (primary analysis unchanged)

70. Metadata (NCBI E-utilities, retrieved 2026-10-26T04:39:41Z UTC):
    - BioSample and BioProject XML for PRJNA565137 → `data/metadata/sources/PRJNA565137_{biosamples,bioproject}.xml`.
    - BioSample attributes are only source_name, cultivar, age, virus infection and tissue. No batch, date, plant or pairing fields.
    - Fukuhara et al. 2020 full text: subscription-only (Springer redirect; no Europe PMC full text), so its methods are UNKNOWN.
71. `bash scripts/bash/fastq_header_audit.sh` (env `stv-rnaseq`; the first attempt from another environment failed because `pigz` was missing).
    - All six libraries: SN679, run 271, flowcell HN3JNBCXY, lanes 1 + 2, 64 tiles each.
    - Output: `results/qc/sample_structure/fastq_header_lanes.tsv`.
72. STV references:
    - LC429302.1 (strain M82, authors' deposit) and KX949574.1 (isolate Florida), fetched with `efetch` (GenBank flat files archived; FASTA in git-ignored `reference/STV/`).
    - The KX949574.1 publication was resolved through Crossref: 10.1128/genomea.01374-16.
73. `bash scripts/bash/stv_signal_check.sh` (env `stv-rnaseq`, Salmon 2.7.0):
    - A short test run first confirmed `-l U` and `--writeMappings` in Salmon 2.7.0.
    - Competitive index (ITAG4.0 cDNA + 2 STV genomes, k 31) at `reference/STV_diagnostic/` (1.1 GB, git-ignored).
    - Free libraries: 0 STV reads. Infected libraries: 386 / 580 / 304 STV reads, 97.0 – 98.2 % genome coverage. Labels supported; no STOP.
    - Expected warnings: `-l U` vs observed SR.
    - `.gitignore` gained `results/qc/stv_signal/salmon/`.
74. `rrcov` for robust PCA:
    - `conda install --dry-run` showed 4 new packages (r-rrcov 1.7-7, r-robustbase 0.99-7, r-pcapp 2.0-5, r-deoptimr 1.2-1) and no updates or removals.
    - Installed into `stv-r` only (strict priority for that command; global `channel_priority` remains `flexible`).
    - `envs/stv-r.yml` gained `r-rrcov=1.7_7`; lock regenerated (152 → 156 packages).
75. `Rscript scripts/R/03_sample_structure_diagnostics.R` (env `stv-r`): exit 0, 0 warnings. Log: `results/qc/sample_structure/sample_structure_run.log`.
    - Built-in checks: Phase 4 VST distances and top-500 PCA reproduced (tolerance 1e-8); 39 Cook's-filtered genes confirmed.
    - The script was re-run twice: once after adding the robust-PCA degeneracy columns, once for figure layout (clipped labels). Only `results/qc/sample_structure/` and `figures/sample_structure/` were rewritten.
76. Robust PCA:
    - `PcaGrid` picked k = 6 = n (orthogonal distances ≈ 1e-13); `PcaHubert` used h = 6 = n.
    - Assessed as degenerate at n = 6 → subtask stopped; no parameter tuning. The `PcaGrid` flag on STV_free_rep1 is not used as evidence.
77. `Rscript scripts/R/04_geo_historical_structure_check.R` (env `stv-r`): exit 0, 0 warnings.
    - GEO file MD5 matched the sidecar.
    - Outputs labelled "Historical processed-data structure check - not used for differential-expression inference".
    - STV_free_rep1 most remote in the FPKM and "count" columns.
78. Literature:
    - Full texts read: Chen et al. 2020 (PMC7324992), Gierliński et al. 2015 (PMC4754627), Schurch et al. 2016 (PMC4878611), Upton et al. 2023 (Frontiers), Khorramdelazad et al. 2018 (PMC5793396), DESeq2 vignette 3.22, rnaseqGene workflow, Badeli et al. 2026 (Research Square HTML).
    - Crossref checks: 10.12688/f1000research.7035.1, 10.18637/jss.v032.i03, 10.1198/004017004000000563, 10.1016/j.chemolab.2007.01.004, 10.21203/rs.3.rs-8823240/v1, 10.1093/bioinformatics/bty895, 10.1128/genomea.01374-16.
79. Docs:
    - New: `sample_structure_investigation.md`, `external_reanalysis_comparison.md`, `results_by_process.md`.
    - Updated: `methods_references.md` (§8 lock sizes, §15 Phase 4, new §17 – §24), `run_log.md`.
80. Owner decision D-022 (2026-10-25): **Option A approved.**
    - All six samples kept; the six-sample DESeq2 analysis remains primary; STV_free_rep1 not removed; leave-one-out not run.
    - The primary-result statement, QC caveat and limitations were added to `README.md`, `docs/results_by_process.md` and `docs/sample_structure_investigation.md`.
    - `docs/methods_references.md` §19 – §24 now link to D-022; `docs/software_versions.md` gained the rrcov rows.
81. Commit 2 (Phase 5) followed commit 1 (Phase 4).
    - Before each commit, staged files were audited: no FASTQ, index, SAM, RDS, credential or helper/staging files; secret-pattern scan run.
    - Both commits were pushed to `origin/main`, with Windows `gh.exe` as a one-off credential helper (nothing stored in the repo config). Local HEAD = `origin/main` was verified with `git ls-remote`; hashes are in the Phase 5 handoff report.

## Phase 5 — not executed (owner restrictions)

Changes to the primary DESeq2 analysis, removal of STV_free_rep1, design changes, pairing or covariate models, leave-one-out sensitivity, threshold changes or nominal-p gene lists, functional enrichment (GO, agriGO, KEGG, MapMan/Mercator, Plant Reactome), overwriting Phase 4 outputs.

## Phase 6 — functional analysis at the pathway level (2026-10-26; primary DESeq2 result unchanged)

All times UTC.

82. **Read-only investigation** (outside the repo, `/tmp/p6`):
    - SGN `ITAG4.0_goterms.txt`: 33,562 rows; 5,845 with GO; only 1,458 GO-bearing IDs are current ITAG4.0 gene IDs; 21,311 IDs carry other versions.
    - SGN ITAG4.1 GO and proteins: 12,410 GO transcripts match ITAG4.0 exactly, with identical proteins.
    - PLAZA 5.0 dicots tomato GO: ITAG4.0, exact IDs.
    - Plant Reactome tomato files.
    - The ITAG4.0 protein FASTA: 34,075 records, 32,522 terminal and 334 internal stops.
    - The 39 Cook's-filtered genes.
83. **References verified.**
    - Crossref:
      - 10.1093/nar/gkab1024, 10.1177/0962280212460441, 10.1093/nar/gks461, 10.1093/bioinformatics/btw763;
      - 10.1101/767764, 10.1111/j.2517-6161.1995.tb02031.x, 10.3390/plants15182842, 10.3389/fpls.2026.1733427;
      - 10.1186/s13100-025-00381-w, 10.1093/nar/gkz996, 10.1007/978-1-0716-1609-3_9, 10.1101/060012;
      - plus the DOIs listed in the earlier Phase 6 Crossref query: Subramanian 2005, Yu 2012, Wu 2021, Xu 2024, Schwacke 2019, Tian 2017, Kim & Volsky 2005, Ashburner 2000, GO Consortium 2023, Gupta 2024.
    - Full texts (Europe PMC): Tomaszewski et al. 2026 (PMC12950737), Li et al. 2026 (PMC13610134), Kelsey et al. 2025 (PMC12659357).
    - Technical documentation:
      - GSEA User Guide; fgsea vignette (Bioc release); fgsea 1.36.2 source (warning texts);
      - GO "Download the ontology" and "Relations" pages;
      - Plant Reactome "Computational Inference" and Release Summary V24;
      - agriGO v2 home/manual/download pages;
      - Mercator4 page (v8.0; maintenance notice; FAQ/version pages 404).
84. **`stv-r` environment.**
    - Dry run: 5 new packages, no updates or removals.
    - Install: `CONDA_CHANNEL_PRIORITY=strict conda install -n stv-r --override-channels -c conda-forge -c bioconda … "bioconductor-fgsea=1.36.2" "r-ontologyindex=2.12"` (plus the existing pins).
    - Added: fgsea 1.36.2, ontologyIndex 2.12, data.table 1.18.6.1, fastmatch 1.1-8, cowplot 1.2.0. Lock regenerated: 156 → 161 packages.
    - Global `channel_priority` is still `flexible`.
85. `bash scripts/bash/download_functional_annotation.sh`, first run 06:13:30Z.
    - Downloaded into git-ignored `reference/functional/` (86 MB): PLAZA `go.sly.csv.gz`; GO `go-basic.obo` release 2026-08-05 (`data-version releases/2026-07-26`); Plant Reactome `gene_ids_by_pathway_and_species.tab`; ITAG4.1 GO, proteins and README.
    - The first run exited 141: SIGPIPE from `zcat | head` under `pipefail`, after all downloads. Fixed.
    - Re-run 06:16:28Z: sizes verified; `results/functional/annotation_sources_manifest.tsv` written with sha256 and Last-Modified.
86. `Rscript scripts/R/05_functional_ranking_annotation.R` (Tasks 1 – 4; no enrichment test).
    - The first run stopped with the STOP check "GO terms with unknown namespace": `get_ontology` also returns 5 `[Typedef]` relation stanzas. The check was restricted to `GO:` IDs.
    - Second run exit 0 (06:16:36 – 06:16:59Z). Third run at 06:22 – 06:24Z to fix a clipped figure title; the ranking, catalog and mapping tables were sha256-identical.
    - Results: 18,238 ranked genes; GO coverage 82.88 %; tested sets BP 952 / MF 459 / CC 237; Plant Reactome 622 ranked genes (3.41 %), 5 testable pathways → collection D stopped (D-028).
87. **Pre-specification.**
    - `docs/functional_analysis_plan.md` §1 – §7 and `docs/decisions.md` D-023 – D-032 were written before any enrichment test.
    - Frozen copies and sha256 are in `results/functional/prespecification/`: plan `e30d60f6b418aaa642b36d01ca93e654271b8342bdcfd021139f4f3e91d66b12`, decisions `c4fcff5c87de917d8cfafff13b405ec8d5980750fffa9a6db1c15412ba1413b7`, frozen 2026-10-26T06:24:15Z.
88. `Rscript scripts/R/06_functional_gsea.R`.
    - Launches at 06:22:46Z and 06:23:27Z stopped at the pre-specification guard before any test: `fread` parsed the timestamp as a date-time. Fixed with `colClasses = "character"`.
    - **Primary run**: 06:24:15 – 06:25:07Z, exit 0 (local log `logs/phase6_06_gsea_run1.log`; `logs/` is git-ignored).
    - Results:
      - significant sets (padj < 0.05): BP 120 (75 positive / 45 negative NES), MF 51 (21 / 30), CC 44 (40 / 4);
      - 0 warnings, 0 NA; same-seed rerun identical;
      - seed-diagnostic changes: BP 15, MF 0, CC 0;
      - 29 enrichment curves by the fixed rule.
89. **Output-only re-run of `06_functional_gsea.R`** (06:26:47 – 06:27:39Z).
    - Changes: added `gsea_go_*_seed_diagnostic.tsv` and shortened curve titles.
    - `sha256sum -c` confirmed that all primary tables were unchanged: `gsea_go_*_all.tsv`, `*_significant_fdr05.tsv`, summary, leading-edge tables, figure selection.
    - The 15 BP calls that differ under seed 20261026 all have primary padj 0.031 – 0.070.
90. `Rscript scripts/R/07_mapman_mercator_prep.R` (06:32:50Z), exit 0.
    - Mercator4 input manifest and FASTA record (sha256 `84c15ec4…f4bb`, unmodified).
    - MapMan expression table (18,277 rows; apeglm, MLE, Wald) and `MERCATOR4_SUBMISSION.md`. The manual web submission is **pending**.
91. `Rscript scripts/R/08_functional_posthoc_descriptives.R` — **post-hoc, not pre-specified**, descriptive only.
    - 15 of 2,858 leading-edge genes are on SL4.0ch00.
    - They form ≥ 10 % of the leading edge in 27 of 215 significant sets.
92. **Docs.**
    - New: `docs/functional_analysis_plan.md` (§8 outcome appended after the run), `docs/agrigo_feasibility.md` (conclusion C).
    - Updated: `docs/decisions.md` (D-023 – D-032), `docs/methods_references.md` (§25 – §33), `docs/results_by_process.md` (rows 21 – 29), `docs/software_versions.md`, `envs/stv-r.yml`, `README.md`, this log.
93. **Git:** nothing committed or pushed (owner review required). The policy for the new tables > 1 MB is proposed as D-032.

## Phase 6 — not executed (by design / owner restrictions)

- Changes to the DESeq2 analysis; nominal-p or relaxed DEG lists; ORA (agriGO SEA, Mercator4/MapMan enrichment, `enrichGO`).
- Plant Reactome GSEA; PLAZA projected GO; KEGG.
- Sensitivity runs with other sizes, metrics or permutation settings.
- agriGO jobs.
- Mercator4 web submission (manual step pending).
- Leave-one-out.
- Commit/push.

## Phase 6B — functional robustness validation (2026-10-26; primary DESeq2 and Phase 6 GSEA results unchanged)

All times UTC. Scripts write logs to the git-ignored `logs/phase6b_*.log`. Ad-hoc helpers (file transfer, runners, look-ups)
live outside the repository and are not committed. Every reported number comes from a tracked script or table.

94. **Owner decisions recorded** (Phase 6B instruction):
    - D-024 approved (PLAZA 5.0 as the GO source; SGN remains the genome/transcriptome reference);
    - D-032 approved (track the five Phase 6 tables > 1 MB);
    - D-029: Mercator4 is not blocking; manual submission pending; no unofficial API or workaround;
    - D-030 approved (agriGO evaluated and rejected);
    - D-023 left to an audit, which is resolved by D-033.
95. **Task 1 — `p6_le.R` incident.**
    - Background. In Phase 6 the leading-edge summary in `functional_analysis_plan.md` §8.3 was computed by untracked ad-hoc
      helpers outside the repository. `p6_le.R` was rewritten on disk with `sed -i` while an editor buffer held an edited
      copy, so the editor's save failed ("newer on disk").
    - Findings:
      - The file on disk kept the original code, without the `as.numeric()` fix.
      - Its two runs (`p6_le.out`, `p6_le2.out`) stopped with a data.table integer/double type error; `p6_le3.out` is empty.
      - The reported numbers came from `p6_le_v2.R`, run read-only from a `/tmp` copy via `p6_le_v3.sh`.
      - No helper wrote into the repository.
      - The primary GSEA tables were written 06:24 – 06:27Z, before the incident, and are covered by the checksums of entry 96.
    - Resolution (D-041): `scripts/R/09_functional_leading_edge_summary.R` (07:44Z, exit 0) reproduces every number from
      tracked tables:
      - 2,858 unique leading-edge genes (1,639 positive / 1,219 negative statistic);
      - median \|stat\| 0.704; maximum single-gene share 0.228;
      - median leading-edge fraction 0.417 / 0.387 / 0.581;
      - unique genes BP 1,695, MF 1,460, CC 863; sets with padj < 0.01: 59 / 23 / 28.
    - Lesson: read typed columns with `colClasses` and never edit one helper from two tools.
96. **Task 2 — checksums before validation.**
    - `results/functional_robustness/phase6_outputs_sha256.txt`: SHA-256 of the 76 Phase 6 primary outputs (ranking,
      annotation, gene-set catalog, GSEA tables, leading edges, figures).
    - The Phase 6 versions of the shared docs, README and env files were saved outside the repository for the two-commit split.
97. **D-023 audit** → [`cooks_gsea_decision.md`](cooks_gsea_decision.md).
    - Sources: DESeq2 vignette, Love et al. 2014, fgsea 1.36.2 source, TE-Seq code, Cook 1977.
    - No method paper was found, so the evidence is ambiguous → primary unchanged, inclusion sensitivity (D-033).
98. **`stv-r` environment.**
    - Dry run: 3 new packages, no updates or removals.
    - Install: `CONDA_CHANNEL_PRIORITY=strict conda install -n stv-r --override-channels -c conda-forge -c bioconda` with the
      existing pins plus `"bioconductor-limma=3.66.0" "bioconductor-edger=4.8.2"` → limma 3.66.0, edgeR 4.8.2, statmod 1.5.2.
    - Lock regenerated (161 → 164 packages); `envs/stv-r.yml` updated.
99. **Pre-specification frozen at 2026-10-26T07:36:28Z**, before any CAMERA, sensitivity or leave-one-out output existed.
    - Frozen: `docs/functional_robustness_plan.md` §1 – §10, `docs/cooks_gsea_decision.md` and `docs/decisions.md` (D-033 – D-041).
    - Snapshots and SHA-256 are in `results/functional_robustness/prespecification/`.
    - Every Phase 6B R script re-checks them (`fr_check_frozen()`).
100. **`10_functional_cooks_sensitivity.R`** (07:44Z, exit 0).
    - Helper validation against Phase 6: max \|Δp\| 5.6 × 10⁻¹⁶, max \|ΔNES\| 5.3 × 10⁻¹⁵, identical leading edges.
    - Rebuilt gene sets are identical to `go_gene_sets.rds`.
    - Inclusion run (18,277 genes): 210 sets retained, 5 supported, 0 lost, 10 sensitivity-only; 0 material warnings.
101. **`11_functional_camera.R`** (07:45Z, exit 0).
    - `quant.sf` SHA-256 matched.
    - tximport `lengthScaledTPM`: 33,976 × 6.
    - TMM factors 0.920 – 1.061; voom weights 0.58 – 12.9; limma prior df 3.56; 0 genes at adj. P < 0.05.
    - Spearman ρ(limma t, Wald) 0.9626.
    - CAMERA `inter.gene.cor = NA`: **0 significant sets** (BP/MF/CC), minimum FDR 0.997 – 0.999.
    - Secondary 0.01: 40 / 7 / 26.
102. **`12_functional_method_comparison.R`** (07:45Z, exit 0).
    - Tier 1 = 0, Tier 2 = 0, Tier 3 = 210, Tier 3C = 0, Tier 4 = 5.
    - Direction agreement 100 % for the fgsea-significant sets; signed-significance ρ 0.923 / 0.914 / 0.893.
    - No STOP.
103. **`13_functional_replicate_consistency.R`** (07:46Z, exit 0).
    - 144 sets broadly consistent, 69 partially consistent, 2 ambiguous, 0 strongly influenced by STV_free_rep1.
    - F2 figures contain diagnostic examples only, because Tier 1 is empty.
104. **`bash scripts/bash/download_organelle_reference.sh`** (07:46Z).
    - Downloaded NC_007898.3 and NC_035963.1 (genome + CDS proteins) into the git-ignored `reference/functional/organelle/`.
    - Manifest with SHA-256 written.
105. **`bash scripts/bash/ch00_organelle_blast.sh`** (07:46Z, BLAST+ 2.17.0+ from `rachis-qiime2-2026.4`).
    - Queries: 513 ch00 genes; blastn hits for 1,061 ITAG4.0 cDNAs; 222 blastp rows.
    - Warning "Examining 5 or more matches is recommended" is harmless (two subjects).
106. **`15_functional_ch00_investigation.R`** (07:47Z, exit 0).
    - 186 of 513 ch00 genes are near-identical to an organelle genome; 22 ranked genes were excluded.
    - ch00-sensitive sets: 15 / 5 / 4 (12.5 / 9.8 / 9.1 %). No STOP.
107. **`14_functional_loo_sensitivity.R`** (07:47 – 07:49Z, exit 0).
    - Six five-sample fits ran with 0 material fgsea warnings.
    - Maximum sample-dependent share 47.7 % (CC, STV_infected_rep1). No STOP.
108. **`16_functional_robustness_figures_themes.R`** (07:55Z, exit 0).
    - Every theme is EXPLORATORY.
    - Figures F1, F3, F4, F5 written.
109. **Integrity checks.**
    - `sha256sum -c results/functional_robustness/phase6_outputs_sha256.txt`: 76/76 OK.
    - `git diff` of tracked Phase 1 – 5 results, figures, config, data and scripts: empty.
    - Freeze snapshots unchanged.
    - New files > 1 MB: only `camera/voom_limma_gene_statistics.tsv` (1.9 MB), now git-ignored (D-042).
110. **Corrections and post-hoc additions** (plan §11.3, D-043).
    - Fact: "45 ranked ch00 genes" corrected to 47.
    - BLAST `-max_hsps 1` limitation: `18_functional_ch00_posthoc_identity_diagnostic.R` (08:11Z, exit 0) excluded 42 genes;
      ch00-sensitive sets 19 / 3 / 4; STOP assessment unchanged.
    - `17_prior_study_locus_lookup.R` (08:09Z, exit 0): prior-study loci in the primary results.
111. **Literature and references** (2026-10-26).
    - Crossref REST checks: 10.1186/gb-2014-15-2-r29, 10.1093/nar/gkv007, 10.1186/gb-2010-11-3-r25,
      10.12688/f1000research.7563.2, 10.12688/f1000research.9005.3, 10.2307/1268249 (→ 10.1080/00401706.1977.10489493),
      10.1093/bioinformatics/btm051, 10.1186/1471-2105-10-421, 10.1093/hr/uhac251, 10.1007/s00705-019-04436-1,
      10.1016/j.virusres.2008.11.018, 10.1098/rstb.2010.0057, 10.1046/j.1365-313x.2003.01625.x, 10.1093/nar/gkaf018,
      10.1007/s00239-005-0254-5.
    - Europe PMC records: 10.3389/fmicb.2019.02764 and 10.3389/fmicb.2016.01905.
    - Full texts read: Liu et al. 2023 (PMC9832971), Khankhum et al. 2016 (PMC5126043), Badeli et al. 2026 (Research Square HTML).
    - Abstracts read: Fukuhara et al. 2020 (PubMed), Sabanadzovic 2009, Roossinck 2010, Takahashi 2019, Whitham 2003,
      Goeman & Bühlmann 2007.
    - Tomato–virus transcriptome search (two Europe PMC queries): no suitable open-access study; none cited.
112. **Docs.**
    - New: `functional_robustness_plan.md` (+ §11 outcome), `cooks_gsea_decision.md`, `ch00_functional_investigation.md`,
      `biological_comparison_prior_studies.md`.
    - Updated: `decisions.md` (D-023/24/29/30/32 approvals, D-033 – D-043), `methods_references.md` (§34 – §41),
      `results_by_process.md` (rows 30 – 38 and the evidence table), `functional_analysis_plan.md` (pointer to 6B),
      `agrigo_feasibility.md` (owner decision), `software_versions.md`, `envs/stv-r.yml`, `.gitignore`, `README.md`, this log.
113. **Pre-commit audit and Git.**
    - The audit was read-only:
      - 76/76 Phase 6 checksums OK and freeze snapshots intact;
      - no diff in tracked Phase 1 – 5 results, figures, config, data, scripts or workflow;
      - the only candidate files > 1 MB are the five tables approved under D-032;
      - no FASTQ/FASTA/RDS/log/helper/staging path; no secret pattern; no CRLF;
      - all relative Markdown links resolve.
    - Absolute paths appear only in R `sessionInfo` BLAS lines and in earlier log entries, as in the committed Phase 4/5 files.
    - Commit 1 "Phase 6: threshold-free functional enrichment analysis": the Phase 6 scripts, results and figures, with the
      shared docs and env files restored to their saved end-of-Phase-6 versions (SHA-256 checked).
    - Commit 2 "Phase 6B: functional robustness validation": everything else, with the final versions of the shared files.
    - Both commits are pushed to `origin/main` with the Windows `gh.exe` one-off credential helper; nothing is stored in the
      repo config.
    - Commit hashes and the `HEAD == origin/main` check (`git ls-remote`) are in the Phase 6B handoff report.

## Phase 6B — not executed (by design / owner restrictions)

- Changes to the primary DESeq2 analysis or to the Phase 6 GSEA results; removal of STV_free_rep1; FDR relaxation; DEG lists;
  nominal-p ORA.
- CAMERA re-runs per dropped sample; any tuning of CAMERA after results.
- Re-quantification with organelle genomes in the Salmon index (a new analysis; owner decision).
- Mercator4 web submission (manual step still pending; no workaround).

## Phase 6B closure — decision audit and organelle-aware quantification assessment (2026-10-26; no result changed)

114. **Task 0 audit (read-only).**
    - Branch `main`; `HEAD` = `origin/main` = `aee2a9ac27d96a64455abfb17554ea19736d8d9e`, parent `928d763`; tree clean;
      no stash.
    - `sha256sum -c results/functional_robustness/phase6_outputs_sha256.txt`: 76/76 OK.
    - Phase 6 / 6B freeze snapshots OK; final `quant.sf` 6/6 OK against `results/quantification/quant_sf_sha256.tsv`.
    - Ignored intermediates: `txi_final.rds` `82d1b43e312a4c475cbe92c237760108133e05ee2736d5bd5fefe17a2a620457`,
      `deseq2_objects.rds` `3975624c0fa6e5ab0f9a3f81694a6bef076d10b526a4b8a1a942b9676482b95f`, `go_gene_sets.rds`
      `81c6a4b93dd94d456af4637e45318ce482e1214bb2dac5e091e22f1c60869282`; voom table
      `a10edfcb0e7080f3c7e20f8627a5491f746a4f344807be5030f6e56692d09893` (1,952,291 B).
115. **Organelle assessment criteria frozen** (§1 of `docs/organelle_requantification_assessment.md`; written
     2026-10-26T15:17:12Z, frozen 15:18:09Z), before any of the research below.
    - Command (ad-hoc helper):
      `cp docs/organelle_requantification_assessment.md results/organelle_assessment/prespecification/organelle_assessment_question_frozen.md`
    - Its `sha256sum` was written with both UTC times to `prespecification_checksum.tsv`. The Phase 6 and Phase 6B freezes
      were made the same way.
116. **Salmon reference reconstruction (read-only).** Sources: `workflow/modules/salmon_index.nf`,
     `workflow/nextflow.config`, index `info.json`, index log, `duplicate_clusters.tsv` and `meta_info.json`.
    - 13 decoys (SL4.0ch00 – ch12); 33,976 targets; 99 duplicates removed.
117. **Provenance research.**
    - SGN release notes and FASTA headers / gap structure; ITAG4.0 GFF sources and README.
    - NCBI Datasets v2 sequence report for GCA_000188115.5.
    - Ensembl REST `info/genomes/taxonomy/4081`, `info/assembly/solanum_lycopersicum_gca000188115v5cm` and
      `info/assembly/arabidopsis_thaliana`.
    - Crossref for 10.1101/767764 and 10.1186/s13059-020-02151-8.
    - Unsuccessful, recorded and not worked around:
      - Hosmani et al. bioRxiv full text (not retrievable);
      - Ensembl REST `info/assembly/solanum_lycopersicum` (HTTP 400; the species is registered under the GCA-based
        name);
      - Crossref request for 10.1038/nmeth.4197 (no content returned; citation verified via the Srivastava et al.
        reference list);
      - rdrr.io page for `?camera` (redirect; the installed limma 3.66.0 help was read instead).
118. **Method documentation.** Salmon "Preparing transcriptome indices"; Srivastava et al. 2020 full text (PMC7487471);
     nf-core/rnaseq 3.27.0 usage. No forum sources.
119. **`scripts/bash/organelle_multihsp_blast.sh`** (BLAST+ 2.17.0+, env `rachis-qiime2-2026.4`, exit 0).
    - 2,512 HSP rows for 1,061 cDNA queries; 835 rows for the 177 ch00 segments.
    - The first draft split ch00 at upper-case N only, which missed the gaps: they are lower-case `n`. Fixed to
      `[Nn]{100,}` before the recorded run.
120. **`scripts/R/19_organelle_scope_diagnostic.R`** (exit 0).
    - First run failed: data.table scoping, where a column `genes` shadowed a table of the same name. Fixed by precomputing
      `total_baseMean`.
    - Later, the gene table was restricted to ch00 genes + anchored genes with any ≥ 98 % HSP (5.3 MB → 709 rows,
      126,527 B). After the re-run, all other outputs were byte-identical, and the restricted table equals the matching rows
      of the full table (SHA-256 compared).
121. **`scripts/bash/ch00_ncbi_provenance.sh`** (exit 0).
    - 164 / 177 ch00 segments are exactly present in GenBank SL4.0.
    - 12 mito-identical segments (832,437 bp) are absent. ch00 segment 15 = AEKE04000081.1 (complete plastome copy).
    - Before the script, an exploratory ad-hoc check had compared two scaffolds and the scaffold lengths; the script
      supersedes it.
122. **Mercator4** official page (<https://www.plabipd.de/mercator_main.html>), checked once on 2026-10-26: Mercator4
     v8.0, "currently undergoing scheduled maintenance". Nothing was submitted.
123. **Docs.**
    - New: `organelle_requantification_assessment.md` §2 – §11, `phase6b_decision_audit.md`,
      `reproducibility_audit.md`, `pipeline_execution_map.md`.
    - Updated: `decisions.md` (D-042 update, D-044), `methods_references.md` (§40 correction note, §42, §43),
      `results_by_process.md` (rows 39 – 42), correction notes in `ch00_functional_investigation.md` and
      `functional_robustness_plan.md` §7, `functional_analysis_plan.md` (pointer), `software_versions.md`, `README.md`,
      this log.
124. **Not executed:** Salmon index / quant, tximport, DESeq2, fgsea, CAMERA, Mercator4 submission, any change to the primary
     analysis. Phase 7 was not started.
125. **Owner decisions recorded** (2026-10-26).
    - D-044 Option A approved.
    - D-042 REGENERATE ON DEMAND approved.
    - No gene-exclusion fgsea.
    - The "organelle-identical feature class" label is allowed, but never as a statement of origin.
    - Updated: `decisions.md`, `phase6b_decision_audit.md`, `organelle_requantification_assessment.md` §8,
      `results_by_process.md`, `functional_analysis_plan.md` §8.6, `reproducibility_audit.md` R13, `README.md`.
    - Pre-commit checks then ran and one commit "Phase 6B closure: decision audit and organelle assessment" was pushed.
      Hashes and remote verification are in the Phase 7 entries below.

## Phase 7 — biological integration (2026-10-26; frozen results only)

126. **Phase 6B closure commit.**
    - Pre-commit audit: 38 staged files, largest 226,202 B; no forbidden path, secret pattern or CRLF; only docs/README
      modified among tracked files.
    - Commit `13e87be8c8c03702c90b85438d8794ab39d4d1ed`, pushed with the one-off `gh.exe` credential helper.
    - `git ls-remote origin refs/heads/main` = local HEAD; working tree clean; no credential stored in the repo config.
127. **Literature search** (`scripts/bash/literature_search_europepmc.sh`; Europe PMC REST; metadata only).
    - First run: HTTP 503. The cause was the `sort=RELEVANCE` parameter (curl test: 200 without it, 503 with it); the
      parameter was removed and retry/backoff added.
    - Hits: Q1 32, Q2 67, Q3 84, Q4 19.
    - Abstracts read for 23 records through the API (ad-hoc reader; nothing saved). One targeted query timed out and was
      repeated successfully.
    - Full text read: Badeli et al. 2026 (Research Square HTML v1). ICTV *Amalgaviridae* report page read.
    - A `curl` query passed through PowerShell lost its quotes (wrong hit count); it was repeated from a script file.
128. **`scripts/R/20_biological_integration.R`** (exit 0; log `logs/phase7_20_biological_integration.log`).
    - First run failed (`setorder` on an expression); fixed.
    - The ethylene description regex was first applied to full descriptions and matched UniProt identifiers inside the
      AHRD tags (e.g. "EIL"). It now matches short descriptions only.
    - U2 was first coded with LOO alone; revised to also require weak replicate consistency (D-045; disclosed in
      `biological_integration.md` §2).
    - Final: 7 EXPLORATORY, 3 UNSUPPORTED.
129. **`scripts/R/21_go_redundancy_modules.R`** (exit 0).
    - First run failed because `data.table(key = ...)` treats `key` as an argument; the column was renamed `set_key`.
    - Result: 73 modules at cut 0.7.
130. **`scripts/R/22_phase7_figures.R`** (exit 0). F1 – F4.
    - Layout fixes: wrapped caveats (F3); shorter labels and aligned panels (F4).
    - One edit was lost to a race between pulling the file and the editor saving it; it was re-applied and verified on
      disk before syncing.
131. **Consistency checks.**
    - Theme-level LOO dependence recomputed from `loo_dependence_by_set.tsv` matches the frozen
      `loo_dependence_summary.tsv` (e.g. CC STV_infected_rep1 21/44 = 47.7 %). No inconsistency between frozen tables and
      Phase 6B summaries (STOP 8 not triggered).
    - Ribosome theme: 26/34 sets depend on STV_infected_rep1. This is not "entirely one sample", and the collection-level
      dependence was documented in Phase 6B (STOP 7 not triggered). ER: 4/4 on STV_free_rep1, already documented in the
      Phase 6B evidence table.
132. **Mercator4** (<https://www.plabipd.de/mercator_main.html>) checked once: v8.0, "currently undergoing scheduled
     maintenance". Not submitted; Phase 7 does not depend on it.
133. **Docs.**
    - New: `biological_integration.md`, `literature_biological_context.md`, `phase7_biological_interpretation.md`.
    - Updated: `decisions.md` (D-045 – D-048), `methods_references.md` (§44 – §46), `results_by_process.md` (rows
      43 – 48), `software_versions.md`, `README.md`, this log.
134. **Not executed:** Salmon, tximport, DESeq2, fgsea, CAMERA, gene-exclusion fgsea (owner: not justified), Mercator4 /
     MapMan, clean-clone reproducibility audit.
135. **Phase 7 commit.** `d29a849f8481ac495cc48666087b394614da587c` "Phase 7: integrated biological interpretation"
     (43 files; largest 529,969 B). Pushed; `git ls-remote` = local HEAD = `origin/main`; working tree clean.

## D-045 revision (2026-10-26; owner decision REVISE; no new analysis)

136. **Classification framework revised; scripts 20 and 22 re-run on unchanged inputs.**
    - Owner decision: neither U2 v1 nor U2 v2 converts a theme to UNSUPPORTED. Evidence status comes from inferential
      evidence; robustness diagnostics give a separate descriptive descriptor (D-045 revision, D-049).
    - History reconstructed from the session record: v1 coded about 16:30 UTC and run once (UNSUPPORTED: ribosome, ATP
      synthase, defence, ER, cell wall); v2 coded about 16:33 UTC after the v1 matrix was viewed (UNSUPPORTED: ribosome, ER,
      ATP synthase); v2 committed in `d29a849`. Both stay in the matrix as flags and superseded classes.
    - `scripts/R/20_biological_integration.R`: new columns `evidence_status`, `evidence_status_basis`,
      `robustness_descriptor`, `robustness_basis`, `loo_dependent_sets_by_sample`, `rule_U1_v1_prespecified_only`,
      `superseded_phase7_class_v1`; `phase7_class` / `phase7_class_reason` renamed `superseded_phase7_class_v2` /
      `_v2_reason` (values unchanged). `curated_theme_literature.tsv`: `primary_caveat` added; U2 wording removed from
      the interpretation / caveat text.
    - Re-run of script 20 (exit 0): only `theme_evidence_matrix.tsv` changed. 52 of its 54 earlier columns are identical
      (the renamed class columns included); the other two are the interpretation and caveat text. All other script 20
      outputs, including `sessionInfo_20.txt`, are byte-identical.
    - Re-run of script 22 (exit 0): F1 and F3 changed as intended; F2 byte-identical. F4 bytes changed only because its
      point jitter and `ggrepel` label placement are unseeded; the committed F4 was restored with `git checkout` (content
      unchanged; reproducibility note in `biological_integration.md` §7).
    - Result: 9 EXPLORATORY, 1 UNSUPPORTED (ATP synthase, explicit review). Robustness: ribosome/translation and ER strongly
      sample-sensitive; photosynthesis/plastid mixed; ATP synthase technically sensitive; defence, heme, glutathione and
      cell wall moderately sample-sensitive; hormone/ethylene and water/abiotic stress relatively stable.
    - An inline `wsl.exe bash -c` loop lost its quoting in PowerShell (`$f` expanded by PowerShell) and created an empty
      junk file in the repo root; it was removed before any staging. Syncs were then run from a script file.
137. **Methodological sources** (`methods_references.md` §47). Europe PMC abstracts and full texts (Schurch et al. 2016;
     Wu & Smyth 2012; Maciejewski 2014; Geistlinger et al. 2021), abstracts (Goeman & Bühlmann 2007; Tamayo et al. 2016;
     Gatti et al. 2010), GSEA User Guide HTML, installed limma 3.66.0 `?camera`; Crossref metadata for all seven DOIs. The
     current-release `camera` help could not be retrieved (code.bioconductor.org TLS certificate expired; not bypassed).
138. **Not executed:** Salmon, tximport, DESeq2, apeglm, fgsea, CAMERA, LOO, Cook's analysis, ch00 exclusion, any
     enrichment, Mercator4, MapMan.
