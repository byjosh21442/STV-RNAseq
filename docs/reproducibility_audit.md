# Reproducibility audit (Phase 6B closure)

Written 2026-10-26. Read-only audit of whether a third party can regenerate the project's results from this repository.
Severity levels:
- **BLOCKING**: a result cannot be regenerated.
- **IMPORTANT**: it can be regenerated, but only with information that was missing from the repo before this audit.
- **MINOR**: an inconvenience, or an expected limitation of public resources.
- **NONE**: checked, no issue.

Every fix listed under "Fixed" is a documentation change. No script logic, no result and no environment was changed.

## Findings

| # | Area | Finding | Severity | Action |
|---|---|---|---|---|
| R1 | Environments | `envs/stv-rnaseq.yml`, `envs/stv-r.yml` + linux-64 lock files (stv-r lock: 164 packages incl. limma / edgeR / fgsea) | NONE | — |
| R2 | BLAST+ (Phase 6B ch00 and D-044 scripts) | `ch00_organelle_blast.sh` and `organelle_multihsp_blast.sh` run BLAST+ through `conda run -n rachis-qiime2-2026.4`. That is a pre-existing environment that is **not defined in the repo**. Installed package: `blast 2.17.0` (bioconda, build `h66d330f_0`) | IMPORTANT | **Fixed (docs):** README and `software_versions.md` now say how to provide an equivalent environment. The script code is unchanged. Making the env name configurable would be a code change: owner decision |
| R3 | Conda location | Bash scripts that need conda source `$HOME/miniconda3/etc/profile.d/conda.sh` | MINOR | Documented in README. Edit the path on other installations |
| R4 | Run order | README listed Phase 1 – 6B. The D-044 scripts (`organelle_multihsp_blast.sh` → `19_organelle_scope_diagnostic.R` → `ch00_ncbi_provenance.sh`) were missing | IMPORTANT | **Fixed (docs):** README run order and structure updated. Full map in [pipeline_execution_map.md](pipeline_execution_map.md) |
| R5 | Ignored R intermediates | Scripts 10, 11, 13, 14, 15 and 18 read `txi_final.rds` (script 01), `deseq2_objects.rds` (02) and `go_gene_sets.rds` (05). All three are git-ignored. SHA-256 at closure: `82d1b43e312a4c475cbe92c237760108133e05ee2736d5bd5fefe17a2a620457`, `3975624c0fa6e5ab0f9a3f81694a6bef076d10b526a4b8a1a942b9676482b95f`, `81c6a4b93dd94d456af4637e45318ce482e1214bb2dac5e091e22f1c60869282` | MINOR | Regenerable by scripts 01 / 02 / 05. An `.rds` file may not be byte-identical across R builds, so verify the tracked TSV outputs instead. Checksums recorded here and in `run_log.md` |
| R6 | Ignored primary inputs | FASTQ (ENA MD5 tracked), SGN reference (SHA-256 manifest), Salmon `quant.sf` (SHA-256 in `results/quantification/quant_sf_sha256.tsv`), organelle FASTA (SHA-256 in `results/functional_robustness/ch00/organelle_reference_manifest.tsv`) | NONE | — |
| R7 | Salmon re-run identity | Whether a fresh `salmon quant` run reproduces the recorded `quant.sf` SHA-256 bit for bit was **not tested** in this audit (no quantification allowed) | MINOR | `scripts/bash/compare_salmon_runs.sh` compares two runs transcript by transcript |
| R8 | Pre-specification freezes | Phase 6, Phase 6B and D-044 freezes are tracked snapshots + SHA-256. Phase 6B scripts and script 19 verify them. The command that created them was an ad-hoc helper: `cp <doc> <snapshot>; sha256sum <snapshot>` written to `prespecification_checksum.tsv` with UTC times | MINOR | **Fixed (docs):** command recorded in `run_log.md` (entry 115) |
| R9 | Remote resources | PLAZA 5.0, GO OBO, Plant Reactome, NCBI RefSeq / GenBank, SGN may change or move. SHA-256 and retrieval dates are recorded; the NCBI sequence report used by D-044 is stored | MINOR | Expected for public resources; manifests allow detection of changes |
| R10 | Manual step | Mercator4 v8.0 web submission (Phase 6, D-029). Still pending: on 2026-10-26 the official page again showed "currently undergoing scheduled maintenance". No result depends on it | MINOR | Documented in `results/functional/mapman/MERCATOR4_SUBMISSION.md` |
| R11 | Untracked helpers (`.stv_tools`, staging, `p6b_run.sh`) | Used to edit files from Windows into WSL and to run scripts with logging into the ignored `logs/`. No tracked script or result references them. No scientific number depends on them (D-041) | NONE | — |
| R12 | Absolute paths | Only in R `sessionInfo` BLAS lines and historical log/audit entries | NONE | — |
| R13 | Voom gene table | Git-ignored (D-042); regenerable by script 11; SHA-256 recorded in [phase6b_decision_audit.md](phase6b_decision_audit.md) | NONE | Resolved 2026-10-26: owner approved REGENERATE ON DEMAND (D-042) |
| R14 | Network-only provenance script | `ch00_ncbi_provenance.sh` needs NCBI Datasets / E-utilities, with system `python3` 3.14.4 and `curl` 8.18.0 | MINOR | Documented in `software_versions.md` |

**No BLOCKING finding.**
