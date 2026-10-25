# Pipeline design

Status: **design proposal; quantification strategy requires human approval.** Nothing was run on real data.

## Documentation consulted (retrieved 2026-10-25)

| Tool | Version seen in docs | URL |
|---|---|---|
| Nextflow | install page (Java 17–26) | <https://docs.seqera.io/nextflow/install> |
| nf-core/rnaseq | 3.27.0 | <https://nf-co.re/rnaseq/3.27.0/docs/usage> |
| FastQC | 0.12.x (changelog) | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/> |
| MultiQC | current docs | <https://docs.seqera.io/multiqc> |
| Salmon | 1.12.0 docs | <https://salmon.readthedocs.io/en/latest/salmon.html>, <https://salmon.readthedocs.io/en/latest/library_type.html> |
| tximport | 1.40.0 (Bioconductor 3.23) | <https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html> |
| DESeq2 | 1.52.0 (Bioconductor 3.23) | <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html> |
| Cufflinks/Cuffnorm | — | <http://cole-trapnell-lab.github.io/cufflinks/cuffnorm/> |

## Dataset constraints that drive the design

- 6 samples, 1 run each, **single-end 75 nt**, 14–20 M reads/sample, 5.7 GB of `fastq.gz` in total (`dataset_audit.md`).
- Strandedness **not documented** → must be measured, not assumed.
- Single-end: Salmon cannot estimate the fragment-length distribution; effective lengths use `--fldMean/--fldSD` (Salmon docs). True fragment length is UNKNOWN → keep Salmon defaults and document it as a limitation.
- Salmon docs: k = 31 "seems to work well for reads of 75bp or longer" → default k is acceptable for 75 nt; consider a smaller k only if mapping rates are unexpectedly low.

## What we take from nf-core/rnaseq (architecture reference only, not run)

- Samplesheet with `strandedness` that may be `auto`; nf-core infers it by sub-sampling 1 M reads and running Salmon, calling *forward/reverse* at ≥ 0.8 of fragments and *unstranded* when forward/reverse differ by < 0.1, otherwise *undetermined* (possible gDNA contamination).
- Salmon with a **decoy-aware index built from the genome** is the default recommendation when a genome FASTA is supplied.
- Salmon `--seqBias`/`--gcBias` are not on by default in nf-core; DESeq2 docs recommend `--gcBias` unless the data are known not to need it. Whether `--gcBias` is appropriate for **single-end** data: TO BE VERIFIED against Salmon behaviour before use.
- Since nf-core/rnaseq v3.0, **featureCounts is no longer used for quantification** (only biotype QC), citing Zhao et al. 2015 and Soneson et al. 2015; HISAT2 runs produce no quantification there.
- STAR needs "typically around 38GB" RAM for human GRCh37 (nf-core docs); HISAT2 is the documented low-memory alternative.
- Trimming (Trim Galore!/fastp) is on by default in nf-core; **we will not adopt this default** (project rule: trim only if QC justifies it).
- Not adopted (overkill for 6 samples): rRNA removal, contamination screening, UMI handling, RSeQC suite beyond strandedness, Qualimap, dupRadar, BAM post-processing.

## Option comparison

Resource numbers marked *est.* are order-of-magnitude expectations, not measurements on this machine.

| Criterion | Option A: FastQC/MultiQC → Salmon → tximport → DESeq2 | Option B: FastQC/MultiQC → STAR or HISAT2 → featureCounts → DESeq2 |
|---|---|---|
| Computational cost | Low: quantification ≈ minutes/sample on 20 threads (*est.*) | Higher: genome alignment + sorting + counting ≈ several × Option A (*est.*) |
| RAM | Transcript-only index: small. Full-genome-decoy index (≈0.8 Gb genome): "considerably more memory" (Salmon docs) — expected to fit in 15 GiB but **must be measured** | STAR index of ≈0.8 Gb genome ≈ 8–10 GB (*est.*, scaled from STAR/nf-core guidance) — tight in 15 GiB WSL; HISAT2 lower |
| Disk (beyond 5.7 GB FASTQ + ≈1 GB reference) | Index ≈ few GB (*est.*); outputs tens of MB | + STAR index ≈ 6–10 GB (*est.*) + sorted BAMs ≈ 1 GB/sample (*est.*) |
| Runtime, whole dataset | ≈ 1–2 h incl. index build (*est.*) | ≈ several hours (*est.*) |
| Reproducibility | High: single tool + pinned index; `quant.sf` + `lib_format_counts.json` record library type | High, but more tools/parameters (aligner, BAM handling, `-s` strand flag) |
| Interpretability / QC | Mapping rate, auto-detected library type, no genome-browser view | BAMs allow IGV inspection, splice/coverage QC, RSeQC strandedness |
| Multi-mapping reads (paralog families) | Probabilistically assigned (EM) — documented advantage in tximport/DESeq2 docs | featureCounts discards or counts multi-mappers crudely |
| Gene-level DESeq2 | `DESeqDataSetFromTximport()` with average-transcript-length offsets (officially recommended in DESeq2 vignette) | `DESeqDataSetFromMatrix()` with integer counts |
| Suitability, 6 tomato samples / 2 weeks | Very good: fast, low disk, few failure points | Feasible but more setup, RAM pressure and disk use |
| Scientific caveats | Depends on annotation completeness (reads from unannotated loci can be mis-assigned → mitigated by genome decoys); SE fragment length assumed | Gene-level only; ignores isoform length changes; strand flag must be correct |

### Recommendation (requires approval)

**Option A — Salmon (selective alignment, decoy-aware index) → tximport → DESeq2**, matching the DESeq2 authors' recommended pipeline and nf-core's default quantifier, at a fraction of Option B's cost. Option B is kept as a fallback only if genome-level inspection (e.g. IGV) becomes a project requirement.

## Proposed minimal workflow (to be built only after approvals)

```text
data/raw/*.fastq.gz  (ENA download script + MD5 check, outside Nextflow)
   │
   ├─ FASTQC (per sample) ──────────────┐
   │                                    ├─ MULTIQC → results/qc/multiqc_report.html
   ├─ SALMON_QUANT (per sample, -l A) ──┘
   │      ▲
   │  SALMON_INDEX (once; transcripts + genome decoys)
   ▼
results/quantification/salmon/<sample>/quant.sf
   │
scripts/R: tximport (tx2gene from GFF) → DESeq2 (~ condition, ref = STV_free) → functional analysis
```

- Nextflow: one `main.nf` (≈5 processes), one `nextflow.config` with `params`, resource labels and two profiles (`conda`, `docker`). R steps stay as plain scripts in `scripts/R/` (simpler to debug and report).
- Strandedness: run Salmon with `-l A` and read the detected type from each `lib_format_counts.json`; if all six agree (expected `SR` given the kit, but unverified), fix `-l` to that value for the final run and record it in `decisions.md`.
- Preprocessing policy: no trimming by default. Salmon's selective alignment uses extension alignment, so reads "need not map end-to-end" (Salmon docs), which tolerates short adapter tails. Trim only if MultiQC shows material adapter content or quality failures; the threshold must be set before looking at DE results.
- Optional QC (proposal): add the STV genome as an extra target to confirm infection status per sample from viral read counts. The STV reference accession has **not** been identified yet — do not use until verified.
