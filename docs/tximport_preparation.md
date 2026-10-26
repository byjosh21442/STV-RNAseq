# tx2gene and tximport preparation (Phase 3)

Status: the import is **prepared and validated**. No differential-expression analysis, filtering, thresholds, contrasts, interpretation or enrichment has been run.

## 1. tx2gene

Script: [`scripts/bash/build_tx2gene.sh`](../scripts/bash/build_tx2gene.sh) → [`reference/ITAG4.0_tx2gene.tsv`](../reference/ITAG4.0_tx2gene.tsv) (tracked; 34,075 rows + header `TXNAME	GENEID`), validation [`reference/ITAG4.0_tx2gene_validation.tsv`](../reference/ITAG4.0_tx2gene_validation.tsv).

- Source: `ITAG4.0_gene_models.gff`, mRNA → `Parent` gene. Columns use the GFF `Name` attributes verbatim (`Solyc00g500001.1.1` → `Solyc00g500001.1`).
- **No transcript or gene version is stripped:** Salmon's `quant.sf` names, the cDNA headers and the GFF Names are identical strings, so there was no documented reason to strip anything (tximport `ignoreTxVersion = FALSE`).

| Check | Result |
|---|---|
| tx2gene rows / unique transcripts / unique genes | 34,075 / 34,075 / 34,075 |
| Transcripts mapped to > 1 gene | 0 |
| cDNA IDs missing from tx2gene / tx2gene IDs missing from cDNA | 0 / 0 |
| `quant.sf` names (33,976) missing from tx2gene | **0** (exact string match) |
| tx2gene IDs absent from `quant.sf` | 99, all of them the duplicates collapsed by `salmon index` (`duplicate_clusters.tsv`); unexplained: **0** |
| Same transcript set and order in all six `quant.sf` | yes |

## 2. R environment `stv-r`

Isolated from the command-line tool env and from every other conda env. Created with `CONDA_CHANNEL_PRIORITY=strict conda env create -f envs/stv-r.yml` (strict priority set for that command only; global `channel_priority` stays `flexible`; `rachis-qiime2-2026.4` untouched).

| Package | Version | Channel |
|---|---|---|
| R | 4.5.3 | conda-forge |
| tximport | 1.38.2 | bioconda (Bioconductor 3.22) |
| DESeq2 | 1.50.2 | bioconda (Bioconductor 3.22), installed, **not run** |
| jsonlite | 2.0.0 | conda-forge |

Exact lock (139 packages): [`envs/stv-r.linux-64.lock.txt`](../envs/stv-r.linux-64.lock.txt); 1.7 GB on disk; creation took 180 s. Bioconductor 3.23 is the current release, but bioconda had no 3.23 builds of DESeq2/tximport on 2026-10-26, so 3.22 (the newest available on conda) is used. The versioned vignettes are cited in `docs/methods_references.md`.

## 3. Import script

[`scripts/R/01_tximport_import.R`](../scripts/R/01_tximport_import.R) (run with `conda activate stv-r; Rscript scripts/R/01_tximport_import.R`). It:

1. reads the workflow's sample index `results/quantification/salmon/samples.csv` and `config/samplesheet.csv` and checks that the sample IDs, conditions and the 3 + 3 design agree;
2. discovers the six `quant.sf` files and checks that their directories match the workflow index;
3. validates tx2gene against the `quant.sf` names;
4. runs `tximport(files, type = "salmon", tx2gene = tx2gene)` with defaults (`countsFromAbundance = "no"`, i.e. estimated counts + abundance-weighted average transcript length, which `DESeqDataSetFromTximport()` converts to an offset, per the tximport and DESeq2 vignettes);
5. writes [`results/quantification/tximport_validation.tsv`](../results/quantification/tximport_validation.tsv) and [`tximport_sessionInfo.txt`](../results/quantification/tximport_sessionInfo.txt); it stops on any failed check.

Result (all checks OK):

| Check | Value |
|---|---|
| Samples / per condition | 6 / STV_free = 3, STV_infected = 3 |
| `quant.sf` found and matching the workflow index | 6/6 |
| Gene matrix dimension | **33,976 genes × 6 samples** (= genes whose transcript is in `quant.sf`) |
| Gene-level totals vs Salmon `NumReads` totals | identical (max relative difference 0) |
| NA counts / non-positive lengths | 0 / 0 |
| Genes with zero counts in all six samples | 9,796 (reported only; no filtering applied) |

The script deliberately stops before `DESeqDataSetFromTximport()`. The DESeq2 design (proposed earlier: `~ condition`, with `STV_free` as reference level), any pre-filtering, the significance threshold, the LFC threshold/shrinkage method and the library-type choice (D-010) all require owner approval.
