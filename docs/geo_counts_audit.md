# GEO processed-data audit — `GSE137303_count.txt.gz`

- File: `data/processed/original_geo/GSE137303_count.txt.gz` (1,344,528 bytes; MD5 `3800ddd08f44eba1f062f034cce4d329`, stored in `.md5` sidecar)
- URL: <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137303/suppl/GSE137303_count.txt.gz> (retrieved 2026-10-25)
- Reproduce: `bash scripts/bash/fetch_geo_counts.sh` (download + structure) and `bash scripts/bash/verify_geo_counts_annotation.sh` (annotation check; streams an 8 MB NCBI GFF, saves nothing).
- **Not used for DESeq2.**

## What GEO says the file contains

From every GSM `data_processing` block: reads mapped with TopHat 2.0.13/Bowtie2 2.2.3 to "Tomato SL2.50 ITAG2.4"; FPKM computed with **Cuffnorm 2.2.1**; supplementary format/content: **"FPKM, Processed Signal, count"**. GEO does not state whether "count" is raw or normalized, nor what "Processed Signal" is.

## Structure (observed)

- Tab-separated, 1 header + **30,389 data rows**, 22 columns, no duplicated row IDs.
- Columns: `tracking_id`, `NCBI ID`, `Name`, `description`, then per sample `<label>-FPKM` and `<label>-Processed Signal` (labels 22, 412, 422, 12, 432, 433), then six count columns `22-count`, `412count`, `422-count`, `12-count`, `432-count`, `433-count`.
- First rows:

```text
tracking_id  NCBI ID           Name          description                       22-FPKM ... 22-count 412count ...
gene0        GeneID:104645797  LOC104645797                                    0.107   ... 1        0
gene1        GeneID:100147723  GAD3          glutamate decarboxylase isoform3  17.903  ... 385      233
gene10       GeneID:101261948  LOC101261948                                    28.912  ... 1087     760
```

- Value types: FPKM and Processed Signal columns are decimals (3 dp) in all rows; **all six count columns are integers in all rows**.
- Processed Signal range ≈ −13.3 to +13.3 (negative values present) — a transformed/relative scale; exact definition **UNKNOWN**.
- 6,551 rows are zero in all six count columns.

## Column → sample mapping (VERIFIED)

Labels match GEO `Sample_description` ("Sample N") and the original submitted SRA file names (`N.fastq`):

| Label | GSM | sample_id |
|---|---|---|
| 12 | GSM4074979 | STV_free_rep1 |
| 432 | GSM4074980 | STV_free_rep2 |
| 433 | GSM4074981 | STV_free_rep3 |
| 22 | GSM4074982 | STV_infected_rep1 |
| 412 | GSM4074983 | STV_infected_rep2 |
| 422 | GSM4074984 | STV_infected_rep3 |

## Identifier namespace — NCBI RefSeq, not ITAG2.4 (VERIFIED)

- `NCBI ID`: `GeneID:<n>` (30,319 rows) or quoted `"GeneID:<n>,miRBase:MI…"` (70 microRNA rows). `Name`: `LOC…` (28,056) or gene symbols (2,333). **No Solyc/ITAG identifiers anywhere** (0 matches).
- The NCBI RefSeq GFF for SL2.50 (`GCF_000188115.3_SL2.50_genomic.gff.gz`, header `#!annotation-source NCBI Solanum lycopersicum Annotation Release 102`) has exactly **30,389** gene/pseudogene features; every `tracking_id` equals a GFF `ID=` value, and **30,292/30,389** (`ID`, `GeneID`) pairs match exactly. The 97 remaining rows (first examples `gene30247`–`gene30251`) carry a GeneID that is offset relative to the GFF; not investigated further.
- Conclusion: the processed table was produced with the **SL2.50 assembly + NCBI RefSeq Annotation Release 102**, even though GEO text says "SL2.50 ITAG2.4". The assembly matches; the annotation/identifiers do not.

## Are the "count" columns raw counts? — evidence says NO

| Sample | count column total | sequenced reads (SRA) | ratio |
|---|---:|---:|---:|
| 22 | 16,123,093 | 19,230,359 | 0.84 |
| 412 | 17,337,359 | 16,858,900 | **1.03** |
| 422 | 14,417,165 | 19,862,460 | 0.73 |
| 12 | 15,418,310 | 18,080,482 | 0.85 |
| 432 | 15,510,050 | 16,288,249 | 0.95 |
| 433 | 15,682,287 | 14,318,362 | **1.10** |

- Two samples have more assigned "counts" than sequenced reads, and totals are compressed relative to library sizes (largest library scaled down, smallest scaled up) — the pattern expected from library-size-normalized counts.
- Cuffnorm documentation (<http://cole-trapnell-lab.github.io/cufflinks/cuffnorm/>, retrieved 2026-10-25): Cuffnorm reports FPKM and *"normalized estimates for the number of fragments … because these counts are already normalized to account for differences in library size, they should not be used with downstream differential expression tools that require raw counts as input."*
- The values being integers means they were rounded or otherwise post-processed; the exact procedure is UNKNOWN.

**Status:** count columns = INFERRED library-size-normalized (Cuffnorm-style) fragment counts, NOT VERIFIED raw counts. DESeq2 documentation requires un-normalized counts, so this table **must not be used as DESeq2 input**. Acceptable later uses: sanity cross-check of our own quantification (e.g. correlation of log expression after ID mapping), never as primary data.
