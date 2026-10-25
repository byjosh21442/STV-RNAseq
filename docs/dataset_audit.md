# Dataset audit — GSE137303 / PRJNA565137

Retrieved: 2026-10-25 (UTC timestamp in `data/metadata/sources/retrieved_utc.txt`).
Reproduce: `bash scripts/bash/fetch_metadata.sh && bash scripts/bash/build_samples_metadata.sh` (run from the project root in WSL).

## Sources (raw records saved under `data/metadata/sources/`)

| Source | URL | Saved as |
|---|---|---|
| NCBI GEO series + samples (SOFT, brief) | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137303&targ=all&form=text&view=brief> | `GSE137303_brief.soft.txt` |
| NCBI SRA run info (E-utilities, `rettype=runinfo`) | <https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&rettype=runinfo> (ids from `esearch db=sra term=PRJNA565137[BioProject]`) | `PRJNA565137_sra_runinfo.csv` |
| NCBI SRA experiment packages (XML) | same E-utilities `efetch db=sra retmode=xml` | `PRJNA565137_sra_experiment_packages.xml` |
| ENA read_run report | <https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA565137&result=read_run&format=tsv> (explicit field list in script) | `PRJNA565137_ena_read_run.tsv` |
| GEO supplementary listing | <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137303/suppl/> | `GSE137303_suppl_listing.html` |

`build_samples_metadata.sh` joins ENA ↔ GEO on the GSM accession and fails if run/experiment/BioSample disagree between ENA, NCBI runinfo and GEO. All three agreed.

## Series-level facts (VERIFIED, GEO)

- Title: *Effect on tomato plants of asymptomatic infection with southern tomato virus*; PubMed 31620899.
- Overall design: three 1-month-old STV-infected and three STV-free *S. lycopersicum* cv. M82 plants.
- Platform GPL19694 (Illumina HiSeq 2500, *S. lycopersicum*); SRA study SRP221436; BioProject PRJNA565137.
- Supplementary file: `GSE137303_count.txt.gz` (1.3 MB) — audited in `geo_counts_audit.md`.
- Sample `data_processing` (identical for all six): Casava 1.8.2 base calling; TopHat 2.0.13 + Bowtie2 2.2.3 + SAMtools 0.1.19 against "Tomato SL2.50 ITAG2.4"; FPKM by Cuffnorm 2.2.1; supplementary content "FPKM, Processed Signal, count".

## Six-sample mapping (all fields VERIFIED; full table: `data/metadata/samples_metadata.tsv`)

| sample_id | Condition (GEO "virus infection") | GSM | BioSample | SRS | SRX | SRR | Original label (GEO description / submitted file) | Reads | Bases | ENA fastq.gz bytes |
|---|---|---|---|---|---|---|---|---:|---:|---:|
| STV_free_rep1 | STV-free | GSM4074979 | SAMN12735328 | SRS5376422 | SRX6834734 | SRR10102791 | Sample 12 / `12.fastq` | 18,080,482 | 1,356,036,150 | 988,076,517 |
| STV_free_rep2 | STV-free | GSM4074980 | SAMN12735327 | SRS5376423 | SRX6834735 | SRR10102792 | Sample 432 / `432.fastq` | 16,288,249 | 1,221,618,675 | 896,430,396 |
| STV_free_rep3 | STV-free | GSM4074981 | SAMN12735326 | SRS5376424 | SRX6834736 | SRR10102793 | Sample 433 / `433.fastq` | 14,318,362 | 1,073,877,150 | 780,172,066 |
| STV_infected_rep1 | STV-infected | GSM4074982 | SAMN12735325 | SRS5376425 | SRX6834737 | SRR10102794 | Sample 22 / `22.fastq` | 19,230,359 | 1,442,276,925 | 1,049,513,467 |
| STV_infected_rep2 | STV-infected | GSM4074983 | SAMN12735324 | SRS5376426 | SRX6834738 | SRR10102795 | Sample 412 / `412.fastq` | 16,858,900 | 1,264,417,500 | 921,773,507 |
| STV_infected_rep3 | STV-infected | GSM4074984 | SAMN12735323 | SRS5376427 | SRX6834739 | SRR10102796 | Sample 422 / `422.fastq` | 19,862,460 | 1,489,684,500 | 1,082,142,559 |

`sample_id` is a project label built from GEO "virus infection" + "repN" in the GEO title. One run per sample (no technical re-sequencing).

Common to all six (VERIFIED in GEO + SRA/ENA): organism *Solanum lycopersicum* (taxid 4081); cultivar M82; tissue leaf; age 1-month-old; GEO molecule "total RNA"; library strategy RNA-Seq; source TRANSCRIPTOMIC; selection cDNA; instrument Illumina HiSeq 2500; protocol text "TruSeq stranded mRNA sample prep kit" (GEO extract protocol and SRA `LIBRARY_CONSTRUCTION_PROTOCOL`).

## Layout and read length

- **Library layout: SINGLE** — VERIFIED in GEO-linked SRA XML (`<LIBRARY_LAYOUT><SINGLE/>`), NCBI runinfo (`LibraryLayout=SINGLE`, `spots_with_mates=0`), and ENA (`library_layout=SINGLE`).
- **Read length: 75 nt, fixed** — VERIFIED: SRA XML `Statistics nreads="1" … average="75" stdev="0"` for every run; runinfo `avgLength=75`; ENA `base_count / read_count = 75`.

## Data volume

| Quantity | Value | Source |
|---|---:|---|
| Total reads | 104,638,812 | ENA/SRA |
| Total bases | 7,847,910,900 | ENA/SRA |
| ENA `fastq.gz` total (direct download, MD5s available) | 5,718,108,512 bytes (≈5.7 GB) | ENA `fastq_bytes` |
| NCBI SRA normalized total | ≈3.9 GB | SRA XML `RUN size` |
| Original submitted (uncompressed) FASTQ total | 21,607,287,754 bytes (≈21.6 GB) | SRA XML `SRAFile supertype="Original"` |

Practical implication: downloading ENA `fastq.gz` needs ≈5.7 GB; `fasterq-dump` from SRA would produce ≈21–22 GB of uncompressed FASTQ before compression.

## Strandedness

- Metadata states the kit (TruSeq stranded mRNA) but **no source states the read orientation**. No GEO/SRA/ENA field records strandedness.
- Status: **UNKNOWN / TO BE VERIFIED**.
- INFERRED (not used as a parameter): the kit name implies a stranded library; for single-end reads the expected Salmon library type would be `SR` (TopHat `fr-firststrand`). This must be confirmed empirically before quantification (Salmon `-l A` auto-detection on a subsample and/or RSeQC `infer_experiment.py`), as nf-core/rnaseq does with `strandedness=auto`.

## VERIFIED / INFERRED / UNKNOWN summary

| Field | Status | Note |
|---|---|---|
| GSM, BioSample, SRS, SRX, SRR, condition, replicate | VERIFIED | three sources agree |
| Organism, cultivar, tissue, age | VERIFIED | GEO characteristics |
| Library strategy/source/selection, platform/instrument | VERIFIED | GEO + SRA + ENA |
| Layout SINGLE, read length 75 nt | VERIFIED | see above |
| Reads, bases, file sizes, ENA MD5 | VERIFIED | ENA/SRA |
| Kit "TruSeq stranded mRNA" | VERIFIED (as reported text) | — |
| Read orientation (SR/SF/U) | UNKNOWN | expected `SR` is INFERRED only |
| Fragment length distribution | UNKNOWN | not estimable from single-end data; relevant to Salmon `--fldMean/--fldSD` |
| Adapter sequences / adapter content | UNKNOWN | to be assessed by FastQC |
| Genome build actually used for the GEO count table | VERIFIED ≠ GEO text | table is keyed by NCBI RefSeq Annotation Release 102 IDs, not ITAG2.4 (see `geo_counts_audit.md`) |

## Notes / minor inconsistencies

- GEO `molecule` is "total RNA" while the library kit is an mRNA kit; this describes the extracted material and is not necessarily contradictory. Recorded as reported.
- SRA `LibrarySelection=cDNA` (not `PolyA`), as submitted.
- The GEO count file header has `412count` (no hyphen) while the other columns use `<label>-count`.
