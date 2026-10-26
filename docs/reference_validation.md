# Reference acquisition and validation — SGN SL4.0 + ITAG4.0

Decision **D-001 (approved)**: SL4.0 genome assembly + ITAG4.0 annotation, both from the Sol Genomics Network (SGN). Downloaded 2026-10-25 local time (UTC−6) = **2026-10-26 UTC** (the date recorded in the manifest).

## 1. Acquisition

Script: [`scripts/bash/download_reference.sh`](../scripts/bash/download_reference.sh). It fetches each file by HTTPS from `https://solgenomics.net/ftp/tomato_genome/`, resumes partial downloads (`curl -C -`, 5 attempts), checks the local size against the server `Content-Length`, runs `gzip -t` on the compressed genome, and records SHA-256. The run took 111 s; 335 MB were stored under `reference/SGN/`.

Provenance: [`reference/reference_manifest.tsv`](../reference/reference_manifest.tsv) with columns file, source_url, genome_version, annotation_version, size_bytes, sha256, upstream_checksum, upstream_last_modified, retrieval_date, purpose.

| File (under `reference/SGN/`) | Bytes | SHA-256 (first 12) | Server Last-Modified | Purpose |
|---|---:|---|---|---|
| `SL4.0/S_lycopersicum_chromosomes.4.00.fa.gz` | 235,901,169 | `365255b77847` | 2019-09-06 | Genome FASTA → Salmon decoys |
| `SL4.0/FileDescription.txt` | 462 | `99864da6b01c` | 2019-09-06 | Assembly provenance |
| `SL4.0/ReleaseNotes.txt` | 15,705 | `857189ab1916` | 2019-09-06 | Assembly provenance |
| `ITAG4.0/ITAG4.0_README.txt` | 1,638 | `04c5b74d8355` | 2019-09-06 | Annotation provenance |
| `ITAG4.0/ITAG4.0_cDNA.fasta` | 45,858,373 | `2956f6f4c76c` | 2019-09-06 | Transcript FASTA → Salmon targets |
| `ITAG4.0/ITAG4.0_gene_models.gff` | 53,605,921 | `61bdf7c0a607` | 2019-09-06 | Gene models → tx2gene, validation |
| `ITAG4.0/ITAG4.0_proteins.fasta` | 12,380,140 | `84c15ec4bb6b` | 2019-09-06 | Later functional annotation |
| `ITAG4.0/ITAG4.0_descriptions.txt` | 2,554,074 | `2ed164a5d301` | 2021-03-14 | Later functional annotation |
| `ITAG4.0/ITAG4.0_goterms.txt` | 710,463 | `b86324a2285a` | 2021-01-21 | Later GO analysis (see §3) |

- **Upstream checksums:** SGN publishes no checksum files for these directories (listing inspected). The manifest therefore says `none published by SGN`, and the SHA-256 values were computed locally after the size check. Server `Last-Modified` dates are stored as the only upstream integrity metadata.
- **Pairing confirmed upstream:** `FileDescription.txt` names ITAG4.0 as the corresponding annotation of build 4.00. The ITAG4.0 README states "ITAG4.0 release of new annotation on the tomato genome assembly SL4.0 … 34,075 gene models", release date Sep 6, 2019.
- **Not downloaded (not needed for Phase 3):** `ITAG4.0_CDS.fasta`, `ITAG4.0_lncRNA.gff`, repeat GFFs, SolCAP/SGN-marker GFFs. lncRNAs are not in the cDNA FASTA, so they are not quantified (see limitations in `docs/salmon_quantification.md`).
- All files are git-ignored. Only the manifest and the validation tables are tracked.

## 2. Compatibility validation

Script: [`scripts/bash/validate_reference.sh`](../scripts/bash/validate_reference.sh), read-only. It exits 1 on any material inconsistency; the full table is [`reference/reference_validation.tsv`](../reference/reference_validation.tsv). Result: **no material inconsistency** (exit 0).

| Check | Result |
|---|---|
| GFF feature types | gene 34,075; mRNA 34,075; exon 161,508; CDS 153,474; five_prime_UTR 21,452; three_prime_UTR 22,835 |
| Genes / transcripts (mRNA) | **34,075 / 34,075**, matching the README |
| Transcripts mapped to a gene (`Parent` exists) | 34,075 (100 %) |
| Unmapped transcripts (missing parent) | **0** |
| Duplicate gene IDs / mRNA IDs / gene Names / mRNA Names | **0 / 0 / 0 / 0** |
| Ambiguous tx→gene (multiple parents) | **0** |
| Transcripts per gene | 1 for every gene (ITAG4.0 has no alternative isoforms) |
| Genes without mRNA | 0 |
| mRNA Name prefixed by its gene Name | 34,075 / 34,075 |
| mRNA on same sequence/strand as gene | all |
| GFF `ID` format | `gene:<Name>` / `mRNA:<Name>` for every feature |
| cDNA records / duplicate cDNA IDs | 34,075 / 0 |
| cDNA IDs ↔ GFF mRNA `Name` | identical sets (0 missing either way) |
| cDNA length = sum of GFF exon lengths | true for all 34,075 transcripts |
| cDNA length min / median / max | 41 / 1,086 / 22,821 nt |
| Genome sequences | 13 (`SL4.0ch00`–`SL4.0ch12`), 782,520,033 bp, no duplicates |
| GFF seqids missing from the genome | 0 |
| `##sequence-region` lengths = genome lengths | all 13 match |
| Features beyond sequence end / start > end | 0 / 0 |
| Protein records | 34,075, IDs = mRNA Names |
| Descriptions file | 34,075 rows, IDs = gene Names |

### Non-material notes (recorded, nothing repaired)

1. **15 gene models carry an extra version level:** e.g. gene `Solyc01g008710.5.1` → mRNA `Solyc01g008710.5.1.1` (also `Solyc01g090260`, `…099150`, `…100030`, `…102380`, `…110440`, `Solyc02g076980`, `…085180`, `Solyc03g093120`, `Solyc06g007130`, `…051630`, `Solyc08g075540`, `Solyc09g010090`, `…010800`, `…011640`). They are unique, parent-consistent and spelled identically in the cDNA FASTA, so they are kept exactly as published. No locus has more than one gene model.
2. **Cosmetic formatting:** 8 genome FASTA headers have trailing whitespace, and the GFF `##sequence-region` pragma lacks a space (`##sequence-regionSL4.0ch00`). Salmon and `decoys.txt` use the first whitespace-delimited token, and the validator parses the pragma correctly.
3. **`ITAG4.0_goterms.txt` is not version-consistent with the GFF.** Of 33,562 rows, only 5,845 carry GO terms. 21,962 IDs are not ITAG4.0 gene Names; 21,311 of those are the same locus with a different version (e.g. `Solyc01g005000.2` vs GFF `.3`). Only 1,458 IDs carrying GO terms match a gene Name exactly. The file is **not used in Phase 3**; any GO analysis must first reconcile these identifiers (owner decision, later phase).

## 3. Identifier policy

- Transcript IDs used everywhere = cDNA header = GFF mRNA `Name` (e.g. `Solyc00g500001.1.1`). Gene IDs = GFF gene `Name` (e.g. `Solyc00g500001.1`).
- The GFF `ID` type prefixes (`gene:`, `mRNA:`) are not used, because the `Name` attribute already equals the FASTA identifier. No string surgery and **no version stripping**.
