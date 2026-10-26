# Reference decision — tomato genome and annotation

Status: **APPROVED 2026-10-25 — SGN SL4.0 genome + ITAG4.0 annotation, full genome as Salmon decoys (D-001).** The research below (retrieved 2026-10-25) recommended ITAG4.1. The owner considered it and selected **ITAG4.0**, the annotation released together with SL4.0 build 4.00 (both 2019-09-06). Download, provenance and validation: `docs/reference_validation.md`.

Original research status line: research complete, decision pending. Nothing large was downloaded at that stage; only directory listings, README files and streamed file headers were read.

## Verified facts

### What the original study actually used

- GEO text (all six GSMs): reads mapped to "Tomato SL2.50 ITAG2.4" with TopHat/Bowtie2; FPKM with Cuffnorm.
- **But** the GEO processed table is keyed by **NCBI RefSeq Annotation Release 102 on SL2.50** (`GCF_000188115.3`): all 30,389 `tracking_id`s equal the GFF `ID=` values and 30,292 (ID, GeneID) pairs match exactly; no Solyc IDs are present (see `geo_counts_audit.md`). So "reproduce the original" is ambiguous: same *assembly* (SL2.50), but the published table uses NCBI LOC/GeneID models, not ITAG2.4 Solyc models.

### Available references

| Resource | Assembly | Annotation | Gene models | IDs | Source (retrieved 2026-10-25) |
|---|---|---|---:|---|---|
| SGN ITAG2.4 | SL2.50 | ITAG2.4 (files dated 2014) | 34,725 | `SolycCCgNNNNNN.V` | <https://solgenomics.net/ftp/tomato_genome/annotation/ITAG2.4_release/> (README) |
| NCBI RefSeq | SL2.50 (`GCF_000188115.3`, status *suppressed*, files still on FTP) | NCBI Annotation Release 102 | 30,389 gene/pseudogene features | `GeneID`, `LOC…`, `XM_/XR_` | <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.3_SL2.50/> (GFF header) |
| SGN ITAG3.2 | SL3.0 | ITAG3.2 (Jun 2017) | 35,768 | Solyc | SGN `ITAG3.2_release/` README |
| SGN ITAG4.0 | SL4.0 (long-read assembly) | ITAG4.0 (Sep 2019) | 34,075 | Solyc | SGN `ITAG4.0_release/` README |
| **SGN ITAG4.1** (latest SGN release found) | SL4.0 | ITAG4.1 (Jan 2020) | 34,688 | Solyc | <https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.1_release/> (README) |
| Ensembl Plants release 63 | SL4.0 (`GCA_000188115.5`) | "Community" gene set, genebuild updated 2025-06 | — | `gene-Solyc01g080690.3`, `mRNA-Solyc01g080690.3.1`; seq names `CM001064.4` | GTF header streamed from `ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/gtf/solanum_lycopersicum_gca000188115v5cm/` |
| NCBI RefSeq (current) | SLM_r2.1 (`GCF_036512215.1`), **cv. Micro-Tom** | RS_2024_09 | 41,476 genes | LOC/GeneID | NCBI Datasets API |
| SGN M82 genome | M82 1.3 (ONT/Canu, 2019) | MAKER annotation | — | non-ITAG | SGN `genomes/Solanum_lycopersicum/M82/README.txt` |

SGN change logs relative to ITAG2.40 (from the READMEs):

| Release | ITAG2.4 genes preserved | dropped | updated | novel |
|---|---:|---:|---:|---:|
| ITAG3.2 | 30,866 | 3,859 | 20,766 | 4,900 |
| ITAG4.0 | 29,281 | 5,444 | 21,375 | 4,794 |
| ITAG4.1 | 29,281 *(to be updated)* | 5,444 *(to be updated)* | 21,375 *(to be updated)* | 6,476 |

The ITAG4.1 README reuses the ITAG4.0 numbers marked "(to be updated)". SGN also ships `ITAG2.40_Solyc.dropped.ID.txt` (in ITAG3.2) and GO files per release (`ITAG2.4.go.csv`, `ITAG3.2_protein_go.tsv`, `ITAG4.0_goterms.txt`, `ITAG4.1_goterms.txt`).

Sizes relevant to Salmon (from listings): SL2.50 chromosomes `.fa.gz` 213 MB; SL4.0 `.fa.gz` 225 MB (759 MB uncompressed); ITAG2.4 cDNA 50 MB; ITAG4.1 cDNA 47 MB, gene models GFF 52 MB; NCBI SL2.50 `rna.fna.gz` 20 MB, GFF 8.1 MB.

## Options compared

**A. SL2.50 + ITAG2.4** (the build named in GEO) · **A2. SL2.50 + NCBI RefSeq R102** (the build behind the GEO table) · **B. SL4.0 + ITAG4.1** (current SGN release)

| Criterion | A: SL2.50 + ITAG2.4 | A2: SL2.50 + NCBI R102 | B: SL4.0 + ITAG4.1 |
|---|---|---|---|
| Reproducibility vs original | Same assembly and the annotation *named* in GEO; does **not** reproduce the GEO table IDs | Same assembly **and** same gene IDs as the GEO table → direct gene-by-gene comparison with authors' FPKM/counts | Re-analysis, not reproduction; comparison with authors only via ID mapping |
| Annotation quality | Oldest (2014); README: release "covering approximately 84% of the genome" | NCBI Gnomon models incl. lncRNA/miRNA; RefSeq record now suppressed | SL4.0 "assembled from long reads" (README); ~21k updated models, ~5k dropped, ~5–6k novel vs ITAG2.4 |
| Solyc ID compatibility | Solyc v.1/v.2 | **None** (LOC/GeneID) — weakest for plant tools | Solyc v.3 for updated models; locus part kept for "preserved" genes |
| Transcript-level quant (Salmon) | cDNA FASTA + GFF available; tx2gene from GFF | `rna.fna` + GFF; many `XM_/XR_` isoforms; tx2gene via GFF `Parent`/`Dbxref` | cDNA FASTA + GFF available; tx2gene from GFF |
| DESeq2 (gene level via tximport) | Straightforward | Straightforward; extra non-coding genes | Straightforward |
| agriGO v2 | Built-in support for ITAG2.x: TO BE VERIFIED; custom GO background possible (`ITAG2.4.go.csv`) | No Solyc → needs ID conversion first | Built-in tomato set announced is **ITAG3.2** (agriGO news 2018-09-03); ITAG4.x via custom background (`ITAG4.1_goterms.txt`) |
| MapMan / Mercator4 | Mercator4 annotates a supplied protein FASTA → usable with any release (UNVERIFIED today: plabipd.de returned "Service Unavailable") | Possible with RefSeq proteins (same caveat) | Same as A |
| Plant Reactome (Gramene) | Older IDs; mapping likely needed — UNVERIFIED | Needs ID conversion — UNVERIFIED | Ensembl Plants 63 is on SL4.0 with Solyc v.3 IDs, so B is the most likely match — UNVERIFIED for Plant Reactome itself |
| Old ↔ new model mapping | ↔ B: locus-ID match for preserved genes; dropped/split/merged genes need sequence- or coordinate-based mapping | ↔ A: coordinate overlap on the same SL2.50 assembly is feasible; ↔ B: sequence-based (e.g. reciprocal BLAST) | ↔ A: as left; ↔ A2: sequence-based |
| Risk of split/merged/removed/renamed genes | Baseline | Different model set altogether (LOC vs Solyc) | ~5.4k ITAG2.4 genes dropped and ~21k updated: some DE genes from older studies will be absent or redefined |
| Effort / risk in 2 weeks | Low | Medium (ID conversion before every plant tool) | Low–medium (modern IDs; possible `gene-`/`mRNA-` prefixes and chromosome-name differences if mixing SGN and Ensembl files) |
| Salmon + tximport practical? | Yes | Yes | Yes |

Out of scope: ITAG3.2/SL3.0 (superseded by 4.x, although agriGO v2 has it built in); NCBI SLM_r2.1 (different cultivar, LOC IDs); SGN M82 assembly (matches the cultivar but non-standard IDs and minimal functional annotation).

Not verified and to be checked after a decision: number of transcripts per gene in each annotation (determines how much isoform-level information Salmon actually adds); whether Ensembl Plants/Gramene "Community" models equal ITAG4.0 or ITAG4.1.

## BLOCKING DECISION: REFERENCE VERSION

**Recommendation (for approval): Option B — SL4.0 genome + ITAG4.1 annotation from SGN, all files from the single SGN source.**

Rationale:

1. The project is a *re-analysis* with plant-specific functional analysis; B gives the most complete assembly and current Solyc models, and the current SGN GO file (`ITAG4.1_goterms.txt`) plus Ensembl Plants alignment on SL4.0.
2. Option A does not actually buy comparability with the published table (that table uses NCBI R102 IDs), so its main argument — reproducibility — is weaker than it looks.
3. Salmon + tximport + DESeq2 work identically for A, A2 and B; the choice only changes identifiers and gene models.
4. Comparison with the authors remains possible as a secondary sanity check by mapping our genes to the GEO table (A2 IDs) at the end.

Fallbacks: choose **A2** if the course requires gene-by-gene reproduction of the authors' table; choose **A** if the instructor requires literally "SL2.50 + ITAG2.4" as stated in GEO.

Researcher must approve one of: **B (recommended) / A / A2**, and also confirm the use of the full SL4.0 genome as Salmon decoys (≈0.8 Gb; see RAM note in `pipeline_design.md`).

## Outcome (2026-10-25)

- Approved: **SL4.0 + ITAG4.0** (a variant of option B with the earlier annotation release), full-genome decoys. ITAG4.1 was not selected.
- Points left open above, now answered for ITAG4.0: **1 transcript per gene** (34,075 genes = 34,075 mRNAs, so no isoform-level information). Identifiers are SGN Solyc `Name`s; GFF `ID`s carry `gene:`/`mRNA:` prefixes, and sequence names are `SL4.0ch00`–`12`.
- Measured RAM for the full-decoy index: 4.7 GB peak (`docs/salmon_quantification.md`).
- Consequences carried forward: `ITAG4.0_goterms.txt` uses gene versions inconsistent with the GFF for 21,311 IDs. This must be reconciled before GO analysis (`docs/reference_validation.md` §2).
