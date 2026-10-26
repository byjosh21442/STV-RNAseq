# Methods references — STV_GSE137303

One entry per process. **Technical documentation** (manuals, vignettes, release notes, file listings) is kept separate from **scientific literature** (peer-reviewed articles, preprints). All DOIs were resolved against Crossref (journal articles, preprint) or DataCite (Zenodo) on 2026-10-25. Access dates are local time (UTC−6). APA 7 references are provisional and must be checked against the final citation style before submission. Complementary source list (design-phase URLs): [`docs/pipeline_sources.md`](pipeline_sources.md).

Reference hierarchy used for technical decisions: nf-core/rnaseq > Nextflow docs/training > Salmon docs > tximport vignette > DESeq2 vignette > SGN.

---

## 1. Dataset and original study

| Field | Content |
|---|---|
| Process | Definition of the dataset re-analysed (GEO GSE137303, BioProject PRJNA565137, 6 SE runs SRR10102791–96) |
| Software + version | — (metadata retrieved with `curl` from GEO/SRA/ENA APIs, `scripts/bash/fetch_metadata.sh`) |
| Parameters | Accessions above; metadata copies in `data/metadata/sources/` |
| Methodological rationale | Public, MD5-verifiable raw reads of M82 leaves, STV-free vs STV-infected, 3 + 3 |
| Official technical documentation | GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137303> |
| Primary scientific reference | Fukuhara et al. (2020), the study that generated the data; Edgar et al. (2002) for GEO |
| DOI | 10.1007/s00705-019-04436-1; 10.1093/nar/30.1.207 |
| Provisional APA 7 reference | Fukuhara, T., Tabara, M., Koiwa, H., & Takahashi, H. (2020). Effect of asymptomatic infection with southern tomato virus on tomato plants. *Archives of Virology, 165*(1), 11–20. https://doi.org/10.1007/s00705-019-04436-1 <br> Edgar, R., Domrachev, M., & Lash, A. E. (2002). Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. *Nucleic Acids Research, 30*(1), 207–210. https://doi.org/10.1093/nar/30.1.207 |
| Decision supported | Dataset choice; D-005 (authors' count table used for audit only) |
| Access info | GEO page and PubMed 31620899 (NCBI E-utilities esummary), 2026-10-25 |

## 2. Raw-data acquisition (ENA)

| Field | Content |
|---|---|
| Process | Download of the six `fastq.gz` files and integrity verification |
| Software + version | curl 8.18.0; md5sum (uutils coreutils 0.8.0); `scripts/bash/download_fastq.sh` |
| Parameters | HTTPS URL = `https://` + ENA `fastq_ftp`; resume `-C -`; size and MD5 compared with ENA `fastq_bytes` / `fastq_md5` |
| Methodological rationale | ENA serves gzip FASTQ (5.7 GB) with published MD5s, smaller than SRA conversion (≈ 21.6 GB) and needing no extra tool |
| Official technical documentation | ENA "How to Download Data Files": <https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html> (MD5 files, curl/wget access); ENA Portal API: <https://www.ebi.ac.uk/ena/portal/api/> |
| Primary scientific reference | Leinonen et al. (2011) |
| DOI | 10.1093/nar/gkq967 |
| Provisional APA 7 reference | Leinonen, R., Akhtar, R., Birney, E., Bower, L., Cerdeno-Tarraga, A., Cheng, Y., Cleland, I., Faruque, N., Goodgame, N., Gibson, R., Hoad, G., Jang, M., Pakseresht, N., Plaister, S., Radhakrishnan, R., Reddy, K., Sobhany, S., Ten Hoopen, P., Vaughan, R., … Cochrane, G. (2011). The European Nucleotide Archive. *Nucleic Acids Research, 39*(Database issue), D28–D31. https://doi.org/10.1093/nar/gkq967 |
| Decision supported | D-007 (ENA FASTQ + MD5 instead of SRA conversion) |
| Access info | ENA docs page, 2026-10-25; ENA API used 2026-10-25 (Phase 1–2) |

## 3. Read sanity checks

| Field | Content |
|---|---|
| Process | Read counts, length distribution and structure of the raw FASTQ |
| Software + version | SeqKit 2.14.0 (bioconda); `scripts/bash/check_raw_fastq.sh` |
| Parameters | `seqkit stats -a`; gzip integrity (`pigz -t`); first 10,000 records structure check |
| Methodological rationale | Independent confirmation that the files are complete (read count = ENA) and single-end 75 nt before QC |
| Official technical documentation | <https://bioinf.shenwei.me/seqkit/> |
| Primary scientific reference | Shen et al. (2016); the SeqKit site now asks users to cite SeqKit2 (Shen, Sipos & Zhao, 2024) |
| DOI | 10.1371/journal.pone.0163962; 10.1002/imt2.191 |
| Provisional APA 7 reference | Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLOS ONE, 11*(10), e0163962. https://doi.org/10.1371/journal.pone.0163962 <br> Shen, W., Sipos, B., & Zhao, L. (2024). SeqKit2: A Swiss army knife for sequence and alignment processing. *iMeta, 3*(3), e191. https://doi.org/10.1002/imt2.191 |
| Decision supported | Raw-data validation (`docs/raw_data_validation.md`) |
| Access info | SeqKit docs, 2026-10-25 |

## 4. Read quality control — FastQC

| Field | Content |
|---|---|
| Process | Per-sample raw read QC |
| Software + version | FastQC 0.13.0 (bioconda; same version as nf-core/modules) |
| Parameters | `--threads 1 --memory <80 % of task memory>`; reports renamed to sample IDs (`workflow/modules/fastqc.nf`) |
| Methodological rationale | Standard read-level QC; module flags interpreted for RNA-seq (`docs/raw_qc_assessment.md`) |
| Official technical documentation | <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/> (module help pages) |
| Primary scientific reference | Software without a journal article; cited as software (Andrews, 2010) |
| DOI | none |
| Provisional APA 7 reference | Andrews, S. (2010). *FastQC: A quality control tool for high throughput sequence data* [Computer software]. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| Decision supported | D-003 / D-009 (no trimming) |
| Access info | FastQC project page, 2026-10-25 |

## 5. QC aggregation — MultiQC

| Field | Content |
|---|---|
| Process | Aggregation of FastQC (Phase 2) and FastQC + Salmon (Phase 3) reports |
| Software + version | MultiQC 1.35 (bioconda; same version as nf-core/modules) |
| Parameters | `--force --title <title> --filename multiqc_report.html` |
| Methodological rationale | Cross-sample view of QC metrics; parses Salmon `meta_info.json`/`lib_format_counts.json` |
| Official technical documentation | <https://docs.seqera.io/multiqc> |
| Primary scientific reference | Ewels et al. (2016) |
| DOI | 10.1093/bioinformatics/btw354 |
| Provisional APA 7 reference | Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. *Bioinformatics, 32*(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354 |
| Decision supported | Raw-QC conclusion A (D-009); Phase 3 QC overview |
| Access info | MultiQC docs, 2026-10-25 |

## 6. Workflow manager — Nextflow

| Field | Content |
|---|---|
| Process | Orchestration of FastQC → Salmon index/quant → MultiQC with caching/resume and provenance reports |
| Software + version | Nextflow 26.10.3 (build 14229), strict DSL2 syntax; OpenJDK 21.0.10 |
| Parameters | `workflow/nextflow.config` (executor 12 CPU / 12 GB; labels; `ext.args`), typed `params` in `workflow/main.nf` |
| Methodological rationale | Reproducible, resumable execution with per-task hashes, trace (runtime/peak RSS) and versions topic |
| Official technical documentation | <https://docs.seqera.io/nextflow/> (install, strict syntax, typed params, workflow outputs, cache/resume); training: <https://training.nextflow.io/latest/nf4_science/rnaseq/> |
| Primary scientific reference | Di Tommaso et al. (2017) |
| DOI | 10.1038/nbt.3820 |
| Provisional APA 7 reference | Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology, 35*(4), 316–319. https://doi.org/10.1038/nbt.3820 |
| Decision supported | D-002 (workflow architecture), D-008 (provisioning) |
| Access info | Nextflow docs/training, 2026-10-25 |

## 7. Architecture reference — nf-core/rnaseq (read, not executed)

| Field | Content |
|---|---|
| Process | Design template for samplesheet, modules (FastQC, Salmon index/quant, MultiQC), strandedness inference and decoy-aware Salmon index |
| Software + version | nf-core/rnaseq 3.27.0; nf-core/modules `fastqc`, `multiqc`, `salmon/index`, `salmon/quant` (master, pins FastQC 0.13.0, MultiQC 1.35, Salmon 2.7.0) |
| Parameters | Not run; differences documented in `docs/nextflow_pipeline.md` |
| Methodological rationale | Community-curated best practice; highest level in the reference hierarchy |
| Official technical documentation | <https://nf-co.re/rnaseq/3.27.0/docs/usage>; <https://github.com/nf-core/modules/tree/master/modules/nf-core/salmon> |
| Primary scientific reference | Ewels et al. (2020) for the nf-core framework; software DOI (Zenodo concept DOI, all versions) for nf-core/rnaseq |
| DOI | 10.1038/s41587-020-0439-x; 10.5281/zenodo.1400710 |
| Provisional APA 7 reference | Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines. *Nature Biotechnology, 38*(3), 276–278. https://doi.org/10.1038/s41587-020-0439-x <br> Patel, H., Manning, J., Ewels, P., Garcia, M. U., et al. (2026). *nf-core/rnaseq: nf-core/rnaseq v3.27.0 – Vanadium Vulture* (Version 3.27.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.1400710 — *provisional: the Zenodo record lists 30 creators; complete the author list before submission.* |
| Decision supported | D-002, D-004 (strandedness inference), decoy-aware index, module versions |
| Access info | nf-co.re usage page and GitHub modules, 2026-10-25; DataCite record, 2026-10-25 |

## 8. Software provisioning — conda / Bioconda

| Field | Content |
|---|---|
| Process | Pinned, isolated environments `stv-rnaseq` and `stv-r` |
| Software + version | Miniconda (conda); channels conda-forge + bioconda, strict priority per command |
| Parameters | `envs/stv-rnaseq.yml` + lock (152 pkgs); `envs/stv-r.yml` + lock (156 pkgs: 139 at Phase 3, + 13 in Phase 4 for apeglm/plotting, + 4 in Phase 5 for rrcov) |
| Methodological rationale | Exact, rebuildable tool versions without touching global config or other envs |
| Official technical documentation | <https://bioconda.github.io/> |
| Primary scientific reference | Grüning et al. (2018) |
| DOI | 10.1038/s41592-018-0046-7 |
| Provisional APA 7 reference | Grüning, B., Dale, R., Sjödin, A., Chapman, B. A., Rowe, J., Tomkins-Tinch, C. H., Valieris, R., Köster, J., & The Bioconda Team. (2018). Bioconda: Sustainable and comprehensive software distribution for the life sciences. *Nature Methods, 15*(7), 475–476. https://doi.org/10.1038/s41592-018-0046-7 |
| Decision supported | D-008 |
| Access info | 2026-10-25 |

## 9. Reference genome — SGN SL4.0

| Field | Content |
|---|---|
| Process | Genome assembly used as Salmon decoys and as coordinate system of the annotation |
| Software + version | SGN build 4.00 (`S_lycopersicum_chromosomes.4.00.fa.gz`, released 2019-09-06); `scripts/bash/download_reference.sh` |
| Parameters | 13 sequences (`SL4.0ch00`–`12`), 782,520,033 bp; SHA-256 in `reference/reference_manifest.tsv` |
| Methodological rationale | Current SGN long-read assembly paired with ITAG4.0 (owner decision D-001) |
| Official technical documentation | <https://solgenomics.net/ftp/tomato_genome/assembly/build_4.00/> (`FileDescription.txt`, `ReleaseNotes.txt`) |
| Primary scientific reference | Hosmani et al. (2019, preprint) for SL4.0/ITAG4.0; Fernandez-Pozo et al. (2015) for SGN; Tomato Genome Consortium (2012) for the original reference genome project |
| DOI | 10.1101/767764; 10.1093/nar/gku1195; 10.1038/nature11119 |
| Provisional APA 7 reference | Hosmani, P. S., Flores-Gonzalez, M., van de Geest, H., Maumus, F., Bakker, L. V., Schijlen, E., van Haarst, J., Cordewener, J., Sanchez-Perez, G., Peters, S., Fei, Z., Giovannoni, J. J., Mueller, L. A., & Saha, S. (2019). *An improved de novo assembly and annotation of the tomato reference genome using single-molecule sequencing, Hi-C proximity ligation and optical maps* [Preprint]. bioRxiv. https://doi.org/10.1101/767764 <br> Fernandez-Pozo, N., Menda, N., Edwards, J. D., Saha, S., Tecle, I. Y., Strickler, S. R., Bombarely, A., Fisher-York, T., Pujar, A., Foerster, H., Yan, A., & Mueller, L. A. (2015). The Sol Genomics Network (SGN)—from genotype to phenotype to breeding. *Nucleic Acids Research, 43*(D1), D1036–D1041. https://doi.org/10.1093/nar/gku1195 <br> The Tomato Genome Consortium. (2012). The tomato genome sequence provides insights into fleshy fruit evolution. *Nature, 485*(7400), 635–641. https://doi.org/10.1038/nature11119 |
| Decision supported | D-001; decoy-aware index |
| Access info | SGN FTP listings and files, 2026-10-25 (download 2026-10-26 UTC) |

Note: Hosmani et al. (2019) is a preprint (not peer reviewed at the time of access); Crossref now lists it under the container "openRxiv".

## 10. Annotation — SGN ITAG4.0

| Field | Content |
|---|---|
| Process | Gene models, transcript (cDNA) sequences and tx2gene source |
| Software + version | ITAG4.0 release (2019-09-06): `ITAG4.0_gene_models.gff`, `ITAG4.0_cDNA.fasta` (+ proteins, descriptions, GO terms for later phases) |
| Parameters | 34,075 genes / 34,075 mRNAs; validation in `docs/reference_validation.md` |
| Methodological rationale | Annotation paired with SL4.0; current Solyc identifiers. ITAG4.1 was considered and not selected by the owner. |
| Official technical documentation | <https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.0_release/> (`ITAG4.0_README.txt`) |
| Primary scientific reference | Hosmani et al. (2019, preprint); Fernandez-Pozo et al. (2015) |
| DOI | 10.1101/767764; 10.1093/nar/gku1195 |
| Provisional APA 7 reference | See §9 (same references) |
| Decision supported | D-001; tx2gene identifier policy (no version stripping) |
| Access info | SGN FTP, 2026-10-25 (download 2026-10-26 UTC) |

## 11. Salmon index (decoy-aware)

| Field | Content |
|---|---|
| Process | Selective-alignment index over ITAG4.0 cDNA with full SL4.0 genome decoys |
| Software + version | Salmon 2.7.0 (bioconda; Rust implementation; index format v2) |
| Parameters | `salmon index -t gentrome.fa --decoys decoys.txt -k 31 --threads 8 --ramLimit 6`; duplicates collapsed (default) |
| Methodological rationale | Decoys absorb reads from unannotated genomic loci that resemble annotated transcripts; k = 31 adequate for 75 nt reads |
| Official technical documentation | <https://combine-lab.github.io/salmon/getting-started/quick-start/>; <https://combine-lab.github.io/salmon/reference/cli/>; <https://github.com/COMBINE-lab/salmon/blob/master/MIGRATION.md>; <https://github.com/COMBINE-lab/salmon/releases/tag/v2.7.0>; 1.x guidance: <https://salmon.readthedocs.io/en/latest/salmon.html> |
| Primary scientific reference | Patro et al. (2017); Srivastava et al. (2020) for selective alignment and decoys |
| DOI | 10.1038/nmeth.4197; 10.1186/s13059-020-02151-8 |
| Provisional APA 7 reference | Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods, 14*(4), 417–419. https://doi.org/10.1038/nmeth.4197 <br> Srivastava, A., Malik, L., Sarkar, H., Zakeri, M., Almodaresi, F., Soneson, C., Love, M. I., Kingsford, C., & Patro, R. (2020). Alignment and mapping methodology influence transcript abundance estimation. *Genome Biology, 21*(1), 239. https://doi.org/10.1186/s13059-020-02151-8 |
| Decision supported | D-002; index construction (`docs/salmon_quantification.md` §2) |
| Access info | Salmon 2.x docs, GitHub README/MIGRATION/release notes, readthedocs 1.x page: 2026-10-25 |

## 12. Salmon quantification

| Field | Content |
|---|---|
| Process | Transcript quantification of the six SE samples. **12a. Exploratory run** (library-type inference, archived) and **12b. Final run** (used downstream). |
| Software + version | Salmon 2.7.0 |
| Parameters | 12a exploratory: `--libType A -r <fastq> --fldMean 250 --fldSD 25 --threads 6`, no bias flags (run `exotic_ardinghelli`). 12b **final**: `--libType SR -r <fastq> --fldMean 250 --fldSD 25 --gcBias --threads 6` (run `sleepy_mcclintock`). Both: selective alignment, decoy-aware index, deterministic mode. |
| Methodological rationale | `-l A` measured the undocumented strandedness (SR in 6/6); the final run fixes `SR` explicitly. `--gcBias` follows the DESeq2 vignette recommendation in the absence of evidence that it is unnecessary; `--seqBias` not used. SE fragment length is not observable, so the explicit default prior is used (SRA `InsertSize = 0`). |
| Official technical documentation | <https://combine-lab.github.io/salmon/guides/library-types/>; <https://combine-lab.github.io/salmon/guides/mapping-modes/>; <https://combine-lab.github.io/salmon/reference/cli/> (`--gcBias`); <https://combine-lab.github.io/salmon/reference/output-formats/> (`lib_format_counts.json`, `meta_info.json`); `--gcBias` recommendation: DESeq2 vignette 3.22 (§15) |
| Primary scientific reference | Patro et al. (2017); Love, Hogenesch & Irizarry (2016) for fragment GC-bias modelling |
| DOI | 10.1038/nmeth.4197; 10.1038/nbt.3682 |
| Provisional APA 7 reference | See §11 for Patro et al. (2017). <br> Love, M. I., Hogenesch, J. B., & Irizarry, R. A. (2016). Modeling of RNA-seq fragment sequence bias reduces systematic errors in transcript abundance estimation. *Nature Biotechnology, 34*(12), 1287–1291. https://doi.org/10.1038/nbt.3682 |
| Decision supported | D-004 (SR in 6/6), D-010 (explicit SR, approved), D-011 (SE fragment length), D-012 (`--gcBias` on, `--seqBias` off, approved) |
| Access info | 2026-10-25 (Crossref check of 10.1038/nbt.3682: 2026-10-25) |

## 13. Transcript-to-gene mapping (tx2gene)

| Field | Content |
|---|---|
| Process | Mapping of ITAG4.0 transcripts to genes for gene-level import |
| Software + version | Python 3.15.1 (stv-rnaseq env); `scripts/bash/build_tx2gene.sh` |
| Parameters | GFF mRNA `Name` → parent gene `Name`; exact string match with `quant.sf`; no version stripping |
| Methodological rationale | tximport requires a two-column TXNAME→GENEID table whose transcript IDs equal those in the quantification files |
| Official technical documentation | tximport vignette (Bioconductor 3.22): <https://bioconductor.org/packages/3.22/bioc/vignettes/tximport/inst/doc/tximport.html> |
| Primary scientific reference | Soneson et al. (2015) |
| DOI | 10.12688/f1000research.7563.1 |
| Provisional APA 7 reference | See §14 |
| Decision supported | Identifier policy (`docs/reference_validation.md` §3) |
| Access info | 2026-10-25 |

## 14. Import — tximport

| Field | Content |
|---|---|
| Process | Gene-level summarization of Salmon estimates with average-transcript-length offset |
| Software + version | tximport 1.38.2 (Bioconductor 3.22, bioconda); R 4.5.3 (env `stv-r`) |
| Parameters | `tximport(files, type = "salmon", tx2gene)` on the **final** Salmon run (`results/quantification/salmon`, SR + `--gcBias`); `countsFromAbundance = "no"`; `ignoreTxVersion = FALSE`. The exploratory `-l A` import check is archived (`tximport_validation_libA_exploratory.tsv`). |
| Methodological rationale | Recommended Salmon→DESeq2 path; corrects for changes in average transcript length across samples (here including the sample-specific GC-bias effective lengths) |
| Official technical documentation | <https://bioconductor.org/packages/3.22/bioc/vignettes/tximport/inst/doc/tximport.html>; package page <https://bioconductor.org/packages/3.22/bioc/html/tximport.html> |
| Primary scientific reference | Soneson et al. (2015) (citation given by the package; version 2 of the article, 10.12688/f1000research.7563.2, also exists) |
| DOI | 10.12688/f1000research.7563.1 |
| Provisional APA 7 reference | Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: Transcript-level estimates improve gene-level inferences. *F1000Research, 4*, 1521. https://doi.org/10.12688/f1000research.7563.1 |
| Decision supported | D-002; import validation (`docs/tximport_preparation.md`) |
| Access info | 2026-10-25 |

## 15. Differential expression — DESeq2 (Phase 4, primary analysis)

| Field | Content |
|---|---|
| Process | Gene-level differential expression STV_infected vs STV_free (`DESeqDataSetFromTximport`, `scripts/R/02_deseq2.R`) |
| Software + version | DESeq2 1.50.2 (Bioconductor 3.22, bioconda) |
| Parameters | Owner-approved D-013 – D-019: design `~ condition`, reference STV_free; pre-filter count ≥ 10 in ≥ 3 samples; Wald test, parametric dispersion fit, `minReplicatesForReplace = 7`; padj < 0.05 (BH), stringent rule adds \|MLE log2FC\| ≥ 1; QC on `vst(blind = TRUE)`, PCA top 500, Euclidean distance, complete linkage |
| Methodological rationale | Negative-binomial GLM on un-normalized estimated counts with tximport offsets; the vignette's recommended upstream pipeline is Salmon + tximport |
| Official technical documentation | <https://bioconductor.org/packages/3.22/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>; <https://bioconductor.org/packages/3.22/bioc/html/DESeq2.html> |
| Primary scientific reference | Love et al. (2014) |
| DOI | 10.1186/s13059-014-0550-8 |
| Provisional APA 7 reference | Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology, 15*(12), 550. https://doi.org/10.1186/s13059-014-0550-8 |
| Decision supported | D-002; D-005 (DESeq2 requires un-normalized counts → GEO table excluded); D-012 (vignette recommends Salmon `--gcBias` unless data are known not to need it → adopted); D-013 – D-019 |
| Access info | 2026-10-25 |

## 16. Statistical environment — R

| Field | Content |
|---|---|
| Process | Runtime for tximport/DESeq2 scripts |
| Software + version | R 4.5.3 (2026-03-11), conda-forge |
| Parameters | `envs/stv-r.yml` |
| Methodological rationale | Version required by Bioconductor 3.22 |
| Official technical documentation | <https://www.r-project.org/> |
| Primary scientific reference | Software citation |
| DOI | none |
| Provisional APA 7 reference | R Core Team. (2026). *R: A language and environment for statistical computing* (Version 4.5.3) [Computer software]. R Foundation for Statistical Computing. https://www.R-project.org/ |
| Decision supported | R environment isolation |
| Access info | 2026-10-25 |

## 17. Effect-size shrinkage — apeglm (Phase 4)

| Field | Content |
|---|---|
| Process | Shrunken log2 fold changes for ranking and visualization (`lfcShrink(type = "apeglm")`) |
| Software + version | apeglm 1.32.0 (Bioconductor 3.22) |
| Parameters | `coef = "condition_STV_infected_vs_STV_free"`; p-values remain the Wald MLE p-values |
| Methodological rationale | The DESeq2 vignette recommends apeglm for shrinkage; it reduces the noise of low-information MLE log2FC without changing the tests |
| Official technical documentation | DESeq2 vignette 3.22, "Log fold change shrinkage for visualization and ranking" |
| Primary scientific reference | Zhu, Ibrahim & Love (2019) |
| DOI | 10.1093/bioinformatics/bty895 (Crossref-verified 2026-10-25) |
| Provisional APA 7 reference | Zhu, A., Ibrahim, J. G., & Love, M. I. (2019). Heavy-tailed prior distributions for sequence count data: Removing the noise and preserving large differences. *Bioinformatics, 35*(12), 2084–2092. https://doi.org/10.1093/bioinformatics/bty895 |
| Decision supported | D-017 |
| Access info | 2026-10-25 |

## 18. Sample-level QC: VST, sample distances, PCA (Phase 4 primary QC; Phase 5 gene-set sensitivity)

| Field | Content |
|---|---|
| Process | Sample distances, PCA, clustering (Phase 4); correlation, per-sample median distances and PCA with top 100/500/1000/all genes (Phase 5, `scripts/R/03_sample_structure_diagnostics.R`) |
| Software + version | DESeq2 1.50.2 (`vst`, `plotPCA`), stats `prcomp`/`dist`/`cor`, pheatmap 1.0.13, ggplot2 4.0.3, ggrepel 0.9.8 (plotting packages approved as visualization-only, D-020) |
| Parameters | `vst(blind = TRUE)` from Phase 4 (not recomputed); `prcomp(center = TRUE, scale. = FALSE)` on the top-n genes ranked by `rowVars`; Pearson and Spearman over all 18,277 retained genes |
| Methodological rationale | The vignette uses blind VST + `plotPCA` (`ntop = 500`) + Euclidean sample distances for data-quality assessment. Varying `ntop` is an owner-requested **sensitivity diagnostic**: no publication prescribes those set sizes, and no PCA is selected from them. Spearman is an owner-requested rank-based complement with no external threshold |
| Official technical documentation | DESeq2 vignette 3.22, "Data quality assessment by sample clustering and visualization"; rnaseqGene workflow <https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html> |
| Primary scientific reference | Love et al. (2014); Love, Anders, Kim & Huber (2015) |
| DOI | 10.1186/s13059-014-0550-8; 10.12688/f1000research.7035.1 (Crossref-verified 2026-10-25; version 2 exists) |
| Provisional APA 7 reference | Love, M. I., Anders, S., Kim, V., & Huber, W. (2015). RNA-Seq workflow: Gene-level exploratory analysis and differential expression. *F1000Research, 4*, 1070. https://doi.org/10.12688/f1000research.7035.1 |
| Decision supported | D-018 (primary QC); Phase 5 Task 3: PCA structure is independent of gene set (`docs/sample_structure_investigation.md` §4) |
| Access info | 2026-10-25 |

## 19. Replicate correlation and replicate-quality criteria (Phase 5)

| Field | Content |
|---|---|
| Process | Median Pearson/Spearman correlation of each replicate with the other replicates of its condition; framing for option C (leave-one-out) |
| Software + version | R stats `cor` |
| Parameters | Statistic reported per sample; **no cut-off applied** |
| Methodological rationale | **Gierliński et al. (2015)**: median Pearson correlation is one of three bad-replicate criteria, derived from 48 yeast replicates per condition, so their cut-offs do not transfer to n = 3. **Upton et al. (2023)** (plant RNA-seq design review): replicates with poor correlation can be identified, and higher replication supports removing them; ≥ 3 replicates suffice for hypothesis generation but 6 – 12 are advised; randomize library prep and sequencing runs. **Schurch et al. (2016)**: 3 replicates detect 20 – 40 % of significant DE genes; ≥ 6 recommended; "significant minority chance that one or more replicates within each condition should be rejected". **Khorramdelazad et al. (2018)** (lentil, 3 replicates per group): VST + sample distances + PCA; sample 4C clustered with another time point and "was identified as an outlier and was excluded" — a plant-RNA-seq precedent of the practice, not a rule |
| Official technical documentation | — |
| Primary scientific reference | Gierliński et al. (2015); Upton et al. (2023); Schurch et al. (2016); Khorramdelazad et al. (2018) |
| DOI | 10.1093/bioinformatics/btv425; 10.3389/fpls.2023.1135455; 10.1261/rna.053959.115; 10.1186/s12864-018-4488-1 |
| Provisional APA 7 reference | Gierliński, M., Cole, C., Schofield, P., Schurch, N. J., Sherstnev, A., Singh, V., Wrobel, N., Gharbi, K., Simpson, G., Owen-Hughes, T., Blaxter, M., & Barton, G. J. (2015). Statistical models for RNA-seq data derived from a two-condition 48-replicate experiment. *Bioinformatics, 31*(22), 3625–3630. https://doi.org/10.1093/bioinformatics/btv425 <br> Upton, R. N., Correr, F. H., Lile, J., Reynolds, G. L., Falaschi, K., Cook, J. P., & Lachowiec, J. (2023). Design, execution, and interpretation of plant RNA-seq analyses. *Frontiers in Plant Science, 14*, 1135455. https://doi.org/10.3389/fpls.2023.1135455 <br> Schurch, N. J., Schofield, P., Gierliński, M., Cole, C., Sherstnev, A., Singh, V., Wrobel, N., Gharbi, K., Simpson, G. G., Owen-Hughes, T., Blaxter, M., & Barton, G. J. (2016). How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA, 22*(6), 839–851. https://doi.org/10.1261/rna.053959.115 <br> Khorramdelazad, M., Bar, I., Whatmore, P., Smetham, G., Bhaaskaria, V., Yang, Y., Bai, S. H., Mantri, N., Zhou, Y., & Ford, R. (2018). Transcriptome profiling of lentil (*Lens culinaris*) through the first 24 hours of *Ascochyta lentis* infection reveals key defence response genes. *BMC Genomics, 19*, 108. https://doi.org/10.1186/s12864-018-4488-1 |
| Decision supported | Phase 5 Task 3 (correlation reported without thresholds); Task 9 option C assessment (n_free = 2, removal tied to higher replication); D-022 (Option A); n = 3 limitation in `docs/results_by_process.md` |
| Access info | Full texts via Europe PMC/PMC (PMC4754627, PMC4878611, PMC5793396) and Frontiers, 2026-10-25 |

## 20. Count-outlier diagnostics — Cook's distance (Phase 5)

| Field | Content |
|---|---|
| Process | Per-sample Cook's distributions, log10 boxplot, the sample with the largest Cook's per gene, counts above the cutoff, table of the 39 Cook's-filtered genes |
| Software + version | DESeq2 1.50.2 (`assays(dds)[["cooks"]]` from the stored Phase 4 fit) |
| Parameters | Cutoff `qf(0.99, p = 2, m − p = 4)` = 18.00 (DESeq2 default); `boxplot(log10(cooks), range = 0, las = 2)`; filtering and `minReplicatesForReplace = 7` unchanged |
| Methodological rationale | Vignette "Approach to count outliers": genes are flagged when a sample's Cook's distance exceeds the cutoff (needs ≥ 3 replicates); replacement only with ≥ 7 replicates; with many outliers, check whether one sample is consistently higher, using a PCA and this boxplot |
| Official technical documentation | DESeq2 vignette 3.22, "Approach to count outliers" and "Data quality" (quality as fitness for purpose) |
| Primary scientific reference | Love et al. (2014) |
| DOI | 10.1186/s13059-014-0550-8 |
| Provisional APA 7 reference | See §15 |
| Decision supported | Phase 5 Task 4; Task 9 / D-022 (Cook's burden weighed; filtering never disabled) |
| Access info | Vignette 3.22 (DESeq2 1.50.2, built 2025-11-13), read 2026-10-25 |

## 21. Robust PCA outlier classification — rrcov (Phase 5; subtask stopped)

| Field | Content |
|---|---|
| Process | Objective outlier classification of samples by robust PCA (score and orthogonal distances) |
| Software + version | rrcov 1.7-7, robustbase 0.99-7, pcaPP 2.0-5 (conda-forge, added to `stv-r` only; lock updated) |
| Parameters | `PcaGrid(X)`, `PcaHubert(X)` with defaults (`crit.pca.distances = 0.975`, `PcaHubert` alpha 0.75, automatic k); X = samples × 18,277 genes; inputs `rlog(dds, blind = TRUE)` (primary, as in Chen et al.) and the stored VST; seed 20261025 |
| Methodological rationale | Chen et al. (2020): `PcaGrid` on rlog with default parameters reached 100 % sensitivity/specificity on simulated and composite data; `PcaHubert` and `vst` gave false positives; smallest tested set 3 + 3 + 1 (n = 7). **At n = 6 the fits are degenerate** (`PcaGrid` k = n; `PcaHubert` h = n), so the subtask was stopped and its flag is not used |
| Official technical documentation | rrcov reference manual (CRAN) |
| Primary scientific reference | Chen et al. (2020); Todorov & Filzmoser (2009); Croux, Filzmoser & Oliveira (2007); Hubert, Rousseeuw & Vanden Branden (2005) |
| DOI | 10.1186/s12859-020-03608-0; 10.18637/jss.v032.i03; 10.1016/j.chemolab.2007.01.004; 10.1198/004017004000000563 (all Crossref-verified 2026-10-25) |
| Provisional APA 7 reference | Chen, X., Zhang, B., Wang, T., Bonni, A., & Zhao, G. (2020). Robust principal component analysis for accurate outlier sample detection in RNA-Seq data. *BMC Bioinformatics, 21*, 269. https://doi.org/10.1186/s12859-020-03608-0 <br> Todorov, V., & Filzmoser, P. (2009). An object-oriented framework for robust multivariate analysis. *Journal of Statistical Software, 32*(3), 1–47. https://doi.org/10.18637/jss.v032.i03 <br> Croux, C., Filzmoser, P., & Oliveira, M. R. (2007). Algorithms for projection–pursuit robust principal component analysis. *Chemometrics and Intelligent Laboratory Systems, 87*(2), 218–225. https://doi.org/10.1016/j.chemolab.2007.01.004 <br> Hubert, M., Rousseeuw, P. J., & Vanden Branden, K. (2005). ROBPCA: A new approach to robust principal component analysis. *Technometrics, 47*(1), 64–79. https://doi.org/10.1198/004017004000000563 |
| Decision supported | Phase 5 Task 5 (feasibility → stop); Chen et al.'s warning that "removing biological outliers will result in underestimation of the natural biological variance" is weighed in Task 9 and D-022 |
| Access info | Chen et al. full text PMC7324992, 2026-10-25 |

## 22. STV label support — competitive Salmon mapping (Phase 5)

| Field | Content |
|---|---|
| Process | Count STV-derived reads per raw library to check the STV-free / STV-infected labels (`scripts/bash/stv_signal_check.sh`) |
| Software + version | Salmon 2.7.0 (separate diagnostic index; host quantification untouched) |
| Parameters | Index: ITAG4.0 cDNA + LC429302.1 + KX949574.1, `-k 31`, no decoys. Quant: `-l U -r <fastq> -p 16 --writeMappings` (STV records kept) |
| Methodological rationale | Competitive mapping assigns reads to STV only if they fit STV better than any host transcript. `-l U` counts both strands of the dsRNA genome. Interpretation is limited to presence/absence; abundance is not interpreted (Upton et al. 2023: mRNA-focused libraries under-detect non-polyadenylated viral RNA) |
| Official technical documentation | Salmon docs (§11 – §12); GenBank records <https://www.ncbi.nlm.nih.gov/nuccore/LC429302.1> and <https://www.ncbi.nlm.nih.gov/nuccore/KX949574.1> (flat files archived in `data/metadata/sources/`) |
| Primary scientific reference | Patro et al. (2017); Alcalá-Briseño et al. (2017) for KX949574.1; LC429302.1 is a direct GenBank submission by Fukuhara, Tabara & Takahashi (2018) |
| DOI | 10.1038/nmeth.4197; 10.1128/genomea.01374-16 (Crossref-verified 2026-10-25) |
| Provisional APA 7 reference | Alcalá-Briseño, R. I., Coşkan, S., Londoño, M. A., & Polston, J. E. (2017). Genome sequence of *Southern tomato virus* in asymptomatic tomato 'Sweet Hearts'. *Genome Announcements, 5*(7), e01374-16. https://doi.org/10.1128/genomea.01374-16 <br> Fukuhara, T., Tabara, M., & Takahashi, H. (2018). *Southern tomato virus M82 genomic RNA, complete genome* (GenBank LC429302.1) [Data set]. NCBI Nucleotide. https://www.ncbi.nlm.nih.gov/nuccore/LC429302.1 |
| Decision supported | Phase 5 Task 2: labels supported; no STOP; D-022 |
| Access info | NCBI E-utilities efetch, 2026-10-26T04:39:41Z (UTC) |

## 23. Historical processed-data structure check (Phase 5)

| Field | Content |
|---|---|
| Process | Sample-structure description of the authors' GEO table — **historical processed-data structure check, not used for differential-expression inference** (`scripts/R/04_geo_historical_structure_check.R`) |
| Software + version | R 4.5.3 base/stats; ggplot2, ggrepel, pheatmap |
| Parameters | FPKM and Cuffnorm "count" columns; rows non-zero in ≥ 1 sample (23,838); `log2(x + 1)`; Pearson/Spearman, Euclidean distance, PCA top 500 |
| Methodological rationale | The values are normalized (D-005), so they are invalid as DESeq2 input but usable to check whether the sample structure exists independently of our pipeline |
| Official technical documentation | Cuffnorm manual <http://cole-trapnell-lab.github.io/cufflinks/cuffnorm/>; [geo_counts_audit.md](geo_counts_audit.md) |
| Primary scientific reference | Fukuhara et al. (2020) (see §1) |
| DOI | 10.1007/s00705-019-04436-1 |
| Provisional APA 7 reference | See §1 |
| Decision supported | Phase 5 Task 7: structure present in the authors' processed data; D-022 |
| Access info | GEO supplementary file retrieved 2026-10-25 (MD5 sidecar) |

## 24. External reanalysis — methods and sample-structure comparison only (Phase 5)

| Field | Content |
|---|---|
| Process | Comparison of methods and of observed sample structure ([external_reanalysis_comparison.md](external_reanalysis_comparison.md)); genes/pathways **not** used |
| Software + version | — |
| Parameters | — |
| Methodological rationale | Independent pipeline (Rsubread/featureCounts, SL3.0, DESeq2 1.42.0) observed the same STV−_rep1 separation |
| Official technical documentation | <https://www.researchsquare.com/article/rs-8823240/v1.html> |
| Primary scientific reference | Badeli, Kaboosi & Nasrollanejad (2026), **preprint, not peer reviewed** |
| DOI | 10.21203/rs.3.rs-8823240/v1 (Crossref-verified 2026-10-25; posted 2026-02-17) |
| Provisional APA 7 reference | Badeli, G., Kaboosi, K., & Nasrollanejad, S. (2026). *Reanalysis of the tomato transcriptome revealed host responses to Southern tomato virus infection* [Preprint]. Research Square. https://doi.org/10.21203/rs.3.rs-8823240/v1 |
| Decision supported | Phase 5 Task 8; supports that the structure is not pipeline-specific; D-022 |
| Access info | 2026-10-25 |

---

# Phase 6 — functional analysis (pathway level)

Sections 25 – 33 cover Phase 6. Pre-specification: [`functional_analysis_plan.md`](functional_analysis_plan.md). Decisions D-023 – D-032. All DOIs were checked against Crossref on 2026-10-26.

## 25. Ranking metric for preranked GSEA — DESeq2 Wald statistic

| Field | Content |
|---|---|
| Process | Freeze the gene ranking for GSEA (`scripts/R/05_functional_ranking_annotation.R`, Task 1) |
| Software + version | DESeq2 1.50.2 (stored primary fit; direction check only); R 4.5.3, data.table 1.18.6.1 |
| Parameters | `stat` (= MLE log2FC / lfcSE) for STV_infected vs STV_free; universe = finite stat with a DESeq2 p-value (18,238); 39 Cook's-filtered genes excluded; no ties |
| Methodological rationale | The Wald statistic combines effect sign and precision. Ranking by padj loses the sign; ranking only significant genes is impossible (0 DEGs); apeglm LFC is reserved for visualization (D-017) |
| Official technical documentation | DESeq2 1.50.2 vignette ("p-values and adjusted p-values", "Approach to count outliers"); fgsea vignette (example ranks are t-like statistics) |
| Primary scientific reference | Love et al. (2014). Workflow precedent: Kelsey, Kalekar & Sedivy (2025), TE-Seq ("ranked by their Wald statistic (produced by DESeq2)"). Plant precedents: Tomaszewski et al. (2026) (*Physalis*, Solanaceae; fgsea; "pre-ranked by the DESeq2 Wald statistic (sign-consistent with the contrast)") and Li et al. (2026) (PSTVd-infected potato; Wald-statistic ranking; 10 – 500 genes; BH q < 0.05) |
| DOI | 10.1186/s13059-014-0550-8; 10.1186/s13100-025-00381-w; 10.3389/fpls.2026.1733427; 10.3390/plants15182842 |
| Provisional APA 7 reference | Kelsey, M. M. G., Kalekar, R. L., & Sedivy, J. M. (2025). TE-Seq: A transposable element annotation and RNA-Seq pipeline. *Mobile DNA, 16*(1), 44. https://doi.org/10.1186/s13100-025-00381-w <br> Tomaszewski, E., Reem, N., Boham, E. T., Swiler, V., & Van Eck, J. (2026). At the breaking point: Developmental and molecular insights into *Physalis grisea* fruit abscission. *Frontiers in Plant Science, 17*, 1733427. https://doi.org/10.3389/fpls.2026.1733427 <br> Li, G., Hussain, W., Pan, T., Chen, Y., Li, J., Zeng, Y., Zhou, Y., & Lv, D. (2026). Integrated transcriptomic and phytohormone metabolomic analyses reveal distinct molecular states associated with PSTVd-induced potato tuber cracking and knob-like protuberances at tuber eyes. *Plants, 15*(18), 2842. https://doi.org/10.3390/plants15182842 |
| Decision supported | D-023 |
| Access info | Full texts PMC12659357, PMC12950737, PMC13610134 (Europe PMC), 2026-10-26 |

## 26. GO annotation source for ITAG4.0 — PLAZA 5.0 (SGN files audited)

| Field | Content |
|---|---|
| Process | Build the gene → GO mapping for the tested genes (Task 2) |
| Software + version | PLAZA 5.0 dicots data file `GO/go.sly.csv.gz` (instance `dicots_05`, generated 2021-09-20; assembly/annotation ITAG4.0; provider solgenomics.net); `scripts/bash/download_functional_annotation.sh` |
| Parameters | `go_source = primary` (UniProt-GOA, InterProScan, GO-Consortium imports); `is_propagated_from_child_term = 0`. Projections (`from_iortho`, `from_family_enrichment`, `from_tree`) excluded. Exact gene IDs; no version stripping |
| Methodological rationale | Same gene models as the quantification. Per-row evidence code and provider are recorded. Coverage is 82.88 % of ranked genes, compared with SGN `ITAG4.0_goterms.txt` (2.52 %; version-mismatched IDs) and SGN `ITAG4.1_goterms.txt` (45.57 %; other release) |
| Official technical documentation | PLAZA download area <https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_05/>; SGN ITAG releases <https://solgenomics.net/ftp/tomato_genome/annotation/> (ITAG4.0, ITAG4.1 README) |
| Primary scientific reference | Van Bel et al. (2022) (PLAZA 5.0); Hosmani et al. (2019) (SL4.0 / ITAG4.0, **preprint**) |
| DOI | 10.1093/nar/gkab1024; 10.1101/767764 |
| Provisional APA 7 reference | Van Bel, M., Silvestri, F., Weitz, E. M., Kreft, L., Botzki, A., Coppens, F., & Vandepoele, K. (2022). PLAZA 5.0: Extending the scope and power of comparative and functional genomics in plants. *Nucleic Acids Research, 50*(D1), D1468–D1474. https://doi.org/10.1093/nar/gkab1024 <br> Hosmani, P. S., Flores-Gonzalez, M., van de Geest, H., Maumus, F., Bakker, L. V., Schijlen, E., van Haarst, J., Cordewener, J., Sanchez-Perez, G., Peters, S., Fei, Z., Giovannoni, J. J., Mueller, L. A., & Saha, S. (2019). *An improved de novo assembly and annotation of the tomato reference genome using single-molecule sequencing, Hi-C proximity ligation and optical maps* [Preprint]. bioRxiv. https://doi.org/10.1101/767764 |
| Decision supported | D-024 |
| Access info | Files retrieved 2026-10-26 (sha256 in `results/functional/annotation_sources_manifest.tsv`) |

## 27. Gene Ontology version handling and propagation

| Field | Content |
|---|---|
| Process | Map PLAZA GO IDs to a dated ontology, propagate, split into BP/MF/CC and build gene sets (Tasks 2 – 3) |
| Software + version | GO `go-basic.obo`, release 2026-08-05 (`data-version releases/2026-07-26`); ontologyIndex 2.12 (`get_ontology(..., propagate_relationships = c("is_a","part_of"))`) |
| Parameters | `alt_id` → primary; obsolete → single `replaced_by`, otherwise dropped; is_a + part_of only; roots excluded; ranked universe; 15 – 500 genes |
| Methodological rationale | GO docs: go-basic is acyclic and "annotations can safely be propagated up the graph"; is_a and part_of are safe for grouping; grouping over `regulates` "changes the relationship between the GO term and the gene product" |
| Official technical documentation | <https://geneontology.org/docs/download-ontology/>; <https://geneontology.org/docs/ontology-relations/>; release archive <https://release.geneontology.org/2026-08-05/> |
| Primary scientific reference | Ashburner et al. (2000); Gene Ontology Consortium (2023); Greene, Richardson & Turro (2017) (ontologyX/ontologyIndex) |
| DOI | 10.1038/75556; 10.1093/genetics/iyad031; 10.1093/bioinformatics/btw763 |
| Provisional APA 7 reference | Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler, H., Cherry, J. M., Davis, A. P., Dolinski, K., Dwight, S. S., Eppig, J. T., Harris, M. A., Hill, D. P., Issel-Tarver, L., Kasarskis, A., Lewis, S., Matese, J. C., Richardson, J. E., Ringwald, M., Rubin, G. M., & Sherlock, G. (2000). Gene Ontology: Tool for the unification of biology. *Nature Genetics, 25*(1), 25–29. https://doi.org/10.1038/75556 <br> Gene Ontology Consortium. (2023). The Gene Ontology knowledgebase in 2023. *Genetics, 224*(1), iyad031. https://doi.org/10.1093/genetics/iyad031 <br> Greene, D., Richardson, S., & Turro, E. (2017). ontologyX: A suite of R packages for working with ontological data. *Bioinformatics, 33*(7), 1104–1106. https://doi.org/10.1093/bioinformatics/btw763 |
| Decision supported | D-025, D-026 |
| Access info | 2026-10-26 |

## 28. Preranked gene-set enrichment — fgsea multilevel

| Field | Content |
|---|---|
| Process | Primary pathway-level test on GO BP, MF, CC (`scripts/R/06_functional_gsea.R`, Tasks 5 – 6) |
| Software + version | fgsea 1.36.2 (Bioconductor 3.22, bioconda), BiocParallel 1.44.0, data.table 1.18.6.1; repository <https://github.com/ctlab/fgsea/> (package DESCRIPTION) and <https://bioconductor.org/packages/fgsea/> |
| Parameters | `fgseaMultilevel(minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000, sampleSize = 101, scoreType = "std", gseaParam = 1, BPPARAM = SerialParam())`; `set.seed(20261025)`; BH within collection; padj < 0.05 |
| Methodological rationale | Threshold-free competitive test on the full ranking. Multilevel Monte-Carlo gives accurate small p-values (`eps = 0`). The 15 – 500 size range and the 5 % FDR follow the GSEA User Guide for gene-set permutation with few samples |
| Official technical documentation | fgsea vignette "Using fgsea package" (Bioconductor release); GSEA User Guide <https://docs.gsea-msigdb.org/#GSEA/GSEA_User_Guide/> (gene-set size defaults, GSEAPreranked best practices, FDR) |
| Primary scientific reference | Subramanian et al. (2005); Korotkevich et al. (2016/2021, **preprint**); Benjamini & Hochberg (1995) |
| DOI | 10.1073/pnas.0506580102; 10.1101/060012; 10.1111/j.2517-6161.1995.tb02031.x |
| Provisional APA 7 reference | Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S., & Mesirov, J. P. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences, 102*(43), 15545–15550. https://doi.org/10.1073/pnas.0506580102 <br> Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M. N., & Sergushichev, A. (2016). *Fast gene set enrichment analysis* [Preprint]. bioRxiv. https://doi.org/10.1101/060012 <br> Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological), 57*(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x |
| Decision supported | D-026, D-027, D-031 |
| Access info | 2026-10-26 |

## 29. Interpretation caveat — competitive tests and inter-gene correlation

| Field | Content |
|---|---|
| Process | Interpretation of preranked GSEA (and of PAGE, §33) |
| Software + version | — |
| Parameters | — |
| Methodological rationale | Gene-permutation (competitive) nulls treat genes as independent. Co-expressed sets therefore receive optimistic p-values. This is reported as a limitation (plan §8.4); no correlation-adjusted test was run in Phase 6 |
| Official technical documentation | GSEA User Guide ("gene_set permutation … less stringent assessment of significance") |
| Primary scientific reference | Tamayo et al. (2016); Wu & Smyth (2012) |
| DOI | 10.1177/0962280212460441; 10.1093/nar/gks461 |
| Provisional APA 7 reference | Tamayo, P., Steinhardt, G., Liberzon, A., & Mesirov, J. P. (2016). The limitations of simple gene set enrichment analysis assuming gene independence. *Statistical Methods in Medical Research, 25*(1), 472–487. https://doi.org/10.1177/0962280212460441 <br> Wu, D., & Smyth, G. K. (2012). Camera: A competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Research, 40*(17), e133. https://doi.org/10.1093/nar/gks461 |
| Decision supported | D-027 (stricter FDR), D-030 (PAGE adds no protection), plan §8.4 |
| Access info | 2026-10-26 |

## 30. clusterProfiler (framework considered, not used)

| Field | Content |
|---|---|
| Process | — (alternative GSEA front end) |
| Software + version | Not installed. Bioconductor <https://bioconductor.org/packages/clusterProfiler/>; repository <https://github.com/YuLab-SMU/clusterProfiler> |
| Parameters | — |
| Methodological rationale | Widely used for GO GSEA with custom gene sets. fgsea is called directly so that every parameter is explicit and versioned |
| Official technical documentation | clusterProfiler book/vignette (YuLab) |
| Primary scientific reference | Yu et al. (2012); Wu et al. (2021); Xu et al. (2024) |
| DOI | 10.1089/omi.2011.0118; 10.1016/j.xinn.2021.100141; 10.1038/s41596-024-01020-z |
| Provisional APA 7 reference | Yu, G., Wang, L.-G., Han, Y., & He, Q.-Y. (2012). clusterProfiler: An R package for comparing biological themes among gene clusters. *OMICS, 16*(5), 284–287. https://doi.org/10.1089/omi.2011.0118 <br> Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L., Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., & Yu, G. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation, 2*(3), 100141. https://doi.org/10.1016/j.xinn.2021.100141 <br> Xu, S., Hu, E., Cai, Y., Xie, Z., Luo, X., Zhan, L., Tang, W., Wang, Q., Liu, B., Wang, R., Xie, W., Wu, T., Xie, L., & Yu, G. (2024). Using clusterProfiler to characterize multiomics data. *Nature Protocols, 19*(11), 3292–3320. https://doi.org/10.1038/s41596-024-01020-z |
| Decision supported | D-027 |
| Access info | 2026-10-26 |

## 31. Plant Reactome (tomato pathways audited; GSEA not run)

| Field | Content |
|---|---|
| Process | Tomato pathway-membership audit and exact mapping (Task 4) |
| Software + version | Plant Reactome current download `gene_ids_by_pathway_and_species.tab` (site: Version 24, Gramene r68/r69; file Last-Modified 2026-07-08; sha256 in manifest) |
| Parameters | `Solanum lycopersicum` rows; exact ITAG4.0 gene IDs only |
| Methodological rationale | Tomato pathways are orthology projections of 351 curated rice pathways (Ensembl Compara; ≥ 30 % reciprocal identity for dicots). The docs call them "not data, but hypotheses". Coverage: 3.41 % of ranked genes; 5 testable pathways → STOP for collection D |
| Official technical documentation | Plant Reactome "Computational Inference" <https://plantreactome.gramene.org/index.php?option=com_content&view=article&id=68&Itemid=270&lang=en>; Release Summary V24; download area <https://plantreactome.gramene.org/download/current/> |
| Primary scientific reference | Gupta et al. (2024); Naithani et al. (2020) |
| DOI | 10.1093/nar/gkad1052; 10.1093/nar/gkz996 |
| Provisional APA 7 reference | Gupta, P., Elser, J., Hooks, E., D'Eustachio, P., Jaiswal, P., & Naithani, S. (2024). Plant Reactome Knowledgebase: Empowering plant pathway exploration and OMICS data analysis. *Nucleic Acids Research, 52*(D1), D1538–D1547. https://doi.org/10.1093/nar/gkad1052 <br> Naithani, S., Gupta, P., Preece, J., D'Eustachio, P., Elser, J. L., Garg, P., Dikeman, D. A., Kiff, J., Cook, J., Olson, A., Wei, S., Tello-Ruiz, M. K., Mundo, A. F., Munoz-Pomer, A., Mohammed, S., Cheng, T., Bolton, E., Papatheodorou, I., Stein, L., … Jaiswal, P. (2020). Plant Reactome: A knowledgebase and resource for comparative pathway analysis. *Nucleic Acids Research* (published online 2019-11-04; volume/pages to be checked). https://doi.org/10.1093/nar/gkz996 |
| Decision supported | D-028 |
| Access info | 2026-10-26 |

## 32. MapMan4 / Mercator4 (preparation only)

| Field | Content |
|---|---|
| Process | Prepare Mercator4 input and a MapMan expression table (`scripts/R/07_mapman_mercator_prep.R`, Task 7) |
| Software + version | Mercator4 **v8.0** web service (version shown on <https://www.plabipd.de/mercator_main.html>, 2026-10-26; "scheduled maintenance" notice); MapMan desktop (version to be recorded at use) |
| Parameters | Unmodified ITAG4.0 protein FASTA (34,075 records); sequence type protein; Swiss-Prot option not selected; display value apeglm log2FC (D-017) |
| Methodological rationale | MapMan4 bins give a plant-specific functional framework. With 0 FDR-significant DEGs, visualization of continuous effect estimates is descriptive; no bin ORA |
| Official technical documentation | Mercator4 web page (job submission: protein or cDNA FASTA; FASTA validation; optional Swiss-Prot annotation); `results/functional/mapman/MERCATOR4_SUBMISSION.md` |
| Primary scientific reference | Schwacke et al. (2019); Bolger, Schwacke & Usadel (2021) |
| DOI | 10.1016/j.molp.2019.01.003; 10.1007/978-1-0716-1609-3_9 |
| Provisional APA 7 reference | Schwacke, R., Ponce-Soto, G. Y., Krause, K., Bolger, A. M., Arsova, B., Hallab, A., Gruden, K., Stitt, M., Bolger, M. E., & Usadel, B. (2019). MapMan4: A refined protein classification and annotation framework applicable to multi-omics data analysis. *Molecular Plant, 12*(6), 879–892. https://doi.org/10.1016/j.molp.2019.01.003 <br> Bolger, M., Schwacke, R., & Usadel, B. (2021). MapMan visualization of RNA-Seq data using Mercator4 functional annotations. In *Solanum tuberosum: Methods and protocols* (Methods in Molecular Biology, pp. 195–212). Humana. https://doi.org/10.1007/978-1-0716-1609-3_9 |
| Decision supported | D-029 |
| Access info | 2026-10-26 (Springer chapter full text not accessible: login redirect) |

## 33. agriGO v2 / PAGE (feasibility only; not used)

| Field | Content |
|---|---|
| Process | Feasibility assessment (Task 8; [`agrigo_feasibility.md`](agrigo_feasibility.md)) |
| Software + version | agriGO v2.0 web server (site "Last maintain: 2021-04-04"); internal GO datasets 2009-10 and 2016-8 |
| Parameters | — (no job submitted) |
| Methodological rationale | SEA needs a DEG list (none). PAGE would run on agriGO's 2016 ontology with an undocumented ITAG4.0 background or a ≤ 4 MB upload, cannot be re-executed reproducibly, and tests the same competitive question as fgsea → conclusion C |
| Official technical documentation | <https://systemsbiology.cau.edu.cn/agriGOv2/manual.php>; <https://systemsbiology.cau.edu.cn/agriGOv2/download.php> |
| Primary scientific reference | Tian et al. (2017); Kim & Volsky (2005) |
| DOI | 10.1093/nar/gkx382; 10.1186/1471-2105-6-144 |
| Provisional APA 7 reference | Tian, T., Liu, Y., Yan, H., You, Q., Yi, X., Du, Z., Xu, W., & Su, Z. (2017). agriGO v2.0: A GO analysis toolkit for the agricultural community, 2017 update. *Nucleic Acids Research, 45*(W1), W122–W129. https://doi.org/10.1093/nar/gkx382 <br> Kim, S.-Y., & Volsky, D. J. (2005). PAGE: Parametric analysis of gene set enrichment. *BMC Bioinformatics, 6*, 144. https://doi.org/10.1186/1471-2105-6-144 |
| Decision supported | D-030 |
| Access info | 2026-10-26 |

---

# Phase 6B — functional robustness validation

Sections 34 – 41 cover Phase 6B. Pre-specification: [`functional_robustness_plan.md`](functional_robustness_plan.md) (frozen
2026-10-26T07:36:28Z; outcome in §11). Decisions D-033 – D-043. All DOIs below were checked against the Crossref REST API on
2026-10-26. Full texts or abstracts were read on that date via Europe PMC or the publisher page, as stated per entry. A source is
cited only for what it actually reports or performs.

## 34. Correlation-adjusted competitive gene-set test — limma `camera`

| Field | Content |
|---|---|
| Process | Validation of the Phase 6 fgsea calls with a test that accounts for inter-gene correlation (`scripts/R/11_functional_camera.R`, `12_functional_method_comparison.R`; Tasks 3 – 6) |
| Software + version | limma 3.66.0 (Bioconductor 3.22, bioconda) |
| Parameters | `camera(v, index, design, contrast = 2, inter.gene.cor = NA, allow.neg.cor = FALSE, use.ranks = FALSE, trend.var = FALSE, sort = FALSE)` (primary); same call with `inter.gene.cor = 0.01` (secondary, descriptive); BH within collection; FDR < 0.05; indices from the Phase 6 GO sets via `ids2indices` |
| Methodological rationale | CAMERA inflates the variance of the set statistic by the variance inflation factor 1 + (m − 1)ρ̄, estimated from the data. A gene-permutation (preranked) null ignores inter-gene correlation. Goeman & Bühlmann (2007) distinguish tests of a set against the rest of the genes from tests of each set on its own. They show that models taking genes as sampling units rely on an "unrealistic independence assumption", so their p-values "can be wildly anti-conservative". `?camera`: with `inter.gene.cor = NA`, camera "gives rigorous error rate control for all sample sizes and all gene sets"; the default 0.01 is "a useful compromise between strict error rate control and interpretable gene set rankings". The strict mode was fixed as primary before any result. Plant precedent (indirect): Liu et al. (2022) used CAMERA inside the GeneBridge G-MAD module-association method on cross-species expression compendia that include tomato. That is not a two-group validation |
| Official technical documentation | limma 3.66.0 reference manual (`?camera`, `?ids2indices`); limma User's Guide (Bioconductor 3.22) |
| Primary scientific reference | Wu & Smyth (2012); Goeman & Bühlmann (2007); Liu et al. (2022) (indirect precedent) |
| DOI | 10.1093/nar/gks461; 10.1093/bioinformatics/btm051; 10.1093/hr/uhac251 |
| Provisional APA 7 reference | Wu, D., & Smyth, G. K. (2012). See §29. <br> Goeman, J. J., & Bühlmann, P. (2007). Analyzing gene expression data in terms of gene sets: Methodological issues. *Bioinformatics, 23*(8), 980–987. https://doi.org/10.1093/bioinformatics/btm051 <br> Liu, M., Guo, C., Xie, K., Chen, K., Chen, J., Wang, Y., & Wang, X. (2023). A cross-species co-functional gene network underlying leaf senescence. *Horticulture Research, 10*(1), uhac251. https://doi.org/10.1093/hr/uhac251 (published online 2022-11-10) |
| Decision supported | D-034, D-035 |
| Access info | Liu et al.: full text PMC9832971 (Europe PMC), 2026-10-26 (Methods, "GeneBridge analysis"). Goeman & Bühlmann: abstract via Europe PMC (PMID 17303618). Wu & Smyth: as §29 |

## 35. Expression model for CAMERA — voom, limma, edgeR TMM

| Field | Content |
|---|---|
| Process | Independent expression model on the same counts (Task 4) |
| Software + version | limma 3.66.0 (`voom`, `lmFit`, `eBayes`); edgeR 4.8.2 (`DGEList`, `calcNormFactors`, `filterByExpr` diagnostic only); statmod 1.5.2 |
| Parameters | 18,238-gene universe (no re-filtering); TMM; `voom(y, design)` with `design = model.matrix(~ condition)`, reference STV_free; coefficient 2; Spearman ρ(limma moderated t, DESeq2 Wald) reported, ≤ 0 = STOP |
| Methodological rationale | voom estimates the mean–variance trend of log-CPM and passes observation weights to limma (Law et al. 2014). TMM corrects composition (Robinson & Oshlack 2010). The RNAseq123 workflow shows `camera()` on a voom object with a design and contrast (Law et al. 2016). edgeR v4 is the current release line (Chen et al. 2025) |
| Official technical documentation | limma 3.66.0 `?voom`; edgeR 4.8.2 `?calcNormFactors`, `?filterByExpr`; RNAseq123 Bioconductor workflow |
| Primary scientific reference | Law et al. (2014); Ritchie et al. (2015); Robinson & Oshlack (2010); Law et al. (2016); Chen et al. (2025) |
| DOI | 10.1186/gb-2014-15-2-r29; 10.1093/nar/gkv007; 10.1186/gb-2010-11-3-r25; 10.12688/f1000research.9005.3; 10.1093/nar/gkaf018 |
| Provisional APA 7 reference | Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology, 15*(2), R29. https://doi.org/10.1186/gb-2014-15-2-r29 <br> Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research, 43*(7), e47. https://doi.org/10.1093/nar/gkv007 <br> Robinson, M. D., & Oshlack, A. (2010). A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology, 11*(3), R25. https://doi.org/10.1186/gb-2010-11-3-r25 <br> Law, C. W., Alhamdoosh, M., Su, S., Dong, X., Tian, L., Smyth, G. K., & Ritchie, M. E. (2018). RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR (Version 3). *F1000Research, 5*, 1408. https://doi.org/10.12688/f1000research.9005.3 <br> Chen, Y., Chen, L., Lun, A. T. L., Baldoni, P. L., & Smyth, G. K. (2025). edgeR v4: Powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets. *Nucleic Acids Research, 53*(2), gkaf018. https://doi.org/10.1093/nar/gkaf018 |
| Decision supported | D-034 |
| Access info | 2026-10-26 (package help pages from the installed versions; Crossref metadata) |

## 36. Counts for limma-voom from Salmon — tximport `lengthScaledTPM`

| Field | Content |
|---|---|
| Process | Re-import of the same six `quant.sf` files for the voom model (`scripts/R/11_functional_camera.R`); the Phase 3 import (`countsFromAbundance = "no"` + length offsets for DESeq2) is unchanged |
| Software + version | tximport 1.38.2 |
| Parameters | `tximport(files, type = "salmon", tx2gene = reference/ITAG4.0_tx2gene.tsv, countsFromAbundance = "lengthScaledTPM")`; `quant.sf` SHA-256 checked against `results/quantification/quant_sf_sha256.tsv` before import |
| Methodological rationale | The tximport vignette (limma-voom section) recommends scaled counts from abundances, because limma-voom does not use the offset matrix. Soneson et al. (2015) show that length-scaled TPM or offsets correct for differential transcript length |
| Official technical documentation | tximport 1.38.2 vignette, "limma-voom" section |
| Primary scientific reference | Soneson, Love & Robinson (2015) |
| DOI | 10.12688/f1000research.7563.2 |
| Provisional APA 7 reference | Soneson, C., Love, M. I., & Robinson, M. D. (2016). Differential analyses for RNA-seq: Transcript-level estimates improve gene-level inferences (Version 2). *F1000Research, 4*, 1521. https://doi.org/10.12688/f1000research.7563.2 |
| Decision supported | D-034 |
| Access info | 2026-10-26 |

## 37. Cook's-filtered genes in the ranking — inclusion sensitivity

| Field | Content |
|---|---|
| Process | Final D-023 rule and its sensitivity analysis (`scripts/R/10_functional_cooks_sensitivity.R`; [`cooks_gsea_decision.md`](cooks_gsea_decision.md)) |
| Software + version | DESeq2 1.50.2 (stored fit), fgsea 1.36.2 |
| Parameters | Primary: the Phase 6 ranking (18,238 genes). Sensitivity: all 18,277 genes with a finite Wald statistic; gene sets rebuilt with the Phase 6 rules and checked for identity; same fgsea call and seed |
| Methodological rationale | DESeq2 keeps the Wald statistic but sets p/padj to NA for Cook's outliers ("Approach to count outliers"). Love et al. (2014) advise excluding the whole gene in small experiments. Cook (1977) defines case-deletion influence. No method paper on Cook's outliers in preranked GSEA was found, so the evidence is ambiguous and the primary choice was kept, with a labelled sensitivity run |
| Official technical documentation | DESeq2 1.50.2 vignette; fgsea 1.36.2 source (`preparePathwaysAndStats`) |
| Primary scientific reference | Cook (1977); Love et al. (2014); Kelsey et al. (2025) (TE-Seq, implicit inclusion in a public workflow) |
| DOI | 10.1080/00401706.1977.10489493 (alias 10.2307/1268249); 10.1186/s13059-014-0550-8; 10.1186/s13100-025-00381-w |
| Provisional APA 7 reference | Cook, R. D. (1977). Detection of influential observation in linear regression. *Technometrics, 19*(1), 15–18. https://doi.org/10.1080/00401706.1977.10489493 <br> Love et al. (2014): see §15. Kelsey et al. (2025): see §25 |
| Decision supported | D-023, D-033 |
| Access info | 2026-10-26 |

## 38. Replicate consistency and leave-one-sample-out at the pathway level (diagnostics)

| Field | Content |
|---|---|
| Process | Per-sample leading-edge scores (`scripts/R/13_functional_replicate_consistency.R`); six reduced DESeq2 fits + fgsea (`scripts/R/14_functional_loo_sensitivity.R`) |
| Software + version | DESeq2 1.50.2 (blind VST from the stored fit; reduced refits), fgsea 1.36.2 |
| Parameters | As in plan §5 – §6 (z-scored blind VST; classes from P and r₁; LOO refits with genes fixed to the 18,277 tested, `minReplicatesForReplace = 7`, same fgsea call) |
| Methodological rationale | Leave-one-out is case deletion applied at the sample level (Cook 1977); DESeq2's per-gene Cook's distance measures the same kind of influence (Love et al. 2014). Schurch et al. (2016) assess the stability of DE results by subsampling replicates (gene level). **No pathway-level leave-one-out method paper was found**, so the analysis is a labelled sensitivity diagnostic. It is never used to remove a sample (D-022) |
| Official technical documentation | DESeq2 1.50.2 vignette ("Approach to count outliers": outlier flagging needs ≥ 3 replicates) |
| Primary scientific reference | Cook (1977); Love et al. (2014); Schurch et al. (2016) |
| DOI | 10.1080/00401706.1977.10489493; 10.1186/s13059-014-0550-8; 10.1261/rna.053959.115 |
| Provisional APA 7 reference | See §37, §15, §19 |
| Decision supported | D-036, D-037 |
| Access info | 2026-10-26 |

## 39. Sequence similarity search — BLAST+

| Field | Content |
|---|---|
| Process | Organelle similarity of ITAG4.0 genes (`scripts/bash/ch00_organelle_blast.sh`; [`ch00_functional_investigation.md`](ch00_functional_investigation.md)) |
| Software + version | BLAST+ 2.17.0+ (`blastn`, `blastp`, `makeblastdb`) from the pre-existing conda env `rachis-qiime2-2026.4` |
| Parameters | blastn megablast, e-value 1e-20, `-max_target_seqs 2 -max_hsps 1`; blastp e-value 1e-10, `-max_target_seqs 1 -max_hsps 1`; tabular output with `pident`, `qcovs`. Known limitation: `-max_hsps 1` under-reports coverage for spliced genes (plan §11.3) |
| Methodological rationale | Direct test of whether ch00 gene models are (near-)identical to organelle sequences |
| Official technical documentation | BLAST+ command-line user manual (NCBI Bookshelf NBK279690) |
| Primary scientific reference | Camacho et al. (2009) |
| DOI | 10.1186/1471-2105-10-421 |
| Provisional APA 7 reference | Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. *BMC Bioinformatics, 10*, 421. https://doi.org/10.1186/1471-2105-10-421 |
| Decision supported | D-038 |
| Access info | 2026-10-26 |

## 40. Organelle references and the SL4.0 unplaced sequence (ch00)

| Field | Content |
|---|---|
| Process | Organelle reference download (`scripts/bash/download_organelle_reference.sh`) and interpretation of SL4.0ch00 |
| Software + version | NCBI E-utilities `efetch` (`db=nuccore`, `rettype=fasta` and `fasta_cds_aa`); records **NC_007898.3** (*S. lycopersicum* chloroplast, complete genome; 87 CDS proteins) and **NC_035963.1** (*S. lycopersicum* TGRC:LA1421 mitochondrion, complete genome; 39 CDS proteins); SGN `build_4.00/ReleaseNotes.txt` |
| Parameters | SHA-256, size, first header and retrieval time per file: `results/functional_robustness/ch00/organelle_reference_manifest.tsv` |
| Methodological rationale | The SL4.0 FASTA contains no organelle genome. The SGN release notes describe ch00 as unplaced contigs, and earlier builds separated or removed organellar contigs (quotes in [`ch00_functional_investigation.md`](ch00_functional_investigation.md) §1). *Correction (D-044): no **labelled** organelle record, but ch00 contains the complete plastome and mitogenome sequence; see §43 and [`organelle_requantification_assessment.md`](organelle_requantification_assessment.md) §4.3* |
| Official technical documentation | NCBI RefSeq records <https://www.ncbi.nlm.nih.gov/nuccore/NC_007898.3>, <https://www.ncbi.nlm.nih.gov/nuccore/NC_035963.1>; SGN <https://solgenomics.net/ftp/tomato_genome/assembly/build_4.00/ReleaseNotes.txt> |
| Primary scientific reference | Kahlau et al. (2006) (tomato plastome sequence). No publication is cited for NC_035963.1; the RefSeq record is cited as a data source |
| DOI | 10.1007/s00239-005-0254-5 |
| Provisional APA 7 reference | Kahlau, S., Aspinall, S., Gray, J. C., & Bock, R. (2006). Sequence of the tomato chloroplast DNA and evolutionary comparison of solanaceous plastid genomes. *Journal of Molecular Evolution, 63*(2), 194–207. https://doi.org/10.1007/s00239-005-0254-5 |
| Decision supported | D-038 |
| Access info | Records and release notes retrieved 2026-10-26 (UTC times in the manifests) |

## 41. Literature used for plausibility only (Task 10 – 11)

| Field | Content |
|---|---|
| Process | Comparison with prior analyses of GSE137303 and relevance table ([`biological_comparison_prior_studies.md`](biological_comparison_prior_studies.md); `results_by_process.md` Phase 6B evidence table) |
| Software + version | — |
| Parameters | Priority: STV studies > tomato–virus transcriptomics > persistent/endornavirus-like viruses > broad plant–virus (plan §8). Each source is used only for what it reports. No source is used to change a rule or class |
| Methodological rationale | STV: Sabanadzovic et al. (2009) describe the virus (dsRNA; seed-transmitted, not mechanically or graft-transmitted). Fukuhara et al. (2020) and Badeli et al. (2026) analyse the same reads (§1, §24). Persistent viruses: Roossinck (2010), review of persistent vs acute lifestyles; Takahashi et al. (2019), review of latent infections; Khankhum et al. (2016), data paper, 3 vs 3 bean leaf RNA-seq with edgeR, 132 DEGs, Blast2GO distribution (not an enrichment test). Broad plant–virus: Whitham et al. (2003), Arabidopsis microarrays after five acute RNA viruses. Tomato–virus transcriptomics: no source cited; two targeted Europe PMC queries found no open-access study with a tomato host, a virus-only design and a leaf transcriptome |
| Official technical documentation | — |
| Primary scientific reference | Sabanadzovic et al. (2009); Roossinck (2010); Takahashi et al. (2019); Khankhum et al. (2016); Whitham et al. (2003) |
| DOI | 10.1016/j.virusres.2008.11.018; 10.1098/rstb.2010.0057; 10.3389/fmicb.2019.02764; 10.3389/fmicb.2016.01905; 10.1046/j.1365-313x.2003.01625.x |
| Provisional APA 7 reference | Sabanadzovic, S., Valverde, R. A., Brown, J. K., Martin, R. R., & Tzanetakis, I. E. (2009). Southern tomato virus: The link between the families *Totiviridae* and *Partitiviridae*. *Virus Research, 140*(1–2), 130–137. https://doi.org/10.1016/j.virusres.2008.11.018 <br> Roossinck, M. J. (2010). Lifestyles of plant viruses. *Philosophical Transactions of the Royal Society B: Biological Sciences, 365*(1548), 1899–1905. https://doi.org/10.1098/rstb.2010.0057 <br> Takahashi, H., Fukuhara, T., Kitazawa, H., & Kormelink, R. (2019). Virus latency and the impact on plants. *Frontiers in Microbiology, 10*, 2764. https://doi.org/10.3389/fmicb.2019.02764 <br> Khankhum, S., Sela, N., Osorno, J. M., & Valverde, R. A. (2016). RNAseq analysis of endornavirus-infected vs. endornavirus-free common bean (*Phaseolus vulgaris*) cultivar Black Turtle Soup. *Frontiers in Microbiology, 7*, 1905. https://doi.org/10.3389/fmicb.2016.01905 <br> Whitham, S. A., Quan, S., Chang, H.-S., Cooper, B., Estes, B., Zhu, T., Wang, X., & Hou, Y.-M. (2003). Diverse RNA viruses elicit the expression of common sets of genes in susceptible *Arabidopsis thaliana* plants. *The Plant Journal, 33*(2), 271–283. https://doi.org/10.1046/j.1365-313x.2003.01625.x |
| Decision supported | Task 10 – 11 documentation; no decision depends on these sources |
| Access info | 2026-10-26: abstracts via Europe PMC (Sabanadzovic, Roossinck, Takahashi, Whitham); full text PMC5126043 (Khankhum); Fukuhara et al. abstract via PubMed 31620899 (full text not accessible); Badeli et al. HTML full text |

## 42. Organelle sequence and decoy-aware Salmon indexing (assessment only; D-044)

| Field | Content |
|---|---|
| Process | Assessment of whether an organelle-aware re-quantification is needed ([`organelle_requantification_assessment.md`](organelle_requantification_assessment.md) §6 – §8). No index was built |
| Software + version | Salmon 2.7.0 (primary index, unchanged); BLAST+ 2.17.0+ for the multi-HSP diagnostic (`scripts/bash/organelle_multihsp_blast.sh`) |
| Parameters | Union-of-HSP coverage at ≥ 98 % identity; near-identical ≥ 80 %, partial 50 – 80 % (frozen §1.4) |
| Methodological rationale | Official guidance: the whole genome is used as decoy ("as is recommended"). Under the selective-alignment decoy rule, ties between transcript and decoy go to the transcript, and fragments that score strictly better on a decoy are discarded. No source found requires separate organelle records for nuclear plant RNA-seq |
| Official technical documentation | Salmon "Preparing transcriptome indices (mapping-based mode)" <https://salmon.readthedocs.io/en/latest/salmon.html>; nf-core/rnaseq 3.27.0 usage, "Alignment options" and "Reference genome options" <https://nf-co.re/rnaseq/latest/docs/usage/> |
| Primary scientific reference | Srivastava et al. (2020); Patro et al. (2017) (see §11) |
| DOI | 10.1186/s13059-020-02151-8; 10.1038/nmeth.4197 |
| Provisional APA 7 reference | Srivastava, A., Malik, L., Sarkar, H., Zakeri, M., Almodaresi, F., Soneson, C., Love, M. I., Kingsford, C., & Patro, R. (2020). Alignment and mapping methodology influence transcript abundance estimation. *Genome Biology, 21*, 239. https://doi.org/10.1186/s13059-020-02151-8 <br> Patro et al. (2017): see §11 |
| Decision supported | D-044 (recommendation A, owner-pending) |
| Access info | 2026-10-26: Salmon docs and nf-core usage pages (HTML); Srivastava et al. full text via Europe PMC (PMC7487471) and Crossref metadata; Patro et al. citation verified through the Crossref reference list of Srivastava et al. (the direct Crossref request failed) |

## 43. SL4.0 assembly provenance: SGN, NCBI GenBank, Ensembl Plants (assessment only; D-044)

| Field | Content |
|---|---|
| Process | Provenance of SL4.0ch00 and its organelle-identical contigs (`scripts/bash/ch00_ncbi_provenance.sh`; [`organelle_requantification_assessment.md`](organelle_requantification_assessment.md) §3) |
| Software + version | NCBI Datasets v2 REST (`genome/accession/GCA_000188115.5/sequence_reports`), NCBI E-utilities `efetch` (`db=nuccore`, FASTA of the 140 unplaced scaffolds, temp dir only), Python 3 (system) |
| Parameters | Exact sequence identity (forward or reverse complement) of SGN ch00 segments split at ≥ 100 N vs GenBank scaffold pieces split the same way |
| Methodological rationale | Establishes whether the organelle-identical ch00 contigs belong to the official SL4.0 assembly or only to the SGN FASTA. Result: plastid-identical contigs are present in both; 12 mito-identical contigs (832,437 bp) are present only in SGN |
| Official technical documentation | NCBI assembly GCA_000188115.5 <https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000188115.5/>; Ensembl REST `info/assembly/solanum_lycopersicum_gca000188115v5cm` and `info/assembly/arabidopsis_thaliana` (precedent: TAIR10 Mt/Pt top-level) <https://rest.ensembl.org>; SGN release notes (§40) |
| Primary scientific reference | Hosmani et al. (2019), bioRxiv preprint (abstract only; full text not retrieved) |
| DOI | 10.1101/767764 |
| Provisional APA 7 reference | Hosmani, P. S., Flores-Gonzalez, M., van de Geest, H., Maumus, F., Bakker, L. V., Schijlen, E., van Haarst, J., Cordewener, J., Sanchez-Perez, G., Peters, S., Fei, Z., Giovannoni, J. J., Mueller, L. A., & Saha, S. (2019). *An improved de novo assembly and annotation of the tomato reference genome using single-molecule sequencing, Hi-C proximity ligation and optical maps* [Preprint]. bioRxiv. https://doi.org/10.1101/767764 |
| Decision supported | D-044 |
| Access info | 2026-10-26 (retrieval date stored in `results/organelle_assessment/provenance/ncbi_GCA_000188115.5_sequence_report.tsv`); Crossref metadata for 10.1101/767764 |

## 44. STV biology and taxonomy (Phase 7, context only)

| Field | Content |
|---|---|
| Process | STV biology summary ([`literature_biological_context.md`](literature_biological_context.md) §2) |
| Software + version | Europe PMC REST API (`scripts/bash/literature_search_europepmc.sh`, Python 3.14.4 standard library); ICTV web report |
| Parameters | Queries Q1 – Q4 (metadata saved in `results/biological_integration/literature/`); abstracts read via the API; nothing used from title-only records |
| Methodological rationale | Establishes what is known about STV independently of this dataset (L5; plausibility only, D-046) |
| Official technical documentation | ICTV Report, family *Amalgaviridae* (interim, October 2024) <https://ictv.global/report/chapter/amalgaviridae/amalgaviridae>; Europe PMC REST API <https://europepmc.org/RestfulWebService> |
| Primary scientific reference | Sabanadzovic et al. (2009); Fukuhara et al. (2020); Elvira-González et al. (2020, *Eur J Plant Pathol*); Elvira-González et al. (2020, *Virus Genes*); Elvira-González et al. (2021); Hao et al. (2023); Hussain et al. (2022); Hussain et al. (2026); Jindal & Ali (2025); Turco et al. (2018); Alcalá-Briseño et al. (2017); Vendrell-Mir et al. (2021); Badeli et al. (2026) [preprint] |
| DOI | 10.1016/j.virusres.2008.11.018; 10.1007/s00705-019-04436-1; 10.1007/s10658-019-01911-y; 10.1007/s11262-020-01785-x; 10.3390/microorganisms9040689; 10.1007/s41348-023-00729-5; 10.3390/v14112481; 10.1094/pdis-04-25-0930-sc; 10.3390/v17091193; 10.1094/mpmi-12-17-0301-r; 10.1128/genomea.01374-16; 10.1111/tpj.15545; 10.21203/rs.3.rs-8823240/v1 |
| Provisional APA 7 reference | Elvira-González, L., Medina, V., Rubio, L., & Galipienso, L. (2020). The persistent southern tomato virus modifies miRNA expression without inducing symptoms and cell ultra-structural changes. *European Journal of Plant Pathology, 156*, 615–622. <br> Elvira-González, L., Rubio, L., & Galipienso, L. (2020). Geographically distant isolates of the persistent southern tomato virus (STV) show very low genetic diversity in the putative coat protein gene. *Virus Genes, 56*, 668–672. <br> Elvira-González, L., Peiró, R., Rubio, L., & Galipienso, L. (2021). Persistent southern tomato virus (STV) interacts with cucumber mosaic and/or pepino mosaic virus in mixed-infections modifying plant symptoms, viral titer and small RNA accumulation. *Microorganisms, 9*(4), 689. <br> Hao, X., Zheng, Y., Cui, B., & Xiang, B. (2023). Localization of southern tomato virus (STV) in tomato tissues. *Journal of Plant Diseases and Protection*. <br> Hussain, M. D., Farooq, T., Chen, X., Jiang, T., et al. (2022). Molecular detection of southern tomato amalgavirus prevalent in tomatoes and its genomic characterization with global evolutionary dynamics. *Viruses, 14*(11), 2481. <br> Hussain, M. D., Farooq, T., Datt, N., Li, J., et al. (2026). Mapping the enigma of cryptic tomato yellow stunt disease: A synergistic viral infection threatening global tomato health. *Plant Disease, 110*, 592–598. <br> Jindal, S., & Ali, A. (2025). First detection and identification of southern tomato virus infecting tomatoes in Oklahoma with complete genome characterization and insights into global genetic diversity. *Viruses, 17*(9), 1193. <br> Turco, S., Golyaev, V., Seguin, J., Gilli, C., et al. (2018). Small RNA-omics for virome reconstruction and antiviral defense characterization in mixed infections of cultivated *Solanum* plants. *Molecular Plant-Microbe Interactions, 31*(7), 707–723. <br> Alcalá-Briseño, R. I., Coşkan, S., Londoño, M. A., & Polston, J. E. (2017). Genome sequence of *Southern tomato virus* in asymptomatic tomato 'Sweet Hearts'. *Genome Announcements, 5*, e01374-16. <br> Vendrell-Mir, P., Perroud, P.-F., Haas, F. B., Meyberg, R., et al. (2021). A vertically transmitted amalgavirus is present in certain accessions of the bryophyte *Physcomitrium patens*. *The Plant Journal, 108*(6), 1786–1797. <br> (Author lists marked "et al." are abbreviated as in the Europe PMC records; complete them before publication.) Sabanadzovic et al. (2009), Fukuhara et al. (2020) and Badeli et al. (2026): see §1, §24 and §41 |
| Decision supported | D-046 |
| Access info | 2026-10-26: Europe PMC abstracts; Badeli et al. full text (Research Square HTML v1); ICTV report page; Fukuhara et al. full text not accessible |

## 45. Tomato and plant virus transcriptomics (Phase 7, plausibility only)

| Field | Content |
|---|---|
| Process | Theme-level plausibility ([`literature_biological_context.md`](literature_biological_context.md) §3; `results/biological_integration/curated_theme_literature.tsv`) |
| Software + version | Europe PMC REST API (as §44) |
| Parameters | Q3 (tomato + virus + RNA-seq/transcriptome + leaf) and Q4 (persistent viruses); targeted queries for ribosomal-protein genes in virus infection and PepMV tomato transcriptome |
| Methodological rationale | Shows what acute-virus responses look like in tomato leaves; all are acute infections and are marked "not directly comparable" |
| Official technical documentation | — |
| Primary scientific reference | Hanssen et al. (2011); Wang et al. (2024); He et al. (2024); Prigigallo et al. (2025); Romero-Rodríguez et al. (2023); Yue et al. (2024); Kamitani et al. (2016); Zorzatto et al. (2015); Whitham et al. (2003); Khankhum et al. (2016) |
| DOI | 10.1104/pp.111.173906; 10.3390/ijms25074012; 10.3390/plants13223198; 10.1038/s41598-025-88028-x; 10.1186/s12870-023-04534-y; 10.1111/pce.14906; 10.1093/femsec/fiw176; 10.1038/nature14171; 10.1046/j.1365-313x.2003.01625.x; 10.3389/fmicb.2016.01905 |
| Provisional APA 7 reference | Hanssen, I. M., van Esse, H. P., Ballester, A.-R., Hogewoning, S. W., et al. (2011). Differential tomato transcriptomic responses induced by pepino mosaic virus isolates with differential aggressiveness. *Plant Physiology, 156*(1), 301–318. <br> Wang, D., Chen, M., Peng, J., Zheng, H., et al. (2024). Transcriptome analysis of tomato leaves reveals candidate genes responsive to *Tomato brown rugose fruit virus* infection. *International Journal of Molecular Sciences, 25*(7), 4012. <br> He, W.-Z., Rong, T., Liu, X.-Y., & Rao, Q. (2024). Transcriptomic profiling unravels the disruption of photosynthesis apparatuses and induction of immune responses by a bipartite begomovirus in tomato plants. *Plants, 13*(22), 3198. <br> Prigigallo, M. I., Picciotti, U., & Bubici, G. (2025). Resistance-breaking strains of tomato spotted wilt virus hamper photosynthesis and protein synthesis pathways in a virus accumulation-dependent manner in Sw5-carrying tomatoes. *Scientific Reports, 15*, 3630. <br> Romero-Rodríguez, B., Petek, M., Jiao, C., Križnik, M., et al. (2023). Transcriptional and epigenetic changes during tomato yellow leaf curl virus infection in tomato. *BMC Plant Biology, 23*, 651. <br> Yue, H., Chen, G., Zhang, Z., Guo, Z., et al. (2024). Single-cell transcriptome landscape elucidates the cellular and developmental responses to tomato chlorosis virus infection in tomato leaf. *Plant, Cell & Environment, 47*(7), 2660–2674. <br> Kamitani, M., Nagano, A. J., Honjo, M. N., & Kudoh, H. (2016). RNA-Seq reveals virus–virus and virus–plant interactions in nature. *FEMS Microbiology Ecology, 92*(11), fiw176. <br> Zorzatto, C., Machado, J. P. B., Lopes, K. V. G., Nascimento, K. J. T., et al. (2015). NIK1-mediated translation suppression functions as a plant antiviral immunity mechanism. *Nature, 520*, 679–682. <br> Whitham et al. (2003) and Khankhum et al. (2016): see §41 |
| Decision supported | D-046 |
| Access info | 2026-10-26: abstracts via Europe PMC |

## 46. Phase 7 evidence integration and GO redundancy (descriptive)

| Field | Content |
|---|---|
| Process | Theme evidence matrix, ethylene case study, theme deep dives (`scripts/R/20_biological_integration.R`); GO redundancy modules (`scripts/R/21_go_redundancy_modules.R`); figures (`scripts/R/22_phase7_figures.R`) |
| Software + version | R 4.5.3; data.table 1.18.6.1; DESeq2 1.50.2 (blind VST read only); ontologyIndex 2.12; pheatmap 1.0.13; ggplot2 4.0.3; ggrepel 0.9.8; cowplot 1.2.0 (env `stv-r`, unchanged) |
| Parameters | Classification: D-045 as revised by the owner (evidence status plus a descriptive robustness descriptor, D-049; the superseded rules U1 – U3 are kept as history, §47); Jaccard on full membership and leading edges; GO ancestors with is_a + part_of (GO release 2026-08-05); average-linkage clustering of 1 − leading-edge Jaccard within direction, cut 0.7 (0.5 / 0.9 sensitivity) |
| Methodological rationale | Integrates existing evidence without a new test; GO redundancy summarized descriptively, because overlapping gene sets are not independent findings |
| Official technical documentation | `?hclust`, `?cutree` (R stats); ontologyIndex `get_ontology`; GO documentation on the ontology structure (<https://geneontology.org/docs/ontology-relations/>) |
| Primary scientific reference | Jaccard (1912), similarity coefficient. Its use for gene-set overlap is standard practice; no method paper is cited for the descriptive clustering |
| DOI | 10.1111/j.1469-8137.1912.tb05611.x |
| Provisional APA 7 reference | Jaccard, P. (1912). The distribution of the flora in the alpine zone. *New Phytologist, 11*(2), 37–50. https://doi.org/10.1111/j.1469-8137.1912.tb05611.x |
| Decision supported | D-045, D-047, D-048, D-049 |
| Access info | 2026-10-26 (Jaccard reference verified with Crossref) |

## 47. D-045 revision: inferential evidence vs robustness diagnostics (methodological support)

| Field | Content |
|---|---|
| Process | Owner revision of D-045: separate inferential evidence (DESeq2, CAMERA) from exploratory discovery (fgsea) and from robustness diagnostics (replicate consistency, LOO, Cook's, ch00); explicit reviews of ATP synthase and ER ([`biological_integration.md`](biological_integration.md) §2) |
| Software + version | No new analysis. Documentation checked: limma 3.66.0 `?camera` (installed; Bioconductor 3.22); Bioconductor landing page (current release limma 3.68.5, Bioconductor 3.23) |
| Parameters | — |
| Methodological rationale | See the source table below. No source was found that validates a fixed leave-one-out share (such as 50 % of a theme's sets) as an inferential boundary for gene-set results; §38 already records that no pathway-level LOO method paper was found. LOO and the other diagnostics are therefore used descriptively (D-049) |
| Official technical documentation | limma `?camera`; GSEA User Guide <https://docs.gsea-msigdb.org/#GSEA/GSEA_User_Guide/> |
| Primary scientific reference | Schurch et al. (2016); Wu & Smyth (2012); Goeman & Bühlmann (2007); Tamayo et al. (2016); Gatti et al. (2010); Maciejewski (2014). Evaluated, not used for a claim: Geistlinger et al. (2021) |
| DOI | 10.1261/rna.053959.115; 10.1093/nar/gks461; 10.1093/bioinformatics/btm051; 10.1177/0962280212460441; 10.1186/1471-2164-11-574; 10.1093/bib/bbt002; 10.1093/bib/bbz158 |
| Provisional APA 7 reference | Schurch, N. J., Schofield, P., Gierliński, M., Cole, C., Sherstnev, A., Singh, V., Wrobel, N., Gharbi, K., Simpson, G. G., Owen-Hughes, T., Blaxter, M., & Barton, G. J. (2016). How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA, 22*(6), 839–851. https://doi.org/10.1261/rna.053959.115 <br> Wu & Smyth (2012), Tamayo et al. (2016): see §29. Goeman & Bühlmann (2007): see §34. <br> Gatti, D. M., Barry, W. T., Nobel, A. B., Rusyn, I., & Wright, F. A. (2010). Heading down the wrong pathway: On the influence of correlation within gene sets. *BMC Genomics, 11*, 574. https://doi.org/10.1186/1471-2164-11-574 <br> Maciejewski, H. (2014). Gene set analysis methods: Statistical models and methodological differences. *Briefings in Bioinformatics, 15*(4), 504–518. https://doi.org/10.1093/bib/bbt002 <br> Geistlinger, L., Csaba, G., Santarelli, M., Ramos, M., Schiffer, L., Turaga, N., Law, C., Davis, S., Carey, V., Morgan, M., Zimmer, R., & Waldron, L. (2021). Toward a gold standard for benchmarking gene set enrichment analysis. *Briefings in Bioinformatics, 22*(1), 545–556. https://doi.org/10.1093/bib/bbz158 |
| Decision supported | D-045 (owner revision), D-049 |
| Access info | 2026-10-26. Bibliographic metadata of all seven DOIs checked with the Crossref REST API. Abstracts via Europe PMC; full texts via Europe PMC for PMC4878611 (Schurch), PMC3458527 (Wu & Smyth), PMC4103537 (Maciejewski) and PMC7820859 (Geistlinger). GSEA User Guide HTML. The current-release `camera` help text could not be retrieved (code.bioconductor.org TLS certificate expired, not bypassed; reference-manual PDF not parsable), so the quoted text is from limma 3.66.0, the version used |

Source table (each source is cited only for the claim shown):

| Source | What was read | Claim supported (quoted or close paraphrase) | Use in D-045 |
|---|---|---|---|
| Schurch et al. (2016), <https://doi.org/10.1261/rna.053959.115> | Abstract and full text | "With three biological replicates, nine of the 11 tools evaluated found only 20%–40% of the significantly differentially expressed (SDE) genes identified with the full set of 42 clean replicates"; > 85 % only for genes changing more than fourfold; at least six replicates recommended (12 to find SDE genes of all fold changes); the leading tools control FDR at all replicate numbers; "there is a significant minority chance that one or more replicates within each condition should be rejected" | n = 3 means low sensitivity: the gene-level null is not evidence of absence, single replicates can be atypical, and removing a replicate (LOO) lowers sensitivity further. Gene-level evidence; it does not address gene-set tests |
| Wu & Smyth (2012), <https://doi.org/10.1093/nar/gks461> | Abstract and full text | Gene-permutation competitive tests are "extremely sensitive to inter-gene correlation", which is "non-ignorable even for experiments on homogeneous cell populations"; the variance of a set statistic is inflated by VIF = 1 + (m − 1)ρ̄, which "can be much larger than one when the test set contains many genes … even if the inter-gene correlation is quite small"; sample permutation "cannot be applied to experiments with small number of biological replicates"; CAMERA "controls the type I error rate correctly even for small sample sizes" | CAMERA is the correlation-aware inferential check of the fgsea patterns with n = 3; fgsea significance alone is exploratory |
| limma 3.66.0 `?camera` (official package documentation) | Installed help page | camera performs a competitive test in the sense of Goeman & Bühlmann and adjusts the variance of test statistics "by an estimated variance inflation factor"; with `inter.gene.cor = NA` it "gives rigorous error rate control for all sample sizes and all gene sets"; the default 0.01 is "a useful compromise between strict error rate control and interpretable gene set rankings" | The primary CAMERA run (estimated correlation) is the strict inferential check; the 0.01 run stays descriptive (D-034) |
| GSEA User Guide, <https://docs.gsea-msigdb.org/#GSEA/GSEA_User_Guide/> | "Permutation type", "False Discovery Rate (FDR)", "GSEAPreranked Page", "Gene Sets and GSEA", "Enrichment in Phenotype" | Phenotype permutation is recommended with at least seven samples per phenotype; gene_set permutation is for fewer samples, and with it "the correlations between the genes in the dataset and the genes in the gene set are not preserved", whereas preserving them gives "a more biologically reasonable (more stringent) assessment of significance"; GSEAPreranked always permutes by gene set; with few samples and gene_set permutation the assessment is "less stringent" and a stricter FDR (e.g. 5 %) is advised; near-duplicate gene sets can skew the FDR; many enriched sets may "represent the same biological signal" or technical artefacts and then give "potentially biased evidence for interesting hypotheses" | fgsea (preranked, gene-set permutation) is exploratory discovery here; its FDR does not account for gene–gene correlation; GO redundancy limits the number of independent findings |
| Goeman & Bühlmann (2007), <https://doi.org/10.1093/bioinformatics/btm051> | Abstract | Defines competitive vs self-contained tests; models that take genes as sampling units rest on an "unrealistic independence assumption", so their p-values "can be wildly anti-conservative" | Competitive gene-sampling significance is not inferential support on its own |
| Tamayo et al. (2016), <https://doi.org/10.1177/0962280212460441> | Abstract | On 50 benchmark datasets, "gene-gene correlations cannot be ignored due to the significant variance inflation they produced on the enrichment scores and should be taken into account when estimating gene set enrichment significance" | As above (view of GSEA authors) |
| Gatti et al. (2010), <https://doi.org/10.1186/1471-2164-11-574> | Abstract | Over 200 GEO datasets: an independence-based gene-set test "produces very high false positive rates" when treatment groups are randomized, and "gene sets with high internal correlation are more likely to be declared significant" | Indirect (over-representation testing of microarray data, not preranked GSEA): highly correlated sets, such as the ribosome theme here (median VIF 36.6), are the most prone to optimistic significance |
| Maciejewski (2014), <https://doi.org/10.1093/bib/bbt002> | Abstract and full text | In simulations, gene-randomization and gene-sampling parametric methods "generate excessive number of false-positive results, with the effect increasing with growing correlation between genes"; it notes that Nam & Kim recommend gene randomization for small sample sizes | Gene-randomization significance is optimistic under correlation, so the fgsea patterns need a correlation-aware check |
| Geistlinger et al. (2021), <https://doi.org/10.1093/bib/bbz158> | Abstract and full text | Benchmark of 10 enrichment methods: the fraction of significant sets depends on the null hypothesis tested, and results should be interpreted by test type | Evaluated only. Not used as support, because it does not test preranked gene-permutation GSEA with three replicates |
