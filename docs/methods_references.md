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
