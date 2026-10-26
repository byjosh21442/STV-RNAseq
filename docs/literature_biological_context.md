# Literature and biological context (Phase 7)

Status: Phase 7, 2026-10-26. The literature is used **after** the internal evidence (L1 – L4) and only for biological
plausibility. It did not choose any threshold, gene set, gene, sample or class.

## 1. Search and selection

- Reproducible Europe PMC searches (`scripts/bash/literature_search_europepmc.sh`; metadata only, saved to
  `results/biological_integration/literature/`). The queries and hit counts on 2026-10-26 were:
  - Q1, "southern tomato virus": 32;
  - Q2, Amalgavirus / Amalgaviridae: 67;
  - Q3, tomato + virus + RNA-seq/transcriptome + leaf: 84;
  - Q4, persistent / endorna- / partiti- / cryptic virus + transcriptome/expression + plant: 19.
- Records were selected manually by the project priority:
  1. STV in tomato;
  2. STV molecular biology;
  3. STV transcriptomics;
  4. directly relevant persistent-virus biology;
  5. tomato leaf transcriptomic responses to viruses;
  6. plant–virus responses touching the observed themes.
- Abstracts were read through the Europe PMC REST API. Full texts were read where stated.
- Two further targeted Europe PMC queries were run: ribosomal-protein genes in virus infection, and PepMV tomato
  transcriptome.
- Unsuccessful attempts are recorded in `run_log.md`. A query with `sort=RELEVANCE` returned HTTP 503 (dropped), and one
  request timed out (repeated).
- Preprints are labelled **[preprint]**. Reviews and taxonomic reports are labelled as such. Records known only by title
  were **not** used.

## 2. STV biology: established vs our hypotheses (Task 5)

### 2.1 Established from primary sources

| Topic | What is established | Source |
|---|---|---|
| Taxonomy | Species *Amalgavirus lycopersici* (virus: southern tomato virus, STV), genus *Amalgavirus*, family *Amalgaviridae*, order *Durnavirales* | ICTV Report, Amalgaviridae (interim, October 2024); Jindal & Ali 2025 |
| Genome | dsRNA, monopartite, about 3.5 kb. Two partially overlapping ORFs: ORF1 encodes a 377-aa protein of unknown function (p42, a putative coat protein); ORF2 encodes the RdRp, likely expressed by +1 ribosomal frameshifting. Virions are not known | Sabanadzovic et al. 2009; ICTV report; Elvira-González et al. 2020 (no particles seen) |
| Transmission | Efficiently seed-transmitted. Not transmitted mechanically or by grafting. Vertical transmission only. About 86 % transmission to progeny in M82 | Sabanadzovic et al. 2009; Elvira-González et al. 2020; Turco et al. 2018; Fukuhara et al. 2020 |
| Persistence and distribution | Persistent lifestyle; systemic in leaves, stems, seeds, shoot apices and root tips; worldwide, with high incidence (up to about 43 % in surveyed regions of China and Pakistan) | Hao et al. 2023; Hussain et al. 2022; Elvira-González et al. 2020 |
| Genetic diversity | Very low diversity; negative (purifying) selection; no or weakly supported recombination | Elvira-González et al. 2020 (Virus Genes); Jindal & Ali 2025; Hussain et al. 2022 |
| Symptoms | First associated with fruit disorders. Single infections show no apparent symptoms and no leaf cell ultrastructural changes. Found in asymptomatic plants | Elvira-González et al. 2021; Elvira-González et al. 2020; Alcalá-Briseño et al. 2017; ICTV ("asymptomatic") |
| Effects on the host | STV-infected M82 plants were taller, produced more fruit and had higher seed germination than STV-free plants. The authors suggest a mutualistic interaction | Fukuhara et al. 2020 (abstract) |
| Gene expression | Ethylene biosynthesis and signalling genes reported down-regulated (Fukuhara et al. 2020, abstract). Five tomato miRNAs differentially expressed in STV-infected plants, with few STV-derived siRNAs (Elvira-González et al. 2020). A reanalysis of the same RNA-seq data found no DEG at FDR < 0.05 (Badeli et al. 2026 **[preprint]**) | as cited |
| Interaction with other viruses | Synergism with CMV and PepMV in co-infection (higher CMV titre and symptoms early); STV abolished the CMV–PepMV antagonism in triple infection. STV is detected in mixed infections associated with tomato yellow stunt disease (association, not causation) | Elvira-González et al. 2021; Hussain et al. 2022, 2026 |
| Defence | Indirect only: miRNA targets *predicted* in defence and stress pathways; low vsiRNA levels in single infection | Elvira-González et al. 2020, 2021 |
| Hormones / ethylene | Only the transcriptomic claim of Fukuhara et al. (abstract); a non-significant ACO1 trend in Badeli et al. **[preprint]**. No hormone measurement found | as cited |
| Photosynthesis, translation, ribosomes, stress | **No primary evidence found for STV.** The +1 frameshift concerns viral, not host, translation | — |
| Relationship to acute viruses | Persistent viruses form long-term, often asymptomatic, vertically transmitted infections; acute viruses cause disease and broad host reprogramming | Roossinck 2010 (review); Takahashi et al. 2019 (review); Vendrell-Mir et al. 2021 (vertically transmitted amalgavirus in a moss) |

### 2.2 Hypotheses generated from this dataset (not established)

These come only from exploratory, fgsea-only evidence (no CAMERA support, no gene-level significance):
- **H1.** In STV-infected M82 leaves, ethylene-related transcripts and GO sets rank lower; the direction agrees with
  Fukuhara et al. and Badeli et al.
- **H2.** Defence / biotic-response and associated stress-related sets (water deprivation, glutathione) rank lower. This is
  one partly shared tendency, not independent signals.
- **H3.** Plastid / photosynthesis-related sets rank higher, independently of the organelle-identical feature class.

The ribosome / translation and ER tendencies are EXPLORATORY but strongly sample-sensitive (each concentrated on one
sample), so they are not listed as hypotheses here. The ATP-synthase signal is UNSUPPORTED (carried primarily by the
organelle-identical feature class). Classification: D-045 as revised by the owner
([biological_integration.md](biological_integration.md) §2 and §4).

## 3. Source table (Task 4)

Relation codes:
- **agrees** — same direction or same conclusion;
- **disagrees** — opposite;
- **NC** — not directly comparable: different virus type, tissue, design or method.

| # | Source | Type | Organism / virus | Tissue | Design | Method, n | Principal finding | Comparable? | Theme(s) | Relation |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Sabanadzovic et al. 2009, *Virus Res* 140:130 | primary | tomato / STV | field + greenhouse plants | virus discovery | cloning, sequencing, transmission tests | dsRNA 3.5 kb, 2 ORFs, +1 frameshift RdRp; seed-transmitted, not mechanically or by graft | no (biology) | STV biology | context |
| 2 | ICTV Report, *Amalgaviridae* (interim, 2024) | taxonomic report | plants / amalgaviruses | — | — | — | *Amalgavirus lycopersici*; asymptomatic, seed-transmitted | no | STV biology | context |
| 3 | Fukuhara et al. 2020, *Arch Virol* 165:11 | primary (abstract only) | tomato M82 / STV | leaves (RNA-seq, 1-month-old per GEO) | STV-free vs STV-infected | phenotypes; RNA-seq 3 vs 3 (GEO) | infected taller, more fruit, higher germination; ethylene biosynthesis and signalling genes down | **same data** | hormone / ethylene | agrees in direction (pathway level); gene-level significance not reproduced |
| 4 | Badeli et al. 2026 **[preprint]**, Research Square | reanalysis | tomato M82 / STV | leaves | same 6 libraries | Rsubread + featureCounts (SL3.0), DESeq2 | 0 DEG; ACO1 −1.47 (padj 0.225); PR-locus −1.12 (FDR 0.083); no GO term at FDR < 0.05 | **same data** | hormone, defence, cell wall | agrees (no DEG; ethylene direction) |
| 5 | Elvira-González et al. 2020, *Eur J Plant Pathol* 156:615 | primary | tomato / STV | leaves | infected vs healthy | small-RNA seq, TEM | 5 miRNAs differentially expressed; few vsiRNAs; no symptoms or ultrastructural changes | NC (sRNA) | defence / stress (predicted targets); photosynthesis (no ultrastructural change) | NC |
| 6 | Elvira-González et al. 2021, *Microorganisms* 9:689 | primary | tomato / STV + CMV, PepMV | plants | single, double, triple infections | symptoms, titre, sRNA | STV is synergistic with CMV and PepMV; abolishes CMV–PepMV antagonism; asymptomatic alone | NC (mixed infections) | STV biology | context |
| 7 | Hao et al. 2023, *J Plant Dis Prot* | brief report | tomato / STV | leaves, stems, seeds, apex, roots | localization | in situ hybridization | systemic distribution, including leaves | no | STV biology | context |
| 8 | Elvira-González et al. 2020, *Virus Genes* 56:668 | primary | tomato / STV | — | isolates worldwide | CP gene sequencing | very low diversity; negative selection | no | STV biology | context |
| 9 | Jindal & Ali 2025, *Viruses* 17:1193 | primary | tomato / STV | — | Oklahoma isolates | complete genomes, phylogeny | two clades; low diversity; purifying selection | no | STV biology | context |
| 10 | Hussain et al. 2022, *Viruses* 14:2481 | primary | tomato / STV | — | survey China/Pakistan | RT-PCR, genomes | infection rates up to 43 %; low diversity | no | STV biology | context |
| 11 | Hussain et al. 2026, *Plant Dis* 110:592 | primary | tomato / 5 viruses incl. STV | leaves, seeds | mixed-infection disease | sRNA-seq, RT-PCR | STV is part of the ToYSD virus complex (association) | NC | STV biology | context |
| 12 | Turco et al. 2018, *MPMI* 31:707 | primary | tomato, potato / mixed viromes | leaves | cross-protection trial | sRNA-omics | STV reconstructed in some tomato plants; mechanically non-transmissible | NC | STV biology | context |
| 13 | Alcalá-Briseño et al. 2017, *Genome Announc* 5:e01374-16 | genome announcement | tomato 'Sweet Hearts' / STV | — | — | sRNA sequencing | STV genome in asymptomatic plants | no | STV biology | context |
| 14 | Vendrell-Mir et al. 2021, *Plant J* 108:1786 | primary | *Physcomitrium patens* / amalgavirus PHPAV1 | whole plant | accessions, crosses | RNA-seq, genetics | vertical transmission through both germlines | NC | persistent biology | context |
| 15 | Roossinck 2010, *Phil Trans R Soc B* 365:1899 | review | plant viruses | — | — | — | persistent vs acute lifestyles | NC | all | context |
| 16 | Takahashi et al. 2019, *Front Microbiol* 10:2764 | review | plant viruses | — | — | — | latency and impact on plants | NC | all | context |
| 17 | Khankhum et al. 2016, *Front Microbiol* 7:1905 | data paper | common bean / endornaviruses | leaves | infected vs free | RNA-seq 3 vs 3 | 132 DEGs; oxidation-reduction is the most frequent GO category (not an enrichment test) | NC (other host and virus) | heme / glutathione | NC |
| 18 | Kamitani et al. 2016, *FEMS Microbiol Ecol* 92:fiw176 | primary | *Arabidopsis halleri* / natural viromes | leaves | natural population | RNA-seq | asymptomatic infection common; AGO2 up with TuMV; novel partitivirus | NC | defence | NC |
| 19 | Hanssen et al. 2011, *Plant Physiol* 156:301 | primary | tomato / PepMV (acute) | seedlings | mild vs aggressive isolate | Affymetrix array (> 22,000 transcripts) | photosynthesis repressed; SA-mediated defence stronger with the aggressive isolate | NC (acute) | photosynthesis, defence, hormone | NC (opposite direction for photosynthesis) |
| 20 | Wang et al. 2024, *Int J Mol Sci* 25:4012 | primary | tomato / ToBRFV (acute) | leaves, 21 dpi | infected vs mock | RNA-seq | 522 DEGs: wounding, stress, protein folding, defence | NC (acute) | defence, ER, water/abiotic | NC |
| 21 | He et al. 2024, *Plants* 13:3198 | primary | tomato / TYLCTHV (acute) | leaves | infected vs mock | RNA-seq | 2,322 up / 1,377 down; photosynthesis (LHC, PSI/II, cyt b6f, F-type ATPase) disrupted; immunity induced | NC (acute) | photosynthesis, ATP synthase, defence | NC (opposite directions) |
| 22 | Prigigallo et al. 2025, *Sci Rep* 15:3630 | primary | tomato (Sw-5b) / TSWV resistance-breaking (acute) | leaves (field) | titre gradient | RNA-seq, 18 samples | DEGs 33 – 44 % of genes, scaling with titre; photosynthesis and protein biosynthesis down; oxidoreductases and defence up | NC (acute) | photosynthesis, ribosome, defence, heme | NC (opposite directions) |
| 23 | Romero-Rodríguez et al. 2023, *BMC Plant Biol* 23:651 | primary | tomato / TYLCV (acute) | leaves, time course | infected vs mock | RNA-seq, sRNA, methylome | substantial reprogramming; gene silencing and immune response induced | NC (acute) | defence | NC |
| 24 | Yue et al. 2024, *Plant Cell Environ* 47:2660 | primary | tomato / ToCV (acute) | leaves | infected vs healthy | scRNA-seq (23,226 cells) | mesophyll function shift; thinner leaves; fewer chloroplast lamellae | NC (acute) | photosynthesis | NC |
| 25 | Whitham et al. 2003, *Plant J* 33:271 | primary | Arabidopsis / 5 acute RNA viruses | inoculated leaves | time course | Affymetrix arrays | common virus-general defence and stress responses | NC (acute) | defence | NC |
| 26 | Zorzatto et al. 2015, *Nature* 520:679 | primary | Arabidopsis / begomovirus | — | genetics | molecular | NIK1–RPL10–LIMYB represses ribosomal-protein genes as antiviral immunity | NC | ribosome | NC |

DOIs and access details: `methods_references.md` §44 – §46.

**What the literature contributes.**
- STV is consistently described as a persistent, seed-transmitted virus that is asymptomatic alone. The absence of
  gene-level DE in this dataset is **consistent with** that description. With n = 3 per group, it does not prove that
  STV has no effect.
- Every tomato leaf transcriptome found for a virus is an **acute** infection, with hundreds to thousands of DEGs and
  defence induction. They show what a strong viral response looks like, but they are not comparable with a persistent
  virus. No published tomato transcriptome of another persistent virus was found.
- No primary source supports or refutes the ribosome, ER or photosynthesis tendencies **for STV**.

## 4. Comparison with Fukuhara et al. 2020 (Task 6)

Sources:
- for the original study: the GEO series record (archived in `data/metadata/sources/GSE137303_brief.soft.txt`) and the
  article abstract. **The full text was not accessible**, so statements about the statistical method are marked
  unverified.

| Question | Fukuhara et al. (original) | Current reanalysis | Agreement |
|---|---|---|---|
| Samples | 3 STV-free + 3 STV-infected M82, 1-month-old leaves (GEO) | same 6 libraries (ENA MD5-verified) | same |
| Library | TruSeq stranded mRNA (GEO) | reverse-stranded (`SR`) measured by Salmon in 6/6 | consistent |
| Reference assembly | SL2.50 (GEO) | SL4.0 (SGN) | different |
| Annotation | ITAG2.4 per GEO metadata. The GEO processed table is keyed to NCBI RefSeq SL2.50 identifiers (Phase 2 check) | ITAG4.0 | different |
| Alignment / quantification | TopHat 2.0.13 + Bowtie2 2.2.3; Cuffnorm 2.2.1 | Salmon 2.7.0 selective alignment with full-genome decoys; tximport | different |
| Normalization | FPKM (Cuffnorm) | DESeq2 median-of-ratios with tximport length offsets | different |
| Statistical approach | not described in the accessible record. Badeli et al. state that no genome-wide FDR was reported (secondary, **unverified**) | DESeq2 negative-binomial Wald test, prespecified filter | not comparable |
| Multiple-testing correction | unknown | Benjamini–Hochberg, FDR 0.05 | unknown |
| Reported expression change | ethylene biosynthesis and signalling genes down-regulated | 0 genes at padj < 0.05. ACO1 lower in 9/9 sample pairs, padj 0.297. Ethylene GO sets lower-ranked (fgsea padj 0.004 – 0.041; CAMERA FDR 0.999) | **direction agrees; gene-level significance not reproduced** |
| Pathway analysis | not described in the abstract | fgsea + CAMERA + robustness diagnostics | — |
| Phenotypes | taller plants, more fruit, higher germination | not measured (RNA-seq only) | not assessable |
| Biological interpretation | mutualistic interaction; STV-infected plants have superior traits | the transcriptome cannot test mutualism; ethylene tendency is exploratory | not assessable |

**Possible reasons for the different gene-level conclusions.** These are listed neutrally; none was tested here, so none
is claimed as the cause.
1. Reference and annotation versions (SL2.50 / ITAG2.4 or RefSeq vs SL4.0 / ITAG4.0). Gene models, transcript
   lengths and gene identities differ. Badeli et al. note that some original transcripts were near zero in their
   pipeline.
2. Quantification: TopHat / Cuffnorm FPKM vs Salmon estimates.
3. Normalization: FPKM vs median-of-ratios.
4. Statistical framework and multiple-testing correction, which is unknown for the original.
5. Small n (3 vs 3). Schurch et al. (2016) report 20 – 40 % sensitivity with three replicates.
6. Sample heterogeneity. STV_free_rep1 is separated in expression space (Phase 5, D-022), and the ACO1 difference partly
   depends on which free sample is compared.

The original study is not "wrong". The two analyses answer different questions with different statistical frameworks,
and both find lower ethylene-related expression in STV-infected leaves.

## 5. Comparison with Badeli et al. 2026 **[preprint]** (Task 7)

Source: full text read directly, Research Square HTML v1, DOI 10.21203/rs.3.rs-8823240/v1, accessed 2026-10-26. It is
still a preprint (not peer reviewed).

| Item | Badeli et al. (as stated) | Current reanalysis | Agreement |
|---|---|---|---|
| Samples | same 6 libraries (GSE137303) | same | same |
| Pre-processing | Trimmomatic 0.39 trimming | no trimming (D-009) | different |
| Reference | Ensembl SL3.0 genome + STV genome (KX949574) as an extra chromosome; annotation `Solanum_lycopersicum.SL3.0.61` (ITAG3.0) | SL4.0 + ITAG4.0, genome decoys | different |
| Quantification | Rsubread alignment + featureCounts; the preprint reports 34,075 gene features, which is the ITAG4.0 gene count (not resolved here) | Salmon + tximport; 33,976 genes | different |
| Filter, tested genes | < 10 total counts removed; about 25,000 tested | ≥ 10 counts in ≥ 3 samples; 18,277 tested | different |
| DE method | DESeq2 1.42.0, FDR < 0.05 and \|log2FC\| ≥ 1; apeglm; relaxed FDR < 0.1 / 0.2 lists for ORA | DESeq2 1.50.2, padj < 0.05; apeglm; no relaxed list | similar core method |
| DEG result | **0** | **0** | **agree** |
| ACO1 | log2FC −1.47 (described as shrunken), raw p 3.7 × 10⁻⁵, padj 0.225 | MLE −1.81, apeglm −1.56, raw p 4.9 × 10⁻⁵, padj 0.297; 9/9 pairs lower | **agree** (direction, n.s.) |
| PR-locus Solyc12g009240 | −1.12, FDR 0.083 | MLE −1.25, apeglm −1.10, padj 0.193. ITAG4.0 names it an ethylene-responsive TF | agree (direction, n.s.); annotation differs |
| ERF017 (Solyc07g007770), NAC1 (Solyc01g095080) | +3.97, FDR 0.157; about −0.4, FDR > 0.5 | not tested (removed by the count filter) | not comparable |
| Ethylene interpretation | slight down-regulation of ethylene-related genes; "plausible explanation" of the growth phenotype | exploratory directional tendency; no gene-level significance; no causal inference | agree on direction; this project makes no mechanistic claim |
| STV_free_rep1 | PC1 (about 58 %) separates STV−_rep1 | PC1 (49.3 %) separates STV_free_rep1 | **agree** |
| STV reads | < 0.01 % of reads | 0 in free; 15.3 – 34.4 per million in infected | consistent |
| GO analysis | ORA on relaxed lists: no term at FDR < 0.05; "ethylene biosynthetic process" and "response to fungus" at p ≈ 0.1 | fgsea: "response to fungus" NES −1.89, padj 0.013 (Tier 3, CAMERA FDR 0.999); ethylene biosynthetic process not tested (4 genes) | different methods; neither survives an FDR-controlled, correlation-aware test |
| Design statement | "the STV-free line was a cured sub-line of M82" (citing Fukuhara et al.) | not verifiable from the accessible original record; documented as a secondary, unverified statement | — |

**Observations.**
- The preprint describes the ACO1 replicate pattern in two ways: "one control plant had much higher ACO1", and "one
  STV-free sample having atypically low ACO1". In this reanalysis the free samples are 9.77 / 9.66 / 8.53 (blind VST),
  and STV_free_rep3 is the low one.
- None of these points contradicts the sample identity or design used here (STOP 2 not triggered).

## 6. Not found or not used

- No primary measurement of ethylene or other hormones in STV-infected plants.
- No STV study on photosynthesis, ribosomes or translation.
- No transcriptome of another persistent virus in tomato.
- The Fukuhara et al. full text (methods and statistics) was not accessible.
- Records known only by title were not cited.
