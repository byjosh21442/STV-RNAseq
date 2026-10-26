# Mercator4 submission record — STV_GSE137303 (Phase 6, D-029)

**Purpose.** Assign MapMan4 functional bins to the ITAG4.0 proteins and visualize the DESeq2 effect estimates in MapMan.
- This is **descriptive only**: no Mercator4/MapMan over-representation analysis is run, and MapMan colours are **not** evidence of statistical significance.
- The primary DESeq2 analysis has **0** FDR-significant genes.

**Why manual.** Mercator4 is a web service (Jülich Research Centre / PlabiPD). On 2026-10-26:
- The page <https://www.plabipd.de/mercator_main.html> showed **Mercator4 v8.0**, a job-submission form (Protein / DNA; optional Swiss-Prot annotations; FASTA upload; job name; e-mail) and the notice "The Mercator4 website is currently undergoing scheduled maintenance".
- The FAQ and version-history pages returned HTTP 404 or 503, so the input limits could not be verified.
- No official programmatic API was found, and none is invented here.

## 1. Input (prepared by `scripts/R/07_mapman_mercator_prep.R`)

| Item | Value |
|---|---|
| FASTA to upload | `reference/SGN/ITAG4.0/ITAG4.0_proteins.fasta`, **unmodified**. SGN ITAG4.0 release; sha256 and record counts in `mercator4_input_fasta_record.tsv` |
| Sequence type | Protein |
| Identifier manifest | `mercator4_input_manifest.tsv` (protein_id → gene_id, length, stop codons, tested / ranked flags) |
| Known features | 32,522 records end with `*` (stop); 334 records contain internal `*` (see manifest) |
| Expression file for MapMan | `mapman_expression_STV_infected_vs_free.txt`. Column 1 `IDENTIFIER` = ITAG4.0 protein ID exactly as in the FASTA. Display column `log2FC_apeglm_infected_vs_free` (positive = higher in STV_infected). Also included: MLE log2FC and Wald statistic. Rows: the 18,277 tested genes; untested genes are omitted (shown as "no data") |

## 2. Manual steps (to be executed by the owner or under owner supervision)

1. Open <https://www.plabipd.de/mercator_main.html> and note the version string shown (expected "Mercator4 v8.0").
2. Optional: run the site's "FASTA validation" tool on the FASTA.
   - If it rejects the terminal `*`, **do not edit the reference file**. Create a derived file and record its sha256:
     `sed '/^>/!s/\*$//' reference/SGN/ITAG4.0/ITAG4.0_proteins.fasta > results/functional/mapman/mercator4_input/ITAG4.0_proteins_no_terminal_stop.fasta`
   - If it also rejects the internal `*` (334 records), stop and ask for an owner decision. Do not delete or alter those sequences silently.
3. Job submission:
   - sequence type **Protein**;
   - **Swiss-Prot annotations: not selected**: the bins are the target, and the extra BLAST layer is not needed;
   - job name `STV_GSE137303_ITAG4.0_proteins`;
   - e-mail optional (do not record personal e-mail addresses in the repository).
4. Download the annotation result unchanged into `results/functional/mapman/mercator4_output/` (keep the original file name). Record its sha256 in the table below.
5. Check identifier compatibility before loading into MapMan. Compare 10 random identifiers of the Mercator4 output with column 1 of the expression file.
   - If Mercator4 changed the identifier case, create a case-matched copy of the expression file and record the command. Never edit identifiers by hand.
6. MapMan desktop: load the Mercator4 mapping and the expression file; display `log2FC_apeglm_infected_vs_free`.
   - Use a symmetric colour scale; record the range used.
   - Use a caption such as: "MapMan visualization of DESeq2 apeglm log2 fold changes (STV_infected vs STV_free). Colours show effect estimates only; no gene is FDR-significant (padj < 0.05: 0 of 18,277)".

## 3. Submission record (fill in when executed)

| Field | Value |
|---|---|
| Date (UTC) | _pending_ |
| Mercator4 version shown | _pending_ (v8.0 displayed on 2026-10-26) |
| Input file + sha256 | _pending_ |
| Options (type / Swiss-Prot) | _pending_ |
| Job ID | _pending_ |
| Output files + sha256 | _pending_ |
| Identifier check result | _pending_ |
| MapMan version and colour range | _pending_ |
| Executed by | _pending_ |

References:
- Schwacke et al. (2019), *Mol Plant* 12(6):879–892, doi:10.1016/j.molp.2019.01.003.
- Bolger, Schwacke & Usadel (2021), *Methods Mol Biol* pp. 195–212, doi:10.1007/978-1-0716-1609-3_9.
