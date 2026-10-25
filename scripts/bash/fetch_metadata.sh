#!/usr/bin/env bash
# Fetch public metadata for GSE137303 / PRJNA565137 from GEO, NCBI SRA (E-utilities) and ENA.
# Downloads small text/XML records only (no FASTQ, no references). Run from the project root.
set -euo pipefail

OUT="data/metadata/sources"
mkdir -p "$OUT"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$OUT/retrieved_utc.txt"

GSE="GSE137303"
PRJ="PRJNA565137"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# GEO series + all sample records (SOFT text, metadata only)
curl -fsSL "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${GSE}&targ=all&form=text&view=brief" \
  -o "$OUT/${GSE}_brief.soft.txt"

# NCBI SRA: run info (CSV) and full experiment packages (XML) for the BioProject
ids="$(curl -fsSL "${EUTILS}/esearch.fcgi?db=sra&term=${PRJ}%5BBioProject%5D&retmax=100" \
  | grep -o '<Id>[0-9]*</Id>' | sed 's/<[^>]*>//g' | paste -sd, -)"
echo "SRA UIDs: ${ids}"
sleep 1
curl -fsSL "${EUTILS}/efetch.fcgi?db=sra&id=${ids}&rettype=runinfo&retmode=text" -o "$OUT/${PRJ}_sra_runinfo.csv"
sleep 1
curl -fsSL "${EUTILS}/efetch.fcgi?db=sra&id=${ids}&retmode=xml" -o "$OUT/${PRJ}_sra_experiment_packages.xml"

# ENA read_run report (independent mirror; includes FASTQ sizes/MD5s)
ENA_FIELDS="study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,sample_alias,sample_title,experiment_title,scientific_name,tax_id,instrument_platform,instrument_model,library_name,library_strategy,library_source,library_selection,library_layout,nominal_length,read_count,base_count,fastq_bytes,fastq_md5,fastq_ftp,sra_bytes,first_public"
curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PRJ}&result=read_run&fields=${ENA_FIELDS}&format=tsv" \
  -o "$OUT/${PRJ}_ena_read_run.tsv"

# GEO supplementary file listing (sizes only)
curl -fsSL "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/${GSE}/suppl/" -o "$OUT/${GSE}_suppl_listing.html"

ls -l "$OUT"
