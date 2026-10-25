#!/usr/bin/env bash
# Verify which annotation the GEO count table is keyed by: stream the NCBI RefSeq SL2.50 GFF
# (Annotation Release 102, ~8 MB, not saved) and compare its gene ID= / GeneID pairs with the table.
set -eu
F="data/processed/original_geo/GSE137303_count.txt.gz"
GFF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.3_SL2.50/GCF_000188115.3_SL2.50_genomic.gff.gz"
tmp="$(mktemp -d)"

curl -fsSL "$GFF_URL" | zcat | awk -F'\t' '$3=="gene" || $3=="pseudogene" {
    id=""; gid="";
    if (match($9, /(^|;)ID=[^;]+/))            { id=substr($9, RSTART, RLENGTH); sub(/.*ID=/,"",id) }
    if (match($9, /Dbxref=GeneID:[0-9]+/))     { gid=substr($9, RSTART+7, RLENGTH-7) }
    print id "\t" gid }' | sort -u > "$tmp/gff_pairs.tsv"

zcat "$F" | awk -F'\t' 'NR>1 { gid=$2; gsub(/"/,"",gid); sub(/,.*/,"",gid); print $1 "\t" gid }' | sort -u > "$tmp/geo_pairs.tsv"

echo "GFF gene/pseudogene features : $(wc -l < "$tmp/gff_pairs.tsv")"
echo "GEO table rows               : $(wc -l < "$tmp/geo_pairs.tsv")"
echo "Exact (ID, GeneID) matches   : $(comm -12 "$tmp/gff_pairs.tsv" "$tmp/geo_pairs.tsv" | wc -l)"
echo "GEO rows not in GFF          : $(comm -13 "$tmp/gff_pairs.tsv" "$tmp/geo_pairs.tsv" | wc -l)"
comm -13 "$tmp/gff_pairs.tsv" "$tmp/geo_pairs.tsv" | head -n 5
echo "ID-only matches (tracking_id == GFF ID=) : $(join -t$'\t' <(cut -f1 "$tmp/gff_pairs.tsv" | sort -u) <(cut -f1 "$tmp/geo_pairs.tsv" | sort -u) | wc -l)"
echo "GFF pairs for the unmatched IDs:"
comm -13 "$tmp/gff_pairs.tsv" "$tmp/geo_pairs.tsv" | cut -f1 | head -n 3 | grep -F -w -f - "$tmp/gff_pairs.tsv" || true
rm -rf "$tmp"
