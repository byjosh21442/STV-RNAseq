#!/usr/bin/env bash
# STV_GSE137303 - Phase 6B closure, organelle assessment (Task 8): POST-HOC sequence-identity diagnostic.
# Corrects the -max_hsps 1 limitation of scripts/bash/ch00_organelle_blast.sh (kept unchanged): all HSPs are reported so
# that coverage can be computed as the union of query intervals (criteria frozen in
# docs/organelle_requantification_assessment.md section 1.4). Queries: (1) every ITAG4.0 cDNA; (2) the 177 segments of
# SL4.0ch00 delimited by its 100-N gaps. Subjects: tomato plastome NC_007898.3 and mitogenome NC_035963.1 (files and SHA-256
# from scripts/bash/download_organelle_reference.sh). BLAST+ from the existing conda env rachis-qiime2-2026.4 (read-only).
# No quantification is run. Usage (project root): bash scripts/bash/organelle_multihsp_blast.sh
set -euo pipefail
cd "$(dirname "$0")/../.."
source "$HOME/miniconda3/etc/profile.d/conda.sh"
env=rachis-qiime2-2026.4
org=reference/functional/organelle
man=results/functional_robustness/ch00/organelle_reference_manifest.tsv
o=results/organelle_assessment/blast
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$o"
bl() { conda run -n "$env" "$@"; }

# organelle genome files must match the recorded manifest
awk -F'\t' 'NR>1 && $3=="fasta"{print $6"  "$1}' "$man" | sha256sum -c --quiet
{ bl blastn -version | head -n 1; bl makeblastdb -version | head -n 1; } > "$o/blast_version.txt"

cat "$org/NC_007898.3.fna" "$org/NC_035963.1.fna" > "$tmp/organelle_genomes.fna"
bl makeblastdb -in "$tmp/organelle_genomes.fna" -dbtype nucl -out "$tmp/orgnt" > /dev/null

# ch00 segments between 100-N gaps (1-based start/end on SL4.0ch00 in the header)
zcat reference/SGN/SL4.0/S_lycopersicum_chromosomes.4.00.fa.gz | awk '/^>/{p=($1==">SL4.0ch00"); next} p{printf "%s",$0} END{print ""}' |
  python3 -c '
import re, sys
s = sys.stdin.read().strip()
pos = 0
# gap runs are 100 N, in upper or lower case (176 gaps -> 177 segments)
for i, seg in enumerate(re.split(r"([Nn]{100,})", s)):
    if seg and seg[0] not in "Nn":
        print(f">ch00seg{i//2+1:03d}_{pos+1}_{pos+len(seg)}\n{seg}")
    pos += len(seg)
' > "$tmp/ch00_segments.fa"
grep '^>' "$tmp/ch00_segments.fa" | sed 's/^>//' | awk -F_ '{print $1"\t"$2"\t"$3"\t"$3-$2+1}' |
  sed '1i segment\tstart\tend\tlength' > "$o/ch00_segments.tsv"
test "$(($(wc -l < "$o/ch00_segments.tsv") - 1))" -eq 177 || { echo "unexpected ch00 segment count" >&2; exit 1; }

fmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
# all HSPs (no -max_hsps); two subjects only, so -max_target_seqs 2 reports every subject
bl blastn -task megablast -query reference/SGN/ITAG4.0/ITAG4.0_cDNA.fasta -db "$tmp/orgnt" -evalue 1e-20 \
  -max_target_seqs 2 -num_threads 4 -outfmt "$fmt" > "$o/blastn_itag40_cdna_vs_organelle_allhsps.tsv"
bl blastn -task megablast -query "$tmp/ch00_segments.fa" -db "$tmp/orgnt" -evalue 1e-20 \
  -max_target_seqs 2 -num_threads 4 -outfmt "$fmt" > "$o/blastn_ch00_segments_vs_organelle_allhsps.tsv"

printf 'item\tvalue\n' > "$o/blast_run_summary.tsv"
printf 'cdna_queries\t%s\n' "$(grep -c '^>' reference/SGN/ITAG4.0/ITAG4.0_cDNA.fasta)" >> "$o/blast_run_summary.tsv"
printf 'cdna_hsp_rows\t%s\n' "$(wc -l < "$o/blastn_itag40_cdna_vs_organelle_allhsps.tsv")" >> "$o/blast_run_summary.tsv"
printf 'cdna_queries_with_hit\t%s\n' "$(cut -f1 "$o/blastn_itag40_cdna_vs_organelle_allhsps.tsv" | sort -u | wc -l)" >> "$o/blast_run_summary.tsv"
printf 'ch00_segments\t177\n' >> "$o/blast_run_summary.tsv"
printf 'ch00_segment_hsp_rows\t%s\n' "$(wc -l < "$o/blastn_ch00_segments_vs_organelle_allhsps.tsv")" >> "$o/blast_run_summary.tsv"
cat "$o/blast_version.txt" "$o/blast_run_summary.tsv"
