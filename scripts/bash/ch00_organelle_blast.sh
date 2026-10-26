#!/usr/bin/env bash
# STV_GSE137303 - Phase 6B (D-038): sequence identity of ITAG4.0 genes against the tomato plastome (NC_007898.3) and
# mitogenome (NC_035963.1). blastn: every ITAG4.0 cDNA (one transcript per gene, reference/ITAG4.0_tx2gene.tsv) against
# both organelle genomes; blastp: proteins of the SL4.0ch00 genes against the organelle-encoded proteins (support only).
# BLAST+ is taken from the existing conda env rachis-qiime2-2026.4 (version recorded); nothing is installed.
# Usage (project root): bash scripts/bash/ch00_organelle_blast.sh
set -euo pipefail
cd "$(dirname "$0")/../.."
source "$HOME/miniconda3/etc/profile.d/conda.sh"
env=rachis-qiime2-2026.4
o=results/functional_robustness/ch00
org=reference/functional/organelle
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$o"
bl() { conda run -n "$env" "$@"; }
bl blastn -version | head -n 1 > "$o/blast_version.txt"
bl makeblastdb -version | head -n 1 >> "$o/blast_version.txt"

# ch00 genes from the GFF; their transcripts via tx2gene
awk -F'\t' '$1=="SL4.0ch00" && $3=="gene"{match($9,/ID=gene:[^;]+/); print substr($9,RSTART+8,RLENGTH-8)"\t"$4"\t"$5"\t"$7}' \
  reference/SGN/ITAG4.0/ITAG4.0_gene_models.gff | sort > "$o/ch00_genes_gff.tsv"
test "$(wc -l < "$o/ch00_genes_gff.tsv")" -eq 513 || { echo "unexpected ch00 gene count" >&2; exit 1; }
awk -F'\t' 'NR==FNR{g[$1]=1;next} FNR>1 && ($2 in g){print $1}' "$o/ch00_genes_gff.tsv" reference/ITAG4.0_tx2gene.tsv > "$tmp/ch00_tx.txt"
awk 'NR==FNR{k[$1]=1;next} /^>/{split(substr($1,2),a," "); keep=(a[1] in k)} keep' "$tmp/ch00_tx.txt" \
  reference/SGN/ITAG4.0/ITAG4.0_proteins.fasta > "$tmp/ch00_prot.faa"

cat "$org/NC_007898.3.fna" "$org/NC_035963.1.fna" > "$tmp/organelle_genomes.fna"
cat "$org/NC_007898.3.cds_aa.faa" "$org/NC_035963.1.cds_aa.faa" > "$tmp/organelle_prot.faa"
bl makeblastdb -in "$tmp/organelle_genomes.fna" -dbtype nucl -out "$tmp/orgnt" > /dev/null
bl makeblastdb -in "$tmp/organelle_prot.faa" -dbtype prot -out "$tmp/orgaa" > /dev/null

fmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs qcovhsp"
bl blastn -task megablast -query reference/SGN/ITAG4.0/ITAG4.0_cDNA.fasta -db "$tmp/orgnt" -evalue 1e-20 \
  -max_target_seqs 2 -max_hsps 1 -num_threads 4 -outfmt "$fmt" > "$tmp/blastn.tsv"
bl blastp -query "$tmp/ch00_prot.faa" -db "$tmp/orgaa" -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -num_threads 4 \
  -outfmt "$fmt slen" > "$tmp/blastp.tsv"
{ echo "$fmt" | tr ' ' '\t' | cut -f2-; cat "$tmp/blastn.tsv"; } > "$o/blastn_itag40_cdna_vs_organelle_genomes.tsv"
{ echo "$fmt slen" | tr ' ' '\t' | cut -f2-; cat "$tmp/blastp.tsv"; } > "$o/blastp_ch00_proteins_vs_organelle_proteins.tsv"
printf 'query_transcripts_ch00\t%s\nproteins_ch00\t%s\nblastn_hit_rows\t%s\nblastn_queries_with_hit\t%s\nblastp_hit_rows\t%s\n' \
  "$(wc -l < "$tmp/ch00_tx.txt")" "$(grep -c '>' "$tmp/ch00_prot.faa")" "$(wc -l < "$tmp/blastn.tsv")" \
  "$(cut -f1 "$tmp/blastn.tsv" | sort -u | wc -l)" "$(wc -l < "$tmp/blastp.tsv")" > "$o/blast_run_summary.tsv"
cat "$o/blast_version.txt" "$o/blast_run_summary.tsv"
