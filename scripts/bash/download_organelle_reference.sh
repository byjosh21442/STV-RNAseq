#!/usr/bin/env bash
# STV_GSE137303 - Phase 6B (D-038): download the tomato organelle RefSeq genomes and their encoded proteins from NCBI
# E-utilities for the SL4.0ch00 investigation. Files are git-ignored (reference/*); the manifest with SHA-256 is tracked.
# Usage (project root): bash scripts/bash/download_organelle_reference.sh
set -euo pipefail
cd "$(dirname "$0")/../.."
d=reference/functional/organelle
m=results/functional_robustness/ch00/organelle_reference_manifest.tsv
mkdir -p "$d" "$(dirname "$m")"
base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&retmode=text"
printf 'file\taccession\tcontent\tsource_url\tsize_bytes\tsha256\tfirst_header\tretrieval_utc\n' > "$m"
for acc in NC_007898.3 NC_035963.1; do
  for kind in fasta fasta_cds_aa; do
    ext=$([ "$kind" = fasta ] && echo fna || echo cds_aa.faa)
    f="$d/${acc}.${ext}"
    url="${base}&id=${acc}&rettype=${kind}"
    curl -fsS --retry 3 "$url" -o "$f"
    test -s "$f" && head -c 1 "$f" | grep -q '>' || { echo "download failed: $url" >&2; exit 1; }
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$f" "$acc" "$kind" "$url" "$(stat -c %s "$f")" \
      "$(sha256sum "$f" | cut -d' ' -f1)" "$(head -n 1 "$f" | tr '\t' ' ')" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$m"
    sleep 1   # NCBI E-utilities: max 3 requests/s without an API key
  done
done
cat "$m"
grep -c '>' "$d"/*.faa
