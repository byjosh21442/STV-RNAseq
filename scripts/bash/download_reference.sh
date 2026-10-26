#!/usr/bin/env bash
# Download the approved SL4.0 + ITAG4.0 reference files from Sol Genomics Network (D-001) into reference/SGN/,
# keep the official file names, verify size against the server, and write reference/reference_manifest.tsv (SHA256).
# Usage (project root): bash scripts/bash/download_reference.sh
set -euo pipefail

BASE="https://solgenomics.net/ftp/tomato_genome"
OUT="reference/SGN"
MANIFEST="reference/reference_manifest.tsv"
TODAY="$(date -u +%Y-%m-%d)"

# subdir | path under $BASE | genome_version | annotation_version | purpose
FILES=(
  "SL4.0|assembly/build_4.00/S_lycopersicum_chromosomes.4.00.fa.gz|SL4.0|NA|Genome FASTA: Salmon decoy sequences"
  "SL4.0|assembly/build_4.00/FileDescription.txt|SL4.0|NA|Release provenance (assembly)"
  "SL4.0|assembly/build_4.00/ReleaseNotes.txt|SL4.0|NA|Release provenance (assembly)"
  "ITAG4.0|annotation/ITAG4.0_release/ITAG4.0_README.txt|SL4.0|ITAG4.0|Release provenance (annotation)"
  "ITAG4.0|annotation/ITAG4.0_release/ITAG4.0_cDNA.fasta|SL4.0|ITAG4.0|Transcript FASTA: Salmon index targets"
  "ITAG4.0|annotation/ITAG4.0_release/ITAG4.0_gene_models.gff|SL4.0|ITAG4.0|Gene models: tx2gene and validation"
  "ITAG4.0|annotation/ITAG4.0_release/ITAG4.0_proteins.fasta|SL4.0|ITAG4.0|Proteins: later functional annotation"
  "ITAG4.0|annotation/ITAG4.0_release/ITAG4.0_descriptions.txt|SL4.0|ITAG4.0|Gene descriptions: later annotation"
  "ITAG4.0|annotation/ITAG4.0_release/ITAG4.0_goterms.txt|SL4.0|ITAG4.0|GO terms: later functional analysis"
)

printf 'file\tsource_url\tgenome_version\tannotation_version\tsize_bytes\tsha256\tupstream_checksum\tupstream_last_modified\tretrieval_date\tpurpose\n' > "$MANIFEST.tmp"

for entry in "${FILES[@]}"; do
  IFS='|' read -r sub path gver aver purpose <<< "$entry"
  url="$BASE/$path"
  dest="$OUT/$sub/$(basename "$path")"
  mkdir -p "$(dirname "$dest")"

  head=$(curl -fsSIL --max-time 60 "$url" | tr -d '\r')
  remote_size=$(awk -F': ' 'tolower($1)=="content-length" { v=$2 } END { print v }' <<< "$head")
  remote_date=$(awk -F': ' 'tolower($1)=="last-modified" { v=$2 } END { print v }' <<< "$head")
  [ -n "$remote_size" ] || { echo "ERROR: no Content-Length for $url" >&2; exit 1; }

  if [ "$(stat -c %s "$dest" 2>/dev/null || echo 0)" != "$remote_size" ]; then
    echo "== downloading $path ($remote_size bytes)"
    for i in 1 2 3 4 5; do
      curl --fail --location --silent --show-error --connect-timeout 30 -C - -o "$dest.part" "$url" && break
      echo "  curl error (attempt $i/5); retrying in 15 s" >&2; sleep 15
    done
    [ "$(stat -c %s "$dest.part")" = "$remote_size" ] || { echo "ERROR: size mismatch for $dest.part" >&2; exit 1; }
    mv "$dest.part" "$dest"
  else
    echo "== present with server size: $dest"
  fi
  case "$dest" in *.gz) gzip -t "$dest" || { echo "ERROR: gzip test failed for $dest" >&2; exit 1; } ;; esac

  sha=$(sha256sum "$dest" | cut -d' ' -f1)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$dest" "$url" "$gver" "$aver" "$remote_size" "$sha" \
    "none published by SGN" "$remote_date" "$TODAY" "$purpose" >> "$MANIFEST.tmp"
done

mv "$MANIFEST.tmp" "$MANIFEST"
column -t -s $'\t' "$MANIFEST" | cut -c1-220
du -sh "$OUT"
