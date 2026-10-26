#!/usr/bin/env bash
# STV_GSE137303 - Phase 6: download the functional-annotation resources used by scripts/R/05_functional_ranking_annotation.R
# into reference/functional/ (git-ignored, re-downloadable), verify size against the server and write a tracked manifest
# with SHA256 and server Last-Modified dates.
#   - PLAZA 5.0 dicots GO annotation for tomato (ITAG4.0 gene IDs)          -> primary GO source (D-024)
#   - Gene Ontology go-basic.obo, dated release 2026-08-05                     -> ontology for propagation (D-025)
#   - Plant Reactome gene_ids_by_pathway_and_species.tab (current release)    -> tomato pathway audit (D-028)
#   - SGN ITAG4.1 GO terms, proteins and README                               -> audit of the SGN GO alternative (D-024)
# Usage (project root): bash scripts/bash/download_functional_annotation.sh
set -euo pipefail

OUT="reference/functional"
MANIFEST="results/functional/annotation_sources_manifest.tsv"
TODAY="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
mkdir -p "$OUT" "$(dirname "$MANIFEST")"

# subdir | URL | resource version | purpose
FILES=(
  "plaza_dicots_05|https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_05/GO/go.sly.csv.gz|PLAZA 5.0 dicots (instance dicots_05); tomato = ITAG4.0|GO annotation of ITAG4.0 genes (primary GO source)"
  "go/release_2026-08-05|https://release.geneontology.org/2026-08-05/ontology/go-basic.obo|GO release 2026-08-05 (go-basic)|Ontology for alt_id/obsolete handling and is_a/part_of propagation"
  "plant_reactome/current|https://plantreactome.gramene.org/download/current/gene_ids_by_pathway_and_species.tab|Plant Reactome current release (site: Version 24)|Tomato pathway-membership audit"
  "sgn_itag4.1|https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_goterms.txt|ITAG4.1 (SL4.0)|Audit of the SGN GO alternative (not used for gene sets)"
  "sgn_itag4.1|https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_proteins.fasta|ITAG4.1 (SL4.0)|Protein identity check for the ITAG4.1 GO audit"
  "sgn_itag4.1|https://solgenomics.net/ftp/tomato_genome/annotation/ITAG4.1_release/ITAG4.1_README|ITAG4.1 (SL4.0)|Release provenance"
)

printf 'file\tsource_url\tresource_version\tsize_bytes\tsha256\tupstream_last_modified\tverified_utc\tpurpose\n' > "$MANIFEST.tmp"

for entry in "${FILES[@]}"; do
  IFS='|' read -r sub url version purpose <<< "$entry"
  dest="$OUT/$sub/$(basename "$url")"
  mkdir -p "$(dirname "$dest")"

  head=$(curl -fsSIL --max-time 60 "$url" | tr -d '\r')
  remote_size=$(awk -F': ' 'tolower($1)=="content-length" { v=$2 } END { print v }' <<< "$head")
  remote_date=$(awk -F': ' 'tolower($1)=="last-modified" { v=$2 } END { print v }' <<< "$head")
  [ -n "$remote_size" ] || { echo "ERROR: no Content-Length for $url" >&2; exit 1; }

  if [ "$(stat -c %s "$dest" 2>/dev/null || echo 0)" != "$remote_size" ]; then
    echo "== downloading $url ($remote_size bytes)"
    rm -f "$dest.part"
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
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$dest" "$url" "$version" "$remote_size" "$sha" "$remote_date" "$TODAY" \
    "$purpose" >> "$MANIFEST.tmp"
done

# Provenance lines written by the providers inside the files
echo "== GO data-version: $(grep -m1 '^data-version:' "$OUT/go/release_2026-08-05/go-basic.obo")"
echo "== PLAZA header:"; zcat "$OUT/plaza_dicots_05/go.sly.csv.gz" 2>/dev/null | head -n 8 || true

mv "$MANIFEST.tmp" "$MANIFEST"
column -t -s $'\t' "$MANIFEST" | cut -c1-200
du -sh "$OUT"
