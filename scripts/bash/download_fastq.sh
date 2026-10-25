#!/usr/bin/env bash
# Download the six ENA FASTQ files listed in data/metadata/samples_metadata.tsv into data/raw/,
# resuming partial downloads, and accept a file only if size AND MD5 match the ENA-published values.
# Usage (project root): bash scripts/bash/download_fastq.sh
set -euo pipefail

META="data/metadata/samples_metadata.tsv"
OUT="data/raw"
LOG="$OUT/download_log.tsv"
SUMMARY="results/qc/raw_fastq_checksums.tsv"
MAX_CYCLES=2          # full download+verify cycles per file (a failed file is quarantined, then re-fetched once)
MAX_RESUMES=5         # curl invocations per cycle, each resuming the .part file
HEADROOM=$((5 * 1024 * 1024 * 1024))

mkdir -p "$OUT" "$(dirname "$SUMMARY")"

mapfile -t ROWS < <(awk -F'\t' 'NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
  { print $h["sample_id"] "\t" $h["sra_run"] "\t" $h["ena_fastq_ftp"] "\t" $h["ena_fastq_md5"] "\t" $h["ena_fastq_bytes"] }' "$META")
[ "${#ROWS[@]}" -eq 6 ] || { echo "ERROR: expected 6 rows in $META, found ${#ROWS[@]}" >&2; exit 1; }

size_of() { [ -f "$1" ] && stat -c %s "$1" || echo 0; }

# --- disk space check (only bytes still missing are counted) ---
need=0
for r in "${ROWS[@]}"; do
  IFS=$'\t' read -r _ _ path _ bytes <<< "$r"
  f="$OUT/$(basename "$path")"
  [ -s "$f" ] && continue
  need=$(( need + bytes - $(size_of "$f.part") ))
done
avail=$(df --output=avail -B1 "$OUT" | tail -n 1)
echo "Bytes still to download: $need; free on $(df --output=target "$OUT" | tail -n 1): $avail; required headroom: $HEADROOM"
[ "$avail" -gt $(( need + HEADROOM )) ] || { echo "ERROR: insufficient free disk space" >&2; exit 1; }

fetch() {  # url part_file expected_bytes
  local url=$1 part=$2 bytes=$3 i
  for i in $(seq 1 "$MAX_RESUMES"); do
    [ "$(size_of "$part")" -eq "$bytes" ] && return 0
    if [ "$(size_of "$part")" -gt "$bytes" ]; then echo "  .part larger than expected; discarding" >&2; rm -f "$part"; fi
    curl --fail --location --silent --show-error --connect-timeout 30 \
         --speed-limit 10240 --speed-time 120 -C - -o "$part" "$url" && continue
    echo "  curl error (invocation $i/$MAX_RESUMES); retrying in 15 s" >&2
    sleep 15
  done
  [ "$(size_of "$part")" -eq "$bytes" ]
}

printf 'timestamp\tsample\trun\tfile\tbytes_expected\tbytes_observed\tmd5_expected\tmd5_observed\tstatus\tcycles\n' > "$LOG"
printf 'sample\trun\tfastq_file\tbytes\tmd5_expected\tmd5_observed\tstatus\n' > "$SUMMARY"
failures=0

for r in "${ROWS[@]}"; do
  IFS=$'\t' read -r sample run path md5 bytes <<< "$r"
  url="https://${path}"
  file="$OUT/$(basename "$path")"
  status="FAIL"; observed="NA"; cycles=0
  echo "== $sample $run -> $file"

  remote=$(curl --fail --location --silent --head "$url" | awk 'tolower($1)=="content-length:" { v=$2 } END { gsub(/\r/,"",v); print v }') || remote=""
  if [ -n "$remote" ] && [ "$remote" != "$bytes" ]; then
    echo "  ERROR: remote size $remote != ENA metadata $bytes; skipping" >&2
  else
    for cycle in $(seq 1 "$MAX_CYCLES"); do
      cycles=$cycle
      if [ ! -s "$file" ]; then
        fetch "$url" "$file.part" "$bytes" || { echo "  download incomplete (cycle $cycle)" >&2; continue; }
        mv "$file.part" "$file"
      fi
      observed=$(md5sum "$file" | cut -d' ' -f1)
      if [ "$(size_of "$file")" -eq "$bytes" ] && [ "$observed" = "$md5" ]; then status="PASS"; break; fi
      echo "  MD5/size mismatch (cycle $cycle): quarantined as $file.corrupt.$cycle" >&2
      mv "$file" "$file.corrupt.$cycle"
    done
  fi

  [ "$status" = "PASS" ] || failures=$(( failures + 1 ))
  echo "  $status md5=$observed"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$sample" "$run" "$file" \
    "$bytes" "$(size_of "$file")" "$md5" "$observed" "$status" "$cycles" >> "$LOG"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$run" "$(basename "$file")" "$(size_of "$file")" "$md5" "$observed" "$status" >> "$SUMMARY"
done

column -t -s $'\t' "$SUMMARY"
df -h "$OUT" | tail -n 1
[ "$failures" -eq 0 ] || { echo "ERROR: $failures file(s) failed validation; do NOT start QC for them" >&2; exit 1; }
echo "All files PASS."
