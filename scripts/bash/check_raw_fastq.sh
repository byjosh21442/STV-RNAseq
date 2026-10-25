#!/usr/bin/env bash
# Read-only sanity checks of the downloaded raw FASTQ (gzip integrity, FASTQ structure, read length, read counts).
# Requires the stv-rnaseq env (pigz, seqkit). Usage (project root): bash scripts/bash/check_raw_fastq.sh
set -euo pipefail

SHEET="config/samplesheet.csv"
META="data/metadata/samples_metadata.tsv"
CHECKSUMS="results/qc/raw_fastq_checksums.tsv"
OUT="results/qc/raw_fastq_sanity.tsv"
N_RECORDS=10000      # records inspected for FASTQ structure
THREADS=4

# Only files whose MD5 already PASSed are checked.
mapfile -t ROWS < <(awk -F',' 'NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next } { print $h["sample"] "\t" $h["run_accession"] "\t" $h["fastq_1"] "\t" $h["read_length"] }' "$SHEET")

printf 'sample\trun\tfile\tgzip_integrity\tstructure_first_%s\theader_example\treads_seqkit\treads_ena\treads_match\tmin_len\tavg_len\tmax_len\tlength_ok\tgc_pct\tq30_pct\tstatus\n' "$N_RECORDS" > "$OUT"
failures=0

for r in "${ROWS[@]}"; do
  IFS=$'\t' read -r sample run fq expected_len <<< "$r"
  echo "== $sample ($fq)"
  if ! awk -F'\t' -v s="$sample" '$1==s && $7=="PASS" { ok=1 } END { exit !ok }' "$CHECKSUMS"; then
    echo "  MD5 not PASS -> skipped"; failures=$((failures + 1)); continue
  fi

  gz="FAIL"; pigz -t -p 2 "$fq" && gz="OK"

  # '|| true': zcat is killed by SIGPIPE once head has enough lines.
  structure=$( { zcat "$fq" || true; } | head -n $(( N_RECORDS * 4 )) | awk -v n="$N_RECORDS" '
    NR%4==1 && substr($0,1,1)!="@" { bad++ }
    NR%4==2 { seq=$0; if (seq !~ /^[ACGTN]+$/) bad++ }
    NR%4==3 && substr($0,1,1)!="+" { bad++ }
    NR%4==0 { if (length($0)!=length(seq)) bad++; recs++ }
    END { print (recs==n && bad==0) ? "OK" : "FAIL(records=" recs ",bad=" bad+0 ")" }')
  header=$( { zcat "$fq" || true; } | head -n 1 | cut -c1-60)

  stats=$(seqkit stats -a -T -j "$THREADS" "$fq" | awk -F'\t' 'NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    { print $h["num_seqs"] "\t" $h["min_len"] "\t" $h["avg_len"] "\t" $h["max_len"] "\t" $h["GC(%)"] "\t" $h["Q30(%)"] }')
  IFS=$'\t' read -r nseq minl avgl maxl gc q30 <<< "$stats"
  ena=$(awk -F'\t' -v s="$sample" 'NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next } $h["sample_id"]==s { print $h["read_count"] }' "$META")

  reads_ok="NO"; [ "$nseq" = "$ena" ] && reads_ok="YES"
  len_ok="NO"; [ "$minl" = "$expected_len" ] && [ "$maxl" = "$expected_len" ] && len_ok="YES"
  status="PASS"
  { [ "$gz" = "OK" ] && [ "$structure" = "OK" ] && [ "$reads_ok" = "YES" ] && [ "$len_ok" = "YES" ]; } || { status="FAIL"; failures=$((failures + 1)); }

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$run" "$(basename "$fq")" "$gz" "$structure" \
    "$header" "$nseq" "$ena" "$reads_ok" "$minl" "$avgl" "$maxl" "$len_ok" "$gc" "$q30" "$status" >> "$OUT"
  echo "  gzip=$gz structure=$structure reads=$nseq (ENA $ena) len=$minl-$maxl -> $status"
done

column -t -s $'\t' "$OUT"
[ "$failures" -eq 0 ] || { echo "ERROR: $failures sample(s) failed sanity checks" >&2; exit 1; }
