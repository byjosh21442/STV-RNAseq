#!/usr/bin/env bash
# STV_GSE137303 - Phase 5 diagnostic: count reads per Illumina instrument:run:flowcell:lane from the original read
# names retained in the ENA FASTQ headers ("@SRRxxxx.n INSTRUMENT:RUN:FLOWCELL:LANE:TILE:X:Y/1").
# Usage (project root, env stv-rnaseq): bash scripts/bash/fastq_header_audit.sh
# Output: results/qc/sample_structure/fastq_header_lanes.tsv
set -euo pipefail
out_dir=results/qc/sample_structure
mkdir -p "$out_dir"
out="$out_dir/fastq_header_lanes.tsv"
printf 'sample\trun_accession\tinstrument\trun_number\tflowcell\tlane\treads\ttiles\n' > "$out"
tail -n +2 config/samplesheet.csv | while IFS=, read -r sample _cond _rep _gsm run _rest; do
  fq="data/raw/${run}.fastq.gz"
  pigz -dc "$fq" | awk -v s="$sample" -v r="$run" 'NR % 4 == 1 {
      split($2, a, ":"); k = a[1] "\t" a[2] "\t" a[3] "\t" a[4]; n[k]++; t[k "\t" a[5]] = 1 }
    END { for (k in n) { c = 0; for (x in t) if (index(x, k "\t") == 1) c++; print s "\t" r "\t" k "\t" n[k] "\t" c } }' \
    | sort >> "$out"
done
column -t -s $'\t' "$out"
