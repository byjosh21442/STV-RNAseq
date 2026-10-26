#!/usr/bin/env bash
# STV_GSE137303 - Phase 5 diagnostic: detect STV-derived reads in each raw library to check the STV-free /
# STV-infected labels. Competitive Salmon selective alignment against host ITAG4.0 cDNA + two STV genomes.
# Separate index and outputs; the host quantification (results/quantification/salmon) is not read or modified.
# Usage (project root, env stv-rnaseq): bash scripts/bash/stv_signal_check.sh
# Outputs: results/qc/stv_signal/{stv_signal_summary.tsv,stv_reference_provenance.tsv,stv_signal_method.txt,
#          stv_mappings/*.stv.sam}; results/qc/stv_signal/salmon/<sample>/ and the index are git-ignored
set -euo pipefail
THREADS="${THREADS:-16}"
host_cdna=reference/SGN/ITAG4.0/ITAG4.0_cDNA.fasta
stv_ids=(LC429302.1 KX949574.1)
stv_refs=(reference/STV/STV_LC429302.1.fa reference/STV/STV_KX949574.1.fa)
idx_dir=reference/STV_diagnostic/salmon_itag4_cdna_plus_stv_k31
out_dir=results/qc/stv_signal
mkdir -p "$out_dir/salmon" "$out_dir/stv_mappings" reference/STV_diagnostic
ids_csv=$(IFS=,; echo "${stv_ids[*]}")

printf 'accession\trole\tlength_nt\tsha256_fasta\tsource\n' > "$out_dir/stv_reference_provenance.tsv"
for i in "${!stv_ids[@]}"; do
  len=$(grep -v '^>' "${stv_refs[$i]}" | tr -d '\n' | wc -c)
  if [ "$i" -eq 0 ]; then role="primary: STV strain M82 deposited by the GSE137303 authors"
  else role="cross-check: STV isolate Florida used by the 2026 reanalysis"; fi
  printf '%s\t%s\t%s\t%s\t%s\n' "${stv_ids[$i]}" "$role" "$len" "$(sha256sum "${stv_refs[$i]}" | cut -d' ' -f1)" \
    "https://www.ncbi.nlm.nih.gov/nuccore/${stv_ids[$i]}" >> "$out_dir/stv_reference_provenance.tsv"
done

if [ ! -f "$idx_dir/.complete" ]; then
  cat "$host_cdna" "${stv_refs[@]}" > reference/STV_diagnostic/itag4_cdna_plus_stv.fa
  salmon index -t reference/STV_diagnostic/itag4_cdna_plus_stv.fa -i "$idx_dir" -k 31 -p "$THREADS" \
    > reference/STV_diagnostic/salmon_index.log 2>&1
  touch "$idx_dir/.complete"
fi

printf 'sample\tcondition\trun_accession\treads_processed\treads_mapped_any_target\tnumreads_LC429302.1\tnumreads_KX949574.1\tstv_numreads_total\tstv_reads_per_million_processed\tstv_percent_of_processed\tstv_primary_alignments_forward\tstv_primary_alignments_reverse\tLC429302.1_positions_covered\tLC429302.1_percent_covered\n' \
  > "$out_dir/stv_signal_summary.tsv"
len1=$(grep -v '^>' "${stv_refs[0]}" | tr -d '\n' | wc -c)
tail -n +2 config/samplesheet.csv | while IFS=, read -r sample cond _rep _gsm run _rest; do
  q="$out_dir/salmon/$sample"
  sam="$out_dir/stv_mappings/${sample}.stv.sam"
  # -l U counts fragments from both strands of the dsRNA genome; only STV mapping records are kept
  salmon quant -i "$idx_dir" -l U -r "data/raw/${run}.fastq.gz" -p "$THREADS" -o "$q" \
    --writeMappings >(awk -v ids="$ids_csv" 'BEGIN { n = split(ids, a, ","); for (i = 1; i <= n; i++) keep[a[i]] = 1 }
                         $1 == "@HD" || $1 == "@PG" || ($1 == "@SQ" && substr($2, 4) in keep) || ($3 in keep)' > "$sam") \
    > "$q.log" 2>&1
  sleep 2
  processed=$(grep -o '"num_processed": [0-9]*' "$q/aux_info/meta_info.json" | grep -o '[0-9]*$')
  mapped=$(grep -o '"num_mapped": [0-9]*' "$q/aux_info/meta_info.json" | grep -o '[0-9]*$')
  n1=$(awk -v id="${stv_ids[0]}" '$1 == id {print $5}' "$q/quant.sf")
  n2=$(awk -v id="${stv_ids[1]}" '$1 == id {print $5}' "$q/quant.sf")
  read -r fwd rev < <(awk '!/^@/ && int($2 / 256) % 2 == 0 { if (int($2 / 16) % 2) r++; else f++ } END { print f + 0, r + 0 }' "$sam")
  cov=$(awk -v id="${stv_ids[0]}" '!/^@/ && $3 == id {
          span = 0; cig = $6
          while (match(cig, /^[0-9]+[MIDNSHP=X]/)) { op = substr(cig, RLENGTH, 1); if (op ~ /[MDN=X]/) span += substr(cig, 1, RLENGTH - 1); cig = substr(cig, RLENGTH + 1) }
          for (p = $4; p < $4 + span; p++) c[p] = 1 }
        END { n = 0; for (p in c) n++; print n }' "$sam")
  awk -v s="$sample" -v c="$cond" -v r="$run" -v p="$processed" -v m="$mapped" -v a="$n1" -v b="$n2" \
      -v f="$fwd" -v v="$rev" -v cv="$cov" -v L="$len1" 'BEGIN { t = a + b; if (cv > L) cv = L
      printf "%s\t%s\t%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.4f\t%.6f\t%d\t%d\t%d\t%.2f\n", s, c, r, p, m, a, b, t, 1e6 * t / p, 100 * t / p, f, v, cv, 100 * cv / L }' \
    >> "$out_dir/stv_signal_summary.tsv"
done
{ echo "salmon $(salmon --version | awk '{print $2}')"
  echo "index: salmon index -t ITAG4.0_cDNA.fasta + ${stv_ids[*]} -k 31 (no decoys)"
  echo "quant: salmon quant -i $idx_dir -l U -r <fastq> -p $THREADS --writeMappings (defaults otherwise; STV records kept)"
  echo "stv_numreads_total = NumReads(LC429302.1) + NumReads(KX949574.1) from quant.sf (EM-assigned)"; } \
  > "$out_dir/stv_signal_method.txt"
column -t -s $'\t' "$out_dir/stv_signal_summary.tsv"
