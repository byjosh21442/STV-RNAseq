#!/usr/bin/env bash
# Build data/metadata/samples_metadata.tsv and config/samplesheet.csv ONLY from the raw records
# saved by fetch_metadata.sh (GEO SOFT + ENA read_run + NCBI SRA runinfo). No values are typed by hand.
set -euo pipefail

SRC="data/metadata/sources"
SOFT="$SRC/GSE137303_brief.soft.txt"
ENA="$SRC/PRJNA565137_ena_read_run.tsv"
RUNINFO="$SRC/PRJNA565137_sra_runinfo.csv"
OUT_TSV="data/metadata/samples_metadata.tsv"
OUT_CSV="config/samplesheet.csv"
mkdir -p config

# 1) GEO sample-level fields, keyed by GSM
geo_tmp="$(mktemp)"
awk -F' = ' '
  { sub(/\r$/, "") }  # GEO SOFT is served with CRLF line endings
  /^\^SAMPLE = / { gsm=$2; order[++n]=gsm; next }
  gsm=="" { next }
  /^!Sample_title = /                 { title[gsm]=$2 }
  /^!Sample_description = /           { desc[gsm]=$2 }
  /^!Sample_molecule_ch1 = /          { mol[gsm]=$2 }
  /^!Sample_characteristics_ch1 = /   { split($2,kv,": "); ch[gsm,kv[1]]=kv[2] }
  /^!Sample_extract_protocol_ch1 = .*TruSeq/ { kit[gsm]="TruSeq stranded mRNA sample prep kit (as reported)" }
  /^!Sample_relation = BioSample: /   { sub(/.*\//,"",$2); bs[gsm]=$2 }
  /^!Sample_relation = SRA: /         { sub(/.*term=/,"",$2); srx[gsm]=$2 }
  /^!Sample_data_processing = Genome_build: / { sub(/^Genome_build: /,"",$2); gb[gsm]=$2 }
  END { for (i=1;i<=n;i++) { g=order[i];
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", g, title[g], desc[g], ch[g,"virus infection"], ch[g,"cultivar"], ch[g,"tissue"], ch[g,"age"], mol[g], kit[g], bs[g], srx[g], gb[g] } }
' "$SOFT" > "$geo_tmp"

# 2) NCBI runinfo size_MB keyed by run
run_tmp="$(mktemp)"
awk -F',' 'NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next } { print $h["Run"] "\t" $h["size_MB"] "\t" $h["avgLength"] }' "$RUNINFO" > "$run_tmp"

# 3) Join on GSM (ENA sample_alias == GEO GSM) and write the metadata table
awk -F'\t' -v OFS='\t' '
  FNR==NR && FILENAME==ARGV[1] { geo[$1]=$0; next }
  FILENAME==ARGV[2] { size[$1]=$2; avglen[$1]=$3; next }
  FNR==1 { for(i=1;i<=NF;i++) h[$i]=i
    print "sample_id","condition","replicate","geo_gsm","geo_title","geo_description","biosample","sra_sample","sra_experiment","sra_run",
          "sra_study","bioproject","organism","tax_id","cultivar","tissue","age","molecule_geo","library_strategy","library_source",
          "library_selection","library_layout","instrument_model","read_count","base_count","avg_read_length","ena_fastq_bytes",
          "ena_fastq_md5","ena_fastq_ftp","ncbi_sra_size_mb","library_kit_reported","strandedness","genome_build_reported","sources"
    next }
  { gsm=$h["sample_alias"]; split(geo[gsm], g, "\t")
    cond=g[4]; gsub(/-/,"_",cond)
    rep=g[2]; sub(/.*rep/,"",rep)
    sid=cond "_rep" rep
    print sid, cond, rep, gsm, g[2], g[3], g[10], $h["secondary_sample_accession"], $h["experiment_accession"], $h["run_accession"],
          $h["secondary_study_accession"], $h["study_accession"], $h["scientific_name"], $h["tax_id"], g[5], g[6], g[7], g[8],
          $h["library_strategy"], $h["library_source"], $h["library_selection"], $h["library_layout"], $h["instrument_model"],
          $h["read_count"], $h["base_count"], avglen[$h["run_accession"]], $h["fastq_bytes"], $h["fastq_md5"], $h["fastq_ftp"],
          size[$h["run_accession"]], g[9], "UNKNOWN_TO_BE_VERIFIED", g[12],
          "GEO:" gsm ";ENA:" $h["run_accession"] ";NCBI_SRA:" $h["run_accession"] }
' "$geo_tmp" "$run_tmp" "$ENA" | { IFS= read -r header; echo "$header"; sort -t$'\t' -k2,2 -k3,3n; } > "$OUT_TSV"

# Consistency checks against the independent NCBI and GEO records
awk -F'\t' 'NR>1 && $11 != "" { n++ } END { if (n!=6) { print "ERROR: expected 6 joined rows, got " n > "/dev/stderr"; exit 1 } }' "$OUT_TSV"
awk -F'\t' 'NR>1 { print $10 "\t" $9 "\t" $7 }' "$OUT_TSV" | sort > "$geo_tmp.chk1"
awk -F',' 'NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next } { print $h["Run"] "\t" $h["Experiment"] "\t" $h["BioSample"] }' "$RUNINFO" | sort > "$geo_tmp.chk2"
diff "$geo_tmp.chk1" "$geo_tmp.chk2" || { echo "ERROR: ENA vs NCBI runinfo mismatch" >&2; exit 1; }
echo "OK: ENA run/experiment/biosample == NCBI runinfo"
awk -F'\t' 'NR>1 { print $4 "\t" $9 "\t" $7 }' "$OUT_TSV" | sort > "$geo_tmp.chk3"
cut -f1,11,10 "$geo_tmp" | awk -F'\t' -v OFS='\t' '{ print $1, $3, $2 }' | sort > "$geo_tmp.chk4"
diff "$geo_tmp.chk3" "$geo_tmp.chk4" || { echo "ERROR: ENA vs GEO SOFT mismatch" >&2; exit 1; }
echo "OK: GSM/SRX/BioSample == GEO SOFT"

# 4) Samplesheet: accession fields + fastq_1 = local ENA file name under data/raw/ (written by download_fastq.sh)
awk -F'\t' -v OFS=',' 'NR==1 { for(i=1;i<=NF;i++) h[$i]=i
    print "sample","condition","replicate","geo_gsm","run_accession","library_layout","read_length","strandedness","ena_fastq_ftp","ena_fastq_md5","fastq_1"; next }
  { n = split($h["ena_fastq_ftp"], parts, "/")
    print $h["sample_id"], $h["condition"], $h["replicate"], $h["geo_gsm"], $h["sra_run"], $h["library_layout"], $h["avg_read_length"],
          "auto", $h["ena_fastq_ftp"], $h["ena_fastq_md5"], "data/raw/" parts[n] }' "$OUT_TSV" > "$OUT_CSV"

rm -f "$geo_tmp" "$geo_tmp".chk* "$run_tmp"
column -t -s$'\t' < <(cut -f1-4,10,22,24,26,27 "$OUT_TSV")
echo; cat "$OUT_CSV"
