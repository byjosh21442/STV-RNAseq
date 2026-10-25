#!/usr/bin/env bash
# Download ONLY the small GEO supplementary processed file for GSE137303 (~1.3 MB) and print a structural summary.
# No pipefail: `zcat | head` would otherwise abort on SIGPIPE.
set -eu

OUT="data/processed/original_geo"
F="$OUT/GSE137303_count.txt.gz"
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137303/suppl/GSE137303_count.txt.gz"
mkdir -p "$OUT"

[ -s "$F" ] || curl -fsSL "$URL" -o "$F"
md5sum "$F" | tee "$F.md5"
ls -l "$F"

echo "--- header (tab-separated fields, one per line)"
zcat "$F" | head -n 1 | tr '\t' '\n' | cat -n
echo "--- first 5 rows"
zcat "$F" | head -n 6
echo "--- dimensions (rows incl. header, fields per row)"
zcat "$F" | awk -F'\t' '{print NF}' | sort | uniq -c
zcat "$F" | wc -l
echo "--- per-column value classes in data rows (int / decimal / other)"
zcat "$F" | awk -F'\t' 'NR>1 { for (i=2;i<=NF;i++) {
    if ($i ~ /^-?[0-9]+$/) c[i,"int"]++; else if ($i ~ /^-?[0-9]*\.[0-9]+([eE][-+]?[0-9]+)?$/) c[i,"dec"]++; else c[i,"other"]++ } ; n=NF }
  END { for (i=2;i<=n;i++) printf "col%d int=%d dec=%d other=%d\n", i, c[i,"int"], c[i,"dec"], c[i,"other"] }'
echo "--- identifier patterns in column 1"
zcat "$F" | awk -F'\t' 'NR>1 {id=$1; gsub(/[0-9]/,"9",id); print id}' | sort | uniq -c | sort -rn | head -n 15
echo "--- duplicated identifiers"
zcat "$F" | awk -F'\t' 'NR>1 {print $1}' | sort | uniq -d | wc -l

echo "--- count column totals (cols 17-22) and all-zero rows"
zcat "$F" | awk -F'\t' 'NR==1 {for(i=17;i<=22;i++) h[i]=$i; next}
  {z=1; for(i=17;i<=22;i++){s[i]+=$i; if($i!=0) z=0}; zr+=z}
  END {for(i=17;i<=22;i++) printf "%s\t%d\n", h[i], s[i]; print "all_zero_rows\t" zr}'
echo "--- FPKM column totals (odd cols 5-15)"
zcat "$F" | awk -F'\t' 'NR==1 {for(i=5;i<=15;i+=2) h[i]=$i; next} {for(i=5;i<=15;i+=2) s[i]+=$i}
  END {for(i=5;i<=15;i+=2) printf "%s\t%.1f\n", h[i], s[i]}'
echo "--- column 2 prefix / column 3 name classes"
zcat "$F" | awk -F'\t' 'NR>1 {split($2,a,":"); p[a[1]]++; if ($3 ~ /^LOC[0-9]+$/) n["LOC"]++; else if ($3=="") n["empty"]++; else n["symbol/other"]++}
  END {for(k in p) print "col2_prefix", k, p[k]; for(k in n) print "col3", k, n[k]}'
echo "--- rows mentioning Solyc / ITAG anywhere"
zcat "$F" | grep -ciE 'solyc|itag' || true
echo "--- sample of non-LOC names"
zcat "$F" | awk -F'\t' 'NR>1 && $3 !~ /^LOC/ {print $1"\t"$2"\t"$3"\t"$4}' | head -n 8
echo "--- min/max of Processed Signal columns (even cols 6-16)"
zcat "$F" | awk -F'\t' 'NR==1{next} NR==2{for(i=6;i<=16;i+=2){mn[i]=$i;mx[i]=$i}} {for(i=6;i<=16;i+=2){if($i+0<mn[i]+0)mn[i]=$i; if($i+0>mx[i]+0)mx[i]=$i}}
  END{for(i=6;i<=16;i+=2) printf "col%d min=%s max=%s\n", i, mn[i], mx[i]}'
