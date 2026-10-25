#!/usr/bin/env bash
# Regenerate the tiny synthetic fixtures used by the -stub wiring test (no biological data).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p reads
seq75() { printf 'ACGT%.0s' $(seq 1 19) | cut -c1-75; }
qual75() { printf 'I%.0s' $(seq 1 75); }
for s in A B; do
  { for i in 1 2 3 4; do printf '@stub_%s_%d\n%s\n+\n%s\n' "$s" "$i" "$(seq75)" "$(qual75)"; done; } | gzip -n > "reads/stub_${s}.fastq.gz"
done
printf '>tx1\n%s\n>tx2\n%s\n' "$(seq75)" "$(seq75 | rev)" > transcripts.fa
printf '>chr_stub\n%s%s\n' "$(seq75)" "$(seq75)" > genome.fa
ls -l reads transcripts.fa genome.fa
