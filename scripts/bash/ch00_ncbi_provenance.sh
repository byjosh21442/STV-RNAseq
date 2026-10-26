#!/usr/bin/env bash
# Post-hoc provenance check (Phase 6B closure, organelle assessment; D-044).
# Compares the SGN SL4.0 ch00 segments (split at runs of >=100 N) with the
# unplaced scaffolds of the NCBI GenBank SL4.0 assembly (GCA_000188115.5) by
# exact sequence identity (forward or reverse complement).
# Network: NCBI Datasets v2 sequence_reports + E-utilities efetch (FASTA kept in a temp dir).
# Outputs (results/organelle_assessment/provenance/):
#   ncbi_GCA_000188115.5_sequence_report.tsv   name, accession, role, length, retrieval date
#   ch00_segment_vs_ncbi_genbank.tsv           per SGN ch00 segment: exact NCBI match or none
#   ch00_vs_ncbi_summary.tsv                   counts / bp by organelle-identity bin
set -euo pipefail
cd "$(dirname "$0")/../.."
OUT=results/organelle_assessment/provenance
SEG=results/organelle_assessment/ch00_segment_organelle_coverage.tsv
GENOME=reference/SGN/SL4.0/S_lycopersicum_chromosomes.4.00.fa.gz
ACC=GCA_000188115.5
mkdir -p "$OUT"
[[ -s "$SEG" ]] || { echo "missing $SEG (run scripts/R/19_organelle_scope_diagnostic.R first)" >&2; exit 1; }
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

curl -fsS "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${ACC}/sequence_reports?page_size=1000" -o "$tmp/report.json"
python3 - "$tmp/report.json" "$OUT/ncbi_${ACC}_sequence_report.tsv" "$tmp/unplaced_ids.txt" <<'EOF'
import json, sys, datetime
rep = json.load(open(sys.argv[1]))["reports"]
day = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d")
with open(sys.argv[2], "w") as o, open(sys.argv[3], "w") as ids:
    o.write("sequence_name\tgenbank_accession\trole\tchr_name\tlength\tretrieved_utc\n")
    for r in rep:
        o.write(f"{r.get('sequence_name','')}\t{r.get('genbank_accession','')}\t{r.get('role','')}\t{r.get('chr_name','')}\t{r.get('length','')}\t{day}\n")
        if r.get("role") == "unplaced-scaffold":
            ids.write(r["genbank_accession"] + "\n")
EOF
curl -fsS -X POST "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
  --data "db=nuccore&rettype=fasta&retmode=text&id=$(paste -sd, "$tmp/unplaced_ids.txt")" -o "$tmp/unplaced.fa"
zcat "$GENOME" | awk '/^>/{p=($1==">SL4.0ch00"); next} p{printf "%s",$0} END{print ""}' > "$tmp/ch00.seq"

python3 - "$tmp" "$SEG" "$OUT" <<'EOF'
import sys, re, csv, hashlib
tmp, segf, out = sys.argv[1:4]
comp = str.maketrans("ACGTN", "TGCAN")
def h(x): return hashlib.sha1(x.encode()).hexdigest()
s = open(f"{tmp}/ch00.seq").read().strip().upper()
pieces, pos = [], 0
for part in re.split(r"(N{100,})", s):
    if part and part[0] != "N":
        pieces.append((pos + 1, part))
    pos += len(part)
# NCBI unplaced scaffolds, split the same way
ncbi, name, buf = {}, None, []
for line in open(f"{tmp}/unplaced.fa"):
    line = line.strip()
    if line.startswith(">"):
        if name: ncbi[name] = "".join(buf).upper()
        name, buf = line[1:].split()[0], []
    elif line:
        buf.append(line)
if name: ncbi[name] = "".join(buf).upper()
idx = {}
for acc, seq in ncbi.items():
    k = 0
    for part in re.split(r"(N{100,})", seq):
        if part and part[0] != "N":
            k += 1
            idx.setdefault(h(part), []).append(f"{acc}:piece{k}:+")
            idx.setdefault(h(part[::-1].translate(comp)), []).append(f"{acc}:piece{k}:-")
segs = list(csv.DictReader(open(segf), delimiter="\t"))
by_start = {int(r["start"]): r for r in segs}
rows = []
for st, part in pieces:
    r = by_start[st]
    hit = idx.get(h(part), [])
    p = float(r["plastid_identical_fraction_of_segment"]); m = float(r["mito_identical_fraction_of_segment"])
    org = "plastid>=90%" if p >= 0.9 else "mito>=90%" if m >= 0.9 else "organelle<10%" if p + m < 0.1 else "intermediate"
    rows.append([r["segment"], st, len(part), f"{p:.4f}", f"{m:.4f}", org, r["genes"], r["genes_tested"],
                 "yes" if hit else "no", ";".join(hit)])
with open(f"{out}/ch00_segment_vs_ncbi_genbank.tsv", "w") as o:
    o.write("segment\tstart\tlength\tplastid_identical_fraction\tmito_identical_fraction\torganelle_bin\tgenes\tgenes_tested\tin_ncbi_genbank_exact\tncbi_match\n")
    for x in rows: o.write("\t".join(map(str, x)) + "\n")
matched_pieces = {m.rsplit(":", 1)[0] for x in rows for m in x[9].split(";") if m}
total_pieces = sum(1 for acc, seq in ncbi.items() for part in re.split(r"(N{100,})", seq) if part and part[0] != "N")
with open(f"{out}/ch00_vs_ncbi_summary.tsv", "w") as o:
    o.write("organelle_bin\tsegments\tbp\tgenes\tgenes_tested\tsegments_in_ncbi\tbp_in_ncbi\tsegments_absent_from_ncbi\tbp_absent_from_ncbi\tgenes_on_absent\tgenes_tested_on_absent\n")
    for b in ["plastid>=90%", "mito>=90%", "intermediate", "organelle<10%", "ALL"]:
        sel = [x for x in rows if b == "ALL" or x[5] == b]
        yes = [x for x in sel if x[8] == "yes"]; no = [x for x in sel if x[8] == "no"]
        o.write("\t".join(map(str, [b, len(sel), sum(x[2] for x in sel), sum(int(x[6]) for x in sel), sum(int(x[7]) for x in sel),
               len(yes), sum(x[2] for x in yes), len(no), sum(x[2] for x in no), sum(int(x[6]) for x in no), sum(int(x[7]) for x in no)])) + "\n")
    o.write(f"#ncbi_unplaced_scaffolds\t{len(ncbi)}\t{sum(len(v) for v in ncbi.values())}\tncbi_pieces\t{total_pieces}\tpieces_matched\t{len(matched_pieces)}\n")
print(open(f"{out}/ch00_vs_ncbi_summary.tsv").read())
EOF
