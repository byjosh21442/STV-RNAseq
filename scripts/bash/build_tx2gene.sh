#!/usr/bin/env bash
# Build reference/ITAG4.0_tx2gene.tsv from the ITAG4.0 GFF (mRNA -> parent gene) and validate it against the
# cDNA FASTA, the Salmon index (collapsed duplicates) and every quant.sf listed in the workflow's samples.csv.
# Identifiers are taken verbatim from the GFF Name attributes (== cDNA headers): no version stripping.
# Usage (project root): bash scripts/bash/build_tx2gene.sh
# Exit 1 on any material mismatch.
set -euo pipefail

python3 - <<'PY'
import csv, json, os, sys
from collections import Counter

gff = "reference/SGN/ITAG4.0/ITAG4.0_gene_models.gff"
cdna = "reference/SGN/ITAG4.0/ITAG4.0_cDNA.fasta"
out = "reference/ITAG4.0_tx2gene.tsv"
val = "reference/ITAG4.0_tx2gene_validation.tsv"
samples_csv = "results/quantification/salmon/samples.csv"
index_dir = "results/genome/index/salmon_index"
rows, fail = [], []

def rec(metric, value, ok=True, note=""):
    rows.append((metric, str(value), "OK" if ok else "FAIL", note))
    if not ok:
        fail.append(metric)

attrs = lambda s: dict(kv.split("=", 1) for kv in s.strip().split(";") if "=" in kv)
gene_name, pairs = {}, []
with open(gff) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        a = attrs(f[8])
        if f[2] == "gene":
            gene_name[a["ID"]] = a["Name"]
        elif f[2] == "mRNA":
            pairs.append((a["Name"], a["Parent"]))

tx2gene = []
for tx, parent in pairs:
    if parent not in gene_name:
        rec("mRNA_parent_missing", tx, False)
        continue
    tx2gene.append((tx, gene_name[parent]))
txc = Counter(t for t, _ in tx2gene)
rec("tx2gene_rows", len(tx2gene))
rec("tx2gene_unique_transcripts", len(txc), len(txc) == len(tx2gene))
rec("tx2gene_unique_genes", len({g for _, g in tx2gene}))
rec("transcripts_mapped_to_more_than_one_gene", sum(1 for v in txc.values() if v > 1), all(v == 1 for v in txc.values()))
rec("transcript_name_starts_with_gene_name", sum(1 for t, g in tx2gene if t.startswith(g + ".")), all(t.startswith(g + ".") for t, g in tx2gene))

cdna_ids = {l[1:].split()[0] for l in open(cdna) if l.startswith(">")}
txset = set(txc)
rec("cDNA_ids_missing_from_tx2gene", len(cdna_ids - txset), not (cdna_ids - txset))
rec("tx2gene_ids_missing_from_cDNA", len(txset - cdna_ids), not (txset - cdna_ids))

with open(out, "w") as fh:
    fh.write("TXNAME\tGENEID\n")
    for t, g in tx2gene:
        fh.write(f"{t}\t{g}\n")

# Salmon index: duplicate sequences collapsed at indexing time (not present in quant.sf).
dup_removed = set()
dup_file = next((os.path.join(index_dir, n) for n in ("duplicate_clusters.tsv",) if os.path.exists(os.path.join(index_dir, n))), None)
if dup_file:
    with open(dup_file) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                dup_removed.add(p[1])
    rec("index_duplicate_clusters_file", dup_file)
    rec("index_collapsed_duplicate_transcripts", len(dup_removed))
else:
    rec("index_duplicate_clusters_file", "absent", True, "no duplicate_clusters.tsv in index dir")

if not os.path.exists(samples_csv):
    rec("quant_sf_check", "PENDING", False, f"{samples_csv} not found; run quantification first")
else:
    ref_names = None
    with open(samples_csv) as fh:
        samples = list(csv.DictReader(fh))
    rec("samples_in_samples_csv", len(samples), len(samples) == 6)
    for s in samples:
        qsf = os.path.join(s["quant_dir"], "quant.sf")
        with open(qsf) as fh:
            next(fh)
            names = [l.split("\t", 1)[0] for l in fh]
        if ref_names is None:
            ref_names = names
            nset = set(names)
            rec("quant_sf_rows", len(names))
            rec("quant_sf_names_not_in_tx2gene", len(nset - txset), not (nset - txset),
                note="exact string match; no version stripping")
            missing = txset - nset
            rec("tx2gene_names_not_in_quant_sf", len(missing), missing <= dup_removed,
                note="must equal transcripts collapsed as duplicates by salmon index")
            rec("unexplained_tx2gene_names_absent_from_quant_sf", len(missing - dup_removed), not (missing - dup_removed),
                note=",".join(sorted(missing - dup_removed)[:5]))
        rec(f"quant_sf_names_identical_to_first:{s['sample']}", names == ref_names, names == ref_names)

with open(val, "w") as fh:
    fh.write("metric\tvalue\tstatus\tnote\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
for r in rows:
    print(f"{r[2]:5s} {r[0]:55s} {r[1][:60]} {r[3][:70]}")
if fail:
    print("FAILED: " + ", ".join(fail)); sys.exit(1)
print(f"tx2gene OK -> {out}")
PY
