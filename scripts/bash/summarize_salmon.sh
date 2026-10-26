#!/usr/bin/env bash
# Summarize Salmon quantification QC for all samples in results/quantification/salmon/samples.csv:
# meta_info.json, lib_format_counts.json and quant.sf -> results/quantification/salmon_qc_summary.tsv
# and results/quantification/salmon_libformat_summary.tsv. Prints consistency checks (no biological claims).
# Usage (project root): bash scripts/bash/summarize_salmon.sh
set -euo pipefail

python3 - <<'PY'
import csv, json, os, statistics, sys

samples_csv = "results/quantification/salmon/samples.csv"
sanity = "results/qc/raw_fastq_sanity.tsv"
sheet = "config/samplesheet.csv"
out_qc = "results/quantification/salmon_qc_summary.tsv"
out_lf = "results/quantification/salmon_libformat_summary.tsv"

raw = {r["sample"]: r for r in csv.DictReader(open(sanity), delimiter="\t")}
srr = {r["sample"]: r["run_accession"] for r in csv.DictReader(open(sheet))}
samples = list(csv.DictReader(open(samples_csv)))

qc_cols = ["sample", "condition", "replicate", "run_accession", "raw_reads", "num_processed", "processed_equals_raw",
           "num_mapped", "percent_mapped", "num_decoy_fragments", "pct_decoy_of_processed", "libtype_requested",
           "detected_library_type", "library_type_applied", "quant_sf_transcripts", "transcripts_numreads_gt0",
           "sum_numreads", "frag_length_mean", "frag_length_sd", "frag_length_source", "mapping_type",
           "num_valid_targets", "num_decoy_targets", "salmon_version", "diagnostics"]
lf_cols = ["sample", "expected_format", "compatible_fragment_ratio", "num_compatible_fragments",
           "num_incompatible_fragments", "num_assigned_fragments", "strand_mapping_bias", "SF", "SR", "U",
           "num_frags_with_concordant_consistent_mappings", "num_frags_with_inconsistent_or_orphan_mappings"]
qc_rows, lf_rows, issues = [], [], []

for s in samples:
    d = s["quant_dir"]
    meta = json.load(open(os.path.join(d, "aux_info", "meta_info.json")))
    lf = json.load(open(os.path.join(d, "lib_format_counts.json")))
    n_tx = n_pos = 0; total = 0.0
    with open(os.path.join(d, "quant.sf")) as fh:
        next(fh)
        for line in fh:
            v = float(line.rstrip("\n").split("\t")[4]); n_tx += 1; total += v; n_pos += v > 0
    rr = int(raw[s["sample"]]["reads_seqkit"])
    proc = int(meta["num_processed"])
    diags = ";".join(f'{x.get("severity")}:{x.get("code")}' for x in meta.get("diagnostics", [])) or "none"
    qc_rows.append({
        "sample": s["sample"], "condition": s["condition"], "replicate": s["replicate"], "run_accession": srr[s["sample"]],
        "raw_reads": rr, "num_processed": proc, "processed_equals_raw": "YES" if proc == rr else "NO",
        "num_mapped": meta["num_mapped"], "percent_mapped": round(float(meta["percent_mapped"]), 2),
        "num_decoy_fragments": meta.get("num_decoy_fragments", "NA"),
        "pct_decoy_of_processed": round(100 * meta.get("num_decoy_fragments", 0) / proc, 2) if proc else "NA",
        "libtype_requested": s["libtype"], "detected_library_type": meta.get("detected_library_type", "NA"),
        "library_type_applied": ",".join(meta.get("library_types", [])), "quant_sf_transcripts": n_tx,
        "transcripts_numreads_gt0": n_pos, "sum_numreads": round(total, 1),
        "frag_length_mean": meta.get("frag_length_mean"), "frag_length_sd": meta.get("frag_length_sd"),
        "frag_length_source": meta.get("frag_length_source", "NA"), "mapping_type": meta.get("mapping_type"),
        "num_valid_targets": meta.get("num_valid_targets"), "num_decoy_targets": meta.get("num_decoy_targets"),
        "salmon_version": meta["salmon_version"], "diagnostics": diags})
    lf_rows.append({k: lf.get(k, "NA") for k in lf_cols if k != "sample"} | {"sample": s["sample"]})

for path, cols, rows in ((out_qc, qc_cols, qc_rows), (out_lf, lf_cols, lf_rows)):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"== {path}")
    for r in rows:
        print("  " + " | ".join(str(r[c]) for c in cols))

# Consistency checks (technical only)
if len(qc_rows) != 6: issues.append(f"expected 6 samples, found {len(qc_rows)}")
for key in ("detected_library_type", "library_type_applied", "salmon_version", "quant_sf_transcripts", "num_valid_targets", "num_decoy_targets"):
    vals = {str(r[key]) for r in qc_rows}
    print(f"check {key}: {sorted(vals)}")
    if len(vals) != 1: issues.append(f"{key} differs between samples: {sorted(vals)}")
for r in qc_rows:
    if r["processed_equals_raw"] != "YES": issues.append(f'{r["sample"]}: num_processed != raw reads')
    if r["diagnostics"] != "none": issues.append(f'{r["sample"]}: diagnostics {r["diagnostics"]}')
    if r["percent_mapped"] < 30: issues.append(f'{r["sample"]}: mapping rate < 30 %')
for r in lf_rows:
    if float(r["compatible_fragment_ratio"]) < 0.95: issues.append(f'{r["sample"]}: compatible_fragment_ratio < 0.95')
by = {}
for r in qc_rows: by.setdefault(r["condition"], []).append(r["percent_mapped"])
for c, v in sorted(by.items()):
    print(f"mapping rate {c}: n={len(v)} mean={statistics.mean(v):.2f} min={min(v):.2f} max={max(v):.2f}")
cfr = [float(r["compatible_fragment_ratio"]) for r in lf_rows]
print(f"compatible_fragment_ratio: min={min(cfr):.4f} max={max(cfr):.4f}")
print("ISSUES: " + ("; ".join(issues) if issues else "none"))
PY
