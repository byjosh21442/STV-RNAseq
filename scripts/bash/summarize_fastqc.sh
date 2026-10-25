#!/usr/bin/env bash
# Summarize raw FastQC results (module flags + key metrics) into small TSV tables for docs/raw_qc_assessment.md.
# Usage (project root, stv-rnaseq env): bash scripts/bash/summarize_fastqc.sh [fastqc_dir]
set -euo pipefail
FASTQC_DIR="${1:-results/qc/fastqc_raw}"
OUT_SUMMARY="results/qc/raw_qc_summary.tsv"
OUT_OVERREP="results/qc/raw_overrepresented.tsv"

python3 - "$FASTQC_DIR" "$OUT_SUMMARY" "$OUT_OVERREP" <<'PY'
import sys, zipfile, glob, os

fastqc_dir, out_summary, out_overrep = sys.argv[1:4]
zips = sorted(glob.glob(os.path.join(fastqc_dir, "*_fastqc.zip")))
if not zips:
    sys.exit(f"no *_fastqc.zip in {fastqc_dir}")

def modules(text):
    mods, headers, name, rows = {}, {}, None, []
    for line in text.splitlines():
        if line.startswith(">>END_MODULE"):
            mods[name] = rows
            name, rows = None, []
        elif line.startswith(">>"):
            name, rows = line[2:].split("\t")[0], []
        elif name is not None and line.startswith("#Total Deduplicated Percentage"):
            rows.append(["__dedup__", line.split("\t")[1]])
        elif name is not None and line.startswith("#"):
            headers[name] = line[1:].split("\t")
        elif name is not None:
            rows.append(line.split("\t"))
    return mods, headers

summary_rows, overrep_rows, flag_cols = [], [], None
for z in zips:
    sample = os.path.basename(z).replace("_fastqc.zip", "")
    with zipfile.ZipFile(z) as zf:
        root = zf.namelist()[0].split("/")[0]
        flags = [l.split("\t") for l in zf.read(f"{root}/summary.txt").decode().splitlines()]
        mods, headers = modules(zf.read(f"{root}/fastqc_data.txt").decode())
    flag_cols = flag_cols or [f[1] for f in flags]
    basic = {r[0]: r[1] for r in mods["Basic Statistics"]}
    pbq = mods["Per base sequence quality"]
    psq = mods["Per sequence quality scores"]
    nct = mods["Per base N content"]
    ada = mods.get("Adapter Content", [])
    dup = dict((r[0], r[1]) for r in mods["Sequence Duplication Levels"] if r[0] == "__dedup__")
    ovr = mods.get("Overrepresented sequences", [])

    ada_max, ada_name = 0.0, "NA"
    if ada:
        header = headers.get("Adapter Content", [])
        for r in ada:
            for i, v in enumerate(r[1:], start=1):
                try:
                    if float(v) > ada_max:
                        ada_max, ada_name = float(v), (header[i] if i < len(header) else f"col{i}")
                except ValueError:
                    pass
    psq_mode = max(psq, key=lambda r: float(r[1]))[0]
    summary_rows.append({
        "sample": sample,
        "total_sequences": basic.get("Total Sequences"),
        "poor_quality": basic.get("Sequences flagged as poor quality", "NA"),
        "sequence_length": basic.get("Sequence length"),
        "gc_pct": basic.get("%GC"),
        "min_mean_base_quality": f"{min(float(r[1]) for r in pbq):.1f}",
        "min_lower_quartile": f"{min(float(r[3]) for r in pbq):.1f}",
        "per_seq_quality_mode": psq_mode,
        "max_n_content_pct": f"{max(float(r[1]) for r in nct):.3f}",
        "max_adapter_pct": f"{ada_max:.3f}",
        "max_adapter_type": ada_name,
        "dedup_pct_remaining": f"{float(dup.get('__dedup__', 'nan')):.2f}",
        "n_overrepresented": str(len(ovr)),
        **{f"flag:{f[1]}": f[0] for f in flags},
    })
    for r in ovr:
        overrep_rows.append([sample] + r)

cols = list(summary_rows[0].keys())
with open(out_summary, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for row in summary_rows:
        fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
with open(out_overrep, "w") as fh:
    fh.write("sample\tsequence\tcount\tpercentage\tpossible_source\n")
    for r in overrep_rows:
        fh.write("\t".join(r) + "\n")
print(f"wrote {out_summary} ({len(summary_rows)} samples) and {out_overrep} ({len(overrep_rows)} rows)")
PY
