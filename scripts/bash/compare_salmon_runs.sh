#!/usr/bin/env bash
# Compare two Salmon result sets (same samples, same index) transcript by transcript.
# Usage (project root): bash scripts/bash/compare_salmon_runs.sh <dir_A> <dir_B> <out.tsv>
#   dir_X = directory containing <sample>/quant.sf and <sample>/aux_info/meta_info.json
# Used for D-010: -l A (results/quantification/salmon) vs explicit -l SR, regenerated with
#   nextflow run workflow/main.nf -resume --run_quant true --salmon_index results/genome/index/salmon_index \
#     --salmon_libtype SR -output-dir /tmp/stv_sr_check
set -euo pipefail
[ $# -eq 3 ] || { echo "usage: $0 <dir_A> <dir_B> <out.tsv>" >&2; exit 2; }

python3 - "$1" "$2" "$3" <<'PY'
import hashlib, json, math, os, sys
da, db, out = sys.argv[1:4]
samples = sorted(d for d in os.listdir(da) if os.path.isfile(os.path.join(da, d, "quant.sf")))
cols = ["sample", "quant_sf_md5_identical", "num_mapped_A", "num_mapped_B", "num_eq_classes_A", "num_eq_classes_B",
        "transcripts", "transcripts_numreads_diff_gt_0.001", "max_abs_numreads_diff", "effective_length_diffs",
        "pearson_log1p_numreads"]
rows = []
def load(p):
    with open(p) as fh:
        next(fh)
        return {f[0]: (float(f[2]), float(f[4])) for f in (l.rstrip("\n").split("\t") for l in fh)}
for s in samples:
    qa, qb = (os.path.join(d, s, "quant.sf") for d in (da, db))
    ma, mb = (json.load(open(os.path.join(d, s, "aux_info", "meta_info.json"))) for d in (da, db))
    a, b = load(qa), load(qb)
    if a.keys() != b.keys():
        sys.exit(f"{s}: transcript sets differ")
    x = [a[k][1] for k in a]; y = [b[k][1] for k in a]
    lx = [math.log1p(v) for v in x]; ly = [math.log1p(v) for v in y]
    mx, my = sum(lx) / len(lx), sum(ly) / len(ly)
    r = sum((p - mx) * (q - my) for p, q in zip(lx, ly)) / math.sqrt(sum((p - mx) ** 2 for p in lx) * sum((q - my) ** 2 for q in ly))
    md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
    rows.append([s, "YES" if md5(qa) == md5(qb) else "NO", ma["num_mapped"], mb["num_mapped"], ma["num_eq_classes"],
                 mb["num_eq_classes"], len(x), sum(1 for p, q in zip(x, y) if abs(p - q) > 1e-3),
                 round(max(abs(p - q) for p, q in zip(x, y)), 3), sum(1 for k in a if abs(a[k][0] - b[k][0]) > 1e-3), f"{r:.6f}"])
with open(out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
for r in rows:
    print("\t".join(map(str, r)))
PY
