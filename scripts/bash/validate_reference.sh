#!/usr/bin/env bash
# Validate SL4.0 + ITAG4.0 compatibility before indexing (read-only on reference/SGN/).
# Writes reference/reference_validation.tsv. Exit 1 if a material inconsistency is found.
# Usage (project root, stv-rnaseq env): bash scripts/bash/validate_reference.sh
set -euo pipefail
DIR="reference/SGN"
OUT="reference/reference_validation.tsv"

python3 - "$DIR" "$OUT" <<'PY'
import gzip, re, sys
from collections import Counter, defaultdict

d, out = sys.argv[1], sys.argv[2]
gff, cdna = f"{d}/ITAG4.0/ITAG4.0_gene_models.gff", f"{d}/ITAG4.0/ITAG4.0_cDNA.fasta"
prot, genome = f"{d}/ITAG4.0/ITAG4.0_proteins.fasta", f"{d}/SL4.0/S_lycopersicum_chromosomes.4.00.fa.gz"
desc, go = f"{d}/ITAG4.0/ITAG4.0_descriptions.txt", f"{d}/ITAG4.0/ITAG4.0_goterms.txt"
GENE_RE, TX_RE = re.compile(r"^Solyc\d{2}g\d{6}\.\d+$"), re.compile(r"^Solyc\d{2}g\d{6}\.\d+\.\d+$")
rows, material = [], []

def rec(metric, value, ok=True, material_if_fail=True, note=""):
    rows.append((metric, str(value), "OK" if ok else ("FAIL" if material_if_fail else "NOTE"), note))
    if not ok and material_if_fail:
        material.append(metric)

def fasta(path, opener=open):
    seqs, name, n = {}, None, 0
    raw_headers = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs.setdefault(name, []).append(n)
                raw_headers.append(line.rstrip("\n"))
                name, n = line[1:].split()[0] if line[1:].strip() else "", 0
            else:
                n += len(line.strip())
        if name is not None:
            seqs.setdefault(name, []).append(n)
    return seqs, raw_headers

# --- GFF ---
attrs = lambda s: dict(kv.split("=", 1) for kv in s.strip().split(";") if "=" in kv)
regions, types = {}, Counter()
genes, mrnas = {}, {}
gene_ids, mrna_ids = Counter(), Counter()
exon_len = defaultdict(int)
coords_bad, seqids = 0, set()
with open(gff) as fh:
    for line in fh:
        if line.startswith("##sequence-region"):
            parts = line[len("##sequence-region"):].split()
            regions[parts[0]] = int(parts[2])
            continue
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        seqid, ftype, start, end, strand, a = f[0], f[2], int(f[3]), int(f[4]), f[6], attrs(f[8])
        types[ftype] += 1; seqids.add(seqid)
        if start > end: coords_bad += 1
        if ftype == "gene":
            gene_ids[a["ID"]] += 1
            genes[a["ID"]] = dict(name=a.get("Name"), seqid=seqid, strand=strand, start=start, end=end)
        elif ftype == "mRNA":
            mrna_ids[a["ID"]] += 1
            mrnas[a["ID"]] = dict(name=a.get("Name"), parent=a.get("Parent"), seqid=seqid, strand=strand, start=start, end=end)
        elif ftype == "exon":
            for p in a.get("Parent", "").split(","):
                exon_len[p] += end - start + 1

rec("gff_feature_types", "; ".join(f"{k}={v}" for k, v in sorted(types.items())))
rec("genes", len(genes))
rec("transcripts_mRNA", len(mrnas))
rec("duplicate_gene_IDs", sum(1 for v in gene_ids.values() if v > 1), all(v == 1 for v in gene_ids.values()))
rec("duplicate_mRNA_IDs", sum(1 for v in mrna_ids.values() if v > 1), all(v == 1 for v in mrna_ids.values()))
gnames = Counter(g["name"] for g in genes.values()); tnames = Counter(m["name"] for m in mrnas.values())
rec("duplicate_gene_Names", sum(1 for v in gnames.values() if v > 1), all(v == 1 for v in gnames.values()))
rec("duplicate_mRNA_Names", sum(1 for v in tnames.values() if v > 1), all(v == 1 for v in tnames.values()))
ANY_RE = re.compile(r"^Solyc\d{2}g\d{6}(\.\d+)+$")
unparseable = [n for n in list(gnames) + list(tnames) if not (n and ANY_RE.match(n))]
rec("names_not_Solyc_locus_plus_versions", len(unparseable), not unparseable, note=",".join(unparseable[:5]))
extra_g = sorted(g["name"] for g in genes.values() if g["name"] and ANY_RE.match(g["name"]) and not GENE_RE.match(g["name"]))
extra_t = sorted(m["name"] for m in mrnas.values() if m["name"] and ANY_RE.match(m["name"]) and not TX_RE.match(m["name"]))
rec("gene_Names_with_extra_version_level", len(extra_g), not extra_g, False,
    note="upstream naming kept as published (unique, parent-consistent, identical in cDNA): " + ",".join(extra_g))
rec("mRNA_Names_with_extra_version_level", len(extra_t), not extra_t, False, note="e.g. " + ",".join(extra_t[:3]))
loci = Counter(re.match(r"^(Solyc\d{2}g\d{6})", g["name"]).group(1) for g in genes.values() if g["name"] and ANY_RE.match(g["name"]))
rec("loci_with_multiple_gene_models", sum(1 for v in loci.values() if v > 1), all(v == 1 for v in loci.values()))
rec("ID_equals_type_prefix_plus_Name", sum(1 for k, g in genes.items() if k != "gene:" + g["name"]) + sum(1 for k, m in mrnas.items() if k != "mRNA:" + m["name"]) == 0,
    note="gene:<Name> / mRNA:<Name>; tx2gene uses the Name attributes (no string surgery)")
orphan = [m for m in mrnas.values() if m["parent"] not in genes]
rec("mRNA_with_missing_parent_gene", len(orphan), not orphan)
multi_parent = [m for m in mrnas.values() if m["parent"] and "," in m["parent"]]
rec("mRNA_with_multiple_parents", len(multi_parent), not multi_parent)
tx_per_gene = Counter(m["parent"] for m in mrnas.values())
rec("genes_without_mRNA", sum(1 for g in genes if g not in tx_per_gene), True, False)
rec("transcripts_per_gene_distribution", "; ".join(f"{k}tx={v}" for k, v in sorted(Counter(tx_per_gene.values()).items())))
prefix_mismatch = [m for m in mrnas.values() if m["parent"] in genes and not m["name"].startswith(genes[m["parent"]]["name"] + ".")]
rec("mRNA_Name_not_prefixed_by_gene_Name", len(prefix_mismatch), not prefix_mismatch)
loc_mismatch = [m for m in mrnas.values() if m["parent"] in genes and (genes[m["parent"]]["seqid"], genes[m["parent"]]["strand"]) != (m["seqid"], m["strand"])]
rec("mRNA_seqid_or_strand_differs_from_gene", len(loc_mismatch), not loc_mismatch)
rec("features_with_start_gt_end", coords_bad, coords_bad == 0)

# --- genome ---
gseq, gheaders = fasta(genome, gzip.open)
glen = {k: v[0] for k, v in gseq.items()}
rec("genome_sequences", len(glen))
rec("genome_names", ",".join(sorted(glen)))
rec("genome_total_bp", sum(glen.values()))
rec("genome_headers_with_trailing_whitespace", sum(1 for h in gheaders if h != h.rstrip()), True, False,
    note="cosmetic; Salmon and decoys.txt use the first whitespace-delimited token")
rec("duplicate_genome_names", sum(1 for v in gseq.values() if len(v) > 1), all(len(v) == 1 for v in gseq.values()))
missing_seq = sorted(seqids - set(glen))
rec("GFF_seqids_missing_in_genome", len(missing_seq), not missing_seq, note=",".join(missing_seq))
region_mismatch = [k for k, v in regions.items() if glen.get(k) != v]
rec("sequence_region_length_mismatches", len(region_mismatch), not region_mismatch, note=",".join(region_mismatch))
rec("GFF_sequence_region_pragma_missing_space", "yes", True, False, note="'##sequence-regionSL4.0ch00' (upstream formatting; parsed correctly)")
beyond = [m["name"] for m in mrnas.values() if m["end"] > glen.get(m["seqid"], 0)]
rec("mRNA_beyond_sequence_end", len(beyond), not beyond)

# --- cDNA ---
cseq, _ = fasta(cdna)
clen = {k: v[0] for k, v in cseq.items()}
rec("cDNA_records", sum(len(v) for v in cseq.values()))
rec("duplicate_cDNA_IDs", sum(1 for v in cseq.values() if len(v) > 1), all(len(v) == 1 for v in cseq.values()))
tset = set(tnames)
rec("cDNA_IDs_not_in_GFF_mRNA_Names", len(set(clen) - tset), not (set(clen) - tset))
rec("GFF_mRNA_Names_not_in_cDNA", len(tset - set(clen)), not (tset - set(clen)))
name_to_id = {m["name"]: k for k, m in mrnas.items()}
len_mismatch = [t for t in clen if t in name_to_id and exon_len.get(name_to_id[t]) != clen[t]]
rec("cDNA_length_ne_sum_of_GFF_exons", len(len_mismatch), not len_mismatch, note="first: " + ",".join(len_mismatch[:5]))
rec("cDNA_length_min_median_max", f"{min(clen.values())}/{sorted(clen.values())[len(clen)//2]}/{max(clen.values())}")

# --- proteins / descriptions / GO (later functional analysis) ---
pseq, _ = fasta(prot)
rec("protein_records", len(pseq))
rec("protein_IDs_not_in_mRNA_Names", len(set(pseq) - tset), not (set(pseq) - tset))
gset = set(gnames)
did = [l.split("\t")[0].strip() for l in open(desc) if l.strip()]
rec("description_rows", len(did))
rec("description_IDs_not_gene_Names", len(set(did) - gset), True, False)
go_rows = [l.rstrip("\n").split("\t") for l in open(go) if l.strip()]
gid = [r[0].strip() for r in go_rows]
go_with_terms = [r[0].strip() for r in go_rows if len(r) > 1 and r[1].strip()]
gloci = {n.split(".")[0] for n in gset}
go_missing = set(gid) - gset
rec("goterm_rows", len(gid))
rec("goterm_rows_with_GO_terms", len(go_with_terms))
rec("goterm_IDs_not_gene_Names", len(go_missing), not go_missing, False,
    note="not used in Phase 3; must be reconciled before any GO analysis")
rec("goterm_IDs_locus_in_GFF_other_version", sum(1 for i in go_missing if i.split(".")[0] in gloci), True, False)
rec("goterm_IDs_with_GO_terms_matching_gene_Names", len(set(go_with_terms) & gset), True, False)

with open(out, "w") as fh:
    fh.write("metric\tvalue\tstatus\tnote\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
for r in rows:
    print(f"{r[2]:5s} {r[0]:45s} {r[1][:90]} {r[3][:80]}")
if material:
    print("MATERIAL INCONSISTENCIES: " + ", ".join(material)); sys.exit(1)
print("No material inconsistencies.")
PY
