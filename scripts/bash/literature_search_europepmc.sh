#!/usr/bin/env bash
# STV_GSE137303 - Phase 7 (Task 4): reproducible Europe PMC literature searches (metadata only; no abstracts stored).
# Each query is run once with the Europe PMC REST API (default relevance order; an explicit sort=RELEVANCE returned HTTP 503);
# hit counts and the returned records are saved.
# Selection of records for reading is manual and documented in docs/literature_biological_context.md.
# Usage (project root): bash scripts/bash/literature_search_europepmc.sh
set -euo pipefail
cd "$(dirname "$0")/../.."
out=results/biological_integration/literature
mkdir -p "$out"
python3 - "$out" <<'EOF'
import json, sys, urllib.parse, urllib.request, datetime, time
out = sys.argv[1]
queries = [
    ("Q1_stv", 'TITLE_ABS:"southern tomato virus"', 100),
    ("Q2_amalgavirus", '(TITLE_ABS:"amalgavirus" OR TITLE_ABS:"amalgaviridae")', 50),
    ("Q3_tomato_virus_transcriptome", '(TITLE:"tomato" AND TITLE_ABS:"virus" AND (TITLE_ABS:"RNA-seq" OR TITLE_ABS:"transcriptome" OR TITLE_ABS:"transcriptomic") AND (TITLE_ABS:"leaf" OR TITLE_ABS:"leaves"))', 40),
    ("Q4_persistent_virus_host_expression", '((TITLE_ABS:"persistent virus" OR TITLE_ABS:"persistent viruses" OR TITLE_ABS:"endornavirus" OR TITLE_ABS:"cryptic virus" OR TITLE_ABS:"partitivirus") AND (TITLE_ABS:"transcriptome" OR TITLE_ABS:"RNA-seq" OR TITLE_ABS:"gene expression") AND TITLE_ABS:"plant")', 30),
]
day = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
log = open(f"{out}/europepmc_search_log.tsv", "w")
log.write("query_id\tquery\thit_count\trecords_saved\tretrieved_utc\n")
for qid, q, n in queries:
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
        {"query": q, "format": "json", "pageSize": n, "resultType": "lite"})
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "STV_GSE137303-literature-search/1.0"})
            d = json.load(urllib.request.urlopen(req, timeout=60))
            break
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504) or attempt == 4:
                raise
            time.sleep(5 * (attempt + 1))
    recs = d.get("resultList", {}).get("result", [])
    with open(f"{out}/europepmc_{qid}.tsv", "w") as o:
        o.write("rank\tsource\tid\tpmid\tdoi\tyear\tjournal\tpub_type\ttitle\n")
        for i, r in enumerate(recs, 1):
            o.write("\t".join(str(x).replace("\t", " ") for x in [i, r.get("source", ""), r.get("id", ""), r.get("pmid", ""),
                    r.get("doi", ""), r.get("pubYear", ""), r.get("journalTitle", r.get("bookOrReportDetails", {}).get("publisher", "") if isinstance(r.get("bookOrReportDetails"), dict) else ""),
                    r.get("pubType", ""), r.get("title", "")]) + "\n")
    log.write(f"{qid}\t{q}\t{d.get('hitCount', '')}\t{len(recs)}\t{day}\n")
    print(qid, d.get("hitCount"), len(recs))
    time.sleep(1)
log.close()
EOF
