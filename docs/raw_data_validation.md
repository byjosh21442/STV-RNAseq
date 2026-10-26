# Raw data acquisition and validation (Phases 2D–2F)

Date: 2026-10-25 (download 16:56–19:34 local time, UTC−6) · project `~/projects/STV_GSE137303`.

## Acquisition

- Script: [`scripts/bash/download_fastq.sh`](../scripts/bash/download_fastq.sh). URLs = `https://` + the ENA `fastq_ftp` value and MD5/size = ENA `fastq_md5`/`fastq_bytes`, all read from `data/metadata/samples_metadata.tsv` (Phase 1 ENA record). No URL was typed or reconstructed by hand.
- Safety: disk check before starting (needed 5,718,108,512 B + 5 GiB headroom; free 939 GB); download to `.part` with `curl -C -` resume; a file is accepted only if **size and MD5** both match; mismatching files would be quarantined as `*.corrupt.N` and re-fetched once; non-zero exit if any sample fails.
- Network: ≈0.6 MB/s on average; total wall time 9,489 s. The ENA HTTPS server repeatedly closed connections (`curl: (56) OpenSSL SSL_read: unexpected eof while reading`). Each interruption resumed from the partial file. For SRR10102795 and SRR10102796 the first cycle exhausted its 5 resume attempts, so the second cycle continued the same `.part` file. All final MD5s match ENA, so the interruptions did not corrupt any file.
- Files keep their ENA names in `data/raw/` (git-ignored); per-attempt log in `data/raw/download_log.tsv` (git-ignored); final table in [`results/qc/raw_fastq_checksums.tsv`](../results/qc/raw_fastq_checksums.tsv).
- Disk: `/` 939 GB free before → 934 GB free after; `data/raw` = 5.4 GiB. No decompressed copies were written.

## Checksum table (all PASS)

| sample | SRR | FASTQ file | size (bytes) | expected MD5 (ENA) | observed MD5 | status |
|---|---|---|---:|---|---|---|
| STV_free_rep1 | SRR10102791 | SRR10102791.fastq.gz | 988,076,517 | b0f6e47d0f0eeb1446f58734031294c1 | b0f6e47d0f0eeb1446f58734031294c1 | PASS |
| STV_free_rep2 | SRR10102792 | SRR10102792.fastq.gz | 896,430,396 | 2963114662b147c393b34eb282cf3519 | 2963114662b147c393b34eb282cf3519 | PASS |
| STV_free_rep3 | SRR10102793 | SRR10102793.fastq.gz | 780,172,066 | 122d5fdae38bc48689b5f724dff17425 | 122d5fdae38bc48689b5f724dff17425 | PASS |
| STV_infected_rep1 | SRR10102794 | SRR10102794.fastq.gz | 1,049,513,467 | 3230fe1c5114ad4ec9190cee89107d89 | 3230fe1c5114ad4ec9190cee89107d89 | PASS |
| STV_infected_rep2 | SRR10102795 | SRR10102795.fastq.gz | 921,773,507 | fbca3d445f6d761e0ad972b40d5f62ff | fbca3d445f6d761e0ad972b40d5f62ff | PASS |
| STV_infected_rep3 | SRR10102796 | SRR10102796.fastq.gz | 1,082,142,559 | 78ccfe8c2bba181343905693ea63efc7 | 78ccfe8c2bba181343905693ea63efc7 | PASS |
| **Total** | | | **5,718,108,512** | | | **6/6 PASS** |

## Sanity checks (read-only)

Script: [`scripts/bash/check_raw_fastq.sh`](../scripts/bash/check_raw_fastq.sh) (runs only on MD5-PASS files; 116 s). Output: [`results/qc/raw_fastq_sanity.tsv`](../results/qc/raw_fastq_sanity.tsv).

- gzip integrity: `pigz -t` over the whole file.
- FASTQ structure: first 10,000 records — `@` header, `+` separator, sequence alphabet `ACGTN`, quality length = sequence length.
- Read count and length: `seqkit stats -a` over the whole file, compared with ENA `read_count` and the verified 75 nt.

| sample | gzip | structure (10,000 rec.) | reads (SeqKit) | reads (ENA) | match | min / avg / max length | GC % | Q30 % | status |
|---|---|---|---:|---:|---|---|---:|---:|---|
| STV_free_rep1 | OK | OK | 18,080,482 | 18,080,482 | YES | 75 / 75.0 / 75 | 42.81 | 97 | PASS |
| STV_free_rep2 | OK | OK | 16,288,249 | 16,288,249 | YES | 75 / 75.0 / 75 | 44.03 | 97 | PASS |
| STV_free_rep3 | OK | OK | 14,318,362 | 14,318,362 | YES | 75 / 75.0 / 75 | 44.14 | 97 | PASS |
| STV_infected_rep1 | OK | OK | 19,230,359 | 19,230,359 | YES | 75 / 75.0 / 75 | 43.80 | 97 | PASS |
| STV_infected_rep2 | OK | OK | 16,858,900 | 16,858,900 | YES | 75 / 75.0 / 75 | 44.16 | 97 | PASS |
| STV_infected_rep3 | OK | OK | 19,862,460 | 19,862,460 | YES | 75 / 75.0 / 75 | 43.78 | 97 | PASS |

Total reads: 104,638,812 (= ENA total).

## Read-header observation

Read names keep the original Illumina header, e.g. `@SRR10102791.1 SN679:271:HN3JNBCXY:1:1103:1455:2134/1`. In the first 200,000 reads of **every** file, 100 % come from instrument `SN679`, flowcell `HN3JNBCXY`, lane 1. As far as these reads show, all six libraries were sequenced in the same lane, so sequencing lane does not line up with condition. Only the first 200,000 reads per file were checked; this is not a full-file verification.
