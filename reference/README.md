# reference/

Reference data for the approved build **SGN SL4.0 + ITAG4.0** (D-001).

| Path | Tracked | Content |
|---|---|---|
| `SGN/SL4.0/`, `SGN/ITAG4.0/` | no (git-ignored) | Files downloaded from <https://solgenomics.net/ftp/tomato_genome/>, recreated with `bash scripts/bash/download_reference.sh` |
| `reference_manifest.tsv` | yes | Source URL, versions, size, SHA-256, server Last-Modified, retrieval date (UTC), purpose of every downloaded file |
| `reference_validation.tsv` | yes | Output of `scripts/bash/validate_reference.sh` (genome ↔ GFF ↔ cDNA compatibility) |
| `ITAG4.0_tx2gene.tsv` | yes | `TXNAME	GENEID` from the GFF (mRNA `Name` → gene `Name`), built by `scripts/bash/build_tx2gene.sh` |
| `ITAG4.0_tx2gene_validation.tsv` | yes | tx2gene checks against cDNA, Salmon index duplicates and all `quant.sf` |

The Salmon index is not stored here. The workflow publishes it to `results/genome/index/` (git-ignored, 2.1 GB). Verify downloaded files with `sha256sum` against the manifest. Details: `docs/reference_validation.md`.
