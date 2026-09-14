# 00_download_pombase_annotation.py — Download PomBase Gene Annotation

---

## Purpose

Downloads a specific monthly release of the PomBase gene annotation file
(`gene_IDs_names_products.tsv`) and caches it locally. The cached file serves
as the reference for the subsequent gene ID mapping step
(`01_update_systematic_ids.py`), which maps legacy gene identifiers from the
Hayles 2013 supplementary table to their current PomBase equivalents.

---

## Data Source

PomBase monthly release directory:

```
https://www.pombase.org/monthly_releases/{year}/pombase-{release}/
    gene_names_and_identifiers/gene_IDs_names_products.tsv
```

Default version: `2026-06-01` (June 1, 2026 release), about 1.27 MB, 12,685 rows.

---

## File Structure (8 columns)

| Column (internal `PomBaseCol` member) | Description | Example |
|---|---|---|
| `GENE_SYSID` (`gene_systematic_id`) | PomBase systematic identifier | `SPAC1002.01` |
| (not used) `gene_systematic_id_with_prefix` | Systematic ID with database prefix | `PomBase:SPAC1002.01` |
| `GENE_NAME` (`gene_name`) | Standard gene symbol | `mrx11` |
| (not used) `chromosome_id` | Chromosome the gene is on | `chromosome_1` |
| (not used) `gene_product` | Description of the gene product | `conserved fungal protein` |
| (not used) `external_id` | External database ID | `Q9UR06` |
| `GENE_TYPE` (`gene_type`) | Feature type | `protein coding gene` |
| `SYNONYMS` (`synonyms`) | Alternative names (comma-separated) | `SPAC1610.05` |

Covers all annotated feature types: protein-coding genes, tRNAs, rRNAs, snRNAs,
snoRNAs, lncRNAs, pseudogenes, etc.

---

## Implementation Notes

### 1. URL Construction

```python
POMBASE_URL_TEMPLATE = (
    "https://www.pombase.org/monthly_releases/{year}/pombase-{release}/"
    "gene_names_and_identifiers/gene_IDs_names_products.tsv"
)
```

The release directory includes a `pombase-` prefix (`pombase-2026-06-01/`).
`{year}` is extracted from the first 4 characters of the release string.

### 2. Caching

- If the file already exists locally, the download is skipped and only file
  metadata is printed.
- `--force` forces a fresh download.
- On download failure, any partial file is cleaned up (`.unlink()`).

### 3. File Inspection

After download (or cache hit), `inspect_file()` prints:
- File size
- Column count and names
- Comment line count
- Total lines and estimated data rows

### 4. Error Handling

- Non-200 HTTP status → `RuntimeError`
- Network errors (`URLError`) → wrapped as `RuntimeError`
- All core functions use `@logger.catch` for automatic exception capture

---

## CLI Usage

```bash
# Default: download 2026-06-01 release to data/references/
mamba run -n bioinformatics python scripts/00_download_pombase_annotation.py

# Specify a different release
mamba run -n bioinformatics python scripts/00_download_pombase_annotation.py --release 2026-05-01

# Force re-download (overwrite cache)
mamba run -n bioinformatics python scripts/00_download_pombase_annotation.py --force

# Debug mode
mamba run -n bioinformatics python scripts/00_download_pombase_annotation.py --verbose

# Custom output directory
mamba run -n bioinformatics python scripts/00_download_pombase_annotation.py --outdir data/references
```

---

## Why `urllib` + `loguru` Instead of `requests`

- `urllib.request` is part of the standard library — zero extra dependencies,
  sufficient for a simple download task.
- `loguru` is the logging library chosen for the project (available in the
  `bioinformatics` conda/mamba environment).
- Avoiding `requests` keeps the script lightweight.

---

## Output

```
data/references/pombase-{release}_gene_IDs_names_products.tsv
```

Example: `data/references/pombase-2026-06-01_gene_IDs_names_products.tsv`
