# 01_update_systematic_ids.py — Update Gene Systematic IDs from PomBase Annotation

---

## Purpose

Maps the 4,843 gene systematic IDs from the Hayles 2013 supplementary table
against a current PomBase annotation release. Since the original paper was
published in 2013, some gene identifiers may have changed due to genome
annotation updates (gene merges, splits, renames, or reclassification to
pseudogene or non-coding RNA).

This script is the second step in the pipeline, after
`00_download_pombase_annotation.py` has cached the annotation file.

---

## Data Flow

```
data/raw/rsob130053supp2.xlsx          Hayles 2013 raw supplementary table
                             ↓
data/references/pombase-{release}_      PomBase annotation (from step 00)
  gene_IDs_names_products.tsv
                             ↓
                   01_update_systematic_ids.py
                             ↓
data/1_formatted/                      Updated phenotypes (Systematic ID +
  Hayles_2013_OB_formatted_              Gene name columns updated)
  phenotypes.xlsx
                             ↓
data/references/                       Change log (every gene's mapping
  gene_id_mapping_changelog.xlsx         outcome)
```

---

## Mapping Strategy

The `update_sysids()` function implements a three-layer lookup, with a
fourth diagnostic fallback:

### Layer 1 — Direct Systematic ID Match (fast path)
Check if the input ID exists in the current PomBase set of protein-coding
systematic IDs. If yes, no change needed.

### Layer 2 — Gene Name Lookup
If the input isn't a current systematic ID, check if it matches a current
gene symbol (e.g., `pom34` → `SPAC1002.02`).

### Layer 3 — Synonym Lookup
If neither of the above match, check the `synonyms` column in PomBase.
Some old IDs are kept as synonyms when a gene is renamed or merged
(e.g., `SPAC823.02` → `SPAC823.17` via synonym).

### Layer 4 — Reclassification Diagnosis (not found in any lookup)
If all three lookups fail within the protein-coding gene set, the script
searches the **full** annotation (all feature types):
1. Check if the ID still exists as a `gene_systematic_id` in another
   feature type (e.g., pseudogene, lncRNA).
2. Check if the ID appears in the `synonyms` column of a non-coding gene,
   indicating the old systematic ID was retired but the locus was
   re-annotated (e.g., `SPBC8E4.02c` → synonym of `SPNCRNA.9001`, a lncRNA).

---

## Mapping Outcomes

| `update_type` | Meaning | Count |
|---|---|---|
| `NO_CHANGE` | Systematic ID is current | 4,838 |
| `VIA_NAME` | Mapped via gene name (no instances in this dataset) | 0 |
| `VIA_SYNONYM` | Mapped via synonym (old ID found in synonyms column) | 1 |
| `RECLASSIFIED` | Gene exists but no longer annotated as protein-coding | 4 |
| `NOT_FOUND` | Not found anywhere in the annotation | 0 |
| `CONFLICT` | Name/synonym maps to multiple systematic IDs | 0 |

---

## Detailed Changes (2026-06-01 vs 2013)

### SPAC823.02 → SPAC823.17
- The old systematic ID `SPAC823.02` was retired. It now exists as a
  synonym for `SPAC823.17` (gene name: `tom6`).
- The phenotype data follows the current ID.
- **Update type:** VIA_SYNONYM

### SPBC1348.14c (ght7) — Reclassified as pseudogene
- Previously annotated as protein-coding in 2013; now a pseudogene.
- The phenotype observation remains valid; the gene simply no longer
  has a protein product annotation.

### SPBC530.06c (clu1) — Reclassified as pseudogene

### SPCC622.17 (apn1) — Reclassified as pseudogene

### SPBC8E4.02c → SPNCRNA.9001 (prt2) — Reclassified as lncRNA
- The old systematic ID `SPBC8E4.02c` was retired entirely and is now
  a synonym for the lncRNA `SPNCRNA.9001`.
- This is the most significant reannotation — what was thought to be a
  protein-coding gene is now annotated as a long non-coding RNA.

---

## Gene Name Updates

The input table had **1,923 genes** (39.7%) with missing gene names (`NaN`).
After looking up current names from PomBase:

| Metric | Before | After |
|---|---|---|
| Missing gene names | 1,923 | 652 |
| Filled by PomBase | — | 1,271 |

The remaining 652 missing names correspond to genes that have no assigned
gene symbol in the current PomBase annotation (systematic-ID-only entries).

---

## Output Files

### `data/1_formatted/Hayles_2013_OB_formatted_phenotypes.xlsx`

Same 12-column structure as the raw input, but with `Systematic ID` and
`Gene name` columns updated. This file is the direct input to the next
pipeline step (`02_format_phenotype_input.py`).

### `data/references/gene_id_mapping_changelog.xlsx`

Columns:
- `original_systematic_id` — ID from Hayles 2013 table
- `original_gene_name` — Gene name from Hayles 2013 table (may be NaN)
- `updated_systematic_id` — Current systematic ID after mapping
- `updated_gene_name` — Current gene name from PomBase
- `note` — Detailed mapping description
- `update_type` — Categorical outcome (NO_CHANGE, VIA_NAME, VIA_SYNONYM,
  RECLASSIFIED, NOT_FOUND, CONFLICT)

---

## Source Code Notes

### `PomBaseCol` StrEnum

Column references use a `StrEnum` (not bare string constants):

```python
class PomBaseCol(StrEnum):
    GENE_SYSID = "gene_systematic_id"
    GENE_NAME = "gene_name"
    GENE_TYPE = "gene_type"
    SYNONYMS = "synonyms"
```

This groups the column names into a single type, providing IDE
autocompletion and catching typos at import time.

### Version Coupling

Scripts 00 and 01 share a single version constant. `00` uses
`DEFAULT_RELEASE` for the download URL; `01` constructs the annotation
path from the same constant so the two scripts stay in sync.

---

## Usage

```bash
# Default: use raw data, 2026-06-01 annotation, write to data/1_formatted/
mamba run -n bioinformatics python src/01_update_systematic_ids.py

# Debug mode
mamba run -n bioinformatics python src/01_update_systematic_ids.py --verbose
```

### Custom paths

```bash
mamba run -n bioinformatics python src/01_update_systematic_ids.py \
    --raw data/raw/rsob130053supp2.xlsx \
    --annotation data/references/pombase-2026-05-01_gene_IDs_names_products.tsv \
    --output data/1_formatted/Hayles_2013_OB_formatted_phenotypes.xlsx \
    --changelog data/references/gene_id_mapping_changelog.xlsx
```

### Non-protein-coding filter

By default, only protein-coding genes are used for the lookup tables. To
include all feature types:

```bash
mamba run -n bioinformatics python src/01_update_systematic_ids.py \
    --gene-filter "gene_type != 'pseudogene'"
```

---

## Dependencies

- pandas, numpy, loguru (all in the `bioinformatics` conda/mamba environment)
- Standard library: argparse, re, sys, pathlib, urllib (for the download script)
