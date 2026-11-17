# 04_categorize_phenotypes.py — Categorise Growth Phenotypes

---

## Purpose

Reads grouped phenotype data (3 sheets from step 03) and assigns each gene a
fine-grained `Sub_category` using the signal-based `classify_growth()` engine.

`classify_growth()` handles modifier-led comma segments internally —
segments starting with a modifier word (e.g. `some`, `often`, `occasionally`)
are treated as secondary descriptions and their signals are preserved with a
modifier prefix (e.g. `some germinated`, `often divided`).

---

## Data Flow

```
data/3_grouped_genes/                    Grouped genes from step 03
  Hayles_2013_OB_grouped_genes.xlsx      (4 sheets: one/multi/inconsistent/all)
                             ↓
data/references/                         Manual annotations for the
  previous_manual_check_of_                131 inconsistent genes
  insistent_phenotypes/
  Inconsistent_phenotypes_at_25_32_
  manual.xlsx
                             ↓
                  04_categorize_phenotypes.py
                             ↓
data/4_categorized_genes/                Categorized output
  Hayles_2013_OB_categorized_            3 data sheets + All genes
  phenotypes.xlsx
  Hayles_2013_OB_inspection_             12 inspection sheets
  phenotypes.xlsx                        (flat inspection + 11 pivots)
```

---

## Processing Logic

### One Basic Phenotype (~4,400 genes)

Applies `classify_growth()` to the `Basic phenotype` column. Modifier-led
comma segments are handled internally by `classify_growth()`.

### Multi Basic Phenotypes (~310 genes)

Same as above — `classify_growth()` on the `Basic phenotype` column.

### Inconsistent Phenotypes (131 genes)

Merges with the manually curated annotation file (keyed on `Systematic ID`
→ `SysID`). Uses `Category_32` from the manual file as the final
`Sub_category` (normalising legacy `WT` → `WT-like`).

---

## Output Files

### `data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx`

**Data sheets (4):**
- `One basic phenotype`
- `Multi basic phenotypes`
- `Inconsistent phenotypes`
- `All genes`

Each row includes the original columns plus:
- `Sub_category` — fine-grained growth category (e.g., `spores, germinated`)

### `data/4_categorized_genes/Hayles_2013_OB_inspection_phenotypes.xlsx`

**Flat inspection sheet (1):** `Flat inspection` — one row per unique
description with Sub_category, Phenotype_count, Consistency_25_32,
and Gene count.

**Inspection pivot sheets (11):** per-branch Phenotypes, Essentiality,
Classification pivots, plus a multi-level pivot
(full description text × 3-level classification) for visual quality review.

---

## Usage

```bash
mamba run -n bioinformatics python scripts/04_categorize_phenotypes.py
mamba run -n bioinformatics python scripts/04_categorize_phenotypes.py --verbose
```
