# 04_categorize_phenotypes.py — Categorise Growth Phenotypes

---

## Purpose

Reads grouped phenotype data (3 sheets from step 03), assigns each gene a
fine-grained `Category` using the signal-based `classify_growth()` engine,
then re-ranks `Growth_tier` by sign-flipped DIT-HAP DR median (most depleted
= tier 1; DR is negative, so more negative = more depleted).

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
data/raw/                                Sign-flipped DIT-HAP DR for tier ranking
  gene_level_fitting_statistics_
  minusDR.tsv
                             ↓
                  04_categorize_phenotypes.py
                             ↓
data/4_categorized_genes/                Categorized output
  Hayles_2013_OB_categorized_            3 data sheets + All genes
  phenotypes.xlsx
  Hayles_2013_OB_inspection_             15 inspection sheets
  phenotypes.xlsx                        (flat inspection + 14 pivots)
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
`Category` (normalising legacy `WT` → `WT-like`), then derives
`Growth_tier` via a reverse lookup table.

### Growth_tier Re-ranking

After all three branches are classified and concatenated into `All genes`,
`Growth_tier` is re-ranked based on the median DIT-HAP DR value per
`Category`, from most depleted (tier 1) to least depleted (tier N). Because
DR is sign-flipped (negative, more negative = more depleted), tier 1 is the
category with the lowest median DR. This ensures tier numbers reflect growth
defect severity as measured by depletion rate, not just signal hierarchy.

---

## Output Files

### `data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx`

**Data sheets (4):**
- `One basic phenotype`
- `Multi basic phenotypes`
- `Inconsistent phenotypes`
- `All genes`

Each row includes the original columns plus:
- `Category` — fine-grained growth category (e.g., `spores, germinated`)
- `Growth_tier` — DR-median-ranked tier (1 = most depleted, i.e. most negative DR; N = least)

### `data/4_categorized_genes/Hayles_2013_OB_inspection_phenotypes.xlsx`

**Flat inspection sheet (1):** `Flat inspection` — one row per unique
description with Category, Phenotype_count, Growth_tier, Consistency_25_32,
and Gene count.

**Inspection pivot sheets (14):** per-branch Phenotypes, Essentiality,
Classification, Growth_tier pivots, plus a multi-level pivot
(full description text × 4-level classification) for visual quality review.

---

## Usage

```bash
mamba run -n bioinformatics python scripts/04_categorize_phenotypes.py
mamba run -n bioinformatics python scripts/04_categorize_phenotypes.py --verbose
```
