# 04_categorize_phenotypes.py — Categorise Growth Phenotypes

---

## Purpose

Reads grouped phenotype data (3 sheets from step 03), assigns each gene a
fine-grained `Category` and a coarse `Growth_tier` (1–5) using the
signal-based `classify_growth()` engine.

Modifier-led comma segments are filtered out before classification via
`filter_primary_segments()` so that secondary descriptions (e.g.
`some germination long`) do not trigger growth signals.

---

## Data Flow

```
data/3_grouped_genes/                    Grouped genes from step 03
  Hayles_2013_OB_grouped_genes.xlsx      (4 sheets: one/multi/inconsistent/all)
                             ↓
data/previous_manual_check_of_           Manual annotations for the
  insistent_phenotypes/                  131 inconsistent genes
  Inconsistent_phenotypes_at_25_32_
  manual.xlsx
                             ↓
                  04_categorize_phenotypes.py
                             ↓
data/4_categorized_genes/                Categorized output
  Hayles_2013_OB_categorized_            3 data sheets
  phenotypes.xlsx
  Hayles_2013_OB_inspection_             14 inspection pivot sheets
  phenotypes.xlsx
```

---

## Processing Logic

### One Basic Phenotype (~4,400 genes)

Applies `filter_primary_segments()` then `classify_growth()` to the
`Basic phenotype` column. Modifier-led comma segments are dropped before
classification.

### Multi Basic Phenotypes (~310 genes)

Same as above — `filter_primary_segments()` + `classify_growth()` on the
`Basic phenotype` column.

### Inconsistent Phenotypes (131 genes)

Merges with the manually curated annotation file (keyed on `Systematic ID`
→ `SysID`). Uses `Category_32` from the manual file as the final
`Category` (normalising legacy `WT` → `WT-like`), then derives
`Growth_tier` via a reverse lookup table.

---

## Output Files

### `data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx`

**Data sheets (3):**
- `One basic phenotype`
- `Multi basic phenotypes`
- `Inconsistent phenotypes`

Each row includes the original columns plus:
- `Category` — fine-grained growth category (e.g., `germinated, spores`)
- `Growth_tier` — coarse tier (1–5)

### `data/4_categorized_genes/Hayles_2013_OB_inspection_phenotypes.xlsx`

**Inspection pivot sheets (14):** per-branch Phenotypes, Essentiality,
Classification, Growth_tier pivots, plus a multi-level pivot
(full description text × 4-level classification) for visual quality review.

---

## Growth Category → Tier Mapping

| Tier | Categories |
|---|---|
| 1 | `spores` |
| 2 | `germinated`, `germinated and divided`, `germinated, spores`, etc. |
| 3 | `microcolonies` |
| 4 | `very small colonies`, `small colonies` |
| 5 | `WT-like` (default) |

---

## Usage

```bash
mamba run -n bioinformatics python src/04_categorize_phenotypes.py
mamba run -n bioinformatics python src/04_categorize_phenotypes.py --verbose
```