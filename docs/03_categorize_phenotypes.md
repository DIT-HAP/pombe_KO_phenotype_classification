# 02_categorize_phenotypes.py — Categorise Growth Phenotypes

---

## Purpose

Reads grouped phenotype data (3 sheets from step 03), assigns each gene a
fine-grained `Category` and a coarse `Growth_tier` (1–5) using the
signal-based `classify_growth()` engine.

Replaces the three original `categorize_genes_with_*_phenotypes.py` scripts.

---

## Data Flow

```
data/2_grouped_genes/                    Grouped genes from step 02
  Hayles_2013_OB_grouped_genes.xlsx      (3 sheets: one/multi/inconsistent)
                             ↓
data/previous_manual_check_of_           Manual annotations for the
  insistent_phenotypes/                  131 inconsistent genes
  Inconsistent_phenotypes_at_25_32_
  manual.xlsx
                             ↓
                  02_categorize_phenotypes.py
                             ↓
data/3_categorized_genes/                Categorized output
  Hayles_2013_OB_categorized_            3 data sheets + 8 pivot tables
  phenotypes.xlsx
```

---

## Processing Logic

### One Basic Phenotype (4,107 genes)

Applies `classify_growth()` directly to the `Basic phenotype` column.
Produces 8 unique categories.

### Multi Basic Phenotypes (605 genes)

Same as above — `classify_growth()` on the `Basic phenotype` column.
Produces 11 unique categories (including composites).

### Inconsistent Phenotypes (131 genes)

Merges with the manually curated annotation file (keyed on `Systematic ID`
→ `SysID`). Uses `Category_32` from the manual file as the final
`Category`, then derives `Growth_tier` via a reverse lookup table.

---

## Output File

`data/3_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx`

**Data sheets (3):**
- `One basic phenotype`
- `Multi basic phenotypes`
- `Inconsistent phenotypes`

Each row includes the original columns plus:
- `Category` — fine-grained growth category (e.g., `germinated, spores`)
- `Growth_tier` — coarse tier (1–5)

**Pivot tables (8):**
- Phenotypes pivot (Category × Hayles phenotype description)
- Essentiality pivot (Category × dispensability)
- Classification pivot (Category × Hayles classification)
- Growth_tier pivot (Growth_tier × Hayles classification)
- (×2 for one/multi branches)

---

## Growth Category → Tier Mapping

| Tier | Categories |
|---|---|
| 1 | `spores` |
| 2 | `germinated`, `germinated and divided`, `germinated, spores`, etc. |
| 3 | `microcolonies` |
| 4 | `very small colonies`, `small colonies` |
| 5 | `WT` (default) |

---

## Usage

```bash
mamba run -n bioinformatics python src/02_categorize_phenotypes.py
mamba run -n bioinformatics python src/02_categorize_phenotypes.py --verbose
```
