# 05_merge_categories.py — Merge Categorized Phenotypes

---

## Purpose

Concatenates the three phenotype branches (one / multi / inconsistent) from
the categorized output into a single table for downstream analysis.
Preserves both the fine-grained `Category` and the coarse `Growth_tier`
columns.

---

## Data Flow

```
data/4_categorized_genes/                Categorized output from step 04
  Hayles_2013_OB_categorized_            3 data sheets
  phenotypes.xlsx
                             ↓
                   05_merge_categories.py
                             ↓
data/5_merged_categories/                Final merged table + summaries
  Hayles_2013_OB_merged_categories.xlsx

results/                                 Copy for source control
  Hayles_2013_OB_merged_categories.xlsx
```

---

## Output

### `data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx`

**Data sheet:** `All genes` — 4,843 rows.

Key columns (new or modified):
- `Category` — fine-grained growth category
- `Growth_tier` — coarse tier (1–5)
- `Phenotype_count` — Single / Multiple / Temp_mismatch
- `Consistency_25_32` — Consistent / Only_32 / Mismatch

**Summary sheets (4):**
- `Consistency at temperatures` — counts per consistency status
- `Phenotype count` — counts per Single/Multiple/Temp_mismatch
- `Category` — counts per category
- `Growth_tier` — counts per tier

### `results/Hayles_2013_OB_merged_categories.xlsx`

Identical copy, saved to the version-controlled `results/` directory.

---

## Usage

```bash
mamba run -n bioinformatics python src/05_merge_categories.py
```