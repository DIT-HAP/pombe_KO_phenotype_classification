# 03_merge_categories.py — Merge Categorized Phenotypes

---

## Purpose

Concatenates the three phenotype branches (one / multi / inconsistent) from
the categorized output into a single table for downstream analysis.
Preserves both the fine-grained `Category` and the coarse `Growth_tier`
columns.

Replaces the original `merge_the_categories_and_assign_essentiality.py`.

---

## Data Flow

```
data/4_categorized_genes/                Categorized output from step 03
  Hayles_2013_OB_categorized_            3 data sheets
  phenotypes.xlsx
                             ↓
                   03_merge_categories.py
                             ↓
data/5_merged_categories/                Final merged table + summaries
  Hayles_2013_OB_merged_categories.xlsx

results/                                 Copy for source control
  Hayles_2013_OB_merged_categories.xlsx
```

---

## Output

### `data/4_merged_categories/Hayles_2013_OB_merged_categories.xlsx`

**Data sheet:** `All genes` — 4,843 rows, 20 columns.

Key columns (new or modified):
- `Category` — fine-grained growth category
- `Growth_tier` — coarse tier (1–5)
- `Category_25`, `Category_32` — from manual annotation (inconsistent group only)

**Summary sheets (4):**
- `Consistency at temperatures` — counts per consistency status
- `One or multi basic phenotypes` — counts per group
- `Category` — counts per category
- `Growth_tier` — counts per tier

### `results/Hayles_2013_OB_merged_categories.xlsx`

Identical copy, saved to the version-controlled `results/` directory.

---

## Growth Tier Distribution (2026-06-09)

| Tier | Count | % |
|---|---|---|
| 1 (Spores) | 156 | 3.2% |
| 2 (Germinated) | 526 | 10.9% |
| 3 (Microcolonies) | 515 | 10.6% |
| 4 (Small colonies) | 355 | 7.3% |
| 5 (WT) | 3,291 | 68.0% |
| **Total** | **4,843** | **100%** |

---

## Usage

```bash
mamba run -n bioinformatics python src/03_merge_categories.py
```
