# 05_merge_categories.py — Merge Categorized Phenotypes

---

## Purpose

Concatenates the three phenotype branches (one / multi / inconsistent) from
the categorized output into a single table, applies manually revised category
merges, and assigns `Growth_tier` based on DIT-HAP DR median per merged
category.

The fine-grained classification from step 04 is preserved as `Sub_category`,
while the revised/merged classification is stored in `Category`.

---

## Data Flow

```
data/4_categorized_genes/                Categorized output from step 04
  Hayles_2013_OB_categorized_
  phenotypes.xlsx                        (3 data sheets + All genes)
                             ↓
data/4_categorized_genes/                Hand-curated category config
  category_revisions.json                (plot order + Sub_category → Category
                                         merges; shared with step 06)
                             ↓
data/references/                         DIT-HAP DR values for tier ranking
  all_coding_genes_with_DIT_HAP_clustering.tsv
                             ↓
                   05_merge_categories.py
                             ↓
data/5_merged_categories/                Final merged table + summaries
  Hayles_2013_OB_merged_categories.xlsx
                             ↓
results/                                 Copy for source control
  Hayles_2013_OB_merged_categories.xlsx
```

---

## Processing Logic

1. **Load and concat** — merge the three branch sheets into one DataFrame
2. **Apply category merges** — rename the step-04 `Category` to
   `Sub_category`, then map it through `category_revisions.json`
   (`Sub_category → Category`). Sub-categories not listed in the JSON keep
   themselves as their merged `Category`. The mapping is loaded by the shared
   `scripts/category_revisions.py` helper.
3. **Re-rank Growth_tier** — compute the median DIT-HAP DR per merged
   `Category`, rank from highest (tier 1) to lowest (tier N).

---

## Output

### `data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx`

**Data sheet:** `All genes` — 4,843 rows.

Key columns:
- `Sub_category` — fine-grained category from step 04 (e.g., `germinated, often divided`)
- `Category` — merged category after manual revision (e.g., `germinated, divided or microcolonies`)
- `Growth_tier` — DR-median-ranked tier (1 = highest DR median, N = lowest)
- `Phenotype_count` — Single / Multiple / Temp_mismatch
- `Consistency_25_32` — Consistent / Only_32 / Mismatch

**Summary sheets (5):**
- `Consistency_25_32` — counts per consistency status
- `Phenotype_count` — counts per Single/Multiple/Temp_mismatch
- `Category` — counts per merged category
- `Sub_category` — counts per fine-grained sub-category
- `Growth_tier` — counts per tier

### `results/Hayles_2013_OB_merged_categories.xlsx`

Identical copy, saved to the version-controlled `results/` directory.

---

## Usage

```bash
mamba run -n bioinformatics python scripts/05_merge_categories.py
```
