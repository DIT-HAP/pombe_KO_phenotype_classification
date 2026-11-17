# 05_merge_categories.py — Merge Categorized Phenotypes

---

## Purpose

Concatenates the three phenotype branches (one / multi / inconsistent) from
the categorized output into a single table and applies manually revised
category merges.

Step 04 already writes the fine-grained classification as `Sub_category`,
while the revised/merged classification is stored in `Category`.

---

## Data Flow

```
data/4_categorized_genes/                Categorized output from step 04
  Hayles_2013_OB_categorized_
  phenotypes.xlsx                        (3 data sheets + All genes)
                             ↓
data/4_categorized_genes/                Hand-curated category config
  category_revisions.json                (sub-category order, merges,
                                         category order; shared with step 06)
                             ↓
                   05_merge_categories.py
                             ↓
data/5_merged_categories/                Final merged table + summaries
  Hayles_2013_OB_merged_categories.xlsx
                             ↓
results/                                 Flat copy of All genes
  Hayles_2013_OB_merged_categories.xlsx
```

---

## Processing Logic

1. **Load and concat** — merge the three branch sheets into one DataFrame
2. **Apply category merges** — map the step-04 `Sub_category` through
   `category_revisions.json` (`Sub_category → Category`). Sub-categories not
   listed in the JSON keep themselves as their merged `Category`. The mapping
   is loaded by the shared `src/category_revisions.py` helper.

---

## Output

### `data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx`

**Data sheet:** `All genes` — 4,843 rows.

Key columns:
- `Sub_category` — fine-grained category from step 04 (e.g., `germinated, often divided`)
- `Category` — merged category after manual revision (e.g., `germinated, divided or microcolonies`)
- `Phenotype_count` — Single / Multiple / Temp_mismatch
- `Consistency_25_32` — Consistent / Only_32 / Mismatch

**Summary sheets (4):**
- `Consistency_25_32` — counts per consistency status
- `Phenotype_count` — counts per Single/Multiple/Temp_mismatch
- `Category` — counts per merged category
- `Sub_category` — counts per fine-grained sub-category

### `results/Hayles_2013_OB_merged_categories.xlsx`

Flat single-sheet copy of the `All genes` sheet, saved to the
version-controlled `results/` directory.

---

## Usage

```bash
mamba run -n bioinformatics python scripts/05_merge_categories.py
```
