# Test Suite — Phenotype Classification Pipeline

## Overview

Seven test files, 117 tests total, all passing:

```
mamba run -n bioinformatics python -m pytest tests/ -v
117 passed in ~3s
```

- `tests/test_growth_signals.py` — 36 tests (end-to-end classification)
- `tests/test_classification.py` — 39 tests (Single/Multiple logic)
- `tests/test_category_revisions.py` — 13 tests (category config loader)
- `tests/test_dataset_merge.py` — 16 tests (loaders + raw 2010/2013 union + Kim Table S3)
- `tests/test_id_mapping.py` — 7 tests (shared PomBase ID update, paralogs, notes)
- `tests/test_phenotype.py` — 2 tests (shared column names + phenotype formatting)
- `tests/test_pipeline_utils.py` — 4 tests (shared helpers)

The archived `01_format_and_update_ids.py` integration test (17 tests) lives with its script in
`tmp/01_format_and_update_ids/tests/`.

The sections below document the classifier cases in detail; the merge-table
tests are summarised in their own section.

---

## `tests/test_growth_signals.py` — End-to-end classification (36 tests)

Each test case is one line: `(description, expected_category, expected_tier)`.
All descriptions are taken verbatim from the Hayles 2013 supplementary table.

**Note:** Tier values in tests are the **signal-based tiers** (1–5) returned
by `classify_growth()`. The final `Growth_tier` in pipeline output is
re-ranked by DIT-HAP DR median (see `04_categorize_phenotypes.py`).

### WT-like (signal tier 5) — no growth signal

| Description | Category | Tier | Note |
|---|---|---|---|
| `VIABLE WT cells at 25,32` | WT-like | 5 | |
| `VIABLE misshapen cells at 25,32` | WT-like | 5 | morphology only |
| `VIABLE slightly long cells at 25,32` | WT-like | 5 | modifier only |
| `VIABLE WT cells, some germination long at 25,32,` | WT-like | 5 | `some` is modifier → `germination` dropped |
| `VIABLE  WT cells, some germination long at 25,32` | WT-like | 5 | double space variant |
| `VIABLE slightly misshapen cells, initially branched septated slightly long, abnormal colony morphology at 25,32` | WT-like | 5 | `initially` is modifier → dropped |

### Spores (signal tier 1)

| Description | Category | Tier | Note |
|---|---|---|---|
| `ESSENTIAL spores  at 25,32` | spores | 1 | double space |
| `ESSENTIAL spores at 25,32` | spores | 1 | |

### Germinated (signal tier 2)

| Description | Category | Tier | Note |
|---|---|---|---|
| `ESSENTIAL germinated spores long at 25,32` | germinated | 2 | |
| `ESSENTIAL germinated spores at 25,32` | germinated | 2 | |
| `ESSENTIAL germinated spores, often divide once at 25,32` | germinated, often divided | 2 | `often divide` → modifier secondary |
| `ESSENTIAL germinated spores, occasionally long and then divide at 25,32` | germinated | 2 | `occasionally` dropped |

### Germinated and divided (signal tier 2)

| Description | Category | Tier | Note |
|---|---|---|---|
| `ESSENTIAL germinated spores slightly misshapen and divide once at 25,32` | germinated and divided | 2 | `divide` in main segment |
| `ESSENTIAL germinated spores slightly misshapen, divide once or twice at 25,32` | germinated and divided | 2 | `divide` not modifier → kept |

### Microcolonies (signal tier 3)

| Description | Category | Tier | Note |
|---|---|---|---|
| `ESSENTIAL microcolonies misshapen cells at 25,32` | microcolonies | 3 | |
| `ESSENTIAL microcolonies skittle cells at 25,32` | microcolonies | 3 | |
| `ESSENTIAL microcolonies slightly misshapen cells, some long cells at 25,32` | microcolonies | 3 | `some long cells` dropped |

### Small colonies (signal tier 4)

| Description | Category | Tier | Note |
|---|---|---|---|
| `VIABLE small colonies slightly misshapen cells at 25,32` | small colonies | 4 | |
| `VIABLE small colonies long cells, possibly diploidising at 25,32` | small colonies | 4 | `possibly` dropped |

### Very small colonies (signal tier 4)

| Description | Category | Tier |
|---|---|---|
| `VIABLE very small colonies rounded cells at 25,32` | very small colonies | 4 |

### Combined — spores + some germinated (signal tier 2)

| Description | Category | Tier | Note |
|---|---|---|---|
| `ESSENTIAL spores, some germinated spores at 25,32` | spores, some germinated | 2 | primary `spores` + secondary `some germinated` |

### Combined — spores + germinated (signal tier 2)

| Description | Category | Tier | Note |
|---|---|---|---|
| `ESSENTIAL spores, germinated spores at 25,32` | spores, germinated | 2 | primary spores before germinated (tier order) |
| `ESSENTIAL spores, germinated spores slightly misshapen may divide once at 25,32` | spores, germinated and divided | 2 | `may` mid-segment → kept |

### Combined — three parallel growth signals (signal tier 3)

| Description | Category | Tier | Note |
|---|---|---|---|
| `ESSENTIAL spores, germinated spores, microcolonies misshapen cells at 25,32` | spores, germinated, microcolonies | 3 | |
| `ESSENTIAL spores, germinated spores, microcolonies long cells, occasionally misshapen branched at 25,32` | spores, germinated, microcolonies | 3 | modifier dropped |
| `ESSENTIAL spores, germinated spores, microcolonies slightly misshapen cells at 25,32` | spores, germinated, microcolonies | 3 | |

### Combined — "occasional" normalised to "occasionally" (signal tier 3)

| Description | Category | Tier | Note |
|---|---|---|---|
| `ESSENTIAL spores, germinated spores, occasional microcolonies of WT/ rounded cells at 25,32` | spores, germinated, occasionally microcolonies | 3 | `occasional` → `occasionally` |
| `ESSENTIAL spores, germinated spores, occasionally microcolonies misshapen cells at 25,32` | spores, germinated, occasionally microcolonies | 3 | |

### Combined — multi-temperature (signal tier 3)

| Description | Category | Tier | Note |
|---|---|---|---|
| `ESSENTIAL spores, germinated spores, microcolonies slightly misshapen cells at 32, spores, germinated spores at 25` | spores, germinated, microcolonies | 3 | |

### Combined — other

| Description | Category | Tier |
|---|---|---|
| `ESSENTIAL misshapen germinated spores and microcolonies misshapen cells at 32` | germinated, microcolonies | 3 |
| `ESSENTIAL spores, germinated spores, small colonies long cells at 25,32` | spores, germinated, small colonies | 4 |
| `ESSENTIAL microcolonies skittle cells, small colonies WT cells at 25,32` | microcolonies, small colonies | 4 |
| `ESSENTIAL spores, microcolonies misshapen cells at 25,32` | spores, microcolonies | 3 |

---

## `tests/test_classification.py` — Single/Multiple logic (39 tests)

Tests `filter_primary_segments()`, `count_growth_segments()`, and
`classify_phenotype_count()` — the functions that decide whether a
description has one or multiple parallel growth phenotypes.

### `filter_primary_segments()` — which comma segments are kept?

| Input | After filtering | Reason |
|---|---|---|
| `VIABLE WT cells, some germination long` | `VIABLE WT cells` | `some` is modifier → dropped |
| `ESSENTIAL spores, germinated spores, some division` | `ESSENTIAL spores, germinated spores` | `some division` dropped |
| `spores, germinated spores, microcolonies long cells, occasionally misshapen branched` | `spores, germinated spores, microcolonies long cells` | modifier dropped |
| `ESSENTIAL germinated spores, slightly misshapen cells` | `ESSENTIAL germinated spores` | modifier dropped |
| `ESSENTIAL germinated spores very large swollen branched, initially germinated spores long` | `ESSENTIAL germinated spores very large swollen branched` | modifier dropped |
| `ESSENTIAL spores, germinated spores, microcolonies slightly misshapen cells` | unchanged | **KEPT** — `microcolonies` starts segment |
| `ESSENTIAL germinated spores, divide once` | unchanged | **KEPT** — `divide` not modifier |
| `VIABLE cells, occasionally long, slightly misshapen` | `VIABLE cells` | all others modifier-led |
| `initially rounded, slightly misshapen` | `initially rounded` | fallback to first segment |
| `VIABLE small colonies long cells, possibly diploidising` | `VIABLE small colonies long cells` | modifier dropped |
| `ESSENTIAL germinated spores, often divide once` | `ESSENTIAL germinated spores` | modifier dropped |

### `count_growth_segments()`

| Description | Count | Note |
|---|---|---|
| `ESSENTIAL germinated spores` | 1 | no comma, 1 signal |
| `VIABLE WT cells` | 0 | no comma, 0 signals |
| `spores, germinated spores` | 2 | two parallel |
| `spores, germinated spores, microcolonies, occasionally misshapen` | 3 | modifier not counted |
| `germinated spores slightly misshapen, some division` | 1 | modifier not counted |
| `germinated spores, divide once` | 2 | `divide` not modifier |
| `np.nan` | 0 | non-string |
| `ESSENTIAL spores, germinated spores, often divide once` | 2 | modifier not counted |

### `classify_phenotype_count()` — Single vs Multiple

| Description | Result | Note |
|---|---|---|
| `ESSENTIAL germinated spores` | Single | 1 signal |
| `spores, germinated spores` | Multiple | 2 signals |
| `spores, germinated spores, microcolonies, occasionally misshapen` | Multiple | 3 signals |
| `ESSENTIAL germinated spores slightly misshapen, some division` | Single | 1 signal |
| `ESSENTIAL spores, germinated spores slightly misshapen, some division` | Multiple | 2 signals |
| `ESSENTIAL microcolonies skittle cells, small colonies WT cells` | Multiple | 2 signals |
| `VIABLE WT cells, some germination long` | Single | 0 signals |
| `ESSENTIAL germinated spores, divide once` | Multiple | 2 signals |
| `VIABLE WT cells` | Single | 0 signals |
| `VIABLE misshapen cells` | Single | 0 signals |
| `ESSENTIAL spores` | Single | 1 signal |
| `ESSENTIAL microcolonies` | Single | 1 signal |
| Pipeline: `VIABLE WT cells, some germination long` → filter → classify | WT-like (5) | integration test |

### Constant integrity checks

| Check | Expected |
|---|---|
| MODIFIER_WORDS has no duplicates | True |
| GROWTH_KEYWORDS has no duplicates | True |
| All MODIFIER_WORDS lowercase | True |
| All GROWTH_KEYWORDS lowercase | True |

---

## `tests/test_category_revisions.py` — category config loader (13 tests)

Tests `load_category_config()` from `src/category_revisions.py`, the shared
loader for the category JSON (`sub_category_order` / `revisions` /
`category_order`).

| Test | Expected |
|---|---|
| Default path points at JSON | `data/4_categorized_genes/category_revisions.json` |
| Loads config | 30 sub-categories, 25 revisions, 7 categories |
| Sub-category order bounds | starts `spores`, ends `WT-like` |
| Category order | exact list; includes merge target `divided` (regression) |
| Identity entries are not revisions | `spores` in order but not in revisions |
| Spot checks | e.g. `germinated and divided` → `divided` |
| No chained merges | a merge target is not itself merged |
| Missing file | raises `FileNotFoundError` |
| Non-object JSON | raises `ValueError` |
| Missing section | raises `ValueError` |
| Non-string values | raises `ValueError` |
| Duplicate order entry | raises `ValueError` |
| Non-idempotent revisions | raises `ValueError` |

---

## `tests/test_dataset_merge.py` — loaders + 2010/2013 union table (16 tests)

Loader/normalisation tests pin the input parsing; synthetic-frame tests cover
every `Presence` value and the consistency/change-class logic; one real-data
test checks the actual tables.

| Test | What it pins down |
|---|---|
| `test_normalize_id_strips_all_whitespace` | whitespace removal and missing-value handling |
| `test_normalize_call_maps_words_and_stray_spaces` | `viable`/`essential` and stray spaces → `V`/`E` |
| `test_reconstruct_gene_name_from_date` | Excel date cells → `Oct-1`, `Sep-11`, `Sep-1` |
| `test_parse_table_s3_repairs_truncated_id` | Table S3 `SPAC13D6.0` → `SPAC13D6.05` |
| `test_merged_presence_column` | `both` / `2010_only` / `2013_only` / `neither` assignment |
| `test_merged_consistency_column` | `consistent` when analysis == study; `inconsistent` otherwise; `not_comparable` for one-sided rows |
| `test_merged_change_classes_and_discrepancy` | `no_change`, `documented_change`, `matches_kim2010_deletion_not_analysis`, `kim2010_deletion_analysis_discrepancy`, `hayles2013_addition`, `no_hayles2013_row`, `no_data_either`; discrepancy text and Table S3 comment |
| `test_merged_gene_name_prefers_2013_then_falls_back_to_2010` | canonical `Gene name` selection |
| `test_merged_keeps_source_columns` | 2010 columns survive under `Kim2010_*` |
| `test_merged_without_table_s3_marks_changed_rows_accordingly` | graceful degradation to `changed_table_s3_not_available` |
| `test_real_inputs_build_expected_union` | raw merge before 01's ID collapse: 4,914 rows; 4,836 constructed; 4,843 Hayles rows; the three repaired Excel-date names; presence 4,825 / 18 / 11 / 60; consistency 4,796 / 29 / 89; change-class counts (Table S3 counts only when the docx is present) |
| `test_kim_table_s3_matches_table_s1` | the 40-row Kim Table S3 TSV: numbering, calls, criteria counts and footnote rows match Table S1 |
| `test_kim_table_s3_validation_rejects_transcription_errors` | a dropped row, a non-disagreement ID, criteria drift and a moved footnote all raise |

---

## `tests/test_id_mapping.py` — shared ID update (7 tests)

| Test | What it pins down |
|---|---|
| `test_normalize_sysid_case_and_suffix` | case normalisation; `SPAC1556.06.1` → `SPAC1556.06`; non-strings pass through |
| `test_build_update_maps_layers` | protein-coding filter; name → ID and synonym → ID maps |
| `test_update_sysids_all_layers` | direct match, name update, synonym update, reclassification, not found |
| `test_collect_mapping_ids_deduplicates_main_and_paralog` | main IDs, Hayles IDs and Kim paralogs are unique and keep first-seen order |
| `test_build_paralog_columns_updates_and_drops_self` | current IDs, current names, self-reference removal, `being curated` |
| `test_combine_mapping_notes_single_and_collapsed` | single note and the `COLLAPSED: ...` format |
| `test_read_pombase_annotation_strips_whitespace` | TSV loading and whitespace stripping |

## `tests/test_phenotype.py` — shared formatting (2 tests)

| Test | What it pins down |
|---|---|
| `test_column_names` | `Phenotype_description`, `Morphological_classification`, `Gene_dispensability` |
| `test_format_hayles_phenotypes_normalises_text` | whitespace, trailing periods, and `25, 32` / `25 32` / `32, 25` → `25,32` |

---

## Running the tests

```bash
# Full suite
mamba run -n bioinformatics python -m pytest tests/ -v

# End-to-end classification only
mamba run -n bioinformatics python -m pytest tests/test_growth_signals.py -v

# Single/Multiple logic only
mamba run -n bioinformatics python -m pytest tests/test_classification.py -v

# Revision mapping loader only
mamba run -n bioinformatics python -m pytest tests/test_category_revisions.py -v

# 2010/2013 union table only
mamba run -n bioinformatics python -m pytest tests/test_dataset_merge.py -v

# Shared ID mapping only
mamba run -n bioinformatics python -m pytest tests/test_id_mapping.py -v

# Shared phenotype formatting only
mamba run -n bioinformatics python -m pytest tests/test_phenotype.py -v

# Stop on first failure
mamba run -n bioinformatics python -m pytest tests/ -x
```
