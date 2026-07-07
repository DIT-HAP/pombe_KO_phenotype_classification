# Test Suite — Phenotype Classification Pipeline

## Overview

Two test files, 75 tests total, all passing:

```
mamba run -n bioinformatics python -m pytest tests/ -v
75 passed in 0.12s
```

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

## Running the tests

```bash
# Full suite
mamba run -n bioinformatics python -m pytest tests/ -v

# End-to-end classification only
mamba run -n bioinformatics python -m pytest tests/test_growth_signals.py -v

# Single/Multiple logic only
mamba run -n bioinformatics python -m pytest tests/test_classification.py -v

# Stop on first failure
mamba run -n bioinformatics python -m pytest tests/ -x
```
