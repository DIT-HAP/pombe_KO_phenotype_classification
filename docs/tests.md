# Test Suite — Phenotype Classification Pipeline

This document describes the test suite in `tests/`. The tests cover the
three critical layers of the pipeline:

1. **Signal engine** — `growth_signals.classify_growth()`
2. **Segment filtering** — `growth_signals.filter_primary_segments()`
3. **Single/Multiple classification** — `classify_phenotype_count()`

All 77 tests pass:

```
mamba run -n bioinformatics python -m pytest tests/ -v
77 passed in 0.21s
```

---

## File: `tests/test_growth_signals.py`

Tests for `classify_growth()` as a pure function — given a description text,
what `(category, tier)` pair is returned?

### No signal → WT-like (tier 5)

| Description | Category | Tier |
|---|---|---|
| `VIABLE WT cells` | WT-like | 5 |
| ` ` (empty string) | WT-like | 5 |
| `VIABLE misshapen cells at 25,32` | WT-like | 5 |
| `VIABLE slightly long cells` | WT-like | 5 |
| `VIABLE cells at 25,32` | WT-like | 5 |

### Single signals

| Description | Category | Tier |
|---|---|---|
| `ESSENTIAL spores` | spores | 1 |
| `ESSENTIAL spores misshapen cells` | spores | 1 |
| `ESSENTIAL spores at 25,32` | spores | 1 |
| `ESSENTIAL germinated spores` | germinated | 2 |
| `ESSENTIAL germinated spores long` | germinated | 2 |
| `VIABLE germination at 25` | germinated | 2 |
| `VIABLE WT cells, some germination long` | germinated | 2 |
| `ESSENTIAL germinated spores divide` | germinated and divided | 2 |
| `ESSENTIAL germinated spores division` | germinated and divided | 2 |
| `germinated spores divides once` | germinated and divided | 2 |
| `germinated spores, division` | germinated and divided | 2 |
| `ESSENTIAL microcolonies` | microcolonies | 3 |
| `ESSENTIAL microcolonies misshapen cells` | microcolonies | 3 |
| `VIABLE small colonies` | small colonies | 4 |
| `VIABLE small colonies long cells` | small colonies | 4 |
| `VIABLE very small colonies` | very small colonies | 4 |
| `VIABLE very small colonies misshapen cells` | very small colonies | 4 |

### Combined / mixed signals

| Description | Category | Tier |
|---|---|---|
| `spores, germinated spores` | germinated, spores | 2 |
| `spores, germinated spores, microcolonies` | germinated, microcolonies, spores | 3 |
| `spores, germinated spores, microcolonies long cells, occasionally misshapen branched` | germinated, microcolonies, spores | 3 |
| `ESSENTIAL spores, germinated spores` | germinated, spores | 2 |
| `germinated spores, microcolonies` | germinated, microcolonies | 3 |

### Standalone-spores logic

The rule: if `spores` appears only as part of "germinated spores" (not as a
standalone word), it is implied and not listed separately.

| Description | Category | Tier | Reason |
|---|---|---|---|
| `germinated spores long` | germinated | 2 | `spores` only inside "germinated spores" → implied |
| `spores, germinated spores` | germinated, spores | 2 | standalone `spores` exists → mixed population |

### Cross-cutting

| Description | Category | Tier |
|---|---|---|
| `VIABLE small colonies long cells` | small colonies | 4 |
| `ESSENTIAL microcolonies misshapen cells` | microcolonies | 3 |

### Case / whitespace robustness

| Description | Category | Tier |
|---|---|---|
| `Essential Spores` (mixed case) | spores | 1 |
| `  ESSENTIAL spores  ` (leading/trailing spaces) | spores | 1 |
| `germinated spores` | germinated | 2 |
| `germinated, microcolonies` | germinated, microcolonies | 3 |
| `spores only here` | spores | 1 |
| `small colonies` | small colonies | 4 |
| `WT-like` | WT-like | 5 |

---

## File: `tests/test_classification.py`

Tests for three functions: `filter_primary_segments()`,
`count_growth_segments()`, and `classify_phenotype_count()`.

### `filter_primary_segments()` — which comma segments are kept?

The rule: a comma segment that **starts with** a modifier word is a secondary
description and is dropped. All other segments are kept.

| Input | After filtering | Reason |
|---|---|---|
| `VIABLE WT cells, some germination long` | `VIABLE WT cells` | `some` is a modifier → segment dropped |
| `ESSENTIAL spores, germinated spores, some division` | `ESSENTIAL spores, germinated spores` | `some division` dropped |
| `spores, germinated spores, microcolonies long cells, occasionally misshapen branched` | `spores, germinated spores, microcolonies long cells` | `occasionally misshapen` dropped |
| `ESSENTIAL germinated spores, slightly misshapen cells` | `ESSENTIAL germinated spores` | `slightly misshapen` dropped |
| `ESSENTIAL germinated spores very large swollen branched, initially germinated spores long` | `ESSENTIAL germinated spores very large swollen branched` | `initially germinated` dropped |
| `ESSENTIAL spores, germinated spores, microcolonies slightly misshapen cells` | `ESSENTIAL spores, germinated spores, microcolonies slightly misshapen cells` | **KEPT** — `microcolonies` starts the segment, not `slightly` |
| `ESSENTIAL germinated spores, divide once` | `ESSENTIAL germinated spores, divide once` | **KEPT** — `divide` is not a modifier |
| `VIABLE cells, occasionally long, slightly misshapen` | `VIABLE cells` | All other segments start with modifiers → fallback to first |
| `initially rounded, slightly misshapen` | `initially rounded` | Both start with modifiers → fallback to first |
| `VIABLE small colonies long cells, possibly diploidising` | `VIABLE small colonies long cells` | `possibly diploidising` dropped |
| `ESSENTIAL germinated spores, often divide once` | `ESSENTIAL germinated spores` | `often divide once` dropped |
| `ESSENTIAL spores` (no comma) | `ESSENTIAL spores` (pass-through) | No comma → unchanged |
| `123` (non-string) | `123` (pass-through) | Non-string → unchanged |

### `count_growth_segments()` — growth signal segments

| Description | Count | Notes |
|---|---|---|
| `ESSENTIAL germinated spores` | 1 | No comma, 1 signal |
| `VIABLE WT cells` | 0 | No comma, 0 signals |
| `spores, germinated spores` | 2 | Two parallel segments |
| `spores, germinated spores, microcolonies, occasionally misshapen` | 3 | `occasionally misshapen` → not counted (modifier-led) |
| `germinated spores slightly misshapen, some division` | 1 | `some division` → not counted (modifier-led) |
| `germinated spores, divide once` | 2 | `divide` is not a modifier → counted |
| `np.nan` | 0 | Non-string → 0 |
| `ESSENTIAL spores, germinated spores, often divide once` | 2 | `often divide once` → not counted (modifier-led) |

### `classify_phenotype_count()` — Single vs Multiple

| Description | Result | Reason |
|---|---|---|
| `ESSENTIAL germinated spores` (no comma) | Single | 1 signal |
| `spores, germinated spores` | Multiple | 2 signals |
| `spores, germinated spores, microcolonies, occasionally misshapen` | Multiple | 3 signals |
| `ESSENTIAL germinated spores slightly misshapen, some division` | Single | 1 signal (`some division` → skipped) |
| `ESSENTIAL spores, germinated spores slightly misshapen, some division` | Multiple | 2 signals (spores + germinated) |
| `ESSENTIAL microcolonies skittle cells, small colonies WT cells` | Multiple | 2 signals |
| `VIABLE WT cells, some germination long` | Single | 0 signals |
| `ESSENTIAL germinated spores, divide once` | Multiple | 2 signals |
| `VIABLE WT cells` (no comma) | Single | 0 signals |
| `VIABLE misshapen cells` (no comma, only morphology) | Single | 0 signals |
| `ESSENTIAL spores` (no comma) | Single | 1 signal |
| `ESSENTIAL microcolonies` (no comma) | Single | 1 signal |

### Constant integrity checks

| Check | Expected |
|---|---|
| MODIFIER_WORDS has no duplicates | True |
| GROWTH_KEYWORDS has no duplicates | True |
| All MODIFIER_WORDS entries are lowercase | True |
| All GROWTH_KEYWORDS entries are lowercase | True |

---

## Running the tests

```bash
# Full suite
mamba run -n bioinformatics python -m pytest tests/ -v

# Single file
mamba run -n bioinformatics python -m pytest tests/test_growth_signals.py -v

# Single test class
mamba run -n bioinformatics python -m pytest tests/ -k "TestClassifyPhenotypeCount" -v

# Stop on first failure
mamba run -n bioinformatics python -m pytest tests/ -x
```