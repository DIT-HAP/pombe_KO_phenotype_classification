# growth_signals.py — Shared Growth Signal Detection Engine

---

## Purpose

Defines the canonical `GrowthSignal` dataclass, the signal table, and the
`classify_growth()` engine that replaces the fragile if-elif chains in
the original categorize scripts.

This module is not a standalone script — it is imported by
`04_categorize_phenotypes.py`, `03_group_genes.py`, and `02_keyword_profile.py`.

---

## Signal-Based Tier Hierarchy

The signal table defines a biological hierarchy of growth defect severity:

| Signal tier | Label | Description | Examples |
|---|---|---|---|
| 1 | Spores | Spores visible, fail to germinate | `spores` |
| 2 | Germinated | Germinated, limited or no division | `germinated`, `divide`, `division` |
| 3 | Microcolonies | Microscopic colonies, severely limited growth | `microcolonies` |
| 4 | Small colonies | Visible but smaller than WT | `small colonies`, `very small colonies` |
| 5 | WT-like | Normal growth, no defect signals | `WT cells` |

**Note:** The signal tier (1–5) is used internally by `classify_growth()` to
determine the biological hierarchy. The final `Growth_tier` in the output
tables is re-ranked by DIT-HAP DR median (see `04_categorize_phenotypes.py`
and `05_merge_categories.py`).

---

## Detection Logic

`classify_growth()` processes a phenotype description as follows:

1. Splits the description into comma-separated segments
2. For each segment, detects all growth signals (case-insensitive substring)
3. Segments starting with a modifier word (e.g. `some`, `often`,
   `occasionally`) are treated as secondary — their signals are preserved
   with a modifier prefix (e.g. `some germinated`, `often divided`)
4. Primary segments (non-modifier-led) contribute plain signal names
5. Sorts categories: primary first (by signal tier ascending), then
   modifier-led. `germinated` sorts before `germinated and divided` within
   the same tier
6. When all segments are modifier-led, the most severe (lowest signal tier)
   is promoted to primary
7. Removes redundant broader categories (e.g. `germinated` is implied by
   `germinated and divided`; `small colonies` by `very small colonies`;
   `spores` by `germinated` when no standalone spores)
8. If no signals are found, returns `("WT-like", 5)`

### Modifier handling

Modifier words are defined in `MODIFIER_STARTS` and include frequency
(`some`, `often`, `occasionally`, `possibly`, `may`, etc.), degree
(`slightly`, `very`, `weak`), and quantity modifiers. The canonical form
normalises variants (e.g. `occasional` → `occasionally`).

### Examples

| Description | Category | Signal tier |
|---|---|---|
| `VIABLE WT cells at 25,32` | WT-like | 5 |
| `ESSENTIAL spores at 25,32` | spores | 1 |
| `ESSENTIAL germinated spores at 25,32` | germinated | 2 |
| `ESSENTIAL spores, germinated spores at 25,32` | spores, germinated | 2 |
| `ESSENTIAL spores, some germinated spores at 25,32` | spores, some germinated | 2 |
| `ESSENTIAL germinated spores, often divide at 25,32` | germinated, often divided | 2 |
| `ESSENTIAL microcolonies misshapen at 32` | microcolonies | 3 |
| `ESSENTIAL germinated spores and then divide` | germinated and divided | 2 |
| `VIABLE very small colonies at 32` | very small colonies | 4 |
| `ESSENTIAL microcolonies, occasionally spores, occasionally germinated` | microcolonies, occasionally spores, occasionally germinated | 3 |

---

## Signal Table (6 entries)

| Keyword | Category | Signal tier | Notes |
|---|---|---|---|
| `very small colon` | very small colonies | 4 | Checked before `small colon` |
| `small colon` | small colonies | 4 | Only if `very small` not present |
| `microcolonies` | microcolonies | 3 | — |
| `germinated` | germinated | 2 | Deduped if `germinated and divided` present |
| `divide` | germinated and divided | 2 | — |
| `division` | germinated and divided | 2 | — |
| `spores` | spores | 1 | — |

**Note:** `germination` (noun) is intentionally NOT a signal — in the data it
only appears in morphological contexts (`germination long`, `poor germination`),
not as a growth state. The adjective `germinated` covers the growth signal.

---

## Other exported functions

- `filter_primary_segments(description)` — strips modifier-led segments,
  returns the remaining primary description. Used by test suite.
- `count_growth_segments(description)` — counts segments containing growth
  signals (after filtering modifiers). Used by `03_group_genes.py` for
  Single/Multiple classification.
- `classify_growth_batch(descriptions)` — batch wrapper around
  `classify_growth()`.

---

## Dependencies

- Python 3.12+ standard library only (`dataclasses`)
