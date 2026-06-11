# growth_signals.py — Shared Growth Signal Detection Engine

---

## Purpose

Defines the canonical `GrowthSignal` dataclass, the signal table, and the
`classify_growth()` engine that replaces the fragile if-elif chains in
the original categorize scripts.

This module is not a standalone script — it is imported by
`02_categorize_phenotypes.py`.

---

## Growth Tier Definitions

| Tier | Label | Description | Examples |
|---|---|---|---|
| 1 | Spores | Spores visible, fail to germinate | `spores` |
| 2 | Germinated | Germinated, limited or no division | `germinated spores`, `germinated and then divide` |
| 3 | Microcolonies | Microscopic colonies, severely limited growth | `microcolonies` |
| 4 | Small colonies | Visible but smaller than WT | `small colonies`, `very small colonies` |
| 5 | WT | Normal growth | `WT cells` |

## Detection Logic

For a phenotype description, `classify_growth()`:

1. Scans for each keyword in `GROWTH_SIGNALS` (case-insensitive substring
   match)
2. Collects all matching signals
3. Composes category names (sorted alphabetically, deduplicated)
4. Removes redundant broader categories when a more specific one is present
   (e.g., `germinated` is implied by `germinated and divided`; `small
   colonies` is implied by `very small colonies`)
5. Assigns the **highest** tier among matched signals — this represents the
   most advanced growth stage the gene can support
6. Returns `(composed_category_name, growth_tier)`
7. If no signals are found, returns `("WT", 5)`

### Examples

| Description | Category | Tier |
|---|---|---|
| `VIABLE WT cells at 25,32` | WT | 5 |
| `ESSENTIAL microcolonies misshapen at 32` | microcolonies | 3 |
| `ESSENTIAL germinated spores at 25,32` | germinated, spores | 2 |
| `ESSENTIAL germinated spores and then divide` | germinated and divided, spores | 2 |
| `spores` | spores | 1 |
| `VIABLE very small colonies at 32` | very small colonies | 4 |

---

## Signal Table (7 entries)

| Keyword | Category | Tier | Notes |
|---|---|---|---|
| `very small colon` | very small colonies | 4 | Checked before `small colon` |
| `small colon` | small colonies | 4 | Only if `very small` not present |
| `microcolonies` | microcolonies | 3 | — |
| `germinated` | germinated | 2 | Deduped if `germinated and divided` present |
| `divide` | germinated and divided | 2 | — |
| `division` | germinated and divided | 2 | — |
| `spores` | spores | 1 | — |

---

## Dependencies

- Python 3.12+ standard library only (`dataclasses`, `re`)
