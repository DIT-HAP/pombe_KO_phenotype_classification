# AGENTS.md

Phenotype classification pipeline for the Hayles 2013 *S. pombe* deletion
mutant screen. Raw natural-language microcopy → structured growth categories →
DIT-HAP/gRNA depletion-rate validation.

## Environment

- Everything runs in the `bioinformatics` conda/mamba env:
  `mamba run -n bioinformatics python ...`
- No `requirements.txt` / `pyproject.toml`. Runtime deps: `pandas`, `numpy`,
  `loguru`, `matplotlib`, `openpyxl` (xlsx I/O), `scipy`, `pytest`.
- Python 3.12+.

## Running scripts

- **Always run from the repo root.** All scripts use hardcoded relative paths
  (`data/...`, `results/...`) as argparse defaults; running from `scripts/` breaks
  them.
- Scripts are numbered, not auto-chained. Real dependency order:

  `00 → 01 → 03 → 04 → 05 → {02, 06}`

  `02_keyword_profile.py` is a **re-analysis of the final merged output** (it
  reads `results/Hayles_2013_OB_merged_categories.xlsx`), so despite its number
  it runs *after* `05`, not before `03`. `06_plot_dr_distribution.py` also runs
  after `05`.
- `05` writes the merged table to both `data/5_merged_categories/` and a
  committed copy in `results/`. `06` writes each figure as both PNG and PDF to
  `results/`.

## Testing

```bash
mamba run -n bioinformatics python -m pytest tests/ -v            # full suite (86 tests, ~0.2s)
mamba run -n bioinformatics python -m pytest tests/test_growth_signals.py -v
mamba run -n bioinformatics python -m pytest tests/test_classification.py -v
mamba run -n bioinformatics python -m pytest tests/ -k <pattern>  # single case
```

- `tests/conftest.py` inserts `scripts/` into `sys.path`, so tests import
  `growth_signals` directly (`from growth_signals import classify_growth`).
- Test files are the spec for the classifier. `docs/tests.md` documents every
  case and its expected category/tier.

## Architecture

- `scripts/growth_signals.py` — shared engine (`classify_growth()`); imported by
  `03`/`04` (`02` imports only its modifier groups).
- `scripts/pipeline_utils.py` — shared `setup_logger`, `load_dit_hap`,
  `category_dr_medians`, and the `DEFAULT_RELEASE` constant.
- `scripts/category_revisions.py` — shared loader for
  `data/4_categorized_genes/category_revisions.json`, a single ordered
  `name → Category` map (key order = draw order; `key != value` are the merges).
  Imported by `05` and `06`.
- **Two different "tiers"** — don't conflate them:
  - *signal tier* (1–5): fixed biological hierarchy, internal to
    `classify_growth()`.
  - `Growth_tier`: output column, re-ranked by median DIT-HAP DR per category
    in `04`/`05`. Test expectations use signal tiers.
- `germination` (noun) is intentionally **not** a growth signal — it appears
  only in morphological contexts; the adjective `germinated` is the signal.
- Plot draw order and the `Sub_category → Category` merge both live in
  `data/4_categorized_genes/category_revisions.json`, a single ordered
  `name → Category` map. Edit that JSON to reorder the figures or change merges;
  there is no hardcoded order constant anymore.

## Data gotchas

- **`data/references/` is gitignored** (`.gitignore` line `references/`) — the
  reference TSVs are required inputs but are *not* in the repo. Missing files
  must be obtained separately:
  - `pombase-2026-06-01_gene_IDs_names_products.tsv` (produced by `00`)
  - `all_coding_genes_with_DIT_HAP_clustering.tsv` (`04`/`05`/`06`)
  - `260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv` (`06`)
  - `DIT_HAP.mplstyle` (`06`)
  The one exception is the hand-curated
  `previous_manual_check_of_insistent_phenotypes/` subfolder, re-included via
  `.gitignore` and tracked.
- Hand-curated, do **not** regenerate or overwrite:
  - `data/references/previous_manual_check_of_insistent_phenotypes/Inconsistent_phenotypes_at_25_32_manual.xlsx`
  - `data/4_categorized_genes/category_revisions.json` — single source of truth for
    the plot draw order and the `Sub_category → Category` merge used by `05`/`06`.
    One ordered `name → Category` map (key order = draw order; `key != value` are
    the merges). Edit this JSON directly. The old
    `Hayles_2013_OB_inspection_phenotypes_category_revised_20260707.xlsx` is
    frozen; no script reads it anymore.
- `data/arc/` and `scripts/arc/` are the **deprecated previous implementation**,
  not the live pipeline. `scripts/arc/script_flow.sh` references deleted filenames
  — ignore it.

## Docs

Per-step design notes live in `docs/` (`00_`, `01_`, `04_`, `05_`,
`growth_signals.md`, `tests.md`). `README.md` has the paper context and concept
glossary (`Sub_category` vs `Category` vs `Growth_tier`).
