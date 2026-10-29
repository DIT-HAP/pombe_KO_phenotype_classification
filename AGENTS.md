# AGENTS.md

Phenotype classification pipeline for the Hayles 2013 *S. pombe* deletion
mutant screen. Raw natural-language microcopy → structured growth categories →
DIT-HAP/gRNA depletion-rate validation.

## Environment

- Steps `00`–`05` run in the `bioinformatics` conda/mamba env:
  `mamba run -n bioinformatics python ...`
- `06`/`07` also need `cnsplots`. On this machine it is installed in the
  separate `cnsplots` env: `mamba run -n cnsplots python scripts/06_...py`.
- No `requirements.txt` / `pyproject.toml`. Runtime deps: `pandas`, `numpy`,
  `loguru`, `matplotlib`, `openpyxl` (xlsx I/O), `xlrd` (xls I/O for the
  Kim 2010 supplementary table), `scipy`, `pytest`.
- Python 3.12+.

## Running scripts

- **Always run from the repo root.** All scripts use hardcoded relative paths
  (`data/...`, `results/...`) as argparse defaults; running from `scripts/` breaks
  them.
- Scripts are numbered, not auto-chained. Real dependency order:

  `00 → 01b → 03 → 04 → 05 → {02, 06, 07}`

  `01b` is the first phenotype step. It updates the Kim and Hayles IDs through
  `src/id_mapping.py`, takes Kim's `Paralog` only, and writes three files: the
  Kim paralog pre-merge (`Kim2010_paralog_premerge.xlsx`, 4,913 rows), the
  Hayles pre-merge (`Hayles2013_phenotypes_premerge.xlsx`, 4,843 rows) and the
  12-column Hayles phenotype table (`Hayles2013_phenotypes_annotated.xlsx`,
  4,843 rows). The last file feeds `03`.

  The old `01_format_and_update_ids.py` and its six 25-column/union outputs
  are archived in `tmp/01_format_and_update_ids/`; it is not part of the live
  pipeline. `src/dataset_merge.py` still contains the union builder used by
  the archived script and its tests.

  `02_keyword_profile.py` is a **re-analysis of the final merged output** (it
  reads `results/Hayles_2013_OB_merged_categories.xlsx`), so despite its number
  it runs *after* `05`, not before `03`. `06_plot_dr_distribution.py` and
  `07_plot_dr_distribution_pure.py` also run after `05`.
- `05` writes the merged table to both `data/5_merged_categories/` and a
  committed copy in `results/`. `06`/`07` write each figure as both PNG and
  PDF to `results/`.

## Testing

```bash
mamba run -n bioinformatics python -m pytest tests/ -v            # full suite (117 tests, ~3s)
mamba run -n bioinformatics python -m pytest tests/test_growth_signals.py -v
mamba run -n bioinformatics python -m pytest tests/test_classification.py -v
mamba run -n bioinformatics python -m pytest tests/test_dataset_merge.py -v
mamba run -n bioinformatics python -m pytest tests/test_id_mapping.py -v
mamba run -n bioinformatics python -m pytest tests/ -k <pattern>  # single case
```

- `tests/conftest.py` inserts `src/` and `scripts/` into `sys.path`, so tests
  import `growth_signals` and `id_mapping` directly
  (`from growth_signals import classify_growth`).
- Test files are the spec for the classifier. `docs/tests.md` documents every
  case and its expected category/tier.

## Architecture

- `scripts/` holds numbered entry points only. Shared library code lives in
  `src/`, sibling to `scripts/`.
- `src/growth_signals.py` — shared engine (`classify_growth()`); imported by
  `03`/`04` (`02` imports only its modifier groups).
- `src/dataset_merge.py` — shared Kim 2010/Hayles 2013 loaders, the Table S3
  parser, `collapse_duplicate_ids()`, and `build_merged_dataset()`. Imported by
  `01b` and the archived `01`. The union's `Presence` /
  `Dispensability_consistency` / `Dispensability_discrepancy` /
  `Dispensability_change_class` columns are defined here.
- `src/id_mapping.py` — shared PomBase ID update: `read_pombase_annotation`,
  `collect_mapping_ids`, `update_sysids`, `build_paralog_columns`,
  `combine_mapping_notes`, plus `PomBaseCol`, `build_update_maps` and
  `normalize_sysid`. Imported by `01b` and the archived `01`.
- `src/phenotype.py` — canonical 01b column names (`Phenotype_description`,
  `Morphological_classification`, `Gene_dispensability`) and
  `format_hayles_phenotypes()`. Imported by `01b`, `02`, `03`, `04` and `07`.
- `src/pipeline_utils.py` — shared `setup_logger`, `load_dit_hap`,
  `category_dr_medians`, and the `DEFAULT_RELEASE` constant.
- `src/category_revisions.py` — shared loader for
  `data/4_categorized_genes/category_revisions.json`, which has three sections:
  `sub_category_order` (original figure), `revisions` (`Sub_category → Category`
  merges), `category_order` (revised figure). Imported by `05` and `06`.
- **Two different "tiers"** — don't conflate them:
  - *signal tier* (1–5): fixed biological hierarchy, internal to
    `classify_growth()`.
  - `Growth_tier`: output column, re-ranked by median DIT-HAP DR per category
    in `04`/`05`. Test expectations use signal tiers.
- `germination` (noun) is intentionally **not** a growth signal — it appears
  only in morphological contexts; the adjective `germinated` is the signal.
- Both figure draw orders and the `Sub_category → Category` merges live in
  `data/4_categorized_genes/category_revisions.json`
  (`sub_category_order` / `revisions` / `category_order`). Edit that JSON to
  reorder the figures or change merges; there is no hardcoded order constant.
  Categories present in the data but missing from an order list are appended,
  not dropped.

## Data gotchas

- **The 2010 raw `.xls` is the full Supplementary Table 1** (4,914 candidate
  genes, 4,836 `Constructed? = yes`). Of the 78 not-constructed candidates,
  18 are the Hayles 2013 additions and 60 are absent from 2013; they have no
  2010 V/E result. **11 constructed Kim strains are absent from 2013** — in the
  archived 01 union they are the 10 rows with `Presence == '2010_only'`, plus
  `SPAC27D7.10c` folded into the `SPAC27D7.09c` row (the IDs are also listed
  in `README.md`); do not treat the 60 not-constructed candidates as
  missing results.
- Both `01b` and the archived `01` update IDs before merging. `SPAC27D7.10c`
  (the 2010-only ts remake) is a PomBase synonym of `SPAC27D7.09c`. In the
  live `01b` Kim pre-merge the two rows collapse: `Kim_ID` =
  `SPAC27D7.09c; SPAC27D7.10c` and `Kim_Note` starts with `COLLAPSED`. In the
  archived 01 union, `Original Systematic ID` = `SPAC27D7.09c; SPAC27D7.10c`,
  `Presence` = `both`, and `note` starts with `COLLAPSED`. The raw
  strain-level reconciliation still has 11 absent strains; the union has 10
  `2010_only` rows. The 4,843-row Hayles table has unique current IDs.
- The live `01b` table has 12 columns: `Systematic ID`, `Gene name`,
  `Gene product`, `Gene type`, `Hayles_ID`, `Kim_ID`, `Hayles_Note`,
  `Kim_Note`, `Phenotype_description`, `Morphological_classification`,
  `Gene_dispensability`, `Kim_Paralog`. `Gene name`/`Gene product`/`Gene type`
  come from PomBase; an empty `Gene name` falls back to `Systematic ID`.
  `Kim_Paralog` comes from Kim's `Paralog` with current IDs.
  The archived `01` union has 25 columns and adds `Presence`,
  `Dispensability_consistency`/`_discrepancy`/`_change_class`, the three
  Hayles Table S3 columns, `Paralog`/`Paralog_gene_name`, the Kim 2010 calls,
  and the original IDs/names. That table is in
  `tmp/01_format_and_update_ids/`; it is not the live input.
- **`data/references/` is gitignored** (`.gitignore` line `references/`) — the
  reference files are required inputs but are *not* in the repo. Missing files
  must be obtained separately:
  - `kim2010/kim2010_supp_table1.xls` (Kim 2010 full Supplementary Table 1;
    `01b` and the archived `01`)
  - `kim2010/kim2010_supp_table3.tsv` (Kim Table S3: 40 literature
    disagreements + judgment criteria; hand-transcribed because the table is
    only published in the supplementary PDF. Tracked via a `.gitignore`
    exception, so it *is* in the repo. Used only by the archived `01`)
  - `hayles2013/hayles2013_supp_table3.xlsx` (Hayles 2013 main table; `01b`
    and the archived `01`)
  - `hayles2013/hayles2013_supp.docx` (Hayles Table S3; used only by the
    archived `01` for the documented-change classes)
  - `pombase-2026-06-01_gene_IDs_names_products.tsv` (produced by `00`)
  The fitted DIT-HAP/gRNA parameters live in `data/raw/` (not gitignored):
  - `gene_level_fitting_statistics_minusDR.tsv` (sign-flipped DR; `04`/`05`/`06`/`07`)
  - `260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv` (`06`/`07`)
  The other exception is the hand-curated
  `previous_manual_check_of_insistent_phenotypes/` subfolder, re-included via
  `.gitignore` and tracked.
- Hand-curated, do **not** regenerate or overwrite:
  - `data/references/previous_manual_check_of_insistent_phenotypes/Inconsistent_phenotypes_at_25_32_manual.xlsx`
  - `data/references/kim2010/kim2010_supp_table3.tsv` — hand-transcribed Kim
    Table S3 (40 published-vs-own disagreements, judgment criteria and the
    four footnotes). The archived `01` validates its `published_data` /
    `our_data` / `final_judgment` columns against Kim Table S1 at runtime.
  - `data/4_categorized_genes/category_revisions.json` — single source of truth for
    both figure draw orders and the `Sub_category → Category` merge used by
    `05`/`06`. Three sections: `sub_category_order`, `revisions`, `category_order`.
    Edit this JSON directly. The old
    `Hayles_2013_OB_inspection_phenotypes_category_revised_20260707.xlsx` is
    frozen; no script reads it anymore.
- `data/arc/` and `scripts/arc/` are the **deprecated previous implementation**,
  not the live pipeline. `scripts/arc/script_flow.sh` references deleted filenames
  — ignore it.

## Docs

Per-step design notes live in `docs/` (`00_`, `01b_`, `04_`, `05_`,
`growth_signals.md`, `tests.md`, `kim2010_hayles2013_reconciliation.md`).
The reconciliation note describes the archived 01 union; the live pipeline
starts at `01b`.
`README.md` has the paper context and concept glossary (`Sub_category` vs
`Category` vs `Growth_tier`).
