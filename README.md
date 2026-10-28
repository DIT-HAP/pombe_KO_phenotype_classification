# pombe_KO_phenotype_classification

Classifying fission yeast knockout phenotypes from the Hayles et al. 2013 genome-wide microscopy screen.

## What this is

The Hayles paper screened 4,843 *S. pombe* deletion mutants under the microscope and wrote down what each one looked like — stuff like "cells are long and branched" or "small round cells at 25°C". That's a lot of natural language descriptions to wade through. This project takes those descriptions and sorts them into clean, structured phenotype categories.

## Key papers

**Kim DU, Hayles J, Kim D, Wood V, Park HO, Won M, Yoo HS, et al.**  
*Analysis of a genome-wide set of gene deletions in the fission yeast Schizosaccharomyces pombe.*  
Nature Biotechnology, 28, 617–623 (2010).  
DOI: [10.1038/nbt.1628](https://doi.org/10.1038/nbt.1628)

Built the knockout library. 4,836 heterozygous deletion mutants, covering 98.4% of the *S. pombe* genome. Essential genes are more likely to be single-copy, broadly conserved, and contain introns. Fission yeast has more essential genes than budding yeast.

**Hayles J, Wood V, Jeffery D, Hoe KL, Kim DU, Park HO, Salas-Pino S, Heichinger C, Nurse P.**  
*A genome-wide resource of cell cycle and cell shape genes of fission yeast.*  
Open Biology, 3(5), 130053 (2013).  
DOI: [10.1098/rsob.130053](https://doi.org/10.1098/rsob.130053)

Microscope screen of 4,843 haploid deletion mutants. Defined 14 phenotype categories: WT, spores, germination, misshapen essential, misshapen viable, misshapen weak viable, long high penetrance, long low penetrance, long branched, rounded, stubby, curved, small, and skittle. Found 513 cell cycle genes (276 newly described) and 333 cell shape genes.

The source tables are in `data/references/kim2010/` and
`data/references/hayles2013/`, downloaded from the Nature Biotechnology and
Open Biology websites.

## Current structure

```
pombe_KO_phenotype_classification/
├── README.md
├── .gitignore              ignores tmp/, references/, etc. (data/ is tracked)
├── data/
│   ├── raw/                read-only. fitted DIT-HAP/gRNA parameters.
│   │   ├── gene_level_fitting_statistics_minusDR.tsv
│   │   └── 260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv
│   ├── references/         external data (PomBase + source papers/tables)
│   │   ├── kim2010/         Kim 2010 PDFs + supplementary tables
│   │   ├── hayles2013/      Hayles 2013 PDFs + supplementary tables
│   │   └── previous_manual_check_.../ (hand-curated)
│   ├── 1_formatted/        ID-updated per-source tables, 2013 subset + 2010/2013 union
│   ├── 2_keyword_profile/  keyword-level analysis per Growth_tier
│   ├── 3_grouped_genes/    grouped by phenotype consistency
│   ├── 4_categorized_genes/ categorized phenotypes, inspection pivots
│   │                        + category_revisions.json (orders + merges)
│   └── 5_merged_categories/ final merged table + summaries
├── results/
│   ├── Hayles_2013_OB_merged_categories.xlsx
│   ├── DR_um_distribution_original.png
│   ├── DR_um_distribution_original.pdf
│   ├── DR_um_distribution_revised.png
│   ├── DR_um_distribution_revised.pdf
│   ├── DR_um_distribution_pure.png
│   └── DR_um_distribution_pure.pdf
└── scripts/
    ├── 00_download_pombase_annotation.py
    ├── 01_format_and_update_ids.py
    ├── 02_keyword_profile.py
    ├── 03_group_genes.py
    ├── 04_categorize_phenotypes.py
    ├── 05_merge_categories.py
    ├── 06_plot_dr_distribution.py       DR/um distribution (violin + stats)
    ├── 07_plot_dr_distribution_pure.py  five pure-category figure
    ├── growth_signals.py             shared signal-detection engine
    ├── category_revisions.py         shared category order + merge mapping
    ├── dataset_merge.py              shared 2010/2013 loaders + union merge
    └── arc/                          deprecated previous implementation
```

## Kim 2010 vs Hayles 2013 reconciliation

`01_format_and_update_ids.py` updates IDs/names in both tables **before**
building a 4,913-row gene-level union (via `scripts/dataset_merge.py`) with
`Presence` (`both` / `2010_only` / `2013_only` / `neither`),
`Dispensability_consistency`, `Dispensability_discrepancy` and
`Dispensability_change_class`, plus `Gene product` (PomBase) and
`Paralog`/`Paralog_gene_name` (Kim data, current IDs/names). Key facts:

- The 2010 `.xls` is the **full Supplementary Table 1**: 4,914 candidate
  genes, of which 4,836 were constructed.
- The 2013 table has 4,843 genes: 4,825 Kim-constructed strains plus 18
  Hayles additions.
- **11 constructed 2010 strains are absent from 2013** (they have a 2010 V/E
  call but no 2013 row): `SPAC19B12.13`, `SPAC1F3.02c`, `SPAC23A1.10`,
  `SPAC25B8.12c`, `SPAC27D7.10c`, `SPAC977.07c`, `SPCC16A11.02`,
  `SPCC16C4.16c`, `SPCC736.05`, `SPCC737.04`, `SPCC970.11c`.
- At the gene level, `SPAC27D7.10c` is a synonym of `SPAC27D7.09c` in current
  PomBase, so the union collapses it into the `.09c` row (both original IDs
  are kept). `Presence == 2010_only` therefore has 10 rows; the 11-strain
  statement holds at the strain level.
- Another 60 Kim rows are `Constructed? = no` candidates that were never
  constructed or screened; they never had a 2010 dispensability result.
- The 18 Hayles additions (Table S2) are Kim candidates that Kim did not
  construct: 6 newly constructed + 12 genes with <40% deletion coverage.
- Table S3 documents 21 dispensability changes; four further V/E changes are
  present in the 2013 main table but not in Table S3.

See `docs/kim2010_hayles2013_reconciliation.md` for the per-gene record and
data-quality notes, and `docs/01_update_systematic_ids.md` for the
merge/output specification.

## Pipeline

Scripts run in this order:

1. **00_download_pombase_annotation.py** — downloads the PomBase gene annotation
2. **01_format_and_update_ids.py** — updates IDs and gene names in both source
   tables before merging, builds the gene-level union with
   presence/dispensability comparison columns, adds `Gene product` (PomBase),
   `Paralog`/`Paralog_gene_name` (Kim data), and the Hayles Table S3
   strain-status / 2010-change annotations, then formats the phenotypes.
   Writes the ID-updated per-source tables
   (`Kim2010_premerge_formatted.xlsx`, `Hayles2013_premerge_formatted.xlsx`),
   the union, and the 4,843-row Hayles subset that feeds step 03. The Kim
   intermediate records the published-vs-own comparison and the 40 Kim
   Table S3 judgment criteria; the Hayles intermediate records the
   2013-vs-Kim verdict.
3. **03_group_genes.py** — groups genes by phenotype consistency at different temperatures
4. **04_categorize_phenotypes.py** — categorises growth phenotypes using signal-based detection;
   re-ranks Growth_tier by sign-flipped DIT-HAP DR median (most depleted = tier 1)
5. **05_merge_categories.py** — merges three branches, applies the hand-curated
   `Sub_category → Category` merges (`data/4_categorized_genes/category_revisions.json`),
   re-ranks Growth_tier by DR median
6. **02_keyword_profile.py** — keyword-level analysis validating word categories
   (growth signals, morphology, modifiers) per Growth_tier. Re-analysis of the
   final merged output, so it runs **after 05**, despite its number.
7. **06_plot_dr_distribution.py** — cnsplots violin/box + statistics-table plots of
   sign-flipped depletion (DR for DIT-HAP, um for gRNA; more negative = more depleted);
   generates two figures (each as PNG + PDF): original (fine-grained Sub_category) and
   revised (merged Category, with Mann-Whitney U p-value annotations)
8. **07_plot_dr_distribution_pure.py** — the same cnsplots layout restricted to five
   pure categories (spores, germinated, microcolonies, small colonies, WT-like)

Final output lands in `results/Hayles_2013_OB_merged_categories.xlsx` and
`results/DR_um_distribution_*.{png,pdf}`.

## Key concepts

- **Sub_category** — fine-grained classification from `classify_growth()` (e.g. `germinated, often divided`)
- **Category** — merged classification after manual revision (e.g. `germinated, divided or microcolonies`)
- **Growth_tier** — ranked by sign-flipped DIT-HAP DR median per category (tier 1 = most depleted = most negative DR, most severe growth defect)
- **Signal tier** — biological hierarchy (1=spores, 2=germinated, 3=microcolonies, 4=small colonies, 5=WT-like); used internally by `classify_growth()` but not in final output

## Test suite

125 tests covering signal detection, modifier handling, segment filtering,
Single/Multiple classification, the Kim-2010/Hayles-2013 union, and the
per-source ID-update intermediates:

```bash
mamba run -n bioinformatics python -m pytest tests/ -v
```

See `docs/tests.md` for detailed test documentation.
