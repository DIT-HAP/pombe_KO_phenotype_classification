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

The raw supplementary table is `rsob130053supp2.xlsx`, downloaded from the Open Biology website.

## Current structure

```
pombe_KO_phenotype_classification/
├── README.md
├── .gitignore              ignores tmp/, references/, etc. (data/ is tracked)
├── data/
│   ├── raw/                read-only. the original supplementary table.
│   │   └── rsob130053supp2.xlsx
│   ├── references/         external data (DIT-HAP, gRNA, PomBase)
│   ├── 1_formatted/        formatted phenotype descriptions
│   ├── 2_keyword_profile/  keyword-level analysis per Growth_tier
│   ├── 3_grouped_genes/    grouped by phenotype consistency
│   ├── 4_categorized_genes/ categorized phenotypes, inspection pivots
│   │                        + category_revisions.json (order + merges)
│   ├── 5_merged_categories/ final merged table + summaries
│   └── previous_manual_check_of_insistent_phenotypes/
├── results/
│   ├── Hayles_2013_OB_merged_categories.xlsx
│   ├── DR_um_distribution_original.png
│   └── DR_um_distribution_revised.png
└── scripts/
    ├── 00_download_pombase_annotation.py
    ├── 01_format_and_update_ids.py
    ├── 02_keyword_profile.py
    ├── 03_group_genes.py
    ├── 04_categorize_phenotypes.py
    ├── 05_merge_categories.py
    ├── 06_plot_dr_distribution.py     DR/um distribution plots
    ├── growth_signals.py             shared signal-detection engine
    ├── category_revisions.py         shared category order + merge mapping
    └── arc/                          deprecated previous implementation
```

## Pipeline

Scripts run in this order:

1. **00_download_pombase_annotation.py** — downloads the PomBase gene annotation
2. **01_format_and_update_ids.py** — formats raw data and updates systematic IDs
3. **02_keyword_profile.py** — keyword-level analysis validating word categories
   (growth signals, morphology, modifiers) per Growth_tier
4. **03_group_genes.py** — groups genes by phenotype consistency at different temperatures
5. **04_categorize_phenotypes.py** — categorises growth phenotypes using signal-based detection;
   re-ranks Growth_tier by DIT-HAP DR median
6. **05_merge_categories.py** — merges three branches, applies the hand-curated
   `Sub_category → Category` merges (`data/4_categorized_genes/category_revisions.json`),
   re-ranks Growth_tier by DR median
7. **06_plot_dr_distribution.py** — plots DR (DIT-HAP) and um (gRNA) distributions per category;
   generates two figures: original (fine-grained Sub_category) and revised (merged Category,
   with Mann-Whitney U p-value annotations)

Final output lands in `results/Hayles_2013_OB_merged_categories.xlsx` and
`results/DR_um_distribution_*.png`.

## Key concepts

- **Sub_category** — fine-grained classification from `classify_growth()` (e.g. `germinated, often divided`)
- **Category** — merged classification after manual revision (e.g. `germinated, divided or microcolonies`)
- **Growth_tier** — ranked by DIT-HAP DR median per category (tier 1 = highest median, most severe growth defect)
- **Signal tier** — biological hierarchy (1=spores, 2=germinated, 3=microcolonies, 4=small colonies, 5=WT-like); used internally by `classify_growth()` but not in final output

## Test suite

75 tests covering signal detection, modifier handling, segment filtering, and
Single/Multiple classification:

```bash
mamba run -n bioinformatics python -m pytest tests/ -v
```

See `docs/tests.md` for detailed test documentation.
