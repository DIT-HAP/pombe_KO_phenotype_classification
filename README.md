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
├── .gitignore              ignores data/, tmp/, references/, etc.
├── data/
│   ├── raw/                read-only. the original supplementary table.
│   │   └── rsob130053supp2.xlsx
│   ├── 1_formatted/        formatted phenotype descriptions
│   ├── 2_phenotype_description_keyword_analysis/
│   ├── 3_grouped_genes/
│   ├── 4_categorized_genes/
│   ├── 5_merged_categories/
│   └── previous_manual_check_of_insistent_phenotypes/
├── results/
│   └── Hayles_2013_OB_merged_categories.xlsx
└── src/
    ├── script_flow.sh                    pipeline entry point
    ├── format_phenotype_input.py
    ├── phenotype_description_keyword_analysis.py
    ├── groups_genes.py
    ├── categorize_genes_with_one_phenotype.py
    ├── categorize_genes_with_multi_phenotypes.py
    ├── categorize_genes_with_inconsistent_phenotypes.py
    └── merge_the_categories_and_assign_essentiality.py
```

## Pipeline

Scripts run in this order via `src/script_flow.sh`:

1. **format_phenotype_input.py** — parses the raw supplementary table into clean, structured format
2. **phenotype_description_keyword_analysis.py** — extracts keywords from the phenotype descriptions
3. **groups_genes.py** — groups genes by phenotype description patterns
4. **categorize_genes_with_one_phenotype.py** — genes with a single, consistent phenotype
   **categorize_genes_with_multi_phenotypes.py** — genes with multiple phenotypes
   **categorize_genes_with_inconsistent_phenotypes.py** — genes with inconsistent observations (e.g. different results at 25°C vs 32°C)
5. **merge_the_categories_and_assign_essentiality.py** — merges the three classification outputs and adds essentiality info

Final output lands in `results/Hayles_2013_OB_merged_categories.xlsx`.

## What's been done so far

- Renamed from `Classify_phenotypes_of_microscopy` to `pombe_KO_phenotype_classification`
- Created this README
- Moved raw data from `data/0_raw/` → `data/raw/` (the `0_` prefix was unnecessary)
- Wrote `.gitignore` covering `data/`, `tmp/`, Python/Jupyter artifacts, macOS system files

## Still to do

- Rename the other data subdirectories (2_phenotype_description_keyword_analysis → keywords, etc.)
- Add numbered prefixes to script filenames so the execution order is obvious
- Decide what to do about the `results/` vs `data/5_merged_categories/` duplication
- Any further restructuring or refactoring of the classification logic
