#!/bin/bash

# This script is used to run the whole pipeline of the project

# 1.format the data
python format_phenotype_input.py

# 2.keyword analysis of the phenotype description
python phenotype_description_keyword_analysis.py

# 3.Split the genes to different groups
python groups_genes.py

# 4.categorize the genes
python categorize_genes_with_one_phenotype.py
python categorize_genes_with_multi_phenotypes.py
python categorize_genes_with_inconsistent_phenotypes.py

# 5.merge the results
python merge_the_categories_and_assign_essentiality.py

# move final results to the results folder for source control
cp ../data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx ../results