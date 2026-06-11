"""
This script will analyze the keywords in the phenotype descriptions for better understanding of the phenotype description.
"""
# import pyforest

import pandas as pd
from pathlib import Path

formatted_phenotypes = pd.read_excel(
    "../data/1_formatted/Hayles_2013_OB_formatted_phenotypes.xlsx"
)

phenotype_description_keywords = (
    formatted_phenotypes["Deletion mutant phenotype description"]
    .str.split(expand=True)
    .stack()
    .rename("Count")
    .droplevel(1, axis=0)
    .str.strip(",")
    .str.strip()
    .value_counts()
    .rename_axis("Keyword")
)

outputFolder = Path("../data/2_phenotype_description_keyword_analysis")
outputFolder.mkdir(exist_ok=True)

phenotype_description_keywords.to_excel(
    outputFolder / "Hayles_2013_OB_phenotype_description_keywords.xlsx"
)
