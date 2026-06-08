"""
This script will merge the categories of all genes.
"""
# %%
# import pyforest

from pathlib import Path
import pandas as pd

# load data
genes_with_one_phenotype = pd.read_excel(
    "../data/4_categorized_genes/Hayles_2013_OB_category_genes_with_one_phenotype.xlsx",
    sheet_name="One basic phenotype",
)
genes_with_multi_phenotypes = pd.read_excel(
    "../data/4_categorized_genes/Hayles_2013_OB_category_genes_with_multi_phenotypes.xlsx",
    sheet_name="Multi basic phenotypes",
)
genes_with_inconsistent_phenotypes = pd.read_excel(
    "../data/4_categorized_genes/Hayles_2013_OB_category_genes_with_inconsistent_phenotypes.xlsx",
    sheet_name="Inconsistent phenotypes",
)

# concat the categories of all genes
all_genes_with_categories = pd.concat(
    [
        genes_with_one_phenotype,
        genes_with_multi_phenotypes,
        genes_with_inconsistent_phenotypes,
    ],
    ignore_index=True,
)

simplified_results = all_genes_with_categories.drop(
    columns=["Basic phenotype", "Additional phenotype", "Category_25", "Category_32"]
)

# save the results
outputFolder = Path("../data/5_merged_categories")
outputFolder.mkdir(exist_ok=True)

with pd.ExcelWriter(outputFolder / "Hayles_2013_OB_merged_categories.xlsx") as writer:
    simplified_results.to_excel(writer, sheet_name="All genes", index=False)

    simplified_results.value_counts("Consistency at temperatures").rename(
        "Count"
    ).to_excel(writer, sheet_name="Consistency at temperatures")

    simplified_results.value_counts("One or multi basic phenotypes").rename(
        "Count"
    ).to_excel(writer, sheet_name="One or multi basic phenotypes")

    simplified_results.value_counts("Category").rename("Count").to_excel(
        writer, sheet_name="Category"
    )
# %%
