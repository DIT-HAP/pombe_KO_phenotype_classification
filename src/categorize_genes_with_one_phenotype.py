"""
This script will category the genes with one phenotype to different levels of growth phenotypes.
"""

# import pyforest
# print(pyforest.active_imports())

import pandas as pd
from pathlib import Path


def main():
    # load the genes with one phenotype
    genes_with_one_phenotype = pd.read_excel(
        "../data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx",
        sheet_name="One basic phenotype",
    )

    genes_with_one_phenotype["Category"] = genes_with_one_phenotype[
        "Basic phenotype"
    ].apply(category_one_phenotype)

    # pivot the table
    Phenotypes_pivot = pivot_analysis_of_the_categories(
        genes_with_one_phenotype, "Deletion mutant phenotype description", "Category"
    )
    Essentiality_pivot = pivot_analysis_of_the_categories(
        genes_with_one_phenotype, "Gene dispensability. This study", "Category"
    )
    Classification_pivot = pivot_analysis_of_the_categories(
        genes_with_one_phenotype,
        "Phenotypic classification used for analysis",
        "Category",
    )

    # reorder the columns
    reordered_columns = [
        "spores",
        "germinated",
        "germinated and divided",
        "microcolonies",
        "very small colonies",
        "small colonies",
        "WT",
    ]
    Phenotypes_pivot = Phenotypes_pivot[reordered_columns].sort_values(
        reordered_columns, ascending=False
    )
    Essentiality_pivot = Essentiality_pivot[reordered_columns].sort_values(
        reordered_columns, ascending=False
    )
    Classification_pivot = Classification_pivot[reordered_columns].sort_values(
        reordered_columns, ascending=False
    )

    # save the result
    outputFolder = Path("../data/4_categorized_genes")
    outputFolder.mkdir(exist_ok=True)

    with pd.ExcelWriter(
        outputFolder / "Hayles_2013_OB_category_genes_with_one_phenotype.xlsx"
    ) as writer:
        genes_with_one_phenotype.to_excel(
            writer, sheet_name="One basic phenotype", index=False
        )
        Phenotypes_pivot.style.background_gradient(cmap="Blues", axis=1).to_excel(
            writer, sheet_name="Phenotypes pivot"
        )
        Essentiality_pivot.style.background_gradient(cmap="Blues", axis=1).to_excel(
            writer, sheet_name="Essentiality pivot"
        )
        Classification_pivot.style.background_gradient(cmap="Blues", axis=1).to_excel(
            writer, sheet_name="Classification pivot"
        )


# function for category the genes with one phenotype to different levels of growth phenotypes
def category_one_phenotype(phenotype):
    if ("spores" in phenotype) and ("germinated" not in phenotype):
        return "spores"
    elif (
        ("germinated" in phenotype)
        and ("spores" in phenotype)
        and (("divide" not in phenotype) and ("division" not in phenotype))
    ):
        return "germinated"
    elif (
        ("germinated" in phenotype)
        and ("spores" in phenotype)
        and (("divide" in phenotype) or ("division" in phenotype))
    ):
        return "germinated and divided"
    elif "microcolonies" in phenotype:
        return "microcolonies"
    elif "very small colonies" in phenotype:
        return "very small colonies"
    elif "small colonies" in phenotype:
        return "small colonies"
    else:
        return "WT"


# pivot function
def pivot_analysis_of_the_categories(df, row_name, col_name):
    pivot_result = (
        df[[row_name, col_name]]
        .value_counts()
        .rename("Count")
        .reset_index()
        .pivot(index=row_name, columns=col_name, values="Count")
        .fillna(0)
        .astype(int)
    )

    return pivot_result


if __name__ == "__main__":
    main()
