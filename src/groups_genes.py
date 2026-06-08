# %%
"""
This script will group the genes based on the consistency of their phenotype descriptions at different temperatures.
"""

# import pyforest
# print(pyforest.active_imports())

import numpy as np
import pandas as pd
from pathlib import Path


def main():
    formatted_phenotypes = pd.read_excel(
        "../data/1_formatted/Hayles_2013_OB_formatted_phenotypes.xlsx"
    )

    # check the consistency of the phenotype descriptions at different temperatures
    index_consistent_at_both_temperatures = formatted_phenotypes[
        "Deletion mutant phenotype description"
    ].str.contains("25,32")
    index_only_at_32 = formatted_phenotypes[
        "Deletion mutant phenotype description"
    ].str.contains("32") & (
        ~formatted_phenotypes["Deletion mutant phenotype description"].str.contains(
            "25"
        )
    )
    index_inconsistent = (~index_consistent_at_both_temperatures) & (~index_only_at_32)

    # create a new column to store the consistency of the phenotype descriptions at different temperatures
    formatted_phenotypes.loc[
        index_consistent_at_both_temperatures, "Consistency at temperatures"
    ] = "Consistent"
    formatted_phenotypes.loc[
        index_only_at_32, "Consistency at temperatures"
    ] = "Only at 32"
    formatted_phenotypes.loc[
        index_inconsistent, "Consistency at temperatures"
    ] = "Inconsistent"

    # split the phenotypes to basic phenotypes and additional phenotypes for consistent phenotypes at both temperatures
    formatted_phenotypes.loc[
        index_consistent_at_both_temperatures, "Basic phenotype"
    ] = (
        formatted_phenotypes.loc[
            index_consistent_at_both_temperatures,
            "Deletion mutant phenotype description",
        ]
        .str.split(" 25,32", expand=True)[0]
        .str.rstrip(" at")
    )
    formatted_phenotypes.loc[
        index_consistent_at_both_temperatures, "Additional phenotype"
    ] = (
        formatted_phenotypes.loc[
            index_consistent_at_both_temperatures,
            "Deletion mutant phenotype description",
        ]
        .str.split(" 25,32", expand=True)[1]
        .str.lstrip(", ")
        .str.strip()
    )

    # split the phenotypes to basic phenotypes and additional phenotypes for phenotypes only at 32
    formatted_phenotypes.loc[index_only_at_32, "Basic phenotype"] = (
        formatted_phenotypes.loc[
            index_only_at_32, "Deletion mutant phenotype description"
        ]
        .str.split(" 32", expand=True)[0]
        .str.rstrip(" at")
    )
    formatted_phenotypes.loc[index_only_at_32, "Additional phenotype"] = (
        formatted_phenotypes.loc[
            index_only_at_32, "Deletion mutant phenotype description"
        ]
        .str.split(" 32", expand=True)[1]
        .str.lstrip(", ")
        .str.strip()
    )

    # Have one or more basic phenotypes
    def one_or_multi_basic_phenotypes(phenotype):
        if type(phenotype) == float:
            return np.nan
        elif "," in phenotype:
            return "Multi phenotypes"
        else:
            return "One phenotype"

    formatted_phenotypes["One or multi basic phenotypes"] = formatted_phenotypes[
        "Basic phenotype"
    ].apply(one_or_multi_basic_phenotypes)

    # check the number of genes
    print(f"Number of genes: {formatted_phenotypes.shape[0]}")
    print("*" * 20)
    print(
        f"Number of genes with consistent phenotypes at both temperatures: {formatted_phenotypes[formatted_phenotypes['Consistency at temperatures']=='Consistent'].shape[0]}"
    )
    print(
        f"Number of genes with phenotypes only at 32: {formatted_phenotypes[formatted_phenotypes['Consistency at temperatures']=='Only at 32'].shape[0]}"
    )
    print(
        f"Number of genes with inconsistent phenotypes: {formatted_phenotypes[formatted_phenotypes['Consistency at temperatures']=='Inconsistent'].shape[0]}"
    )
    print("*" * 20)
    print(
        f"Number of genes with one basic phenotype: {formatted_phenotypes[formatted_phenotypes['One or multi basic phenotypes']=='One phenotype'].shape[0]}"
    )
    print(
        f"Number of genes with multi basic phenotypes: {formatted_phenotypes[formatted_phenotypes['One or multi basic phenotypes']=='Multi phenotypes'].shape[0]}"
    )
    print(
        f"Number of genes with inconsistent phenotypes: {formatted_phenotypes[formatted_phenotypes['One or multi basic phenotypes'].isna()].shape[0]}"
    )

    # output the results into three sheets in one excel file
    outputFolder = Path("../data/3_grouped_genes")
    outputFolder.mkdir(exist_ok=True)

    with pd.ExcelWriter(outputFolder / "Hayles_2013_OB_grouped_genes.xlsx") as writer:
        formatted_phenotypes[
            formatted_phenotypes["One or multi basic phenotypes"] == "One phenotype"
        ].to_excel(writer, sheet_name="One basic phenotype", index=False)
        formatted_phenotypes[
            formatted_phenotypes["One or multi basic phenotypes"] == "Multi phenotypes"
        ].to_excel(writer, sheet_name="Multi basic phenotypes", index=False)

        # remove the unuseful columns and output the inconsistent phenotypes into a new sheet
        formatted_phenotypes[
            formatted_phenotypes["One or multi basic phenotypes"].isna()
        ].drop(
            columns=[
                "Basic phenotype",
                "Additional phenotype",
                "One or multi basic phenotypes",
            ]
        ).to_excel(
            writer, sheet_name="Inconsistent phenotypes", index=False
        )


if __name__ == "__main__":
    main()
