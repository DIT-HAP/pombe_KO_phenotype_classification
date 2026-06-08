# %%
"""
This script formats the phenotype input file to be used in the classification.
"""

# import pyforest
# print(pyforest.active_imports())

from pathlib import Path
import pandas as pd


def main():
    # Read in the phenotype input file
    phenotype_input = pd.read_excel("../data/0_raw/rsob130053supp2.xlsx")

    # remove the blank spaces
    phenotype_input["Systematic ID"] = phenotype_input["Systematic ID"].str.strip()
    phenotype_input["Deletion mutant phenotype description"] = (
        phenotype_input["Deletion mutant phenotype description"]
        .str.strip()
        .str.strip(".")
        .str.strip()
    )
    phenotype_input["Phenotypic classification used for analysis"] = phenotype_input[
        "Phenotypic classification used for analysis"
    ].str.strip()
    phenotype_input["Gene dispensability. This study"] = phenotype_input[
        "Gene dispensability. This study"
    ].str.strip()

    # replace the representations of the same phenotype in 25,32
    phenotype_input.replace(
        to_replace="25, 32", value="25,32", regex=True, inplace=True
    )
    phenotype_input.replace(to_replace="25 32", value="25,32", regex=True, inplace=True)
    phenotype_input.replace(
        to_replace="32, 25", value="25,32", regex=True, inplace=True
    )

    # drop the final column
    phenotype_input.drop(columns=["Unnamed: 12"], inplace=True)

    # check the classification
    print(
        "The classification categories for check: ",
        phenotype_input["Phenotypic classification used for analysis"].value_counts(),
        sep="\n",
    )

    # check the dispensability
    print(
        "The dispensability categories for check: ",
        phenotype_input["Gene dispensability. This study"].value_counts(),
        sep="\n",
    )

    # save the formatted phenotype input file
    outputFolder = Path("../data/1_formatted")
    outputFolder.mkdir(exist_ok=True)
    phenotype_input.to_excel(
        outputFolder / "Hayles_2013_OB_formatted_phenotypes.xlsx", index=False
    )


if __name__ == "__main__":
    main()
