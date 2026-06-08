"""
This script will category the genes with multi phenotypes to different levels of growth phenotypes.
"""

# import pyforest
# print(pyforest.active_imports())

import re
from pathlib import Path

import pandas as pd


def main():
    # load the genes with multi phenotypes
    genes_with_multi_phenotypes = pd.read_excel(
        "../data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx",
        sheet_name="Multi basic phenotypes",
    )

    genes_with_multi_phenotypes["Category"] = genes_with_multi_phenotypes[
        "Basic phenotype"
    ].apply(category_multi_phenotype)

    # pivot the table
    Phenotypes_pivot = pivot_analysis_of_the_categories(
        genes_with_multi_phenotypes, "Deletion mutant phenotype description", "Category"
    )
    Essentiality_pivot = pivot_analysis_of_the_categories(
        genes_with_multi_phenotypes, "Gene dispensability. This study", "Category"
    )
    Classification_pivot = pivot_analysis_of_the_categories(
        genes_with_multi_phenotypes,
        "Phenotypic classification used for analysis",
        "Category",
    )

    # reorder the columns
    reordered_columns = [
        "spores",
        "spores, germinated spores",
        "spores, germinated spores and divided",
        "spores, germinated spores and colonies",
        "spores, colonies",
        "germinated",
        "germinated and divided",
        "germinated and colonies",
        "microcolonies",
        "microcolonies, small colonies",
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
        outputFolder / "Hayles_2013_OB_category_genes_with_multi_phenotypes.xlsx"
    ) as writer:
        genes_with_multi_phenotypes.to_excel(
            writer, sheet_name="Multi basic phenotypes", index=False
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


# function for category the genes with multi phenotypes to different levels of growth phenotypes
def category_multi_phenotype(x):
    low_frequency_words = ("occasionally", "often", "occasional", "may", "some")

    modfication_words = "|( (slightly|very|possibly|larger|better|highly|slight|large|YES|high|slightlylarger|some|once|occasionally|often|more|twice|many|lots|one|multi|a|no|less|almost|occasional|several|mostly|wee|barely))"
    shape_words = "misshapen|long|skittle|rounded|curved|stubby|branched|swollen|longer|wide|T-shaped|shaped|skittles|wider|septated|septum|aberrantly"
    additional_phenotype_words = "dead|lysis|vacuolated|dark|misplaced|suppressors|diploidising|poor|dumbbell|stationary|sporulating|phase|morphology|diploids|edge|reverting|chains|die|edged|misplace|centre|lyses|colour|supressors|wavy|edges|piled|up|diploidises|inviable"

    x_groups = pd.Series(x.split(","))

    filtered_strs = x_groups[
        ~x_groups.str.strip()
        .str.strip(",")
        .str.strip(".")
        .str.startswith(low_frequency_words)
    ].values

    x = ",".join(filtered_strs)
    x = re.sub(
        f"(({modfication_words}) ({shape_words})|({additional_phenotype_words}))*",
        "",
        x,
    )
    x = re.sub("  ", " ", x)
    x = x.strip()

    if ("spores" in x) and ("germinated" not in x) and ("colon" not in x):
        return "spores"
    elif "spores, some germinated spores" in x:
        return "spores"
    elif (
        ("spores, germinated spores" in x)
        and ("microcolonies" not in x)
        and ("small colonies" not in x)
        and (("divide" not in x) and ("division" not in x))
    ):
        return "spores, germinated spores"
    elif (
        ("spores, germinated spores" in x)
        and ("microcolonies" not in x)
        and ("small colonies" not in x)
        and (("divide" in x) or ("division" in x))
        and (
            ("occasionally" in x)
            or ("often" in x)
            or ("occasional" in x)
            or ("may" in x)
            or ("some" in x)
        )
    ):
        return "spores, germinated spores"
    elif (
        ("spores, germinated spores" in x)
        and ("microcolonies" not in x)
        and ("small colonies" not in x)
        and (("divide" in x) or ("division" in x))
    ):
        return "spores, germinated spores and divided"
    elif (
        ("spores" in x)
        and ("germinated" not in x)
        and ("colon" in x)
        and (
            ("occasionally" in x)
            or ("often" in x)
            or ("occasional" in x)
            or ("may" in x)
            or ("some" in x)
        )
    ):
        return "spores, germinated spores and divided"
    elif ("spores" in x) and ("germinated" not in x) and ("colon" in x):
        return "spores, colonies"
    elif ("spores, germinated spores" in x) and ("colon" in x):
        return "spores, germinated spores and colonies"
    elif (
        ("germinated" in x)
        and ("microcolonies" not in x)
        and ("small colonies" not in x)
        and (("divide" not in x) and ("division" not in x))
    ):
        return "germinated"
    elif (
        ("germinated" in x)
        and ("microcolonies" not in x)
        and ("small colonies" not in x)
        and (("divide" in x) or ("division" in x))
    ):
        return "germinated and divided"
    elif ("germinated" in x) and ("colon" in x):
        return "germinated and colonies"
    elif ("microcolonies" in x) and ("small colonies" not in x):
        return "microcolonies"
    elif (
        ("small colonies" in x)
        and ("microcolonies" not in x)
        and ("very small colonies" in x)
    ):
        return "very small colonies"
    elif (
        ("small colonies" in x)
        and ("microcolonies" not in x)
        and ("very small colonies" not in x)
    ):
        return "small colonies"
    elif ("microcolonies" in x) and ("small colonies" in x):
        return "microcolonies, small colonies"
    else:
        return "WT"


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
