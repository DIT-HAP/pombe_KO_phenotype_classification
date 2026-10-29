"""Tests for shared Hayles phenotype formatting and column names."""

from __future__ import annotations

import pandas as pd

from phenotype import (
    DISPENSABILITY_COL,
    MORPHOLOGY_COL,
    PHENOTYPE_DESC_COL,
    format_hayles_phenotypes,
)


def test_column_names():
    assert PHENOTYPE_DESC_COL == "Phenotype_description"
    assert MORPHOLOGY_COL == "Morphological_classification"
    assert DISPENSABILITY_COL == "Gene_dispensability"


def test_format_hayles_phenotypes_normalises_text():
    df = pd.DataFrame(
        {
            "phenotype_description": [
                "small colonies at 25, 32.",
                "elongated at 25 32",
                "round at 32, 25.",
                "WT cells.",
            ],
            "phenotype_classification": [" small ", "long high penetrance", "rounded", "WT"],
            "hayles2013_study": [" V ", "E", "V", " V "],
        }
    )

    out = format_hayles_phenotypes(df)

    assert out["phenotype_description"].tolist() == [
        "small colonies at 25,32",
        "elongated at 25,32",
        "round at 25,32",
        "WT cells",
    ]
    assert out["phenotype_classification"].tolist() == [
        "small",
        "long high penetrance",
        "rounded",
        "WT",
    ]
    assert out["hayles2013_study"].tolist() == ["V", "E", "V", "V"]
