"""Tests for the Hayles phenotype-table formatter."""

from __future__ import annotations

import importlib

import pandas as pd

format_hayles = importlib.import_module("01_format_hayles_phenotypes")


def test_add_annotation_direct_and_synonym_fallback():
    annotation = pd.DataFrame(
        {
            "gene_systematic_id": ["SPAC1.01", "SPAC2.02", "SPAC3.03"],
            "gene_name": ["alpha", "beta", ""],
            "gene_product": ["p1", "p2", "p3"],
            "gene_type": ["protein-coding"] * 3,
            "synonyms": [None, "old2, alias2", None],
        }
    )
    merged = pd.DataFrame({"Systematic ID": ["SPAC1.01", "alias2", "SPAC3.03"]})

    result = format_hayles.add_annotation(merged, annotation)

    assert result["Gene name"].tolist() == ["alpha", "beta", "SPAC3.03"]
    assert result["Gene product"].tolist() == ["p1", "p2", "p3"]
    assert result["Gene type"].tolist() == ["protein-coding"] * 3


def test_add_annotation_prefers_direct_id_over_synonym():
    annotation = pd.DataFrame(
        {
            "gene_systematic_id": ["SPAC1.01", "SPAC2.02"],
            "gene_name": ["alpha", "beta"],
            "gene_product": ["p1", "p2"],
            "gene_type": ["protein-coding", "protein-coding"],
            "synonyms": [None, "SPAC1.01"],
        }
    )
    merged = pd.DataFrame({"Systematic ID": ["SPAC1.01"]})

    result = format_hayles.add_annotation(merged, annotation)

    assert result.loc[0, "Gene name"] == "alpha"
    assert result.loc[0, "Gene product"] == "p1"


def test_add_annotation_falls_back_to_id_for_blank_name_and_missing_id():
    annotation = pd.DataFrame(
        {
            "gene_systematic_id": ["SPAC1.01"],
            "gene_name": ["   "],
            "gene_product": ["p1"],
            "gene_type": ["protein-coding"],
            "synonyms": [None],
        }
    )
    merged = pd.DataFrame({"Systematic ID": ["SPAC1.01", "unknown"]})

    result = format_hayles.add_annotation(merged, annotation)

    assert result["Gene name"].tolist() == ["SPAC1.01", "unknown"]
    assert result.loc[0, "Gene product"] == "p1"
    assert pd.isna(result.loc[1, "Gene product"])
    assert pd.isna(result.loc[1, "Gene type"])
