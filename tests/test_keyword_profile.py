"""Tests for keyword-profile category grouping."""

from __future__ import annotations

from collections import Counter
import importlib

import pandas as pd

keyword_profile = importlib.import_module("02_keyword_profile")


def test_count_words_per_category_groups_by_merged_category():
    df = pd.DataFrame(
        {
            "Category": ["spores", "spores", "divided"],
            "Phenotype_description": ["spores", "spores", "elongated"],
        }
    )

    counts = keyword_profile.count_words_per_category(df)

    assert counts["spores"]["spores"] == 2
    assert counts["divided"]["elongated"] == 1


def test_build_term_table_uses_category_columns():
    table = keyword_profile.build_term_table(
        {"spores": Counter({"spores": 2}), "divided": Counter({"long": 1})},
    )

    assert list(table.columns) == ["Word category", "divided", "spores", "Total"]
