"""Tests for keyword-profile tier labels."""

from __future__ import annotations

from collections import Counter
import importlib

import pandas as pd

keyword_profile = importlib.import_module("02_keyword_profile")


def test_tier_labels_from_frame_uses_categories():
    df = pd.DataFrame(
        {
            "Growth_tier": [1, 1, 3, 6],
            "Category": ["spores", "spores", "divided", "small colonies"],
        }
    )

    assert keyword_profile.tier_labels_from_frame(df) == {
        1: "spores",
        3: "divided",
        6: "small colonies",
    }


def test_tier_label_falls_back_for_unknown_tier():
    assert keyword_profile.tier_label(9, {1: "spores"}) == "Tier 9"


def test_build_term_table_uses_data_derived_labels():
    table = keyword_profile.build_term_table(
        {1: Counter({"spores": 2}), 3: Counter({"long": 1})},
        tier_labels={1: "spores", 3: "divided"},
    )

    assert list(table.columns) == [
        "Category",
        "Tier_1_spores",
        "Tier_3_divided",
        "Total",
    ]
