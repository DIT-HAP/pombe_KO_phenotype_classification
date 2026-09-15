"""Tests for shared pipeline utilities (DR ranking + reference table loading).

DR is sign-flipped: negative and more negative = more depleted. These tests
pin the tier-1 semantics (most depleted first) and the gRNA sign flip.
"""

from __future__ import annotations

import pandas as pd

from pipeline_utils import category_dr_medians, load_dit_hap, load_grna


def test_category_dr_medians_sorted_most_depleted_first():
    df = pd.DataFrame(
        {
            "Systematic ID": ["g1", "g2", "g3", "g4", "g5", "g6"],
            "Category": ["least", "least", "middle", "middle", "most", "most"],
        }
    )
    dit_hap = pd.DataFrame(
        {
            "Systematic ID": ["g1", "g2", "g3", "g4", "g5", "g6"],
            # Sign-flipped: "most" has the most negative median.
            "DR": [0.2, 0.3, -0.4, -0.6, -1.5, -1.9],
        }
    )

    medians = category_dr_medians(df, dit_hap)

    assert list(medians.index) == ["most", "middle", "least"]
    assert medians.iloc[0] < medians.iloc[-1]


def test_category_dr_medians_tier_one_is_most_negative():
    df = pd.DataFrame(
        {"Systematic ID": ["a", "b"], "Category": ["X", "Y"]}
    )
    dit_hap = pd.DataFrame(
        {"Systematic ID": ["a", "b"], "DR": [-2.0, 0.5]}
    )

    medians = category_dr_medians(df, dit_hap)
    tier_map = {cat: i + 1 for i, cat in enumerate(medians.index)}

    assert tier_map == {"X": 1, "Y": 2}


def test_load_dit_hap_drops_null_dr(tmp_path):
    path = tmp_path / "dit_hap.tsv"
    path.write_text(
        "Systematic ID\tDR\n"
        "SPAC1.01\t-0.5\n"
        "SPAC1.02\t\n"
        "SPAC1.03\t-1.2\n",
        encoding="utf-8",
    )

    df = load_dit_hap(path)

    assert list(df.columns) == ["Systematic ID", "DR"]
    assert list(df["Systematic ID"]) == ["SPAC1.01", "SPAC1.03"]


def test_load_grna_negates_um_and_drops_nulls(tmp_path):
    path = tmp_path / "grna.tsv"
    path.write_text(
        "Systematic ID\tum\n"
        "SPAC1.01\t0.8\n"
        "SPAC1.02\t\n"
        "SPAC1.03\t-0.2\n",
        encoding="utf-8",
    )

    df = load_grna(path)

    assert list(df.columns) == ["Systematic ID", "um"]
    assert list(df["Systematic ID"]) == ["SPAC1.01", "SPAC1.03"]
    # Positive source um becomes negative (more negative = more depleted).
    assert df.loc[df["Systematic ID"] == "SPAC1.01", "um"].iloc[0] == -0.8
    assert df.loc[df["Systematic ID"] == "SPAC1.03", "um"].iloc[0] == 0.2
