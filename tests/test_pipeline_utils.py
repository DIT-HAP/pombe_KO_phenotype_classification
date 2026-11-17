"""Tests for shared pipeline utilities (reference table loading + validators).

DR is sign-flipped: negative and more negative = more depleted. These tests
pin the reference-table loading semantics, including the gRNA sign flip.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from pipeline_utils import (
    load_dit_hap,
    load_grna,
    require_columns,
    require_sheets,
)


def test_require_columns_rejects_missing_columns():
    df = pd.DataFrame({"present": [1]})

    with pytest.raises(ValueError, match="missing required columns"):
        require_columns(df, ("present", "absent"), Path("input.tsv"))


def test_require_sheets_rejects_missing_sheets():
    sheets = {"Present": pd.DataFrame()}

    with pytest.raises(ValueError, match="missing required sheets"):
        require_sheets(sheets, ("Present", "Absent"), Path("input.xlsx"))


def test_load_dit_hap_drops_null_dr(tmp_path):
    path = tmp_path / "dit_hap.tsv"
    path.write_text(
        "Systematic ID\tDR\nSPAC1.01\t-0.5\nSPAC1.02\t\nSPAC1.03\t-1.2\n",
        encoding="utf-8",
    )

    df = load_dit_hap(path)

    assert list(df.columns) == ["Systematic ID", "DR"]
    assert list(df["Systematic ID"]) == ["SPAC1.01", "SPAC1.03"]


def test_load_grna_negates_um_and_drops_nulls(tmp_path):
    path = tmp_path / "grna.tsv"
    path.write_text(
        "Systematic ID\tum\nSPAC1.01\t0.8\nSPAC1.02\t\nSPAC1.03\t-0.2\n",
        encoding="utf-8",
    )

    df = load_grna(path)

    assert list(df.columns) == ["Systematic ID", "um"]
    assert list(df["Systematic ID"]) == ["SPAC1.01", "SPAC1.03"]
    # Positive source um becomes negative (more negative = more depleted).
    assert df.loc[df["Systematic ID"] == "SPAC1.01", "um"].iloc[0] == -0.8
    assert df.loc[df["Systematic ID"] == "SPAC1.03", "um"].iloc[0] == 0.2
