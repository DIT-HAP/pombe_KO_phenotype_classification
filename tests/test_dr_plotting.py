"""Tests for the shared DR/um plotting helpers."""

from __future__ import annotations

import pandas as pd

from dr_plotting import (
    GROUP_COL,
    ROW_COL,
    VALUE_COL,
    build_long_frame,
    format_pvalue,
)


def test_format_pvalue_scientific_notation():
    assert format_pvalue(1.2e-9) == "p=1.2E-9"
    assert format_pvalue(8.7e-128) == "p=8.7E-128"


def test_build_long_frame_without_group():
    merged = pd.DataFrame(
        {
            "Systematic ID": ["a", "b"],
            "Row label": ["x", "y"],
        }
    )
    values = pd.DataFrame(
        {
            "Systematic ID": ["a", "b"],
            "DR": [-1.0, -2.0],
        }
    )

    frame = build_long_frame(merged, values, "DR", row_col="Row label")

    assert list(frame.columns) == [VALUE_COL, ROW_COL]
    assert frame[ROW_COL].tolist() == ["x", "y"]
    assert frame[VALUE_COL].tolist() == [-1.0, -2.0]


def test_build_long_frame_with_group():
    merged = pd.DataFrame(
        {
            "Systematic ID": ["a"],
            "Row label": ["x"],
            "Category": ["g"],
        }
    )
    values = pd.DataFrame(
        {
            "Systematic ID": ["a"],
            "DR": [-1.0],
        }
    )

    frame = build_long_frame(
        merged, values, "DR", row_col="Row label", group_col="Category",
    )

    assert list(frame.columns) == [VALUE_COL, ROW_COL, GROUP_COL]
    assert frame[GROUP_COL].tolist() == ["g"]
