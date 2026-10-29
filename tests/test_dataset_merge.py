"""Tests for the Kim 2010 and Hayles 2013 source loaders."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

from dataset_merge import (
    DEFAULT_HAYLES2013,
    DEFAULT_KIM2010,
    load_hayles2013,
    load_kim2010,
    normalize_call,
    normalize_id,
    reconstruct_gene_name_from_date,
)

ROOT = Path(__file__).resolve().parent.parent


def test_normalize_id_strips_all_whitespace():
    assert normalize_id("  SPAC1F3.02c ") == "SPAC1F3.02c"
    assert normalize_id("SPAC 1002.01") == "SPAC1002.01"
    assert normalize_id("") is None
    assert normalize_id(pd.NA) is None


def test_normalize_call_maps_words_and_stray_spaces():
    assert normalize_call("viable") == "V"
    assert normalize_call("essential") == "E"
    assert normalize_call(" E ") == "E"
    assert normalize_call("V ") == "V"
    assert normalize_call("unknown") is None


@pytest.mark.parametrize(
    ("date", "expected"),
    [
        (datetime(2001, 10, 1), "Oct-1"),
        (datetime(2001, 9, 1), "Sep-1"),
        (datetime(2011, 9, 1), "Sep-11"),
        (datetime(2001, 10, 5), "Oct-5"),
    ],
)
def test_reconstruct_gene_name_from_date(date: datetime, expected: str):
    assert reconstruct_gene_name_from_date(date) == expected


def test_real_inputs_load_expected_rows():
    kim = load_kim2010(ROOT / DEFAULT_KIM2010)
    hayles = load_hayles2013(ROOT / DEFAULT_HAYLES2013)

    assert len(kim) == 4914
    assert int(kim["constructed"].sum()) == 4836
    assert len(hayles) == 4843

    # Excel converted these three gene names into dates; the loader repairs them.
    names = hayles.set_index("id")["gene_name_2013"]
    assert names["SPAC1F3.10c"] == "Oct-1"
    assert names["SPAC5D6.05"] == "Sep-11"
    assert names["SPBC4C3.12"] == "Sep-1"
