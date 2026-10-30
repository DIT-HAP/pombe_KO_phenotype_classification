"""Tests for the Kim 2010 and Hayles 2013 source loaders."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from dataset_merge import (
    DEFAULT_HAYLES2013,
    DEFAULT_KIM2010,
    load_hayles2013,
    load_kim2010,
    normalize_call,
    normalize_id,
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


def test_real_inputs_load_expected_rows():
    kim = load_kim2010(ROOT / DEFAULT_KIM2010)
    hayles = load_hayles2013(ROOT / DEFAULT_HAYLES2013)

    assert list(kim.columns) == ["id", "paralogues_2010"]
    assert list(hayles.columns) == [
        "id",
        "hayles2013_study",
        "phenotype_description",
        "phenotype_classification",
    ]
    assert len(kim) == 4914
    assert len(hayles) == 4843
    assert kim["id"].is_unique
    assert hayles["id"].is_unique
