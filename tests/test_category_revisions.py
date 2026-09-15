"""Tests for the shared category config loader (order + revisions)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from category_revisions import DEFAULT_CONFIG, load_category_config

ROOT = Path(__file__).resolve().parent.parent

CONFIG = ROOT / DEFAULT_CONFIG

# Non-identity entries: Sub_category -> merged Category.
EXPECTED_REVISIONS = {
    "germinated and divided": "germinated, divided or microcolonies",
    "germinated, microcolonies": "germinated, divided or microcolonies",
    "spores, some germinated": "spores, germinated",
    "spores, some germinated, some divided": "spores, germinated",
    "small colonies, some microcolonies": "small colonies",
    "spores, microcolonies": "spores, miscellaneous",
}


def test_default_path_points_at_json():
    assert DEFAULT_CONFIG == Path("data/4_categorized_genes/category_revisions.json")


def test_loads_config():
    config = load_category_config(CONFIG)
    assert len(config.order) == 33
    assert len(config.revisions) == 22
    assert config.order[0] == "spores"
    assert config.order[-1] == "WT-like"


def test_order_covers_both_plot_rows():
    config = load_category_config(CONFIG)
    # Merged-only categories must be in the order for the revised figure.
    for name in ("spores, miscellaneous", "germinated, divided or microcolonies"):
        assert name in config.order
    # Previously-missing sub-categories must be present now.
    for name in ("germinated and divided", "small colonies, some microcolonies"):
        assert name in config.order


def test_identity_entries_are_not_revisions():
    config = load_category_config(CONFIG)
    assert "spores" in config.order
    assert "spores" not in config.revisions


def test_order_is_unique_strings():
    config = load_category_config(CONFIG)
    assert all(isinstance(c, str) and c for c in config.order)
    assert len(config.order) == len(set(config.order))


def test_spot_checks():
    config = load_category_config(CONFIG)
    for sub_category, category in EXPECTED_REVISIONS.items():
        assert config.revisions[sub_category] == category


def test_no_chained_merges():
    config = load_category_config(CONFIG)
    # A merge target must not itself be merged further.
    assert set(config.revisions.values()) & set(config.revisions) == set()


def test_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_category_config(ROOT / "data" / "does_not_exist.json")


def test_non_object_json_raises(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(json.dumps(["not", "an", "object"]), encoding="utf-8")
    with pytest.raises(ValueError):
        load_category_config(path)


def test_non_string_values_raise(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(json.dumps({"a": 1}), encoding="utf-8")
    with pytest.raises(ValueError):
        load_category_config(path)


def test_non_idempotent_mapping_raises(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(json.dumps({"a": "b", "b": "c"}), encoding="utf-8")
    with pytest.raises(ValueError):
        load_category_config(path)
