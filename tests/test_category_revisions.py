"""Tests for the shared category config loader (order + revisions)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from category_revisions import DEFAULT_CONFIG, load_category_config

ROOT = Path(__file__).resolve().parent.parent
CONFIG = ROOT / DEFAULT_CONFIG

EXPECTED_CATEGORY_ORDER = [
    "spores",
    "germinated",
    "divided",
    "microcolonies",
    "very small colonies",
    "small colonies",
    "WT-like",
]

# Non-identity entries: Sub_category -> merged Category.
EXPECTED_REVISIONS = {
    "spores, some germinated": "germinated",
    "spores, germinated and divided": "divided",
    "germinated and divided": "divided",
    "germinated, microcolonies": "microcolonies",
    "microcolonies, small colonies": "small colonies",
    "small colonies, some microcolonies": "small colonies",
}


def test_default_path_points_at_json():
    assert DEFAULT_CONFIG == Path("data/4_categorized_genes/category_revisions.json")


def test_loads_config():
    config = load_category_config(CONFIG)
    assert len(config.sub_category_order) == 30
    assert len(config.revisions) == 25
    assert len(config.category_order) == 7


def test_sub_category_order_bounds():
    config = load_category_config(CONFIG)
    assert config.sub_category_order[0] == "spores"
    assert config.sub_category_order[-1] == "WT-like"


def test_category_order_is_explicit_and_includes_merge_targets():
    config = load_category_config(CONFIG)
    assert config.category_order == EXPECTED_CATEGORY_ORDER
    # Regression: "divided" is a merge target and must be drawn.
    assert "divided" in config.category_order
    assert "divided" in set(config.revisions.values())


def test_identity_entries_are_not_revisions():
    config = load_category_config(CONFIG)
    assert "spores" in config.sub_category_order
    assert "spores" not in config.revisions


def test_spot_checks():
    config = load_category_config(CONFIG)
    for sub_category, category in EXPECTED_REVISIONS.items():
        assert config.revisions[sub_category] == category


def test_no_chained_merges():
    config = load_category_config(CONFIG)
    assert set(config.revisions.values()) & set(config.revisions) == set()


def test_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_category_config(ROOT / "data" / "does_not_exist.json")


def test_non_object_json_raises(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(json.dumps(["not", "an", "object"]), encoding="utf-8")
    with pytest.raises(ValueError):
        load_category_config(path)


def test_missing_section_raises(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(
        json.dumps({"sub_category_order": [], "revisions": {}}), encoding="utf-8"
    )
    with pytest.raises(ValueError):
        load_category_config(path)


def test_non_string_revision_values_raise(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(
        json.dumps(
            {
                "sub_category_order": ["a"],
                "revisions": {"a": 1},
                "category_order": ["b"],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        load_category_config(path)


def test_duplicate_order_entry_raises(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(
        json.dumps(
            {
                "sub_category_order": ["a", "a"],
                "revisions": {},
                "category_order": ["b"],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        load_category_config(path)


def test_non_idempotent_revisions_raise(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(
        json.dumps(
            {
                "sub_category_order": ["a", "b"],
                "revisions": {"a": "b", "b": "c"},
                "category_order": ["c"],
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        load_category_config(path)
