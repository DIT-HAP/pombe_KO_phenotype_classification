"""Tests for the shared category config loader (order + revisions)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from category_revisions import DEFAULT_CONFIG, load_category_config

ROOT = Path(__file__).resolve().parent.parent

CONFIG = ROOT / DEFAULT_CONFIG

# Every revision key is a Sub_category, every value a merged Category.
EXPECTED_SPOT_CHECKS = {
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
    assert len(config.revisions) == 22
    assert config.order[0] == "spores"
    assert config.order[-1] == "septated"
    # Order must include entries only present in the code-side constant.
    for sentinel in ("some microcolonies, small colonies", "?", "septated"):
        assert sentinel in config.order


def test_mappings_are_flat_strings():
    config = load_category_config(CONFIG)
    for key, value in config.revisions.items():
        assert isinstance(key, str) and key
        assert isinstance(value, str) and value


def test_order_is_unique_strings():
    config = load_category_config(CONFIG)
    assert all(isinstance(c, str) and c for c in config.order)
    assert len(config.order) == len(set(config.order))


def test_spot_checks():
    config = load_category_config(CONFIG)
    for sub_category, category in EXPECTED_SPOT_CHECKS.items():
        assert config.revisions[sub_category] == category


def test_no_chained_mappings():
    config = load_category_config(CONFIG)
    assert set(config.revisions) & set(config.revisions.values()) == set()


def test_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_category_config(ROOT / "data" / "does_not_exist.json")


def test_non_object_json_raises(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(json.dumps(["not", "an", "object"]), encoding="utf-8")
    with pytest.raises(ValueError):
        load_category_config(path)


def test_non_string_revision_values_raise(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(json.dumps({"order": [], "revisions": {"a": 1}}), encoding="utf-8")
    with pytest.raises(ValueError):
        load_category_config(path)


def test_invalid_order_raises(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(json.dumps({"order": [1, 2], "revisions": {}}), encoding="utf-8")
    with pytest.raises(ValueError):
        load_category_config(path)


def test_duplicate_order_raises(tmp_path):
    path = tmp_path / "revisions.json"
    path.write_text(
        json.dumps({"order": ["a", "a"], "revisions": {}}), encoding="utf-8"
    )
    with pytest.raises(ValueError):
        load_category_config(path)
