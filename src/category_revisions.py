"""
Category Configuration — Order and Merges
==========================================

Single source of truth for how the fine-grained ``Sub_category`` labels (step
04) map to merged ``Category`` labels, and for the draw order of both figures
(step 06).

The file has three explicit sections::

    {
      "sub_category_order": ["spores", "spores, some germinated", ...],
      "revisions": {"spores, some germinated": "germinated", ...},
      "category_order": ["spores", "germinated", "divided", ...]
    }

- ``sub_category_order`` — draw order for the original figure (one row per
  ``Sub_category``). Categories present in the data but not listed are appended.
- ``revisions`` — the ``Sub_category -> Category`` merges used by step 05 and
  for grouping in step 06. Only changed sub-categories are listed.
- ``category_order`` — draw order for the revised figure (one row per merged
  ``Category``). Categories present in the data but not listed are appended, so
  a merge target can never be silently dropped.

Editing
-------
Edit ``data/4_categorized_genes/category_revisions.json`` directly. The legacy
``Hayles_2013_OB_inspection_phenotypes_category_revised_20260707.xlsx`` is
frozen: no script reads it anymore.

Input
-----
- ``data/4_categorized_genes/category_revisions.json``

Output
------
- ``CategoryConfig`` — ``sub_category_order``, ``revisions``, ``category_order``.

Usage
-----
    from category_revisions import load_category_config

    config = load_category_config()
    category = config.revisions.get(sub_category, sub_category)

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-09-15
Version:  4.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import json
from dataclasses import dataclass
from pathlib import Path

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_CONFIG = Path("data/4_categorized_genes/category_revisions.json")


# =============================================================================
# DATA MODEL
# =============================================================================


@dataclass(frozen=True)
class CategoryConfig:
    """Draw order for both figures plus the Sub_category -> Category merges."""

    sub_category_order: list[str]
    revisions: dict[str, str]
    category_order: list[str]


# =============================================================================
# LOADING
# =============================================================================


def _string_list(data: dict, key: str, path: Path) -> list[str]:
    """Validate that ``data[key]`` is a list of unique non-empty strings."""
    value = data.get(key)
    if not isinstance(value, list) or not all(isinstance(x, str) and x for x in value):
        raise ValueError(f"Category config '{key}' must be a list of strings: {path}")
    if len(value) != len(set(value)):
        raise ValueError(f"Category config '{key}' contains duplicates: {path}")
    return value


def load_category_config(path: Path = DEFAULT_CONFIG) -> CategoryConfig:
    """Load the category configuration JSON.

    Args:
        path: Path to the category configuration JSON file.

    Returns:
        ``CategoryConfig`` with ``sub_category_order``, ``revisions`` and
        ``category_order``.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        ValueError: If the structure or value types are invalid, or if the
            revision mapping is not idempotent (a target is itself revised).
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Category config file not found: {path}. "
            f"Expected {DEFAULT_CONFIG}."
        )

    with path.open(encoding="utf-8") as fh:
        data = json.load(fh)

    if not isinstance(data, dict):
        raise ValueError(
            f"Category config must be a JSON object, got {type(data).__name__}: {path}"
        )

    sub_category_order = _string_list(data, "sub_category_order", path)
    category_order = _string_list(data, "category_order", path)

    revisions = data.get("revisions")
    if not isinstance(revisions, dict):
        raise ValueError(f"Category config 'revisions' must be an object: {path}")
    for key, value in revisions.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ValueError(
                f"Category revisions must map strings to strings; "
                f"offending entry: {key!r} -> {value!r}"
            )
    if set(revisions.values()) & set(revisions):
        raise ValueError(
            f"Category revisions are not idempotent (a merge target is itself "
            f"revised): {sorted(set(revisions.values()) & set(revisions))}"
        )

    return CategoryConfig(
        sub_category_order=sub_category_order,
        revisions=revisions,
        category_order=category_order,
    )
