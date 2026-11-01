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
    category = config.apply(sub_category)

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-09-15
Version:  4.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_CONFIG = Path("data/4_categorized_genes/category_revisions.json")


# =============================================================================
# DATA MODEL
# =============================================================================


@dataclass(kw_only=True, frozen=True)
class CategoryConfig:
    """Draw order for both figures plus the Sub_category -> Category merges."""

    sub_category_order: list[str]
    revisions: dict[str, str]
    category_order: list[str]

    def apply(self, sub_category: str) -> str:
        """Return the merged Category for a Sub_category."""
        return self.revisions.get(sub_category, sub_category)


# =============================================================================
# LOADING
# =============================================================================


def _string_list(data: dict[str, Any], key: str, path: Path) -> list[str]:
    """Validate that ``data[key]`` is a list of unique non-empty strings."""
    value = data.get(key)
    if not isinstance(value, list) or not all(isinstance(x, str) and x for x in value):
        raise ValueError(f"Category config '{key}' must be a list of strings: {path}")
    if len(value) != len(set(value)):
        raise ValueError(f"Category config '{key}' contains duplicates: {path}")
    return value


def load_category_config(path: Path = DEFAULT_CONFIG) -> CategoryConfig:
    """Load and validate the category order and merge configuration JSON."""
    if not path.exists():
        raise FileNotFoundError(
            f"Category config file not found: {path}. Expected {DEFAULT_CONFIG}."
        )

    with path.open(encoding="utf-8") as fh:
        data = json.load(fh)

    if not isinstance(data, dict):
        raise TypeError(f"Category config must be a JSON object, got {type(data).__name__}: {path}")

    sub_category_order = _string_list(data, "sub_category_order", path)
    category_order = _string_list(data, "category_order", path)

    revisions = data.get("revisions")
    if not isinstance(revisions, dict):
        raise TypeError(f"Category config 'revisions' must be an object: {path}")
    for key, value in revisions.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise TypeError(
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
