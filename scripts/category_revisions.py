"""
Category Configuration — Plot Order and Sub_category to Category Revisions
==========================================================================

Single source of truth for the two hand-curated pieces of category metadata:

- ``order`` — the display order of ``Sub_category`` labels on the DR/um
  distribution figures (step 06). Categories absent from the data are ignored;
  categories present in the data but missing from ``order`` are appended at the
  end.
- ``revisions`` — the ``Sub_category -> Category`` merge mapping used by step 05
  and for grouping in step 06. Only revised ``Sub_category`` values appear as
  keys; any ``Sub_category`` not listed keeps itself as its ``Category``
  (equivalent to ``fillna`` on the original sub-category).

Both live in one tracked JSON file (``data/4_categorized_genes/category_revisions.json``).

Editing the mapping
-------------------
Edit ``data/4_categorized_genes/category_revisions.json`` directly. The legacy
``Hayles_2013_OB_inspection_phenotypes_category_revised_20260707.xlsx`` is
frozen: no script reads it anymore.

Input
-----
- ``data/4_categorized_genes/category_revisions.json``

Output
------
- ``CategoryConfig`` — ``order`` (``list[str]``) and ``revisions``
  (``dict[str, str]``).

Usage
-----
    from category_revisions import load_category_config

    config = load_category_config()
    category = config.revisions.get(sub_category, sub_category)

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-09-14
Version:  2.0.0
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
    """Hand-curated category metadata: display order and merge revisions."""

    order: list[str]
    revisions: dict[str, str]


# =============================================================================
# LOADING
# =============================================================================


def load_category_config(path: Path = DEFAULT_CONFIG) -> CategoryConfig:
    """Load the category order and revision mapping from JSON.

    Args:
        path: Path to the category configuration JSON file.

    Returns:
        ``CategoryConfig`` with the display ``order`` and ``revisions`` mapping.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        ValueError: If the JSON structure or value types are invalid.
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

    order = data.get("order")
    if not isinstance(order, list) or not all(isinstance(c, str) for c in order):
        raise ValueError(f"Category config 'order' must be a list of strings: {path}")
    if len(order) != len(set(order)):
        raise ValueError(f"Category config 'order' contains duplicates: {path}")

    revisions = data.get("revisions")
    if not isinstance(revisions, dict):
        raise ValueError(f"Category config 'revisions' must be an object: {path}")
    for key, value in revisions.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ValueError(
                f"Category revisions must map strings to strings; "
                f"offending entry: {key!r} -> {value!r}"
            )

    return CategoryConfig(order=order, revisions=revisions)
