"""
Category Configuration — Sub_category to Category and Plot Order
================================================================

Single source of truth for how fine-grained ``Sub_category`` labels (step 04)
map to the merged ``Category`` labels, and in what order they are drawn on the
DR/um distribution figures (step 06).

The file is one ordered JSON object. Each key is a name that appears in a
figure; the value is its final category (identical to the key when unchanged):

    {
      "spores": "spores",
      "germinated and divided": "germinated, divided or microcolonies",
      ...
    }

- **Order** — the key order is the display order. The object must contain
  every ``Sub_category`` (for the original figure) plus the merged ``Category``
  names that are not themselves sub-categories (so the revised figure can show
  them). Names absent from the data are ignored; names present in the data but
  missing from the file are dropped from the figure.
- **Merges** — entries where key and value differ are the ``Sub_category ->
  Category`` revisions used by step 05 and for grouping in step 06. Identity
  entries (key == value) are only there to record draw order.

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
- ``CategoryConfig`` — ``order`` (``list[str]``) and ``revisions``
  (``dict[str, str]``, non-identity entries only).

Usage
-----
    from category_revisions import load_category_config

    config = load_category_config()
    category = config.revisions.get(sub_category, sub_category)

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-09-14
Version:  3.0.0
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
    """Display order plus the non-identity merge mapping."""

    order: list[str]
    revisions: dict[str, str]


# =============================================================================
# LOADING
# =============================================================================


def load_category_config(path: Path = DEFAULT_CONFIG) -> CategoryConfig:
    """Load draw order and merge revisions from the category JSON.

    Args:
        path: Path to the category configuration JSON file.

    Returns:
        ``CategoryConfig`` with the display ``order`` (all keys) and the
        ``revisions`` mapping (key != value only).

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        ValueError: If the JSON is not a flat ``str -> str`` object, or if the
            mapping is not idempotent (a merged category is itself merged).
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

    for key, value in data.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ValueError(
                f"Category config must map strings to strings; "
                f"offending entry: {key!r} -> {value!r}"
            )
        if data.get(value, value) != value:
            raise ValueError(
                f"Category config is not idempotent: {key!r} -> {value!r}, "
                f"but {value!r} is itself merged to {data[value]!r}"
            )

    revisions = {k: v for k, v in data.items() if k != v}
    return CategoryConfig(order=list(data), revisions=revisions)
