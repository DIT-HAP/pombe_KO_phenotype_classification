"""
Hayles 2013 Phenotype Formatting
================================

Shared column names and text formatting for the Hayles phenotype table.
`01_format_hayles_phenotypes.py` writes the columns with these names,
and the downstream classification scripts read them.

Formatting rules match the original `01_format_and_update_ids.py` step:

- strip whitespace;
- strip a trailing period from the phenotype description;
- normalise `25, 32` / `25 32` / `32, 25` to `25,32`.

Usage
-----
    from phenotype import format_hayles_phenotypes

    hayles = format_hayles_phenotypes(hayles)

Author:   Yusheng Yang (guidance) + Claude Code (implementation)
Date:     2026-10-29
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
# (none)

# 2. Data Processing Imports
import pandas as pd

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

# Canonical columns in the Hayles 2013 phenotype table.
PHENOTYPE_DESC_COL = "Phenotype_description"
MORPHOLOGY_COL = "Morphological_classification"
DISPENSABILITY_COL = "Gene_dispensability"

# Canonical loader columns in dataset_merge.load_hayles2013().
RAW_PHENOTYPE_DESC_COL = "phenotype_description"
RAW_MORPHOLOGY_COL = "phenotype_classification"
RAW_DISPENSABILITY_COL = "hayles2013_study"

# Temperature text variants that mean both temperatures.
TEMPERATURE_VARIANTS = (
    ("25, 32", "25,32"),
    ("25 32", "25,32"),
    ("32, 25", "25,32"),
)

# =============================================================================
# CORE LOGIC
# =============================================================================


def format_hayles_phenotypes(df: pd.DataFrame) -> pd.DataFrame:
    """Clean Hayles phenotype text in the canonical loader columns."""
    out = df.copy()
    out[RAW_PHENOTYPE_DESC_COL] = (
        out[RAW_PHENOTYPE_DESC_COL]
        .str.strip()
        .str.strip(".")
        .str.strip()
    )
    out[RAW_MORPHOLOGY_COL] = out[RAW_MORPHOLOGY_COL].str.strip()
    out[RAW_DISPENSABILITY_COL] = out[RAW_DISPENSABILITY_COL].str.strip()

    for old, new in TEMPERATURE_VARIANTS:
        out.replace(to_replace=old, value=new, regex=True, inplace=True)

    return out
