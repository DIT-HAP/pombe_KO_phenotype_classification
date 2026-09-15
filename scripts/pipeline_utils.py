"""
Shared Pipeline Utilities
=========================

Small helpers used by more than one pipeline step: logging setup, the
PomBase release constant, and reference-table loading (DIT-HAP / gRNA) with
DR-median ranking.

DIT-HAP ``DR`` is sign-flipped: it is negative, and a more negative value
means more severe depletion. ``category_dr_medians`` therefore sorts
ascending so that tier 1 is the most depleted (most negative) category.

Keeping these here removes the copy-pasted ``setup_logger`` blocks and the
duplicated DR-ranking logic in ``04_categorize_phenotypes.py`` and
``05_merge_categories.py``.

Usage
-----
    from pipeline_utils import setup_logger, load_dit_hap, category_dr_medians

    setup_logger()
    medians = category_dr_medians(df, load_dit_hap(dit_hap_path))

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-09-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import sys
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

# PomBase monthly release, shared by 00 (download) and 01 (ID mapping).
DEFAULT_RELEASE = "2026-06-01"


# =============================================================================
# LOGGING
# =============================================================================


def setup_logger(log_level: str = "INFO") -> None:
    """Configure Loguru to log to stdout with a single consistent format."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level:<8} | {message}",
        level=log_level,
    )


# =============================================================================
# REFERENCE TABLES (DIT-HAP / gRNA)
# =============================================================================


def load_dit_hap(path: Path) -> pd.DataFrame:
    """Load DIT-HAP values, keeping systematic IDs with a non-null DR.

    ``DR`` is sign-flipped (negative = depleted); the source file already
    stores the flipped values, so no negation happens here.
    """
    df = pd.read_csv(path, sep="\t")
    return df[["Systematic ID", "DR"]].dropna(subset=["DR"]).copy()


def category_dr_medians(
    df: pd.DataFrame,
    dit_hap: pd.DataFrame,
    category_col: str = "Category",
) -> pd.Series:
    """Median DIT-HAP DR per category, sorted lowest first (tier 1 = lowest).

    Because DR is sign-flipped (more negative = more depleted), sorting
    ascending puts the most depleted category first, i.e. tier 1.
    """
    joined = df.merge(dit_hap, on="Systematic ID", how="left")
    return joined.groupby(category_col)["DR"].median().sort_values(ascending=True)


def load_grna(path: Path) -> pd.DataFrame:
    """Load gRNA values, keeping systematic IDs with a non-null um.

    The source ``um`` is a positive depletion rate; it is negated on load so
    that, like the DIT-HAP DR, more negative means more depleted.
    """
    df = pd.read_csv(path, sep="\t")
    df = df[["Systematic ID", "um"]].dropna(subset=["um"]).copy()
    df["um"] = -df["um"]
    return df
