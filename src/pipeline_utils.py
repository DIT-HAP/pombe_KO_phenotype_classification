"""
Shared Pipeline Utilities
=========================

Small helpers used by more than one pipeline step: logging setup, the
PomBase release constant, and reference-table loading (DIT-HAP / gRNA).

DIT-HAP ``DR`` is sign-flipped: it is negative, and a more negative value
means more severe depletion.

Keeping these here removes the copy-pasted ``setup_logger`` blocks.

Usage
-----
    from pipeline_utils import setup_logger, load_dit_hap

    setup_logger()
    dit_hap = load_dit_hap(dit_hap_path)

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-09-14
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path
import sys

# 2. Third-party Imports
from loguru import logger

# 3. Data Processing Imports
import pandas as pd

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
# BOUNDARY VALIDATION
# =============================================================================


def require_columns(
    df: pd.DataFrame,
    required: tuple[str, ...],
    source: Path,
) -> None:
    """Raise ValueError when *source* lacks any required column."""
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"{source}: missing required columns {sorted(missing)}")


def require_sheets(
    sheets: dict[str, pd.DataFrame],
    required: tuple[str, ...],
    source: Path,
) -> None:
    """Raise ValueError when *source* lacks any required worksheet."""
    missing = set(required) - set(sheets)
    if missing:
        raise ValueError(f"{source}: missing required sheets {sorted(missing)}")


# =============================================================================
# REFERENCE TABLES (DIT-HAP / gRNA)
# =============================================================================


def _load_metric(path: Path, value_col: str, sign: int = 1) -> pd.DataFrame:
    """Load one metric by systematic ID, dropping missing values."""
    df = pd.read_csv(path, sep="\t")
    columns = ("Systematic ID", value_col)
    require_columns(df, columns, path)
    df = df[list(columns)].dropna(subset=[value_col]).copy()
    if sign != 1:
        df[value_col] = sign * df[value_col]
    return df


def load_dit_hap(path: Path) -> pd.DataFrame:
    """Load DIT-HAP values, keeping systematic IDs with a non-null DR.

    ``DR`` is sign-flipped (negative = depleted); the source file already
    stores the flipped values, so no negation happens here.
    """
    return _load_metric(path, "DR")


def load_grna(path: Path) -> pd.DataFrame:
    """Load gRNA values, keeping systematic IDs with a non-null um.

    The source ``um`` is a positive depletion rate; it is negated on load so
    that, like the DIT-HAP DR, more negative means more depleted.
    """
    return _load_metric(path, "um", sign=-1)
