"""
Shared Kim 2010 / Hayles 2013 loaders
=====================================

Library module used by ``01_format_hayles_phenotypes.py`` for the two source
loaders and duplicate-ID collapse. It canonicalises exactly the columns the
live pipeline consumes: Kim's systematic ID and Paralog, and Hayles' ID,
dispensability, phenotype description and classification. It does not build
the retired 25-column union.

The archived ``01_format_and_update_ids.py`` vendors its own frozen copy of
this module under ``tmp/01_format_and_update_ids/src/`` (gitignored and
untracked); it does not import this file.

Input
-----
- Kim 2010 full Supplementary Table 1 ``.xls``.
- Hayles 2013 main ``.xlsx``.

Output
------
- Canonical DataFrames for the two sources.

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-10-30
Version: 2.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from pathlib import Path
import re

# 2. Third-party Imports
from loguru import logger

# 3. Data Processing Imports
import pandas as pd

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_KIM2010 = Path("data/references/kim2010/kim2010_supp_table1.xls")
DEFAULT_HAYLES2013 = Path("data/references/hayles2013/hayles2013_supp_table3.xlsx")

KIM_ID = "Systemic ID"  # sic: the published table spells it "Systemic"
HAYLES_ID = "Systematic ID"
HAYLES_STUDY = "Gene dispensability. This study"
HAYLES_PHENOTYPE = "Phenotypic classification used for analysis"
HAYLES_DESCRIPTION = "Deletion mutant phenotype description"

# Accepted spellings of the two dispensability calls.
DISPENSABILITY_CALLS: dict[str, str] = {
    "V": "V",
    "VIABLE": "V",
    "E": "E",
    "ESSENTIAL": "E",
}


# =============================================================================
# VALUE NORMALIZATION
# =============================================================================


def normalize_id(value: object) -> str | None:
    """Strip all whitespace from a systematic ID and drop missing values."""
    if value is None or value is pd.NA:
        return None
    text = re.sub(r"\s+", "", str(value))
    if not text or text.lower() == "nan":
        return None
    return text


def normalize_call(value: object) -> str | None:
    """Map viable/essential (and stray-spaced V/E) to ``V`` or ``E``."""
    if value is None or value is pd.NA:
        return None
    text = str(value).strip().upper()
    if not text or text == "NAN":
        return None
    call = DISPENSABILITY_CALLS.get(text)
    if call is None:
        logger.warning(f"Unrecognised dispensability call {value!r}; treating as missing")
    return call


def clean_text_column(series: pd.Series) -> pd.Series:
    """Strip a text column, turning empty strings into missing values."""
    return series.astype("string").str.strip().replace("", pd.NA)


# =============================================================================
# LOADERS
# =============================================================================


def load_kim2010(path: Path) -> pd.DataFrame:
    """Load Kim 2010 full Supplementary Table 1 into the live contract.

    The live table contributes the updated systematic ID and Kim's Paralog
    field; Kim phenotype calls are not merged.
    """
    df = pd.read_excel(path, sheet_name=0)
    df.columns = [str(c).strip() for c in df.columns]

    required = {KIM_ID, "Paralog"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")

    out = pd.DataFrame(
        {
            "id": df[KIM_ID].map(normalize_id),
            "paralogues_2010": clean_text_column(df["Paralog"]),
        }
    )
    if out["id"].isna().any():
        raise ValueError(f"{path}: {out['id'].isna().sum()} rows have no systematic ID")
    if out["id"].duplicated().any():
        dup = out.loc[out["id"].duplicated(), "id"].tolist()
        raise ValueError(f"{path}: duplicate systematic IDs {dup[:10]}")

    logger.info(f"Kim 2010 full Table S1: {len(out):,} rows")
    return out


def load_hayles2013(path: Path) -> pd.DataFrame:
    """Load the Hayles 2013 main supplementary table into the live contract."""
    df = pd.read_excel(path, sheet_name=0)
    df.columns = [str(c).strip() for c in df.columns]
    required = {
        HAYLES_ID,
        HAYLES_STUDY,
        HAYLES_PHENOTYPE,
        HAYLES_DESCRIPTION,
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")

    out = pd.DataFrame(
        {
            "id": df[HAYLES_ID].map(normalize_id),
            "hayles2013_study": df[HAYLES_STUDY].map(normalize_call),
            "phenotype_description": clean_text_column(df[HAYLES_DESCRIPTION]),
            "phenotype_classification": clean_text_column(df[HAYLES_PHENOTYPE]),
        }
    )
    if out["id"].isna().any():
        raise ValueError(f"{path}: {out['id'].isna().sum()} rows have no systematic ID")
    if out["id"].duplicated().any():
        dup = out.loc[out["id"].duplicated(), "id"].tolist()
        raise ValueError(f"{path}: duplicate systematic IDs {dup[:10]}")
    logger.info(f"Hayles 2013: {len(out):,} genes")
    return out


# =============================================================================
# DUPLICATE ID COLLAPSE
# =============================================================================


def collapse_duplicate_ids(
    df: pd.DataFrame,
    id_col: str,
    original_id_col: str,
) -> pd.DataFrame:
    """Collapse rows sharing a current ID, combining their original IDs."""
    if df[id_col].is_unique:
        return df

    rows = []
    for gene_id, group in df.groupby(id_col, sort=False):
        if len(group) == 1:
            rows.append(group.iloc[0])
            continue

        # Prefer the row whose original ID is already the current ID.
        primary_mask = group[original_id_col].astype("string").str.upper() == str(gene_id).upper()
        primary = group.loc[primary_mask].iloc[0] if primary_mask.any() else group.iloc[0]
        row = primary.copy()

        orig_ids = list(dict.fromkeys(str(value) for value in group[original_id_col]))
        row[original_id_col] = "; ".join(orig_ids)

        logger.warning(f"Collapsed {len(group)} rows onto {gene_id}: {orig_ids}")
        rows.append(row)

    return pd.DataFrame(rows).reset_index(drop=True)
