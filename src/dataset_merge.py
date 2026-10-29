"""
Shared Kim 2010 / Hayles 2013 loaders
=====================================

Library module used by ``01_format_hayles_phenotypes.py`` for the two source
loaders and duplicate-ID collapse. It canonicalises the raw column names and
values the live pipeline consumes; it does not build the retired 25-column
union.

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
Date:     2026-10-29
Version: 2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import re
from datetime import datetime
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_KIM2010 = Path("data/references/kim2010/kim2010_supp_table1.xls")
DEFAULT_HAYLES2013 = Path(
    "data/references/hayles2013/hayles2013_supp_table3.xlsx"
)

KIM_ID = "Systemic ID"  # sic: the published table spells it "Systemic"
KIM_DELETION = "Deletion dataset"
KIM_ANALYSIS = "Analysis dataset"
HAYLES_ID = "Systematic ID"
HAYLES_STUDY = "Gene dispensability. This study"
HAYLES_PHENOTYPE = "Phenotypic classification used for analysis"


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
    if text in {"V", "VIABLE"}:
        return "V"
    if text in {"E", "ESSENTIAL"}:
        return "E"
    logger.warning(f"Unrecognised dispensability call {value!r}; treating as missing")
    return None


def clean_text_column(series: pd.Series) -> pd.Series:
    """Strip a text column, turning empty strings into missing values."""
    return series.astype("string").str.strip().replace("", pd.NA)


def optional_column(df: pd.DataFrame, name: str) -> pd.Series:
    """Return a column, or an all-missing column of the same length."""
    if name in df.columns:
        return df[name]
    return pd.Series(pd.NA, index=df.index, dtype="string")


def reconstruct_gene_name_from_date(value: datetime) -> str:
    """Rebuild an Excel auto-converted gene name such as ``Oct-1`` from its date.

    Excel reads ``Oct-1`` as 1 Oct of the current year, but reads ``Sep-11``
    as 1 Sep 2011 (the 11 becomes a year). When the day is 1 and the year is
    after 2001, the intended day is therefore ``year - 2000``.
    """
    day = value.year - 2000 if value.day == 1 and value.year > 2001 else value.day
    return f"{value.strftime('%b')}-{day}"


# =============================================================================
# LOADERS
# =============================================================================


@logger.catch(reraise=True)
def load_kim2010(path: Path) -> pd.DataFrame:
    """Load Kim 2010 full Supplementary Table 1 into canonical columns.

    Every constructed gene has a ``kim2010_analysis`` call; the 78
    ``Constructed? = no`` rows mostly have no call (two have a deletion-only
    partial call).
    """
    df = pd.read_excel(path, sheet_name=0)
    df.columns = [str(c).strip() for c in df.columns]

    if "Constructed?" not in df.columns:
        raise ValueError(
            f"{path}: this is not the full Kim 2010 Supplementary Table 1 "
            "(missing 'Constructed?'). The 1:1-ortholog Table 11 subset is "
            "not sufficient for the complete reconciliation."
        )
    return _load_kim_full(df, path)


def _load_kim_full(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Canonicalize Kim 2010 full Supplementary Table 1."""
    required = [KIM_ID, "Constructed?", KIM_DELETION, KIM_ANALYSIS, "Paralog"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    out = pd.DataFrame(
        {
            "id": df[KIM_ID].map(normalize_id),
            "gene_name_2010": clean_text_column(df["Gene name"]),
            "pct_ko_2010": pd.to_numeric(df["KO%"], errors="coerce"),
            "kim2010_deletion": df[KIM_DELETION].map(normalize_call),
            "kim2010_analysis": df[KIM_ANALYSIS].map(normalize_call),
            "published_phenotype": clean_text_column(df["Published phenotype"]),
            "constructed": clean_text_column(df["Constructed?"])
            .str.lower()
            .eq("yes"),
            "construction_method": clean_text_column(
                optional_column(df, "Construction method")
            ),
            "gene_product_2010": clean_text_column(
                optional_column(df, "Gene function")
            ),
            "paralogues_2010": clean_text_column(df["Paralog"]),
            "sc_orthologue_2010": clean_text_column(
                optional_column(df, "Sc ortholog")
            ),
            "ts_remake": clean_text_column(
                optional_column(
                    df,
                    "Heterozygous diploid containing ts mutation (being remade)",
                )
            ),
            "source_format": "full_table_s1",
        }
    )
    if out["id"].isna().any():
        raise ValueError(f"{path}: {out['id'].isna().sum()} rows have no systematic ID")
    if out["id"].duplicated().any():
        dup = out.loc[out["id"].duplicated(), "id"].tolist()
        raise ValueError(f"{path}: duplicate systematic IDs {dup[:10]}")

    out["kim2010_call"] = out["kim2010_analysis"].where(
        out["kim2010_analysis"].notna(), out["kim2010_deletion"]
    )
    logger.info(
        f"Kim 2010 full Table S1: {len(out):,} rows; "
        f"{int(out['constructed'].sum()):,} constructed"
    )
    return out


@logger.catch(reraise=True)
def load_hayles2013(path: Path) -> pd.DataFrame:
    """Load the Hayles 2013 main supplementary table."""
    df = pd.read_excel(path, sheet_name=0)
    df.columns = [str(c).strip() for c in df.columns]
    required = [
        HAYLES_ID,
        HAYLES_STUDY,
        HAYLES_PHENOTYPE,
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    # Excel auto-converted three gene names into dates (Oct-1, Sep-11, Sep-1);
    # reconstruct them before the column is cast to strings.
    gene_names = df["Gene name"].map(
        lambda value: reconstruct_gene_name_from_date(value)
        if isinstance(value, datetime)
        else value
    )

    out = pd.DataFrame(
        {
            "id": df[HAYLES_ID].map(normalize_id),
            "gene_name_2013": clean_text_column(gene_names),
            "working_code": clean_text_column(df["Working code"]),
            "pct_ko_2013": pd.to_numeric(df["% KO"], errors="coerce"),
            "hayles2013_study": df[HAYLES_STUDY].map(normalize_call),
            "phenotype_description": clean_text_column(
                df["Deletion mutant phenotype description"]
            ),
            "phenotype_classification": clean_text_column(df[HAYLES_PHENOTYPE]),
            "gene_product_2013": clean_text_column(df["Gene product description"]),
            "paralogues_2013": clean_text_column(df["Paralogues"]),
            "cell_length": clean_text_column(
                optional_column(
                    df, "Cell length >15.6 (st dev) for 34 novel long genes"
                )
            ),
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
    id_col: str = "id",
    original_id_col: str = "original_id",
) -> pd.DataFrame:
    """Collapse rows sharing a current ID, combining their original IDs."""
    if original_id_col not in df.columns:
        df = df.assign(**{original_id_col: df[id_col]})
    if df[id_col].is_unique:
        return df

    rows = []
    for gene_id, group in df.groupby(id_col, sort=False):
        if len(group) == 1:
            rows.append(group.iloc[0])
            continue

        # Prefer the row whose original ID is already the current ID.
        primary_mask = (
            group[original_id_col].astype("string").str.upper()
            == str(gene_id).upper()
        )
        primary = (
            group.loc[primary_mask].iloc[0] if primary_mask.any() else group.iloc[0]
        )
        row = primary.copy()

        orig_ids = list(dict.fromkeys(str(value) for value in group[original_id_col]))
        row[original_id_col] = "; ".join(orig_ids)

        if "constructed" in group.columns:
            row["constructed"] = bool(group["constructed"].fillna(False).any())

        # Take a non-null gene name / ts-remake flag if the primary row lacks it.
        for col in ("gene_name_2010", "updated_gene_name_2010", "ts_remake"):
            if col in group.columns and pd.isna(row[col]):
                non_null = group[col].dropna()
                if not non_null.empty:
                    row[col] = non_null.iloc[0]

        for col in ("kim2010_deletion", "kim2010_analysis"):
            if col in group.columns:
                values = group[col].dropna().unique().tolist()
                if len(values) > 1:
                    logger.warning(f"Conflicting {col} for {gene_id}: {values}")

        logger.warning(f"Collapsed {len(group)} rows onto {gene_id}: {orig_ids}")
        rows.append(row)

    return pd.DataFrame(rows).reset_index(drop=True)
