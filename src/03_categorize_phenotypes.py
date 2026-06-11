"""
Categorise Growth Phenotypes
=============================

Reads grouped phenotype data (3 sheets from step 03), assigns each gene a
fine-grained ``Category`` and a coarse ``Growth_tier`` (1–5) using the
signal-based ``classify_growth()`` engine.

Replaces the three old ``categorize_genes_with_*_phenotypes.py`` scripts.

Input
-----
- ``data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx``
  3 sheets: One basic phenotype, Multi basic phenotypes,
  Inconsistent phenotypes.

- ``data/previous_manual_check_of_insistent_phenotypes/
   Inconsistent_phenotypes_at_25_32_manual.xlsx``
  Manual annotations for the 131 inconsistent genes.

Output
------
- ``data/3_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx``
  4 data sheets + pivot tables. Each data sheet includes ``Category``
  and ``Growth_tier`` columns.
    - ``One basic phenotype`` — single-phenotype, consistent at 25/32°C
    - ``Multi basic phenotypes`` — multi-phenotype, consistent at 25/32°C
    - ``Inconsistent phenotypes`` — temperature-inconsistent, manual annotation
    - ``All genes`` — all 4,843 genes concatenated
  Pivot tables use the full ``Deletion mutant phenotype description``
  as row labels (not ``Basic phenotype``), covering each branch plus
  an all‑genes combined view.

Usage
-----
    mamba run -n bioinformatics python src/02_categorize_phenotypes.py
    mamba run -n bioinformatics python src/02_categorize_phenotypes.py --verbose

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-06-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from pathlib import Path
from typing import Any

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 3b. Local module
from growth_signals import GROWTH_SIGNALS, classify_growth

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

# Paths relative to project root
DEFAULT_GROUPED = Path("data/2_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx")
DEFAULT_MANUAL = Path(
    "data/previous_manual_check_of_insistent_phenotypes"
    "/Inconsistent_phenotypes_at_25_32_manual.xlsx"
)
DEFAULT_OUTPUT = Path(
    "data/3_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx"
)

# Reverse lookup: known category name → growth tier.
# Derived from the canonical signal table; fallback to tier 5 (WT).
CATEGORY_TO_TIER: dict[str, int] = {}
for sig in GROWTH_SIGNALS:
    CATEGORY_TO_TIER[sig.category] = sig.tier

# =============================================================================
# LOGGING SETUP
# =============================================================================


def setup_logger(log_level: str = "INFO") -> None:
    """Configure the Loguru logger."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level:<8} | {message}",
        level=log_level,
    )


setup_logger()

# =============================================================================
# CORE LOGIC
# =============================================================================


def classify_one_phenotype(df: pd.DataFrame) -> pd.DataFrame:
    """Assign Category and Growth_tier for one-phenotype genes."""
    results = df["Basic phenotype"].apply(classify_growth)
    df = df.copy()
    df[["Category", "Growth_tier"]] = pd.DataFrame(
        results.tolist(), index=df.index
    )
    return df


def classify_multi_phenotype(df: pd.DataFrame) -> pd.DataFrame:
    """Assign Category and Growth_tier for multi-phenotype genes."""
    results = df["Basic phenotype"].apply(classify_growth)
    df = df.copy()
    df[["Category", "Growth_tier"]] = pd.DataFrame(
        results.tolist(), index=df.index
    )
    return df


def classify_inconsistent_phenotype(
    df: pd.DataFrame,
    manual_path: Path,
) -> pd.DataFrame:
    """Merge manual annotations and derive Category / Growth_tier.

    Uses the pre-annotated ``Category_32`` from the manual file.  If the
    merge produces an unrecognised category name, falls back to ``WT`` /
    tier 5.
    """
    logger.info(f"Loading manual annotations: {manual_path}")
    manual = pd.read_excel(manual_path)

    df = df.merge(
        manual[["SysID", "Category_25", "Category_32"]],
        left_on="Systematic ID",
        right_on="SysID",
        how="left",
    ).drop(columns="SysID")

    # Use Category_32 as the final category
    df["Category"] = df["Category_32"].fillna("WT")

    # Derive Growth_tier from category name
    df["Growth_tier"] = df["Category"].map(CATEGORY_TO_TIER).fillna(5).astype(int)

    return df


def build_pivot(df: pd.DataFrame, row_col: str, col_col: str) -> pd.DataFrame:
    """Build a pivot table counting rows by two categorical columns."""
    return (
        df[[row_col, col_col]]
        .value_counts()
        .rename("Count")
        .reset_index()
        .pivot(index=row_col, columns=col_col, values="Count")
        .fillna(0)
        .astype(int)
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Categorise growth phenotypes using signal-based detection.",
    )
    parser.add_argument(
        "--grouped",
        type=Path,
        default=DEFAULT_GROUPED,
        help=f"Grouped genes xlsx (default: {DEFAULT_GROUPED})",
    )
    parser.add_argument(
        "--manual",
        type=Path,
        default=DEFAULT_MANUAL,
        help=f"Manual annotations for inconsistent genes (default: {DEFAULT_MANUAL})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output xlsx path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable DEBUG level logging",
    )
    return parser.parse_args()


def main() -> int:
    """Main orchestrator."""
    args = parse_args()
    if args.verbose:
        setup_logger(log_level="DEBUG")

    grouped_path = args.grouped.resolve()
    manual_path = args.manual.resolve()
    output_path = args.output.resolve()

    logger.info(f"Reading grouped phenotypes: {grouped_path}")
    grouped = pd.read_excel(grouped_path, sheet_name=None)

    # ------------------------------------------------------------------
    # 1. One basic phenotype
    # ------------------------------------------------------------------
    logger.info(
        f"Processing 'One basic phenotype' ({len(grouped['One basic phenotype'])} genes)…"
    )
    one = classify_one_phenotype(grouped["One basic phenotype"])

    # ------------------------------------------------------------------
    # 2. Multi basic phenotypes
    # ------------------------------------------------------------------
    logger.info(
        f"Processing 'Multi basic phenotypes' ({len(grouped['Multi basic phenotypes'])} genes)…"
    )
    multi = classify_multi_phenotype(grouped["Multi basic phenotypes"])

    # ------------------------------------------------------------------
    # 3. Inconsistent phenotypes
    # ------------------------------------------------------------------
    logger.info(
        f"Processing 'Inconsistent phenotypes' ({len(grouped['Inconsistent phenotypes'])} genes)…"
    )
    inconsistent = classify_inconsistent_phenotype(
        grouped["Inconsistent phenotypes"],
        manual_path,
    )

    # ------------------------------------------------------------------
    # 4. Build pivot tables
    # ------------------------------------------------------------------
    pivots: dict[str, pd.DataFrame] = {}

    for label, df in [
        ("One basic phenotype", one),
        ("Multi basic phenotypes", multi),
        ("Inconsistent phenotypes", inconsistent),
    ]:
        pivots[f"Phenotypes pivot ({label})"] = build_pivot(
            df, "Deletion mutant phenotype description", "Category"
        )
        pivots[f"Essentiality pivot ({label})"] = build_pivot(
            df, "Gene dispensability. This study", "Category"
        )
        pivots[f"Classification pivot ({label})"] = build_pivot(
            df, "Phenotypic classification used for analysis", "Category"
        )
        pivots[f"Growth_tier pivot ({label})"] = build_pivot(
            df, "Phenotypic classification used for analysis", "Growth_tier"
        )

    # Also build a combined "All genes" pivot for the full description
    all_genes = pd.concat([one, multi, inconsistent], ignore_index=True)
    pivots["Phenotypes pivot (All genes)"] = build_pivot(
        all_genes, "Deletion mutant phenotype description", "Category"
    )

    # ------------------------------------------------------------------
    # 5. Save
    # ------------------------------------------------------------------
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_path) as writer:
        # Data sheets (3 branches + all genes)
        one.to_excel(writer, sheet_name="One basic phenotype", index=False)
        multi.to_excel(writer, sheet_name="Multi basic phenotypes", index=False)
        inconsistent.to_excel(writer, sheet_name="Inconsistent phenotypes", index=False)
        all_genes.to_excel(writer, sheet_name="All genes", index=False)
        for sheet_name, pivot_df in pivots.items():
            # Truncate Excel sheet names to 31 chars
            safe_name = sheet_name[:31]
            pivot_df.to_excel(writer, sheet_name=safe_name)

    logger.success(
        f"Categorized phenotypes saved: {output_path}\n"
        f"  One basic phenotype:    {len(one)} genes\n"
        f"  Multi basic phenotypes: {len(multi)} genes\n"
        f"  Inconsistent phenotypes: {len(inconsistent)} genes\n"
        f"  Pivot tables: {len(pivots)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
