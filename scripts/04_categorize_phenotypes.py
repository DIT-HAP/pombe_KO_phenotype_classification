"""
Categorise Growth Phenotypes
=============================

Reads grouped phenotype data (3 sheets from step 03), assigns each gene a
fine-grained ``Category`` using the signal-based ``classify_growth()``
engine, then re-ranks ``Growth_tier`` by DIT-HAP DR median (highest median
= tier 1).

Input
-----
- ``data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx``
  3 sheets: One basic phenotype, Multi basic phenotypes,
  Inconsistent phenotypes.

- ``data/previous_manual_check_of_insistent_phenotypes/
   Inconsistent_phenotypes_at_25_32_manual.xlsx``
  Manual annotations for the 131 inconsistent genes.

- ``data/references/all_coding_genes_with_DIT_HAP_clustering.tsv``
  DIT-HAP DR values used for Growth_tier re-ranking.

Output
------
- ``data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx``
  Main data output — 4 sheets with ``Category`` and ``Growth_tier`` columns:
    - ``One basic phenotype`` — single-phenotype, consistent at 25/32°C
    - ``Multi basic phenotypes`` — multi-phenotype, consistent at 25/32°C
    - ``Inconsistent phenotypes`` — temperature-inconsistent, manual annotation
    - ``All genes`` — all 4,843 genes concatenated

- ``data/4_categorized_genes/Hayles_2013_OB_inspection_phenotypes.xlsx``
  Inspection pivot tables — separate from the main output for easy review:
    - ``Flat inspection`` — one row per unique description with gene counts
    - Per‑branch Phenotypes, Essentiality, Classification, Growth_tier pivots
    - ``Multi-level pivot (All genes)`` — full description text × 4‑level
      column hierarchy with conditional formatting and binary‑signature sorting

Usage
-----
    mamba run -n bioinformatics python src/04_categorize_phenotypes.py
    mamba run -n bioinformatics python src/04_categorize_phenotypes.py --verbose

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-06-09
Version:  1.2.0
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
DEFAULT_GROUPED = Path("data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx")
DEFAULT_MANUAL = Path(
    "data/previous_manual_check_of_insistent_phenotypes"
    "/Inconsistent_phenotypes_at_25_32_manual.xlsx"
)
DEFAULT_OUTPUT = Path(
    "data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx"
)
DEFAULT_DIT_HAP = Path(
    "data/references/all_coding_genes_with_DIT_HAP_clustering.tsv"
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


def rerank_growth_tier_by_dr(df: pd.DataFrame, dit_hap_path: Path) -> pd.DataFrame:
    """Reassign Growth_tier based on DR median per Category, ranked high→low."""
    if not dit_hap_path.exists():
        logger.warning(f"DIT-HAP file not found: {dit_hap_path}, keeping signal-based tiers")
        return df
    dit_hap = pd.read_csv(dit_hap_path, sep="\t")
    dit_hap = dit_hap[["Systematic ID", "DR"]].dropna(subset=["DR"])
    joined = df.merge(dit_hap, on="Systematic ID", how="left")
    cat_dr_median = (
        joined.groupby("Category")["DR"]
        .median()
        .sort_values(ascending=False)
    )
    tier_map = {cat: i + 1 for i, cat in enumerate(cat_dr_median.index)}
    df["Growth_tier"] = df["Category"].map(tier_map).fillna(len(tier_map) + 1).astype(int)
    logger.info(f"  Growth_tier re-ranked by DR median ({len(tier_map)} categories)")
    return df


def classify_one_phenotype(df: pd.DataFrame) -> pd.DataFrame:
    """Assign Category and Growth_tier for one-phenotype genes.

    ``classify_growth`` handles modifier‑led segments internally —
    primary signals get plain names, modifier‑led signals get prefixed
    names (e.g. ``often divided``).
    """
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
    merge produces an unrecognised category name, falls back to ``WT-like``.
    Growth_tier is later re-ranked by DIT-HAP DR median.
    """
    logger.info(f"Loading manual annotations: {manual_path}")
    manual = pd.read_excel(manual_path)

    df = df.merge(
        manual[["SysID", "Category_25", "Category_32"]],
        left_on="Systematic ID",
        right_on="SysID",
        how="left",
    ).drop(columns="SysID")

    # Use Category_32 as the final category, normalising legacy "WT" → "WT-like"
    df["Category"] = df["Category_32"].fillna("WT-like").replace({"WT": "WT-like"})

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

    # Re-rank Growth_tier by DR median (unified with 05_merge_categories)
    logger.info("Re-ranking Growth_tier by DIT-HAP DR median…")
    all_genes = rerank_growth_tier_by_dr(all_genes, DEFAULT_DIT_HAP)
    # Apply same re-ranking to individual branches
    n_one = len(one)
    n_multi = len(multi)
    one = all_genes.iloc[:n_one].copy()
    multi = all_genes.iloc[n_one:n_one + n_multi].copy()
    inconsistent = all_genes.iloc[n_one + n_multi:].copy()
    pivots["Phenotypes pivot (All genes)"] = build_pivot(
        all_genes, "Deletion mutant phenotype description", "Category"
    )

    # ------------------------------------------------------------------
    # 4b. Multi-level pivot for quality inspection
    #     Rows: full description text.
    #     Columns: (Growth_tier, Category, Phenotype_count, Consistency_25_32).
    #     Cells: gene count. Conditional format: colour only non-zero cells.
    # ------------------------------------------------------------------
    logger.info("Building multi-level inspection pivot…")

    # Group counts by all four dimensions
    grouped = (
        all_genes
        .groupby([
            "Deletion mutant phenotype description",
            "Growth_tier",
            "Category",
            "Phenotype_count",
            "Consistency_25_32",
        ])
        .size()
        .rename("Count")
        .reset_index()
    )

    # Pivot with four-level column index
    multi_pivot = grouped.pivot_table(
        index="Deletion mutant phenotype description",
        columns=["Growth_tier", "Category", "Phenotype_count", "Consistency_25_32"],
        values="Count",
        fill_value=0,
    ).astype(int)

    # Sort columns: Growth_tier ascending (1→5), then alphabetical for lower levels
    multi_pivot = multi_pivot.sort_index(axis=1)

    # Sort rows so that within each column's non-zero block, values are
    # sorted descending (largest / most representative at the top).
    # Columns are prioritised left-to-right (Growth_tier 1 then 2 then 3…):
    # column 0 is the primary sort key, column 1 the secondary, etc.
    # All-zero rows sink to the bottom.
    sort_cols = list(multi_pivot.columns)
    multi_pivot = multi_pivot.sort_values(
        by=sort_cols,
        ascending=[False] * len(sort_cols),
    )

    styled_pivot = multi_pivot.style.map(
        lambda v: "background-color: #DCE6F1" if v > 0 else "",
    )

    # ------------------------------------------------------------------
    # 4c. Flat inspection table — one row per unique description + gene count
    # ------------------------------------------------------------------
    logger.info("Building flat inspection table…")
    flat_full = all_genes[
        [
            "Deletion mutant phenotype description",
            "Category",
            "Phenotype_count",
            "Growth_tier",
            "Consistency_25_32",
        ]
    ].copy()

    # Count genes per (description, Category, Phenotype_count, Growth_tier, Consistency)
    gene_counts = (
        flat_full.groupby(
            [
                "Deletion mutant phenotype description",
                "Category",
                "Phenotype_count",
                "Growth_tier",
                "Consistency_25_32",
            ]
        )
        .size()
        .rename("Gene count")
        .reset_index()
    )

    flat_inspection = gene_counts.rename(
        columns={
            "Deletion mutant phenotype description": "Phenotype description",
            "Phenotype_count": "Single/Multiple/Temp_mismatch",
        }
    )
    # Column order: description, Category, S/M/T, Consistency, Growth_tier, Gene count
    flat_inspection = flat_inspection[
        [
            "Phenotype description",
            "Category",
            "Single/Multiple/Temp_mismatch",
            "Consistency_25_32",
            "Growth_tier",
            "Gene count",
        ]
    ]
    flat_inspection.sort_values(
        by=["Category", "Single/Multiple/Temp_mismatch", "Growth_tier", "Gene count"],
        ascending=[True, True, True, False],
        inplace=True,
    )
    flat_inspection.reset_index(drop=True, inplace=True)

    # ------------------------------------------------------------------
    # 5. Save
    # ------------------------------------------------------------------
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Main data output: data sheets only
    with pd.ExcelWriter(output_path) as writer:
        one.to_excel(writer, sheet_name="One basic phenotype", index=False)
        multi.to_excel(writer, sheet_name="Multi basic phenotypes", index=False)
        inconsistent.to_excel(writer, sheet_name="Inconsistent phenotypes", index=False)
        all_genes.to_excel(writer, sheet_name="All genes", index=False)

    logger.success(
        f"Categorized phenotypes saved: {output_path}\n"
        f"  One basic phenotype:    {len(one)} genes\n"
        f"  Multi basic phenotypes: {len(multi)} genes\n"
        f"  Inconsistent phenotypes: {len(inconsistent)} genes\n"
        f"  All genes:              {len(all_genes)} genes",
    )

    # Inspection output: pivot tables + multi-level pivot
    inspection_path = output_path.with_name(
        output_path.stem.replace("_categorized", "_inspection") + output_path.suffix,
    )
    with pd.ExcelWriter(inspection_path) as writer:
        # Flat inspection table first — most useful for manual review
        flat_inspection.to_excel(
            writer, sheet_name="Flat inspection", index=False,
        )

        for sheet_name, pivot_df in pivots.items():
            safe_name = sheet_name[:31]
            pivot_df.to_excel(writer, sheet_name=safe_name)

        # Multi-level inspection pivot (styled)
        styled_pivot.to_excel(
            writer, sheet_name="Multi-level pivot (All genes)",
        )

    n_sheets = 1 + len(pivots) + 1
    logger.success(f"Inspection pivots saved: {inspection_path}  ({n_sheets} sheets)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
