"""
Categorise Growth Phenotypes
=============================

Reads grouped phenotype data (3 sheets from step 02), assigns each gene a
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
- ``data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx``
  Main data output — 4 sheets with ``Category`` and ``Growth_tier`` columns:
    - ``One basic phenotype`` — single-phenotype, consistent at 25/32°C
    - ``Multi basic phenotypes`` — multi-phenotype, consistent at 25/32°C
    - ``Inconsistent phenotypes`` — temperature-inconsistent, manual annotation
    - ``All genes`` — all 4,843 genes concatenated

- ``data/4_categorized_genes/Hayles_2013_OB_inspection_pivot.xlsx``
  Inspection pivot tables — separate from the main output for easy review:
    - Per‑branch Phenotypes, Essentiality, Classification, Growth_tier pivots
    - ``Multi-level pivot (All genes)`` — full description text × 4‑level
      column hierarchy with conditional formatting and binary‑signature sorting

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
DEFAULT_GROUPED = Path("data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx")
DEFAULT_MANUAL = Path(
    "data/previous_manual_check_of_insistent_phenotypes"
    "/Inconsistent_phenotypes_at_25_32_manual.xlsx"
)
DEFAULT_OUTPUT = Path(
    "data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx"
)

# Modifier words — when a comma‑separated segment starts with one of these,
# it is treated as a secondary description and excluded from classification.
# Must match the list in 03_group_genes.py.
MODIFIER_WORDS = (
    "occasionally", "often", "occasional", "sometimes",
    "mostly", "rarely", "frequently", "frequency", "rare",
    "possible", "may", "possibly", "rapidly", "initially",
    "slightly", "very", "highly", "barely", "slight", "high", "weak",
    "some", "many", "few", "lots", "several", "multiple",
    "once", "twice", "more", "multi",
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
    """Assign Category and Growth_tier for one-phenotype genes.

    When a comma‑separated segment starts with a modifier word (e.g.
    ``some``, ``occasionally``, ``slightly``), it is treated as a
    secondary supplement and excluded from growth‑signal detection.
    All other comma segments are parallel growth phenotypes and are
    kept for classification.
    """
    def _filter_primary(text: str) -> str:
        """Keep only comma segments that do NOT start with a modifier word."""
        if not isinstance(text, str) or "," not in text:
            return text
        segments = [seg.strip() for seg in text.split(",")]
        primary = [seg for seg in segments
                   if not any(seg.lower().startswith(m) for m in MODIFIER_WORDS)]
        return ", ".join(primary) if primary else segments[0]

    basic = df["Basic phenotype"].apply(_filter_primary)
    results = basic.apply(classify_growth)
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
    df["Category"] = df["Category_32"].fillna("WT-like")

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
        for sheet_name, pivot_df in pivots.items():
            safe_name = sheet_name[:31]
            pivot_df.to_excel(writer, sheet_name=safe_name)

        # Multi-level inspection pivot (styled)
        styled_pivot.to_excel(
            writer, sheet_name="Multi-level pivot (All genes)",
        )

    logger.success(f"Inspection pivots saved: {inspection_path}  ({len(pivots) + 1} sheets)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
