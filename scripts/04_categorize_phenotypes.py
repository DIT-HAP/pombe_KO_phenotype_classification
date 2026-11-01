"""
Categorise Growth Phenotypes
=============================

Reads grouped phenotype data (3 sheets from step 03), assigns each gene a
fine-grained ``Category`` using the signal-based ``classify_growth()``
engine, then re-ranks ``Growth_tier`` by DIT-HAP DR median (most depleted
= tier 1; DR is sign-flipped, so more negative = more depleted).

Input
-----
- ``data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx``
  3 sheets: One basic phenotype, Multi basic phenotypes,
  Inconsistent phenotypes.

- ``data/references/previous_manual_check_of_insistent_phenotypes/
   Inconsistent_phenotypes_at_25_32_manual.xlsx``
  Manual annotations for the 131 inconsistent genes.

- ``data/raw/gene_level_fitting_statistics_minusDR.tsv``
  Sign-flipped DIT-HAP DR values used for Growth_tier re-ranking
  (negative = more depleted).

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
    - Per-branch Phenotypes, Essentiality, Classification, Growth_tier pivots
    - ``Multi-level pivot (All genes)`` — full description text x 4-level
      column hierarchy with conditional formatting and binary-signature sorting

Usage
-----
    mamba run -n bioinformatics python scripts/04_categorize_phenotypes.py
    mamba run -n bioinformatics python scripts/04_categorize_phenotypes.py --verbose

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-10-28
Version:  1.2.1
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
from pathlib import Path
import sys

# 2. Third-party Imports
from loguru import logger

# 3. Data Processing Imports
import pandas as pd

# 4. Project Path Setup
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# 5. Local imports
from growth_signals import SIGNAL_TIERS, classify_growth  # noqa: E402
from phenotype import (  # noqa: E402
    DISPENSABILITY_COL,
    MORPHOLOGY_COL,
    PHENOTYPE_DESC_COL,
)
from pipeline_utils import (  # noqa: E402
    dr_tier_map,
    load_dit_hap,
    require_columns,
    require_sheets,
    setup_logger,
)

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

# Paths relative to project root
DEFAULT_GROUPED = Path("data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx")
DEFAULT_MANUAL = Path(
    "data/references/previous_manual_check_of_insistent_phenotypes"
    "/Inconsistent_phenotypes_at_25_32_manual.xlsx"
)
DEFAULT_OUTPUT = Path("data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx")
DEFAULT_DIT_HAP = Path("data/raw/gene_level_fitting_statistics_minusDR.tsv")

DATA_SHEETS: tuple[str, ...] = (
    "One basic phenotype",
    "Multi basic phenotypes",
    "Inconsistent phenotypes",
)
GROUPED_COLUMNS: tuple[str, ...] = (
    "Systematic ID",
    "Phenotype_count",
    "Consistency_25_32",
    DISPENSABILITY_COL,
    MORPHOLOGY_COL,
    PHENOTYPE_DESC_COL,
)
BASIC_COLUMNS: tuple[str, ...] = ("Basic phenotype",)
MANUAL_COLUMNS: tuple[str, ...] = ("SysID", "Category_25", "Category_32")

DESC_COL = PHENOTYPE_DESC_COL

# =============================================================================
# CORE LOGIC
# =============================================================================


def rerank_growth_tier_by_dr(df: pd.DataFrame, dit_hap: pd.DataFrame) -> pd.DataFrame:
    """Reassign Growth_tier based on DR median per Category, most depleted first.

    DR is sign-flipped (more negative = more depleted), so the lowest median
    becomes tier 1.
    """
    tier_map, _ = dr_tier_map(df, dit_hap)
    tiers = df["Category"].map(tier_map).astype(int)
    logger.info(f"  Growth_tier re-ranked by DR median ({len(tier_map)} categories)")
    return df.assign(Growth_tier=tiers)


def classify_basic_phenotype(df: pd.DataFrame) -> pd.DataFrame:
    """Assign Category and Growth_tier from the ``Basic phenotype`` column.

    Used for both the one-phenotype and multi-phenotype branches —
    ``classify_growth`` handles modifier-led segments internally.
    """
    results = df["Basic phenotype"].map(classify_growth)
    df = df.copy()
    df[["Category", "Growth_tier"]] = pd.DataFrame(results.tolist(), index=df.index)
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
    require_columns(manual, MANUAL_COLUMNS, manual_path)

    df = df.merge(
        manual[["SysID", "Category_25", "Category_32"]],
        left_on="Systematic ID",
        right_on="SysID",
        how="left",
    ).drop(columns="SysID")

    # Use Category_32 as the final category, normalising legacy "WT" → "WT-like"
    df["Category"] = df["Category_32"].fillna("WT-like").replace({"WT": "WT-like"})

    # Derive Growth_tier from category name
    df["Growth_tier"] = df["Category"].map(SIGNAL_TIERS).fillna(5).astype(int)

    return df


def build_pivot(df: pd.DataFrame, row_col: str, col_col: str) -> pd.DataFrame:
    """Build a pivot table counting rows by two categorical columns."""
    return pd.crosstab(df[row_col], df[col_col])


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


def load_grouped(path: Path) -> dict[str, pd.DataFrame]:
    """Load the grouped workbook and validate its sheets and required columns."""
    grouped = pd.read_excel(path, sheet_name=None)
    require_sheets(grouped, DATA_SHEETS, path)
    for sheet_name in DATA_SHEETS:
        require_columns(grouped[sheet_name], GROUPED_COLUMNS, path)
    for sheet_name in DATA_SHEETS[:2]:
        require_columns(grouped[sheet_name], BASIC_COLUMNS, path)
    return grouped


def classify_branches(
    grouped: dict[str, pd.DataFrame],
    manual_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Classify the one-, multi-, and inconsistent-phenotype branches."""
    logger.info(f"Processing 'One basic phenotype' ({len(grouped['One basic phenotype'])} genes)…")
    one = classify_basic_phenotype(grouped["One basic phenotype"])

    logger.info(
        f"Processing 'Multi basic phenotypes' ({len(grouped['Multi basic phenotypes'])} genes)…"
    )
    multi = classify_basic_phenotype(grouped["Multi basic phenotypes"])

    logger.info(
        f"Processing 'Inconsistent phenotypes' ({len(grouped['Inconsistent phenotypes'])} genes)…"
    )
    inconsistent = classify_inconsistent_phenotype(
        grouped["Inconsistent phenotypes"],
        manual_path,
    )
    return one, multi, inconsistent


def build_pivots(
    branches: dict[str, pd.DataFrame],
    all_genes: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Build one pivot sheet per branch and dimension, plus the combined branch."""
    pivots: dict[str, pd.DataFrame] = {}
    for label, frame in branches.items():
        pivots[f"Phenotypes pivot ({label})"] = build_pivot(
            frame,
            DESC_COL,
            "Category",
        )
        pivots[f"Essentiality pivot ({label})"] = build_pivot(
            frame,
            DISPENSABILITY_COL,
            "Category",
        )
        pivots[f"Classification pivot ({label})"] = build_pivot(
            frame,
            MORPHOLOGY_COL,
            "Category",
        )
        pivots[f"Growth_tier pivot ({label})"] = build_pivot(
            frame,
            MORPHOLOGY_COL,
            "Growth_tier",
        )
    pivots["Phenotypes pivot (All genes)"] = build_pivot(
        all_genes,
        DESC_COL,
        "Category",
    )
    return pivots


def split_reranked_branches(
    all_genes: pd.DataFrame,
    n_one: int,
    n_multi: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Slice the concatenated, re-ranked table back into its three branches."""
    return (
        all_genes.iloc[:n_one].copy(),
        all_genes.iloc[n_one : n_one + n_multi].copy(),
        all_genes.iloc[n_one + n_multi :].copy(),
    )


def build_inspection_pivot(all_genes: pd.DataFrame) -> pd.DataFrame:
    """Pivot gene counts by description and the four classification dimensions."""
    dimensions = [
        DESC_COL,
        "Growth_tier",
        "Category",
        "Phenotype_count",
        "Consistency_25_32",
    ]
    counts = all_genes.groupby(dimensions).agg(Count=(DESC_COL, "size")).reset_index()
    multi_pivot = counts.pivot_table(
        index=DESC_COL,
        columns=["Growth_tier", "Category", "Phenotype_count", "Consistency_25_32"],
        values="Count",
        fill_value=0,
    ).astype(int)

    # Growth_tier ascending (1 to 5), then alphabetical for lower levels.
    multi_pivot = multi_pivot.sort_index(axis=1)
    # Within each column's non-zero block, the largest count sorts first; the
    # leftmost column is the primary key and all-zero rows sink to the bottom.
    sort_cols = list(multi_pivot.columns)
    return multi_pivot.sort_values(
        by=sort_cols,
        ascending=[False] * len(sort_cols),
    )


def build_flat_inspection(all_genes: pd.DataFrame) -> pd.DataFrame:
    """Build one row per unique description and classification combination."""
    dimensions = [
        DESC_COL,
        "Category",
        "Phenotype_count",
        "Growth_tier",
        "Consistency_25_32",
    ]
    gene_counts = (
        all_genes.groupby(dimensions)
        .agg(**{"Gene count": (DESC_COL, "size")})
        .reset_index()
        .rename(
            columns={
                DESC_COL: "Phenotype description",
                "Phenotype_count": "Single/Multiple/Temp_mismatch",
            }
        )
    )
    ordered = gene_counts[
        [
            "Phenotype description",
            "Category",
            "Single/Multiple/Temp_mismatch",
            "Consistency_25_32",
            "Growth_tier",
            "Gene count",
        ]
    ]
    return ordered.sort_values(
        by=["Category", "Single/Multiple/Temp_mismatch", "Growth_tier", "Gene count"],
        ascending=[True, True, True, False],
    ).reset_index(drop=True)


def write_data_workbook(
    output_path: Path,
    one: pd.DataFrame,
    multi: pd.DataFrame,
    inconsistent: pd.DataFrame,
    all_genes: pd.DataFrame,
) -> None:
    """Write the four data sheets of the categorized-phenotype workbook."""
    with pd.ExcelWriter(output_path) as writer:
        one.to_excel(writer, sheet_name="One basic phenotype", index=False)
        multi.to_excel(writer, sheet_name="Multi basic phenotypes", index=False)
        inconsistent.to_excel(writer, sheet_name="Inconsistent phenotypes", index=False)
        all_genes.to_excel(writer, sheet_name="All genes", index=False)


def write_inspection_workbook(
    inspection_path: Path,
    flat_inspection: pd.DataFrame,
    pivots: dict[str, pd.DataFrame],
    multi_pivot: pd.DataFrame,
) -> None:
    """Write the flat inspection table, branch pivots, and multi-level pivot."""
    styled_pivot = multi_pivot.style.map(
        lambda value: "background-color: #DCE6F1" if value > 0 else "",
    )
    with pd.ExcelWriter(inspection_path) as writer:
        flat_inspection.to_excel(writer, sheet_name="Flat inspection", index=False)
        for sheet_name, pivot_df in pivots.items():
            pivot_df.to_excel(writer, sheet_name=sheet_name[:31])
        styled_pivot.to_excel(writer, sheet_name="Multi-level pivot (All genes)")


def build_categorized(
    grouped_path: Path,
    manual_path: Path,
    output_path: Path,
) -> None:
    """Categorise the three phenotype branches and write both workbooks."""
    logger.info(f"Reading grouped phenotypes: {grouped_path}")
    grouped = load_grouped(grouped_path)

    one, multi, inconsistent = classify_branches(grouped, manual_path)
    branches = {
        "One basic phenotype": one,
        "Multi basic phenotypes": multi,
        "Inconsistent phenotypes": inconsistent,
    }
    all_genes = pd.concat([one, multi, inconsistent], ignore_index=True)
    pivots = build_pivots(branches, all_genes)

    logger.info("Re-ranking Growth_tier by DIT-HAP DR median…")
    all_genes = rerank_growth_tier_by_dr(all_genes, load_dit_hap(DEFAULT_DIT_HAP))
    one, multi, inconsistent = split_reranked_branches(
        all_genes,
        len(one),
        len(multi),
    )

    logger.info("Building multi-level inspection pivot…")
    multi_pivot = build_inspection_pivot(all_genes)
    logger.info("Building flat inspection table…")
    flat_inspection = build_flat_inspection(all_genes)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_data_workbook(output_path, one, multi, inconsistent, all_genes)
    logger.success(
        f"Categorized phenotypes saved: {output_path}\n"
        f"  One basic phenotype:    {len(one)} genes\n"
        f"  Multi basic phenotypes: {len(multi)} genes\n"
        f"  Inconsistent phenotypes: {len(inconsistent)} genes\n"
        f"  All genes:              {len(all_genes)} genes",
    )

    inspection_path = output_path.with_name(
        output_path.stem.replace("_categorized", "_inspection") + output_path.suffix,
    )
    write_inspection_workbook(inspection_path, flat_inspection, pivots, multi_pivot)
    n_sheets = 1 + len(pivots) + 1
    logger.success(f"Inspection pivots saved: {inspection_path}  ({n_sheets} sheets)")


def main() -> int:
    """Run the categorization and propagate fatal errors as exit code 1."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        build_categorized(
            args.grouped.resolve(),
            args.manual.resolve(),
            args.output.resolve(),
        )
    except Exception as error:
        logger.exception(f"Phenotype categorization failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
