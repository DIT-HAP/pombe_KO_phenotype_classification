"""
Merge Categorized Phenotypes
=============================

Concatenates the three phenotype branches (one / multi / inconsistent) from
the categorized output into a single table, applies manually revised category
merges (Sub_category → Category), and re-ranks ``Growth_tier`` by DIT-HAP
DR median per merged category (DR is sign-flipped: most depleted = tier 1).

Input
-----
- ``data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx``
  (3 data sheets + All genes, output of 04_categorize_phenotypes.py)

- ``data/4_categorized_genes/category_revisions.json``
  Hand-curated category config (plot order + ``Sub_category -> Category``
  merge mapping), shared with step 06.

- ``data/raw/gene_level_fitting_statistics_minusDR.tsv``
  Sign-flipped DIT-HAP DR values used for Growth_tier re-ranking
  (negative = more depleted).

Output
------
- ``data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx``
  Single merged table + summary sheets (Category, Sub_category,
  Consistency_25_32, Phenotype_count, Growth_tier).
- ``results/Hayles_2013_OB_merged_categories.xlsx`` (flat single-sheet copy
  of ``All genes`` for source control)

Usage
-----
    mamba run -n bioinformatics python scripts/05_merge_categories.py

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-06-09
Version:  1.2.0
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

# 5. Local Imports
from category_revisions import (  # noqa: E402
    DEFAULT_CONFIG,
    CategoryConfig,
    load_category_config,
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

DEFAULT_INPUT = Path("data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx")
DEFAULT_OUTPUT = Path("data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx")
DEFAULT_RESULTS = Path("results/Hayles_2013_OB_merged_categories.xlsx")
DEFAULT_DIT_HAP = Path("data/raw/gene_level_fitting_statistics_minusDR.tsv")

DATA_SHEETS: tuple[str, ...] = (
    "One basic phenotype",
    "Multi basic phenotypes",
    "Inconsistent phenotypes",
)
REQUIRED_COLUMNS: tuple[str, ...] = (
    "Systematic ID",
    "Category",
    "Consistency_25_32",
    "Phenotype_count",
)

# =============================================================================
# CORE LOGIC
# =============================================================================


def load_and_concat(input_path: Path) -> pd.DataFrame:
    """Load all 3 data sheets, validate their columns, and concatenate them."""
    sheets = pd.read_excel(input_path, sheet_name=None)
    require_sheets(sheets, DATA_SHEETS, input_path)

    frames: list[pd.DataFrame] = []
    for sheet_name in DATA_SHEETS:
        df = sheets[sheet_name]
        require_columns(df, REQUIRED_COLUMNS, input_path)
        logger.info(f"  {sheet_name}: {len(df)} genes")
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)
    logger.info(f"Total: {len(merged)} genes")
    return merged


def apply_category_merge(
    merged: pd.DataFrame,
    config: CategoryConfig,
    dit_hap: pd.DataFrame,
) -> pd.DataFrame:
    """Rename fine-grained Category to Sub_category and apply revised merges.

    The revisions file maps ``Sub_category -> Category``. Sub-categories not
    listed in the mapping keep themselves as their merged Category.

    Growth_tier is reassigned based on the DR median of each merged Category,
    ranked from the most depleted (lowest, i.e. most negative DR) as tier 1
    to the least depleted as tier N.
    """
    # 1. Rename old fine-grained Category → Sub_category
    merged = merged.rename(columns={"Category": "Sub_category"})

    # 2. Apply: if sub-category is in the mapping, use revised; else keep it
    merged["Category"] = merged["Sub_category"].map(config.apply)

    logger.info(
        f"  {len(config.revisions)} sub-categories revised, "
        f"{merged['Category'].nunique()} merged categories"
    )

    # 3. Compute median DR per merged Category (1 = most depleted)
    tier_map, cat_dr_median = dr_tier_map(merged, dit_hap)
    logger.info(f"  DR medians computed for {len(cat_dr_median)} categories")

    merged["Growth_tier"] = merged["Category"].map(tier_map)

    # Log the ranking
    for cat, med in zip(cat_dr_median.index, cat_dr_median, strict=True):
        logger.info(f"    Tier {tier_map[cat]:2d}: {cat:50s} med DR = {med:.4f}")

    return merged


def build_summaries(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Build summary sheets for the output workbook."""
    summaries: dict[str, pd.DataFrame] = {}
    for column in (
        "Consistency_25_32",
        "Phenotype_count",
        "Category",
        "Sub_category",
        "Growth_tier",
    ):
        counts = df.value_counts(column)
        if column == "Growth_tier":
            counts = counts.sort_index()
        summaries[column] = counts.rename("Count").to_frame()
    return summaries


# =============================================================================
# MAIN EXECUTION
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge categorized phenotypes and re-rank Growth_tier by DR median.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Categorized phenotypes xlsx (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output xlsx path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=DEFAULT_RESULTS,
        help=f"Results-copy xlsx path (default: {DEFAULT_RESULTS})",
    )
    parser.add_argument(
        "--revisions",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"Category config JSON (default: {DEFAULT_CONFIG})",
    )
    parser.add_argument(
        "--dit-hap",
        type=Path,
        default=DEFAULT_DIT_HAP,
        help=f"DIT-HAP TSV (default: {DEFAULT_DIT_HAP})",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable DEBUG level logging")
    return parser.parse_args()


def build_merged(
    input_path: Path,
    output_path: Path,
    results_path: Path,
    revisions_path: Path,
    dit_hap_path: Path,
) -> None:
    """Merge the categorized branches, re-rank Growth_tier, and write both copies."""
    logger.info(f"Loading categorized phenotypes: {input_path}")
    merged = load_and_concat(input_path)

    logger.info("Applying revised category merges and re-ranking Growth_tier…")
    config = load_category_config(revisions_path)
    merged = apply_category_merge(merged, config, load_dit_hap(dit_hap_path))
    summaries = build_summaries(merged)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    results_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path) as writer:
        merged.to_excel(writer, sheet_name="All genes", index=False)
        for sheet_name, summary_df in summaries.items():
            summary_df.to_excel(writer, sheet_name=sheet_name)
    logger.success(f"Merged categories saved: {output_path}")

    merged.to_excel(results_path, index=False)
    logger.success(f"Results copy saved: {results_path}")


def main() -> int:
    """Run the category merge and propagate fatal errors as exit code 1."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    try:
        build_merged(
            args.input.resolve(),
            args.output.resolve(),
            args.results.resolve(),
            args.revisions.resolve(),
            args.dit_hap.resolve(),
        )
    except Exception as error:
        logger.exception(f"Category merge failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
