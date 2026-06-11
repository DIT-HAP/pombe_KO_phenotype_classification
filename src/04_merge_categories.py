"""
Merge Categorized Phenotypes and Assign Essentiality
=====================================================

Concatenates the three phenotype branches (one / multi / inconsistent) from
the categorized output into a single table, preserving both the fine-grained
``Category`` and the coarse ``Growth_tier`` columns for downstream analysis.

Replaces the original script that read three separate categorized files.

Input
-----
- ``data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx``
  (3 data sheets, output of 02_categorize_phenotypes.py)

Output
------
- ``data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx``
  Single merged table + summary sheets.
- ``results/Hayles_2013_OB_merged_categories.xlsx`` (copy for source control)

Usage
-----
    mamba run -n bioinformatics python src/03_merge_categories.py

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-06-09
Version:  1.1.0
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

DEFAULT_INPUT = Path(
    "data/3_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx"
)
DEFAULT_OUTPUT = Path(
    "data/4_merged_categories/Hayles_2013_OB_merged_categories.xlsx"
)
DEFAULT_RESULTS = Path(
    "results/Hayles_2013_OB_merged_categories.xlsx"
)

# =============================================================================
# LOGGING SETUP
# =============================================================================

logger.remove()
logger.add(
    sys.stdout,
    format="{time:YYYY-MM-DD HH:mm:ss} | {level:<8} | {message}",
)

# =============================================================================
# CORE LOGIC
# =============================================================================


def load_and_concat(input_path: Path) -> pd.DataFrame:
    """Load all 3 data sheets and concatenate into one DataFrame."""
    sheets = pd.read_excel(input_path, sheet_name=None)

    data_sheets = [
        "One basic phenotype",
        "Multi basic phenotypes",
        "Inconsistent phenotypes",
    ]

    frames: list[pd.DataFrame] = []
    for sn in data_sheets:
        df = sheets[sn]
        logger.info(f"  {sn}: {len(df)} genes")
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)
    logger.info(f"Total: {len(merged)} genes")
    return merged


def build_summaries(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Build summary sheets for the output workbook."""
    summaries: dict[str, pd.DataFrame] = {}

    summaries["Consistency_25_32"] = (
        df.value_counts("Consistency_25_32")
        .rename("Count")
        .to_frame()
    )

    summaries["Phenotype_count"] = (
        df.value_counts("Phenotype_count")
        .rename("Count")
        .to_frame()
    )

    summaries["Category"] = (
        df.value_counts("Category")
        .rename("Count")
        .to_frame()
    )

    summaries["Growth_tier"] = (
        df.value_counts("Growth_tier")
        .sort_index()
        .rename("Count")
        .to_frame()
    )

    return summaries


# =============================================================================
# MAIN EXECUTION
# =============================================================================


def main() -> int:
    """Main orchestrator."""
    input_path = DEFAULT_INPUT.resolve()
    output_path = DEFAULT_OUTPUT.resolve()
    results_path = DEFAULT_RESULTS.resolve()

    # 1. Load and merge
    logger.info(f"Loading categorized phenotypes: {input_path}")
    merged = load_and_concat(input_path)

    # 2. Build summaries
    summaries = build_summaries(merged)

    # 3. Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    results_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_path) as writer:
        merged.to_excel(writer, sheet_name="All genes", index=False)
        for sheet_name, summary_df in summaries.items():
            summary_df.to_excel(writer, sheet_name=sheet_name)

    logger.success(f"Merged categories saved: {output_path}")

    # 4. Copy to results/ for source control
    merged.to_excel(results_path, index=False)
    logger.success(f"Results copy saved: {results_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
