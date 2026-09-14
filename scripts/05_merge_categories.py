"""
Merge Categorized Phenotypes
=============================

Concatenates the three phenotype branches (one / multi / inconsistent) from
the categorized output into a single table, applies manually revised category
merges (Sub_category → Category), and re-ranks ``Growth_tier`` by DIT-HAP
DR median per merged category.

Input
-----
- ``data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx``
  (3 data sheets + All genes, output of 04_categorize_phenotypes.py)

- ``data/4_categorized_genes/Hayles_2013_OB_inspection_phenotypes_category_revised_20260707.xlsx``
  Manually revised category mappings (Revised column in Flat inspection sheet).

- ``data/references/all_coding_genes_with_DIT_HAP_clustering.tsv``
  DIT-HAP DR values used for Growth_tier re-ranking.

Output
------
- ``data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx``
  Single merged table + summary sheets (Category, Sub_category,
  Consistency_25_32, Phenotype_count, Growth_tier).
- ``results/Hayles_2013_OB_merged_categories.xlsx`` (copy for source control)

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
import sys
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_INPUT = Path(
    "data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx"
)
DEFAULT_OUTPUT = Path(
    "data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx"
)
DEFAULT_RESULTS = Path(
    "results/Hayles_2013_OB_merged_categories.xlsx"
)
DEFAULT_REVISED = Path(
    "data/4_categorized_genes/Hayles_2013_OB_inspection_phenotypes_category_revised_20260707.xlsx"
)
DEFAULT_DIT_HAP = Path(
    "data/references/all_coding_genes_with_DIT_HAP_clustering.tsv"
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


def apply_category_merge(
    merged: pd.DataFrame,
    revised_path: Path,
    dit_hap_path: Path,
) -> pd.DataFrame:
    """Rename fine-grained Category to Sub_category and apply revised merges.

    The revised file maps original description → revised Category.
    Genes whose description is not in the revised file keep their original
    Category as the merged Category.

    Growth_tier is reassigned based on the DR median of each merged Category,
    ranked from highest median (tier 1) to lowest (tier N).
    """
    # 1. Rename old fine-grained Category → Sub_category
    merged = merged.rename(columns={"Category": "Sub_category"})

    # 2. Load revised mapping: description → revised category
    rev = pd.read_excel(revised_path, sheet_name="Flat inspection")
    rev_map: dict[str, str] = {}
    for _, r in rev[rev["Revised"].notna()].iterrows():
        rev_map[r["Phenotype description"]] = str(r["Revised"])

    # 3. Apply: if description is in rev_map, use revised; else use Sub_category
    desc_col = "Deletion mutant phenotype description"
    merged["Category"] = merged[desc_col].map(rev_map).fillna(merged["Sub_category"])

    logger.info(f"  {len(rev_map)} descriptions revised, "
                f"{merged['Category'].nunique()} merged categories")

    # 4. Load DIT-HAP DR values
    dit_hap = pd.read_csv(dit_hap_path, sep="\t")
    dit_hap = dit_hap[["Systematic ID", "DR"]].dropna(subset=["DR"])

    # 5. Compute median DR per merged Category
    joined = merged.merge(dit_hap, on="Systematic ID", how="left")
    cat_dr_median = (
        joined.groupby("Category")["DR"]
        .median()
        .sort_values(ascending=False)
    )
    logger.info(f"  DR medians computed for {len(cat_dr_median)} categories")

    # 6. Assign Growth_tier: 1=highest median, 2=second, ...
    tier_map = {cat: i + 1 for i, cat in enumerate(cat_dr_median.index)}
    merged["Growth_tier"] = merged["Category"].map(tier_map)

    # Log the ranking
    for cat, med in cat_dr_median.items():
        logger.info(f"    Tier {tier_map[cat]:2d}: {cat:50s} med DR = {med:.4f}")

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

    summaries["Sub_category"] = (
        df.value_counts("Sub_category")
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

    # 2. Apply revised category merges and re-rank Growth_tier by DR median
    logger.info("Applying revised category merges and re-ranking Growth_tier…")
    merged = apply_category_merge(merged, DEFAULT_REVISED, DEFAULT_DIT_HAP)

    # 3. Build summaries
    summaries = build_summaries(merged)

    # 4. Save
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
