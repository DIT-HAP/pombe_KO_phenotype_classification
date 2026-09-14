"""
Group Genes by Phenotype Consistency at Different Temperatures
===============================================================

Splits the formatted phenotype table into three groups based on whether
the phenotype description is consistent across temperatures, reported
only at 32°C, or inconsistent between 25°C and 32°C.

For consistent entries, the phenotype description is further split into
a ``Basic phenotype`` (before the temperature marker) and an
``Additional phenotype`` (after the marker), which downstream scripts
use for growth classification.

The output adds two classification columns:
- ``Consistency_25_32`` — ``Consistent`` (description spans both
  25°C and 32°C), ``Only_32`` (described only at 32°C), or
  ``Mismatch`` (descriptions differ between 25°C and 32°C).
- ``Phenotype_count`` — ``Single`` (single phenotype), ``Multiple``
  (comma-separated phenotypes), or ``Temp_mismatch`` (cannot be split
  due to inconsistent temperature descriptions).

This grouping strategy avoids the need to manually resolve temperature
differences for the 97% of genes whose descriptions are consistent or
only available at 32°C. Only the 3% inconsistent entries require manual
curation.

Input
-----
- ``data/1_formatted/Hayles_2013_OB_formatted_phenotypes.xlsx``
  Formatted phenotype data (output of 01_format_and_update_ids.py).

Output
------
- ``data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx``
  4 sheets:
    - ``All genes`` — all 4,843 genes with consistency and phenotype split
      columns added
    - ``One basic phenotype`` — genes with a single, consistent phenotype
      at both temperatures
    - ``Multi basic phenotypes`` — genes with multiple comma-separated
      phenotypes at both temperatures
    - ``Inconsistent phenotypes`` — genes with different descriptions at
      25°C and 32°C

Usage
-----
    mamba run -n bioinformatics python scripts/03_group_genes.py
    mamba run -n bioinformatics python scripts/03_group_genes.py --verbose

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-06-09
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Project Imports
from growth_signals import count_growth_segments
from pipeline_utils import setup_logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_INPUT = Path("data/1_formatted/Hayles_2013_OB_formatted_phenotypes.xlsx")
DEFAULT_OUTPUT = Path("data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx")

# Temperature marker used to detect consistency in phenotype descriptions
TEMP_BOTH = "25,32"
TEMP_32 = "32"

# =============================================================================
# LOGGING SETUP
# =============================================================================


setup_logger()

# =============================================================================
# CORE LOGIC
# =============================================================================


@logger.catch
def identify_consistency_groups(df: pd.DataFrame) -> dict[str, np.ndarray]:
    """Return boolean masks for the three consistency groups.

    Returns
    -------
    dict with keys 'both', 'only_32', 'inconsistent', each a boolean array
    of length ``len(df)``.
    """
    desc = df["Deletion mutant phenotype description"]

    both = desc.str.contains(TEMP_BOTH, na=False)
    only_32 = desc.str.contains(TEMP_32, na=False) & (~desc.str.contains("25", na=False))
    inconsistent = (~both) & (~only_32)

    return {"both": both.values, "only_32": only_32.values, "inconsistent": inconsistent.values}


@logger.catch
def split_basic_additional(desc_series: pd.Series, marker: str) -> pd.DataFrame:
    """Split a phenotype description at *marker* into basic and additional parts.

    The basic phenotype is the text before the marker (with trailing " at"
    stripped); the additional phenotype is the text after the marker.
    """
    parts = desc_series.str.split(marker, expand=True)

    basic = parts[0].str.rstrip(" at")
    additional = parts[1].str.lstrip(", ").str.strip() if parts.shape[1] > 1 else pd.Series([None] * len(parts))

    return pd.DataFrame({"Basic phenotype": basic, "Additional phenotype": additional})


@logger.catch
def classify_phenotype_count(phenotype: object) -> str | float:
    """Classify as ``Single`` or ``Multiple``, or NaN for non‑string input.

    Counts how many comma‑separated segments contain a growth signal
    (excluding modifier‑led segments).  More than one → ``Multiple``.
    """
    if not isinstance(phenotype, str):
        return np.nan
    return "Multiple" if count_growth_segments(phenotype) > 1 else "Single"


# =============================================================================
# MAIN EXECUTION
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Group genes by phenotype temperature consistency.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Formatted phenotypes xlsx (default: {DEFAULT_INPUT})",
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

    input_path = args.input.resolve()
    output_path = args.output.resolve()

    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------
    logger.info(f"Loading formatted phenotypes: {input_path}")
    df = pd.read_excel(input_path)
    logger.info(f"Loaded {len(df):,} genes")

    # ------------------------------------------------------------------
    # 2. Identify consistency groups
    # ------------------------------------------------------------------
    groups = identify_consistency_groups(df)
    df["Consistency_25_32"] = pd.Series(index=df.index, dtype="object")
    df.loc[groups["both"], "Consistency_25_32"] = "Consistent"
    df.loc[groups["only_32"], "Consistency_25_32"] = "Only_32"
    df.loc[groups["inconsistent"], "Consistency_25_32"] = "Mismatch"

    # ------------------------------------------------------------------
    # 3. Split phenotype descriptions
    # ------------------------------------------------------------------
    desc_col = "Deletion mutant phenotype description"

    # Consistent at both temperatures — split on "25,32"
    both_idx = df[groups["both"]].index
    both_split = split_basic_additional(df.loc[both_idx, desc_col], f" {TEMP_BOTH}")
    df.loc[both_idx, "Basic phenotype"] = both_split["Basic phenotype"]
    df.loc[both_idx, "Additional phenotype"] = both_split["Additional phenotype"]

    # Only at 32°C — split on " 32"
    only32_idx = df[groups["only_32"]].index
    only32_split = split_basic_additional(df.loc[only32_idx, desc_col], f" {TEMP_32}")
    df.loc[only32_idx, "Basic phenotype"] = only32_split["Basic phenotype"]
    df.loc[only32_idx, "Additional phenotype"] = only32_split["Additional phenotype"]

    # ------------------------------------------------------------------
    # 4. Classify one vs multi basic phenotypes
    # ------------------------------------------------------------------
    df["Phenotype_count"] = (
        df["Basic phenotype"].apply(classify_phenotype_count)
    )

    # For temperature‑inconsistent genes, the phenotype cannot be split into
    # basic/additional, so assign an explicit label instead of NaN.
    df.loc[groups["inconsistent"], "Phenotype_count"] = (
        "Temp_mismatch"
    )

    # ------------------------------------------------------------------
    # 5. Summary statistics
    # ------------------------------------------------------------------
    logger.info(f"Total genes: {len(df):,}")
    logger.info(f"  Consistent at both temp:  {groups['both'].sum():,}")
    logger.info(f"  Only at 32°C:             {groups['only_32'].sum():,}")
    logger.info(f"  Inconsistent:             {groups['inconsistent'].sum():,}")

    n_one = (df["Phenotype_count"] == "Single").sum()
    n_multi = (df["Phenotype_count"] == "Multiple").sum()
    n_temp_inconsistent = (df["Phenotype_count"] == "Temp_mismatch").sum()
    logger.info(f"  One basic phenotype:      {n_one:,}")
    logger.info(f"  Multi basic phenotypes:   {n_multi:,}")
    logger.info(f"  Temperature mismatch:     {n_temp_inconsistent:,}")

    # ------------------------------------------------------------------
    # 6. Save
    # ------------------------------------------------------------------
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_path) as writer:
        # Sheet 1: All genes (total overview)
        df.to_excel(writer, sheet_name="All genes", index=False)

        # Sheet 2: One basic phenotype
        df[df["Phenotype_count"] == "Single"].to_excel(
            writer, sheet_name="One basic phenotype", index=False,
        )
        # Sheet 3: Multi basic phenotypes
        df[df["Phenotype_count"] == "Multiple"].to_excel(
            writer, sheet_name="Multi basic phenotypes", index=False,
        )
        # Sheet 4: Inconsistent phenotypes (drop split columns)
        inconsistent_df = df[df["Phenotype_count"] == "Temp_mismatch"].drop(
            columns=["Basic phenotype", "Additional phenotype"],
            errors="ignore",
        )
        inconsistent_df.to_excel(writer, sheet_name="Inconsistent phenotypes", index=False)

    logger.success(f"Grouped genes saved: {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
