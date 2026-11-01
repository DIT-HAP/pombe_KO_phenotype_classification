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
- ``data/1_formatted/Hayles2013_phenotypes_annotated.xlsx``
  Formatted phenotype data (output of 01_format_hayles_phenotypes.py).

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
from pathlib import Path
import sys

# 2. Third-party Imports
from loguru import logger

# 3. Data Processing Imports
import numpy as np
import pandas as pd

# 4. Project Path Setup
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# 5. Project Imports
from growth_signals import count_growth_segments  # noqa: E402
from phenotype import PHENOTYPE_DESC_COL  # noqa: E402
from pipeline_utils import require_columns, setup_logger  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_INPUT = Path("data/1_formatted/Hayles2013_phenotypes_annotated.xlsx")
DEFAULT_OUTPUT = Path("data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx")

# Temperature marker used to detect consistency in phenotype descriptions
TEMP_BOTH = "25,32"
TEMP_32 = "32"

REQUIRED_COLUMNS: tuple[str, ...] = (PHENOTYPE_DESC_COL,)

# =============================================================================
# CORE LOGIC
# =============================================================================


def identify_consistency_groups(df: pd.DataFrame) -> dict[str, np.ndarray]:
    """Return boolean masks for the three consistency groups."""
    desc = df[PHENOTYPE_DESC_COL]

    both = desc.str.contains(TEMP_BOTH, na=False)
    only_32 = desc.str.contains(TEMP_32, na=False) & (~desc.str.contains("25", na=False))
    inconsistent = (~both) & (~only_32)

    return {
        "both": both.to_numpy(),
        "only_32": only_32.to_numpy(),
        "inconsistent": inconsistent.to_numpy(),
    }


def split_basic_additional(desc_series: pd.Series, marker: str) -> pd.DataFrame:
    """Split a phenotype description at *marker* into basic and additional parts.

    The basic phenotype is the text before the marker (with trailing " at"
    stripped); the additional phenotype is the text after the marker.
    """
    parts = desc_series.str.split(marker, expand=True)

    basic = parts[0].str.removesuffix(" at").str.strip()
    if parts.shape[1] > 1:
        additional = parts[1].str.lstrip(", ").str.strip()
    else:
        additional = pd.Series([None] * len(parts))

    return pd.DataFrame({"Basic phenotype": basic, "Additional phenotype": additional})


def classify_phenotype_count(phenotype: object) -> str | float:
    """Classify as ``Single`` or ``Multiple``, or NaN for non-string input.

    Counts how many comma-separated segments contain a growth signal
    (excluding modifier-led segments).  More than one → ``Multiple``.
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


def add_consistency_column(
    df: pd.DataFrame,
    groups: dict[str, np.ndarray],
) -> pd.DataFrame:
    """Label each gene's temperature consistency."""
    consistency = pd.Series(index=df.index, dtype="object")
    consistency.loc[groups["both"]] = "Consistent"
    consistency.loc[groups["only_32"]] = "Only_32"
    consistency.loc[groups["inconsistent"]] = "Mismatch"
    return df.assign(Consistency_25_32=consistency)


def add_phenotype_splits(
    df: pd.DataFrame,
    groups: dict[str, np.ndarray],
) -> pd.DataFrame:
    """Split descriptions into basic and additional phenotypes per group."""
    out = df.copy()

    both_index = out[groups["both"]].index
    both_split = split_basic_additional(
        out.loc[both_index, PHENOTYPE_DESC_COL],
        f" {TEMP_BOTH}",
    )
    out.loc[both_index, "Basic phenotype"] = both_split["Basic phenotype"]
    out.loc[both_index, "Additional phenotype"] = both_split["Additional phenotype"]

    only32_index = out[groups["only_32"]].index
    only32_split = split_basic_additional(
        out.loc[only32_index, PHENOTYPE_DESC_COL],
        f" {TEMP_32}",
    )
    out.loc[only32_index, "Basic phenotype"] = only32_split["Basic phenotype"]
    out.loc[only32_index, "Additional phenotype"] = only32_split["Additional phenotype"]
    return out


def add_phenotype_count(
    df: pd.DataFrame,
    groups: dict[str, np.ndarray],
) -> pd.DataFrame:
    """Label single vs multiple growth-signal segments, or a temperature mismatch."""
    counts = df["Basic phenotype"].map(classify_phenotype_count)
    counts.loc[groups["inconsistent"]] = "Temp_mismatch"
    return df.assign(Phenotype_count=counts)


def log_summary(df: pd.DataFrame, groups: dict[str, np.ndarray]) -> None:
    """Log group sizes and phenotype-count totals."""
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


def write_grouped_workbook(output_path: Path, df: pd.DataFrame) -> None:
    """Write the four grouped-gene sheets."""
    with pd.ExcelWriter(output_path) as writer:
        df.to_excel(writer, sheet_name="All genes", index=False)
        df[df["Phenotype_count"] == "Single"].to_excel(
            writer,
            sheet_name="One basic phenotype",
            index=False,
        )
        df[df["Phenotype_count"] == "Multiple"].to_excel(
            writer,
            sheet_name="Multi basic phenotypes",
            index=False,
        )
        inconsistent_df = df[df["Phenotype_count"] == "Temp_mismatch"].drop(
            columns=["Basic phenotype", "Additional phenotype"],
        )
        inconsistent_df.to_excel(
            writer,
            sheet_name="Inconsistent phenotypes",
            index=False,
        )


def build_grouped_table(input_path: Path, output_path: Path) -> None:
    """Load the phenotype table, group it by temperature consistency, and write it."""
    logger.info(f"Loading formatted phenotypes: {input_path}")
    df = pd.read_excel(input_path)
    require_columns(df, REQUIRED_COLUMNS, input_path)
    logger.info(f"Loaded {len(df):,} genes")

    groups = identify_consistency_groups(df)
    df = add_consistency_column(df, groups)
    df = add_phenotype_splits(df, groups)
    df = add_phenotype_count(df, groups)
    log_summary(df, groups)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_grouped_workbook(output_path, df)
    logger.success(f"Grouped genes saved: {output_path}")


def main() -> int:
    """Run the grouping and propagate fatal errors as exit code 1."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        build_grouped_table(args.input.resolve(), args.output.resolve())
    except Exception as error:
        logger.exception(f"Gene grouping failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
