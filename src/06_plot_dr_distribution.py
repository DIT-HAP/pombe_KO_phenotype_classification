"""
DR / um Distribution by Phenotype Category
==========================================

Reads the merged phenotype classification and joins it with two external
datasets — DIT-HAP (DR/DL) and gRNA (um/lam) — then plots the distribution
of depletion-rate proxies (DR for DIT-HAP, um for gRNA) across phenotype
categories using boxplot + violinplot with statistics annotations.

Input
-----
- ``data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx``
  (sheet ``All genes``) — must contain ``Systematic ID`` and ``Category``.
- ``data/references/all_coding_genes_with_DIT_HAP_clustering.tsv``
  — must contain ``Systematic ID`` and ``DR``.
- ``data/references/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv``
  — must contain ``Systematic ID`` and ``um``.

Output
------
- ``results/DR_distribution_DIT-HAP.png``
- ``results/um_distribution_gRNA.png``

Usage
-----
    mamba run -n bioinformatics python src/06_plot_dr_distribution.py
    mamba run -n bioinformatics python src/06_plot_dr_distribution.py --verbose

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-07-02
Version:  1.0.0
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from loguru import logger
from scipy.stats import mannwhitneyu

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_MERGED = Path("data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx")
DEFAULT_DIT_HAP = Path("data/references/all_coding_genes_with_DIT_HAP_clustering.tsv")
DEFAULT_GRNA = Path("data/references/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv")
DEFAULT_REVISED = Path("data/4_categorized_genes/Hayles_2013_OB_inspection_phenotypes_category_revised_20260707.xlsx")
DEFAULT_OUTPUT_DIR = Path("results")

AX_WIDTH = 5
AX_HEIGHT = 8

VALUE_RANGE = (-0.3, 2.0)

CATEGORY_COLOR_MAP: dict[str, str] = {
    "WT-like": "#4DAF4A",
    "spores": "#377EB8",
    "germinated": "#984EA3",
    "germinated, spores": "#7B68A6",
    "spores, germinated, divided or microcolonies": "#6B8DB5",
    "some germinated, spores": "#9E7BB5",
    "germinated and divided": "#F781BF",
    "germinated, divided or microcolonies": "#E7298A",
    "germinated and divided, spores": "#E7298A",
    "germinated, some divided": "#FF7F0E",
    "germinated, often divided": "#FFA040",
    "germinated, occasionally divided": "#FFB870",
    "germinated, some divided, spores": "#E67E22",
    "germinated, occasionally divided, spores": "#D35400",
    "germinated, often divided, spores": "#C0392B",
    "germinated, occasionally microcolonies, spores": "#8E44AD",
    "germinated, some microcolonies, spores": "#9B59B6",
    "germinated, occasionally microcolonies": "#A569BD",
    "germinated, some microcolonies": "#BB8FCE",
    "germinated, microcolonies": "#A67800",
    "germinated, microcolonies, spores": "#BF9B30",
    "microcolonies": "#FFD92F",
    "microcolonies, small colonies": "#E6AB02",
    "microcolonies, spores": "#D4B020",
    "some microcolonies, spores": "#C8B860",
    "microcolonies, occasionally germinated, occasionally spores": "#D0B850",
    "small colonies": "#A65628",
    "very small colonies": "#7A3D1B",
    "germinated, small colonies, spores": "#B07843",
    "septated": "#999999",
    "?": "#666666",
    "spores, germinated spores": "#5B7C99",
    "spores, germinated spores and colonies": "#4A6B85",
    "spores, germinated spores and divided": "#3A5A75",
    "some divided, some germinated, spores": "#6B8DB5",
    "germinated, occasional microcolonies, spores": "#8B7500",
    "small colonies, some microcolonies": "#8B4513",
}

CATEGORY_ORDER: list[str] = [
    "spores",
    "spores, some germinated",
    "spores, some germinated, some divided",
    "spores, germinated",
    "spores, germinated, some divided",
    "spores, germinated, occasionally divided",
    "spores, germinated, often divided",
    "spores, germinated, divided or microcolonies",
    "spores, germinated and divided",
    "spores, germinated, some microcolonies",
    "spores, germinated, occasionally microcolonies",
    "spores, germinated, microcolonies",
    "spores, miscellaneous",
    "spores, germinated, small colonies",
    "spores, microcolonies",
    "spores, some microcolonies",
    "germinated",
    "germinated, some divided",
    "germinated, occasionally divided",
    "germinated, often divided",
    "germinated, divided or microcolonies",
    "germinated, some microcolonies",
    "germinated, occasionally microcolonies",
    "germinated, microcolonies",
    "microcolonies",
    "microcolonies, occasionally spores, occasionally germinated",
    "microcolonies, some spores, some germinated",
    "microcolonies, small colonies",
    "some microcolonies, small colonies",
    "very small colonies",
    "small colonies",
    "WT-like",
    "?",
    "septated",
]

# =============================================================================
# LOGGING SETUP
# =============================================================================


def setup_logger(log_level: str = "INFO") -> None:
    """Configure loguru with the given level."""
    logger.remove()
    logger.add(sys.stderr, level=log_level, format="<level>{level: <8}</level> | {message}")


# =============================================================================
# CORE LOGIC
# =============================================================================


def load_merged(path: Path) -> pd.DataFrame:
    """Load merged categories with Systematic ID, description, Sub_category, and Category."""
    df = pd.read_excel(path, sheet_name="All genes")
    cols = ["Systematic ID", "Deletion mutant phenotype description"]
    if "Sub_category" in df.columns:
        cols.append("Sub_category")
    cols.append("Category")
    return df[cols].copy()


def apply_revised(merged: pd.DataFrame, revised_path: Path) -> pd.DataFrame:
    """Override Category with manually revised values from the inspection file."""
    if not revised_path.exists():
        logger.warning(f"Revised file not found: {revised_path}")
        return merged
    rev = pd.read_excel(revised_path, sheet_name="Flat inspection")
    rev = rev[rev["Revised"].notna()][["Phenotype description", "Revised"]].copy()
    rev.columns = ["Deletion mutant phenotype description", "Revised"]
    merged = merged.merge(rev, on="Deletion mutant phenotype description", how="left")
    n_changed = merged["Revised"].notna().sum()
    merged["Category"] = merged["Revised"].fillna(merged["Category"])
    merged = merged.drop(columns=["Revised"])
    logger.info(f"  Applied {n_changed} revised categories")
    return merged


def load_dit_hap(path: Path) -> pd.DataFrame:
    """Load DIT-HAP TSV, keeping Systematic ID and DR."""
    df = pd.read_csv(path, sep="\t")
    return df[["Systematic ID", "DR"]].dropna(subset=["DR"]).copy()


def load_grna(path: Path) -> pd.DataFrame:
    """Load gRNA TSV, keeping Systematic ID and um."""
    df = pd.read_csv(path, sep="\t")
    return df[["Systematic ID", "um"]].dropna(subset=["um"]).copy()


def build_value_dict(merged: pd.DataFrame, values: pd.DataFrame, value_col: str) -> dict[str, list[float]]:
    """Join merged categories with a value column and group by Category."""
    joined = merged.merge(values, on="Systematic ID", how="inner")
    result: dict[str, list[float]] = {}
    # Use CATEGORY_ORDER first, then append any new categories not in the list
    all_cats = list(CATEGORY_ORDER) + [
        c for c in sorted(joined["Category"].unique())
        if c not in CATEGORY_ORDER
    ]
    for cat in all_cats:
        vals = joined.loc[joined["Category"] == cat, value_col].tolist()
        if vals:
            result[cat] = vals
    return result


def horizontal_violin_box(categories: list[str], data: list[list[float]], ax: plt.Axes, colors: list[str]) -> None:
    """Draw horizontal violin + boxplot; x=value, y=category index."""
    positions = range(len(categories))

    # KDE requires ≥2 points; for single-element lists, skip violin
    violin_data: list[list[float]] = []
    violin_pos: list[int] = []
    for i, vals in enumerate(data):
        clean = [v for v in vals if not np.isnan(v)]
        if len(clean) >= 2:
            violin_data.append(clean)
            violin_pos.append(i)

    if violin_data:
        parts = ax.violinplot(
            violin_data,
            positions=violin_pos,
            showmeans=False,
            showmedians=False,
            showextrema=False,
            vert=False,
        )
        for j, pc in enumerate(parts["bodies"]):
            idx = violin_pos[j]
            pc.set_facecolor(colors[idx])
            pc.set_alpha(0.3)
            pc.set_edgecolor(colors[idx])

    bp = ax.boxplot(
        data,
        positions=positions,
        vert=False,
        widths=0.15,
        patch_artist=True,
        showfliers=True,
        flierprops={"marker": ".", "markersize": 2, "alpha": 0.5},
    )
    for i, box in enumerate(bp["boxes"]):
        box.set_facecolor(colors[i])
        box.set_alpha(0.7)
    for element in ["whiskers", "caps", "medians"]:
        for item in bp[element]:
            item.set_color("black")
            item.set_linewidth(0.8)

    ax.set_yticks(positions)
    ax.set_yticklabels(categories, fontsize=7)
    ax.set_xlim(*VALUE_RANGE)


def build_grouped_order(
    dr_dict: dict[str, list[float]],
    um_dict: dict[str, list[float]],
    revised_path: Path | None = None,
) -> tuple[list[str], list[int]]:
    """Build category order, optionally grouped by revised category.

    Returns (ordered_categories, boundary_positions) where boundary_positions
    are the y-index positions after which a horizontal divider should be drawn.
    """
    present = [c for c in CATEGORY_ORDER if c in dr_dict or c in um_dict]
    if revised_path is None or not revised_path.exists():
        return present, []

    # Load revised mapping: original → revised
    rev = pd.read_excel(revised_path, sheet_name="Flat inspection")
    rev_map: dict[str, str] = {}
    for _, r in rev[rev["Revised"].notna()].iterrows():
        rev_map[r["Category"]] = r["Revised"]

    # Group present categories by their revised category (or self if no revision)
    groups: dict[str, list[str]] = {}
    for c in present:
        key = rev_map.get(c, c)
        groups.setdefault(key, []).append(c)

    # Order groups by the position of their first member in CATEGORY_ORDER
    group_keys = sorted(groups.keys(), key=lambda g: min(CATEGORY_ORDER.index(c) for c in groups[g]))

    ordered: list[str] = []
    boundaries: list[int] = []
    for g in group_keys:
        ordered.extend(groups[g])
        if g != group_keys[-1]:
            boundaries.append(len(ordered) - 1)  # draw line after this position

    return ordered, boundaries


def plot_combined(
    dr_dict: dict[str, list[float]],
    um_dict: dict[str, list[float]],
    output_path: Path,
    revised_path: Path | None = None,
    show_pvalues: bool = False,
) -> None:
    """Plot DR (DIT-HAP) and um (gRNA) side by side with statistics panel."""
    all_cats, boundaries = build_grouped_order(dr_dict, um_dict, revised_path)
    colors = [CATEGORY_COLOR_MAP.get(c, "gray") for c in all_cats]

    n_cats = len(all_cats)
    fig_height = max(AX_HEIGHT, n_cats * 0.35)
    fig, axes = plt.subplots(
        1, 3, figsize=(2 * AX_WIDTH + 3, fig_height),
        sharey=True, gridspec_kw={"width_ratios": [4, 4, 3]},
    )

    dr_data = [dr_dict.get(c, [np.nan]) for c in all_cats]
    horizontal_violin_box(all_cats, dr_data, axes[0], colors)
    axes[0].set_title("DR (DIT-HAP)")
    axes[0].set_xlabel("DR")

    um_data = [um_dict.get(c, [np.nan]) for c in all_cats]
    horizontal_violin_box(all_cats, um_data, axes[1], colors)
    axes[1].set_title("um (gRNA)")
    axes[1].set_xlabel("um")

    for row, cat in enumerate(all_cats):
        dr_vals = dr_dict.get(cat, [])
        um_vals = um_dict.get(cat, [])
        dr_str = (
            f"DR:  n={len(dr_vals):<4}, med={np.median(dr_vals):<6.3f}"
            if dr_vals else
            "DR:  n=   -, med=  -   "
        )
        um_str = (
            f"gRNA: n={len(um_vals):<4}, med={np.median(um_vals):<6.3f}"
            if um_vals else
            "gRNA: n=   -, med=  -   "
        )
        axes[2].text(0.05, row, f"{dr_str}  {um_str}", va="center", ha="left", fontsize=5.5, fontweight="bold", family="monospace")
    axes[2].set_ylim(-0.5, n_cats - 0.5)
    axes[2].invert_yaxis()
    axes[2].axis("off")

    # Draw horizontal divider lines between revised groups
    for b in boundaries:
        y = b + 0.5
        for ax in axes[:2]:
            ax.axhline(y=y, color="gray", linewidth=0.8, linestyle="--")

    # P-value annotations between adjacent categories
    if show_pvalues:
        _draw_pvalue_annotations(all_cats, dr_dict, axes[0], n_cats)
        _draw_pvalue_annotations(all_cats, um_dict, axes[1], n_cats)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.success(f"Saved: {output_path}")


def _draw_pvalue_annotations(
    cats: list[str],
    data_dict: dict[str, list[float]],
    ax: plt.Axes,
    n_cats: int,
) -> None:
    """Draw Mann-Whitney U p-values between adjacent category pairs."""
    xlim = ax.get_xlim()
    x_bracket = xlim[1] + 0.02  # just right of the data area
    for i in range(len(cats) - 1):
        vals_a = [v for v in data_dict.get(cats[i], []) if not np.isnan(v)]
        vals_b = [v for v in data_dict.get(cats[i + 1], []) if not np.isnan(v)]
        if len(vals_a) < 2 or len(vals_b) < 2:
            continue
        _, pval = mannwhitneyu(vals_a, vals_b, alternative="two-sided")
        bold = pval < 0.05
        y = i + 0.5
        # Vertical bracket connecting the two categories
        ax.plot([x_bracket, x_bracket], [i, i + 1], color="black", linewidth=0.5, clip_on=False)
        # Small horizontal ticks at each end
        tick_w = 0.01
        ax.plot([x_bracket, x_bracket + tick_w], [i, i], color="black", linewidth=0.5, clip_on=False)
        ax.plot([x_bracket, x_bracket + tick_w], [i + 1, i + 1], color="black", linewidth=0.5, clip_on=False)
        # P-value text to the right of the bracket
        pval_str = f"p={pval:.3f}" if pval >= 0.001 else "p<0.001"
        ax.text(
            x_bracket + tick_w + 0.01, y, pval_str,
            va="center", ha="left", fontsize=4.5,
            fontweight="bold" if bold else "normal",
            color="red" if bold else "gray",
            clip_on=False,
        )


# =============================================================================
# MAIN EXECUTION
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Plot DR/um distribution by phenotype category")
    parser.add_argument("--merged", type=Path, default=DEFAULT_MERGED, help=f"Merged categories xlsx (default: {DEFAULT_MERGED})")
    parser.add_argument("--dit-hap", type=Path, default=DEFAULT_DIT_HAP, help=f"DIT-HAP TSV (default: {DEFAULT_DIT_HAP})")
    parser.add_argument("--grna", type=Path, default=DEFAULT_GRNA, help=f"gRNA TSV (default: {DEFAULT_GRNA})")
    parser.add_argument("--revised", type=Path, default=DEFAULT_REVISED, help=f"Revised categories xlsx (default: {DEFAULT_REVISED})")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("--verbose", action="store_true", help="Enable DEBUG level logging")
    return parser.parse_args()


def main() -> int:
    """Main orchestrator."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")

    logger.info("Loading merged categories…")
    merged = load_merged(args.merged)
    logger.info(f"  {len(merged)} genes, {merged['Category'].nunique()} categories")

    logger.info("Loading DIT-HAP…")
    dit_hap = load_dit_hap(args.dit_hap)
    logger.info(f"  {len(dit_hap)} genes with DR values")

    logger.info("Loading gRNA…")
    grna = load_grna(args.grna)
    logger.info(f"  {len(grna)} genes with um values")

    # Plot 1: original fine-grained categories (Sub_category)
    logger.info("=== Plot 1: original (Sub_category) ===")
    merged_orig = merged.copy()
    merged_orig["Category"] = merged_orig["Sub_category"]
    dr_dict = build_value_dict(merged_orig, dit_hap, "DR")
    um_dict = build_value_dict(merged_orig, grna, "um")
    logger.info(f"  {len(dr_dict)} DR categories, {len(um_dict)} um categories")
    plot_combined(dr_dict, um_dict, args.output_dir / "DR_um_distribution_original.png", revised_path=args.revised)

    # Plot 2: merged categories (Category column already has revised merges)
    logger.info("=== Plot 2: revised (merged Category) ===")
    logger.info(f"  {merged['Category'].nunique()} merged categories")
    dr_dict_rev = build_value_dict(merged, dit_hap, "DR")
    um_dict_rev = build_value_dict(merged, grna, "um")
    logger.info(f"  {len(dr_dict_rev)} DR categories, {len(um_dict_rev)} um categories")
    plot_combined(dr_dict_rev, um_dict_rev, args.output_dir / "DR_um_distribution_revised.png", show_pvalues=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())