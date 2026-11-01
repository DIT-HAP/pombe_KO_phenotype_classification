#!/usr/bin/env python3
"""
DR / um Distribution — Pure Phenotype Categories
==================================================

``cnsplots`` figure restricted to five *pure* phenotype categories — rows
whose ``Sub_category`` is exactly one growth stage rather than a mixed
description. It follows the revised figure of step 06: two horizontal violin
panels (DIT-HAP ``DR`` and gRNA ``um``, both sign-flipped so more negative =
more depleted) plus a table panel with sample size and median per row.

Rows, in draw order:

1. ``spores`` — ``Sub_category == "spores"``
2. ``germinated`` — ``Sub_category == "germinated"``
3. ``microcolonies`` — ``Sub_category == "microcolonies"``
4. ``small colonies`` — ``Sub_category == "small colonies"``
5. ``WT-like (WT)`` — ``Sub_category == "WT-like"`` and
   ``Morphological_classification == "WT"``

Every adjacent pair is annotated with a two-sided Mann-Whitney U p-value
(raw, no multiple-comparison correction — the same choice as step 06).
The bracket column, label format and row layout match the revised cnsplots
figure. Each violin axes is 100 px wide at the 288 dpi export (25 pt); height
grows with the row count and text keeps cnsplots' default sizes. Colours come
from the Cell journal qualitative palette: significant
p-values (p < 0.05) in the palette terracotta and the rest in the palette
grey. The five row colours are the first five Cell colours in draw order,
which keeps them maximally distinct. Both violin panels repeat the row
labels, each with that panel's own sample size in parentheses on the same
line, e.g. ``spores (n = 143)``.

Input
-----
- ``data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx``
  (sheet ``All genes``) — ``Systematic ID``, ``Sub_category``, and
  ``Morphological_classification``.
- ``data/raw/gene_level_fitting_statistics_minusDR.tsv``
  — ``Systematic ID`` and sign-flipped ``DR`` (negative = depleted).
- ``data/raw/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv``
  — ``Systematic ID`` and ``um`` (negated on load).

Output
------
- ``results/DR_um_distribution_pure.{png,pdf}``

Usage
-----
    mamba run -n cnsplots python scripts/07_plot_dr_distribution_pure.py
    mamba run -n cnsplots python scripts/07_plot_dr_distribution_pure.py --verbose

Author:   Yusheng Yang (guidance) + Claude Code (implementation)
Date:     2026-10-27
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================

# 1. Standard Library Imports
import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

# 2. Third-party Imports
import cnsplots as cns
from cnsplots._utils import apply_unicode_font as _cns_apply_unicode_font
from loguru import logger
import matplotlib as mpl
from matplotlib.axes import Axes
import matplotlib.colors as mcolors

# 3. Data Processing Imports
import pandas as pd

# 4. Project Path Setup
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# 5. Local Imports
from dr_plotting import (  # noqa: E402
    CANVAS_SLACK_PT,
    DEFAULT_DIT_HAP,
    DEFAULT_GRNA,
    DEFAULT_MERGED,
    DEFAULT_OUTPUT_DIR,
    FONT_FAMILY,
    ID_COL,
    LABEL_CHAR_PT,
    MARGIN_RIGHT_PT,
    MERGED_SHEET,
    PAD_LEFT_PT,
    PAD_TOP_PT,
    PANEL_WIDTH_PT,
    PVALUE_MARGIN_PT,
    PVALUE_NOTE,
    ROW_COL,
    STATS_PAD_PT,
    STATS_SLACK_PT,
    SUB_CATEGORY_COL,
    VALUE_COL,
    VALUE_RANGE,
    FigureSpec,
    Panel,
    adjacent_pairs,
    annotate_pairs,
    build_long_frame,
    draw_stats_panel,
    log_category_stats,
    log_pvalues,
    stats_columns,
)
from phenotype import MORPHOLOGY_COL  # noqa: E402
from pipeline_utils import (  # noqa: E402
    load_dit_hap,
    load_grna,
    require_columns,
    setup_logger,
)

# Select the headless backend before any figure is created.
mpl.use("Agg")

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================

# Cell journal qualitative palette. The five rows take the first five colours
# in draw order, which keeps adjacent rows maximally distinct.
PALETTE = [mcolors.to_hex(tuple(color)) for color in cns.palettes("Cell")]
PVALUE_SIGNIFICANT_COLOR = PALETTE[0]
PVALUE_NEUTRAL_COLOR = PALETTE[9]

# Panel geometry. The violin axes are 100 px wide in the 288 dpi export
# (PANEL_WIDTH_PT = 25 pt); height grows with the number of rows. Text keeps
# cnsplots' default sizes, so a row only has to clear one 7 pt tick label.
ROW_HEIGHT_PT = 16.89

REQUIRED_COLUMNS: tuple[str, ...] = (
    ID_COL,
    SUB_CATEGORY_COL,
    MORPHOLOGY_COL,
)


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================


@dataclass(kw_only=True, frozen=True)
class PureRow:
    """One drawn row: label, exact Sub_category, and optional morphology filter."""

    label: str
    sub_category: str
    morphology: str | None = None


PURE_ROWS: tuple[PureRow, ...] = (
    PureRow(label="spores", sub_category="spores"),
    PureRow(label="germinated", sub_category="germinated"),
    PureRow(label="microcolonies", sub_category="microcolonies"),
    PureRow(label="small colonies", sub_category="small colonies"),
    PureRow(
        label="WT-like (WT)",
        sub_category="WT-like",
        morphology="WT",
    ),
)


# =============================================================================
# DATA PREPARATION
# =============================================================================


def load_merged(path: Path) -> pd.DataFrame:
    """Load the merged classification with sub-category and morphology columns."""
    df = pd.read_excel(path, sheet_name=MERGED_SHEET)
    require_columns(df, REQUIRED_COLUMNS, path)
    return df[[ID_COL, SUB_CATEGORY_COL, MORPHOLOGY_COL]].copy()


def select_pure_rows(merged: pd.DataFrame) -> pd.DataFrame:
    """Tag every gene in the five pure categories with its draw label."""
    frames = []
    for row in PURE_ROWS:
        mask = merged[SUB_CATEGORY_COL] == row.sub_category
        if row.morphology is not None:
            mask &= merged[MORPHOLOGY_COL] == row.morphology
        frames.append(merged.loc[mask, [ID_COL]].assign(**{ROW_COL: row.label}))

    selected = pd.concat(frames, ignore_index=True)
    counts = selected[ROW_COL].value_counts()
    for row in PURE_ROWS:
        count = int(counts.get(row.label, 0))
        logger.info(f"  {row.label:<20} n={count}")
        if count < 2:
            logger.warning(f"  {row.label!r} has fewer than two genes")
    return selected


def build_row_colors() -> dict[str, str]:
    """Assign the first Cell-palette colours to the rows in draw order."""
    if len(PURE_ROWS) > len(PALETTE):
        logger.warning(
            f"{len(PURE_ROWS)} rows exceed the {len(PALETTE)}-colour Cell palette; colours repeat."
        )
    return {row.label: PALETTE[index % len(PALETTE)] for index, row in enumerate(PURE_ROWS)}


def row_tick_labels(order: list[str], frame: pd.DataFrame) -> list[str]:
    """Tick labels with this panel's own sample size, e.g. ``spores (n = 143)``."""
    counts = frame[ROW_COL].value_counts()
    return [f"{row} (n = {int(counts.get(row, 0))})" for row in order]


# =============================================================================
# PLOTTING
# =============================================================================


def _apply_unicode_font(ax: Axes | None = None, font: str = FONT_FAMILY) -> None:
    """Apply the figure face as the fallback for characters the face lacks."""
    _cns_apply_unicode_font(ax, font=font)


def use_single_face() -> None:
    """Set the figure in ``FONT_FAMILY`` — one face for every element.

    cnsplots rewrites matplotlib's sans-serif list from its own settings each
    time it builds a figure, so the face is set there rather than in rcParams.
    At save time it also re-fonts any text holding a non-ASCII character to
    DejaVu Sans — the em dash in the axis labels and the U+2212 minus signs in
    the DR ticks are enough to trigger it — which would split a figure across
    two faces, so that fallback is pointed at the figure font as well.
    """
    cns.settings.font_sans_serif = (FONT_FAMILY, "Helvetica", "DejaVu Sans")
    # cnsplots re-fonts non-ASCII text to its default fallback at save time;
    # rebind that module-level hook so the figure face is the fallback instead.
    vars(cns._utils)["apply_unicode_font"] = _apply_unicode_font


def draw_figure(spec: FigureSpec, colors: dict[str, str]) -> None:
    """Draw the five-row violin figure plus its stats table and save it."""
    plot_height = spec.row_height_pt * len(spec.row_order)
    tick_labels = [row_tick_labels(spec.row_order, panel.frame) for panel in spec.panels]
    label_width = max(len(label) for labels in tick_labels for label in labels) * LABEL_CHAR_PT
    columns, stats_width = stats_columns(spec)
    stats_width += STATS_SLACK_PT
    margin_right = PVALUE_MARGIN_PT if spec.pairs else MARGIN_RIGHT_PT
    canvas_width = (
        2 * (PANEL_WIDTH_PT + PAD_LEFT_PT + label_width + margin_right)
        + stats_width
        + STATS_PAD_PT
        + CANVAS_SLACK_PT
    )

    mp = cns.multipanel(max_width=canvas_width)
    axes = []
    for index, panel in enumerate(spec.panels):
        ax = mp.panel(
            chr(ord("A") + index),
            width=PANEL_WIDTH_PT,
            height=plot_height,
            pad_left=PAD_LEFT_PT,
            pad_top=PAD_TOP_PT,
            margin_right=margin_right,
        )
        cns.violinplot(
            data=panel.frame,
            x=VALUE_COL,
            y=ROW_COL,
            order=spec.row_order,
            # hue is the row itself, never the merged category: cnsplots seaborn
            # normalises each violin against the densest peak of its own hue
            # group, so sharing a hue would shrink rows whose peak is lower.
            hue=ROW_COL,
            hue_order=spec.row_order,
            palette=[colors[row] for row in spec.row_order],
            legend=False,
            ax=ax,
        )
        ax.set_xlabel(panel.xlabel)
        ax.set_ylabel("")
        ax.set_xlim(*VALUE_RANGE)
        # Repeat the row labels on every violin panel; each panel shows the
        # sample size it actually plotted (DR and um coverage differ).
        ax.set_yticks(
            range(len(spec.row_order)),
            labels=tick_labels[index],
        )
        axes.append(ax)

    for index, panel in enumerate(spec.panels):
        if spec.pairs:
            annotate_pairs(
                axes[index],
                spec.row_order,
                panel.frame,
                spec.pairs,
                PANEL_WIDTH_PT,
                spec.row_height_pt,
                significant_color=PVALUE_SIGNIFICANT_COLOR,
                neutral_color=PVALUE_NEUTRAL_COLOR,
            )

    # Sample sizes and medians live in their own trailing panel, so the tick
    # labels stay one line. It carries a label like the others: an empty one
    # would reserve no label height and multipanel would place this axes higher.
    stats_ax = mp.panel(
        chr(ord("A") + len(spec.panels)),
        width=stats_width,
        height=plot_height,
        pad_left=STATS_PAD_PT,
        pad_top=PAD_TOP_PT,
        margin_right=MARGIN_RIGHT_PT,
    )
    draw_stats_panel(stats_ax, spec, axes[0], columns, stats_width)

    output_path = spec.output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf"):
        cns.savefig(output_path.with_suffix(suffix))
    logger.success(
        f"Saved: {output_path.with_suffix('.png')} and {output_path.with_suffix('.pdf')}"
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot DR/um distributions for five pure phenotype categories with cnsplots"
    )
    parser.add_argument(
        "--merged",
        type=Path,
        default=DEFAULT_MERGED,
        help=f"Merged categories xlsx (default: {DEFAULT_MERGED})",
    )
    parser.add_argument(
        "--dit-hap",
        type=Path,
        default=DEFAULT_DIT_HAP,
        help=f"DIT-HAP TSV (default: {DEFAULT_DIT_HAP})",
    )
    parser.add_argument(
        "--grna", type=Path, default=DEFAULT_GRNA, help=f"gRNA TSV (default: {DEFAULT_GRNA})"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable DEBUG level logging")
    return parser.parse_args()


def build_spec(
    selected: pd.DataFrame,
    dit_hap: pd.DataFrame,
    grna: pd.DataFrame,
    output_dir: Path,
) -> FigureSpec:
    """Build the five-row figure spec from the selected genes."""
    panels = [
        Panel(
            frame=build_long_frame(selected, dit_hap, "DR", row_col=ROW_COL),
            xlabel="DR — DIT-HAP",
            short_label="DR",
        ),
        Panel(
            frame=build_long_frame(selected, grna, "um", row_col=ROW_COL),
            xlabel="um — gRNA",
            short_label="um",
        ),
    ]
    row_order = [row.label for row in PURE_ROWS]
    present = set().union(*(set(panel.frame[ROW_COL]) for panel in panels))
    missing = [row for row in row_order if row not in present]
    if missing:
        raise ValueError(f"Rows missing from at least one proxy: {missing}")

    return FigureSpec(
        name="Pure categories",
        row_order=row_order,
        panels=panels,
        pairs=adjacent_pairs(row_order, panels),
        output_path=output_dir / "DR_um_distribution_pure.png",
        pvalue_note=PVALUE_NOTE,
        row_height_pt=ROW_HEIGHT_PT,
    )


def main() -> int:
    """Main orchestrator."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    use_single_face()

    try:
        logger.info("Loading merged categories…")
        merged = load_merged(args.merged)
        logger.info(f"  {len(merged)} genes")

        logger.info("Selecting pure phenotype rows…")
        selected = select_pure_rows(merged)
        logger.info(f"  {len(selected)} genes in {len(PURE_ROWS)} rows")

        logger.info("Loading DIT-HAP and gRNA…")
        dit_hap = load_dit_hap(args.dit_hap)
        grna = load_grna(args.grna)
        logger.info(f"  {len(dit_hap)} genes with DR, {len(grna)} genes with um")

        colors = build_row_colors()
        spec = build_spec(selected, dit_hap, grna, args.output_dir)
        logger.info(f"=== {spec.name} figure — {len(spec.row_order)} rows ===")
        log_category_stats(spec.row_order, spec.panels)
        log_pvalues(spec.pairs, spec.panels)
        draw_figure(spec, colors)
    except Exception as error:
        logger.exception(f"An unexpected error occurred: {error}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
