"""
Shared DR / um plotting primitives
==================================

Helpers shared by ``06_plot_dr_distribution.py`` and
``07_plot_dr_distribution_pure.py``: long-form data preparation, the
stats-table columns, and the p-value brackets. The figure-specific pieces —
``draw_figure``, colour palettes, row selection and spec builders — stay in
each script.

This module needs matplotlib and scipy but not cnsplots, so it is importable
from the regular test environment.

Usage
-----
    from dr_plotting import FigureSpec, Panel, build_long_frame

Author:   Yusheng Yang (guidance) + Claude Code (implementation)
Date:     2026-10-29
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from dataclasses import dataclass, field
from itertools import pairwise
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger
from matplotlib.axes import Axes
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from scipy.stats import mannwhitneyu

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

# The face every element of both figures is drawn in.
FONT_FAMILY = "Arial"

DEFAULT_MERGED = Path("data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx")
DEFAULT_DIT_HAP = Path("data/raw/gene_level_fitting_statistics_minusDR.tsv")
DEFAULT_GRNA = Path("data/raw/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv")
DEFAULT_OUTPUT_DIR = Path("results")

MERGED_SHEET = "All genes"
ID_COL = "Systematic ID"
SUB_CATEGORY_COL = "Sub_category"

# Internal long-form column names handed to cnsplots.
VALUE_COL = "Value"
ROW_COL = "Row label"
GROUP_COL = "Merged category"

# Shared panel geometry in points (1/72 inch).
# 25 pt is 100 px in the 288 dpi export (1 pt = 4 px), the target width of one
# violin axes area. Height is set per figure from the row count.
PANEL_WIDTH_PT = 25
PAD_LEFT_PT = 12
PAD_TOP_PT = 12
MARGIN_RIGHT_PT = 32
STATS_PAD_PT = 12
STATS_SLACK_PT = 10
CANVAS_SLACK_PT = 190
LABEL_CHAR_PT = 4.4
PVALUE_MARGIN_PT = 60

# Matches the sign-flipped range covered by both proxies.
VALUE_RANGE = (-2.2, 0.3)

# Panel C's table: one set of probe / count / median columns per proxy.
STATS_ALIGN = ("left", "right", "right")
STATS_COLUMN_GAP_PT = 5
STATS_GROUP_GAP_PT = 13
# Custom text matches cnsplots' default tick/legend size (7 pt), so it needs
# no per-figure scaling even when the row count is large.
STATS_FONTSIZE = 7

# P-value column offsets, in points from the axes' right edge.
PVALUE_GAP_PT = 6
PVALUE_TICK_PT = 3
PVALUE_TEXT_PT = 4
PVALUE_INSET_PT = 3
PVALUE_FONTSIZE = 7
PVALUE_ALPHA = 0.05

# Footnote closing the stats panel in the annotated figures.
PVALUE_NOTE = "p: two-sided Mann-Whitney U (uncorrected)"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================


@dataclass(kw_only=True, slots=True, frozen=True)
class Panel:
    """One proxy panel: its long-form frame, value-axis label, and short name."""

    frame: pd.DataFrame
    xlabel: str
    short_label: str


@dataclass(kw_only=True, slots=True, frozen=True)
class FigureSpec:
    """Everything one output figure needs besides the loaded config."""

    name: str
    row_order: list[str]
    panels: list[Panel]
    pairs: list[tuple[str, str]]
    output_path: Path
    pvalue_note: str | None
    row_height_pt: float
    row_groups: dict[str, str] = field(default_factory=dict)
    boundaries: list[float] = field(default_factory=list)
    show_legend: bool = False


# =============================================================================
# DATA PREPARATION
# =============================================================================


def build_long_frame(
    merged: pd.DataFrame,
    values: pd.DataFrame,
    value_col: str,
    *,
    row_col: str,
    group_col: str | None = None,
) -> pd.DataFrame:
    """Join one proxy's values onto the merged table and reshape to long form."""
    joined = merged.merge(values, on=ID_COL, how="inner")
    data = {
        VALUE_COL: joined[value_col].to_numpy(),
        ROW_COL: joined[row_col].to_numpy(),
    }
    if group_col is not None:
        data[GROUP_COL] = joined[group_col].to_numpy()
    # The loaders already drop missing values; the guard keeps the per-row
    # sample sizes honest if a different reference table is passed in.
    return pd.DataFrame(data).dropna(subset=[VALUE_COL])


def build_stats_rows(spec: FigureSpec) -> list[list[str]]:
    """Panel C's cells: one list of that row's column strings per category.

    Each proxy contributes a probe label, a sample size and a median, in that
    order, so the columns repeat every ``len(STATS_ALIGN)`` entries.
    """
    tables = [
        (panel.short_label, panel.frame.groupby(ROW_COL)[VALUE_COL].agg(["size", "median"]))
        for panel in spec.panels
    ]
    rows = []
    for row in spec.row_order:
        cells = []
        for label, table in tables:
            if row in table.index:
                cells += [
                    f"{label}:",
                    f"n={int(table.loc[row, 'size'])}",
                    f"med={table.loc[row, 'median']:.3f}",
                ]
            else:
                cells += [f"{label}:", "n=-", "med=-"]
        rows.append(cells)
    return rows


def text_width_pt(text: str, fontsize: float, weight: str) -> float:
    """Width of one string in points, as the figure's own face draws it."""
    prop = FontProperties(family=FONT_FAMILY, weight=weight, size=fontsize)
    return TextPath((0, 0), text, prop=prop, size=fontsize).get_extents().width


def stats_columns(spec: FigureSpec) -> tuple[list[tuple[float, float]], float]:
    """Panel C's column edges in points — ``(left, right)`` each — and its width.

    Every column is as wide as its widest cell at the size this figure draws its
    text at. The panel is then as wide as either the table or the footnote under
    it, whichever is longer.
    """
    rows = build_stats_rows(spec)
    fontsize = STATS_FONTSIZE
    widths = [
        max(text_width_pt(row[column], fontsize, "bold") for row in rows)
        for column in range(len(rows[0]))
    ]
    edges: list[tuple[float, float]] = []
    cursor = 0.0
    for column, width in enumerate(widths):
        if column:
            gap = STATS_GROUP_GAP_PT if column % len(STATS_ALIGN) == 0 else STATS_COLUMN_GAP_PT
            cursor += gap
        edges.append((cursor, cursor + width))
        cursor += width
    if spec.pvalue_note:
        cursor = max(cursor, text_width_pt(spec.pvalue_note, fontsize, "normal"))
    return edges, cursor


def adjacent_pairs(order: list[str], panels: list[Panel]) -> list[tuple[str, str]]:
    """Adjacent category pairs holding at least two observations in every panel."""
    counts = [panel.frame[ROW_COL].value_counts() for panel in panels]
    pairs = []
    for first, second in pairwise(order):
        sizes_present = [
            int(table.get(row, 0)) for table in counts for row in (first, second)
        ]
        if min(sizes_present) >= 2:
            pairs.append((first, second))
    return pairs


# =============================================================================
# STATISTICS LOGGING
# =============================================================================


def pair_pvalue(frame: pd.DataFrame, first: str, second: str) -> float:
    """Two-sided Mann-Whitney U p-value between two rows of one panel."""
    left = frame.loc[frame[ROW_COL] == first, VALUE_COL]
    right = frame.loc[frame[ROW_COL] == second, VALUE_COL]
    return float(mannwhitneyu(left, right, alternative="two-sided").pvalue)


def format_pvalue(pvalue: float) -> str:
    """Format a p-value in scientific notation, e.g. ``p=1.2E-9``."""
    mantissa, exponent = f"{pvalue:.1e}".split("e")
    return f"p={mantissa.removesuffix('.0')}E{int(exponent)}"


def log_category_stats(order: list[str], panels: list[Panel]) -> None:
    """Log per-category sample size and median for every proxy."""
    for panel in panels:
        stats = panel.frame.groupby(ROW_COL)[VALUE_COL].agg(["size", "median"])
        logger.info(f"  {panel.xlabel} — n / median per category")
        for row in order:
            if row in stats.index:
                logger.info(
                    f"    {row:<60} n={int(stats.loc[row, 'size']):<5} "
                    f"med={stats.loc[row, 'median']:.3f}"
                )


def log_pvalues(pairs: list[tuple[str, str]], panels: list[Panel]) -> None:
    """Log the two-sided Mann-Whitney U p-value of every annotated pair and panel."""
    for panel in panels:
        logger.info(f"  {panel.xlabel} — two-sided Mann-Whitney U (no correction)")
        for first, second in pairs:
            logger.info(
                f"    {first:<32} vs {second:<32} "
                f"p={pair_pvalue(panel.frame, first, second):.3e}"
            )


# =============================================================================
# PLOTTING
# =============================================================================


def annotate_pairs(
    ax: Axes,
    order: list[str],
    frame: pd.DataFrame,
    pairs: list[tuple[str, str]],
    plot_width_pt: float,
    row_height_pt: float,
    *,
    significant_color: str,
    neutral_color: str,
) -> None:
    """Draw every adjacent-pair p-value on one bracket column right of the axes."""
    positions = {row: index for index, row in enumerate(order)}
    # The category axis puts one row every row_height_pt points, so the point
    # inset converts into the data units the bracket ends are measured in.
    inset = PVALUE_INSET_PT / row_height_pt
    transform = ax.get_yaxis_transform()
    inner_x = 1 + PVALUE_GAP_PT / plot_width_pt
    outer_x = 1 + (PVALUE_GAP_PT + PVALUE_TICK_PT) / plot_width_pt
    text_x = 1 + (PVALUE_GAP_PT + PVALUE_TICK_PT + PVALUE_TEXT_PT) / plot_width_pt
    for first, second in pairs:
        pvalue = pair_pvalue(frame, first, second)
        significant = pvalue < PVALUE_ALPHA
        color = significant_color if significant else neutral_color
        top, bottom = positions[first] + inset, positions[second] - inset
        for row in (top, bottom):
            ax.plot(
                [inner_x, outer_x], [row, row],
                transform=transform, color=color, linewidth=0.6, clip_on=False,
            )
        ax.plot(
            [outer_x, outer_x], [top, bottom],
            transform=transform, color=color, linewidth=0.6, clip_on=False,
        )
        ax.text(
            text_x, (top + bottom) / 2, format_pvalue(pvalue),
            transform=transform, ha="left", va="center",
            color=color, fontsize=PVALUE_FONTSIZE,
            fontweight="bold" if significant else "normal", clip_on=False,
        )


def draw_stats_panel(
    ax: Axes,
    spec: FigureSpec,
    row_axis: Axes,
    columns: list[tuple[float, float]],
    width_pt: float,
) -> None:
    """Draw the per-row sample sizes and medians as its own axes-less panel.

    Each cell is placed at its own column edge, in points converted to axes
    fractions, which is what keeps the table aligned in a proportional face.
    The footnote that closes the panel sits one row pitch under the last row —
    outside the copied row limits, which ``ax.text`` renders unclipped — and is
    set in the light weight of the same face to read as a note, not a row.
    """
    fontsize = STATS_FONTSIZE
    rows = build_stats_rows(spec)
    aligns = STATS_ALIGN * len(spec.panels)
    for index, cells in enumerate(rows):
        for (left, right), align, cell in zip(columns, aligns, cells):
            ax.text(
                (right if align == "right" else left) / width_pt, index, cell,
                ha=align, va="center", family=FONT_FAMILY,
                fontsize=fontsize, fontweight="bold",
            )
    if spec.pvalue_note:
        ax.text(
            0, len(rows), spec.pvalue_note, ha="left", va="center",
            family=FONT_FAMILY, fontsize=fontsize, fontweight="normal",
        )
    # Copy the categorical limits (including their inversion) so the text rows
    # line up with the violins; this panel carries no axes of its own.
    ax.set_xlim(0, 1)
    ax.set_ylim(*row_axis.get_ylim())
    ax.axis("off")
