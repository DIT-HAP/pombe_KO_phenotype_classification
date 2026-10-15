"""
DR / um Distribution by Phenotype Category (CNSPlots Style)
===========================================================

``cnsplots`` re-rendering of the violin plots from
``06_plot_dr_distribution.py``. The data selection, category order, grouping
and statistics are identical — only the rendering engine differs. Step ``06``
and its outputs are left untouched so both styles stay available.

Each figure holds three panels that share one category axis:

- **A** — DIT-HAP depletion ``DR`` (sign-flipped: more negative = depleted).
- **B** — gRNA depletion ``um`` (negated on load, same convention).
- **C** — sample size and median per row for both proxies, one table row each,
  in the same layout step ``06`` used. Keeping them off the tick labels leaves
  those single-line, which lets the rows sit closer together. The figure
  carrying p-value brackets closes the panel with one further line naming the
  test behind them.

Rows are phenotype categories; the violin shows the density and the embedded
white box the median/IQR. The two proxies do not cover the same genes (gRNA has
524 genes without a DIT-HAP value), so the per-row sample sizes legitimately
differ between panels A and B — panel C reports them separately rather than as
one shared ``(n = … / …)``.

Two figures are produced:

- **original** — fine-grained ``Sub_category`` rows, filled by their merged
  group. Rows of one merged category are kept together and the groups follow
  the merged-category draw order, with dashed dividers between them.
- **revised** — merged ``Category`` rows with two-sided Mann-Whitney U
  brackets between adjacent categories (raw p-values, no multiple-comparison
  correction — same choice as step ``06``). Every bracket sits on one vertical
  column right of its panel, labelled with the exact p-value: red and bold
  when significant (p < 0.05), grey otherwise. Because every adjacent pair is
  annotated, each bracket stops short of its two rows so neighbouring brackets
  stay separate instead of fusing into one continuous bar.
  ``cnsplots``'s own ``pairs=``
  annotation is not used because it staircases the brackets outwards and
  cannot colour them by significance.

Input
-----
- ``data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx``
  (sheet ``All genes``) — ``Systematic ID``, ``Sub_category``, ``Category``.
- ``data/references/gene_level_fitting_statistics_minusDR.tsv``
  — ``Systematic ID`` and sign-flipped ``DR``.
- ``data/references/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv``
  — ``Systematic ID`` and ``um`` (negated on load).
- ``data/4_categorized_genes/category_revisions.json``
  — draw orders plus the ``Sub_category -> Category`` merges.

Output
------
- ``results/DR_um_distribution_original_cnsplots.{png,pdf}`` — 2616 x 1197 px
- ``results/DR_um_distribution_revised_cnsplots.{png,pdf}`` — 2171 x 594 px

Both are sized to those pixels at the 288 dpi savefig default by per-figure
geometry constants near the top of the file (calibrated by rendering and
measuring; see the per-figure zoom block), so the two figures differ in row
height, panel width and — for the tall one — text scale. Width, unlike height,
is content-tight: the export crops to the drawn extent, so narrowing any panel
narrows the image by the same amount.

Every element of both figures — axis labels, tick labels and panel C's table
alike — is set in one face, Arial (``FONT_FAMILY``).

Usage
-----
    mamba run -n cnsplots python scripts/07_plot_dr_distribution_cns.py
    mamba run -n cnsplots python scripts/07_plot_dr_distribution_cns.py --verbose

Author:   Yusheng Yang (guidance) + Claude Code (implementation)
Date:     2026-10-13
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================

# 1. Standard Library Imports
import argparse
import sys
from dataclasses import dataclass
from functools import partial
from itertools import pairwise
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
import matplotlib
matplotlib.use("Agg")
import cnsplots as cns
import matplotlib.colors as mcolors
from loguru import logger
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Patch
from matplotlib.text import Text
from matplotlib.textpath import TextPath
from scipy.stats import mannwhitneyu

# 4. Local Imports
from category_revisions import DEFAULT_CONFIG, CategoryConfig, load_category_config
from pipeline_utils import load_dit_hap, load_grna, setup_logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================

# The face every element of both figures is drawn in. The fallbacks only come
# into play on a machine without Arial.
FONT_FAMILY = "Arial"

DEFAULT_MERGED = Path("data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx")
DEFAULT_DIT_HAP = Path("data/references/gene_level_fitting_statistics_minusDR.tsv")
DEFAULT_GRNA = Path("data/references/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv")
DEFAULT_OUTPUT_DIR = Path("results")

MERGED_SHEET = "All genes"
ID_COL = "Systematic ID"
SUB_CATEGORY_COL = "Sub_category"
CATEGORY_COL = "Category"

# Internal long-form column names handed to cnsplots.
VALUE_COL = "Value"
ROW_COL = "Row label"
GROUP_COL = "Merged category"

# Panel geometry in points (1/72 inch). Row height must fit one 8 pt tick label
# (≈ 11 pt tall) and one 7 pt stats line, so it bounds how low the figure can go.
# The two figures carry very different row counts — 50 vs 7 — and are limited by
# different dimensions, so each spec scales the base geometry for its own figure
# (see the spec builders): the tall one trades row height for height, the wide
# one trades panel width for width.
PLOT_WIDTH_PT = 185
ROW_HEIGHT_PT = 18
# --- per-figure zoom ---------------------------------------------------------
# Each figure independently scales the design geometry above to reach its target
# pixel size at the 288 dpi savefig default. The values are calibrated by
# rendering and measuring the saved PNG — multipanel's panel units do not map
# linearly to output pixels, so they cannot be computed from the target alone.

# 50 rows, height-limited (≈1200 px tall). Row height the layout and its text
# sizes are designed at: there, the one-line labels of adjacent rows just fit.
DESIGN_ROW_HEIGHT_PT = 14.2
ORIGINAL_ROW_HEIGHT_PT = 9.115
ORIGINAL_PLOT_WIDTH_PT = 142.4
# Text zooms with the rows. Leaving it full size while rows tighten is what makes
# adjacent labels collide, so this is derived, never set by hand.
ORIGINAL_FONT_SCALE = round(ORIGINAL_ROW_HEIGHT_PT / DESIGN_ROW_HEIGHT_PT, 3)

# 7 rows, width-limited (≈600 px tall). Rows stay near the design height, so the
# text keeps its natural size; the x axis keeps three ticks at this panel width
# (matplotlib drops the 0.5 ones).
REVISED_ROW_HEIGHT_PT = 16.89
REVISED_PLOT_WIDTH_PT = 69.25
PAD_LEFT_PT = 12
PAD_TOP_PT = 12
MARGIN_RIGHT_PT = 16
# The stats panel is a table: one set of probe / count / median columns per
# proxy, on a row per category. The figure's face is proportional, so a column
# cannot be aligned by padding spaces inside one string — each cell is placed at
# its own x instead, and the columns are sized from their widest cell (measured
# at STATS_FONTSIZE, i.e. before the per-figure font zoom; see stats_columns).
# Counts and medians are right-aligned, which lines their digits up.
STATS_ALIGN = ("left", "right", "right")
STATS_COLUMN_GAP_PT = 5
STATS_GROUP_GAP_PT = 13
STATS_PAD_PT = 12
STATS_SLACK_PT = 10
STATS_FONTSIZE = 7
# Annotated panels reserve room for the bracket column and its labels
# ("p=8.7E-128" is the longest label at 36 pt in Arial Bold at PVALUE_FONTSIZE,
# so the column needs PVALUE_GAP_PT + PVALUE_TICK_PT + PVALUE_TEXT_PT + 36 ≈
# 49 pt, plus a little air before the next panel's tick labels).
PVALUE_MARGIN_PT = 60
CANVAS_SLACK_PT = 190
LABEL_CHAR_PT = 4.4

# Matches the sign-flipped range covered by both proxies (step 06).
VALUE_RANGE = (-2.2, 0.3)

# P-value column: one column of brackets right of the axes (never staircased,
# which would overlap on a dense category axis). Offsets are points from the
# axes' right edge, converted to axes fractions at draw time; y stays in
# category (data) units. Ticks run outward from each row to the vertical bar,
# so the bracket closes on the outer side and points back at the pair it spans,
# matching statannotations' layout.
PVALUE_GAP_PT = 6
PVALUE_TICK_PT = 3
PVALUE_TEXT_PT = 4
# Every adjacent pair is annotated, so the brackets tile the axis and their bars
# would fuse into one continuous line. Each bracket ends PVALUE_INSET_PT short
# of its two row centres (converted with the row pitch at draw time), leaving a
# gap to the bracket above and below it.
PVALUE_INSET_PT = 3
PVALUE_FONTSIZE = 7
PVALUE_ALPHA = 0.05
PVALUE_SIGNIFICANT_COLOR = "red"
PVALUE_NEUTRAL_COLOR = "0.45"
# Footnote closing panel C in the annotated figure: what the brackets measure.
# It is one line under the last row, so its length also sets the panel width.
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
    row_groups: dict[str, str]
    panels: list[Panel]
    boundaries: list[float]
    pairs: list[tuple[str, str]]
    output_path: Path
    show_legend: bool
    # Line appended under panel C's rows, or ``None`` when the figure carries no
    # p-value brackets to explain.
    pvalue_note: str | None
    row_height_pt: float
    plot_width_pt: float
    font_scale: float


# =============================================================================
# DATA PREPARATION
# =============================================================================


def load_merged(path: Path) -> pd.DataFrame:
    """Load the merged classification with systematic ID, sub-category and category."""
    df = pd.read_excel(path, sheet_name=MERGED_SHEET)
    return df[[ID_COL, SUB_CATEGORY_COL, CATEGORY_COL]].copy()


def build_long_frame(
    merged: pd.DataFrame,
    values: pd.DataFrame,
    value_col: str,
    *,
    row_col: str,
    group_col: str,
) -> pd.DataFrame:
    """Join one proxy's values onto the merged table and reshape to long form."""
    joined = merged.merge(values, on=ID_COL, how="inner")
    frame = pd.DataFrame(
        {
            VALUE_COL: joined[value_col].to_numpy(),
            ROW_COL: joined[row_col].to_numpy(),
            GROUP_COL: joined[group_col].to_numpy(),
        }
    )
    # The loaders already drop missing values; the guard keeps the per-row
    # sample sizes honest if a different reference table is passed in.
    return frame.dropna(subset=[VALUE_COL])


def build_group_order(frames: list[pd.DataFrame], category_order: list[str]) -> list[str]:
    """Order merged categories as configured, appending any unlisted ones."""
    present = set().union(*(set(frame[GROUP_COL]) for frame in frames))
    ordered = [cat for cat in category_order if cat in present]
    return ordered + sorted(present - set(ordered))


def build_group_colors(group_order: list[str]) -> dict[str, str]:
    """Map each merged category to a colour from the cnsplots qualitative palette."""
    palette = [mcolors.to_hex(tuple(color)) for color in cns.palettes(cns.settings.palette_qual)]
    if len(group_order) > len(palette):
        logger.warning(
            f"{len(group_order)} groups exceed the {len(palette)}-colour cnsplots "
            f"palette '{cns.settings.palette_qual}'; colours repeat."
        )
    return {group: palette[i % len(palette)] for i, group in enumerate(group_order)}


def build_sub_category_rows(
    panels: list[Panel],
    sub_category_order: list[str],
    category_order: list[str],
    row_groups: dict[str, str],
) -> tuple[list[str], list[float]]:
    """Group sub-categories under their merged category and return row order plus dividers."""
    # Both proxies must be represented: gRNA covers sub-categories DIT-HAP does not.
    present = set().union(*(set(panel.frame[ROW_COL]) for panel in panels))
    listed = [sub for sub in sub_category_order if sub in present]
    listed += sorted(present - set(listed))

    # Sub-categories keep their configured relative order inside each merged
    # category; the groups themselves follow the merged-category draw order.
    grouped: dict[str, list[str]] = {}
    for sub in listed:
        grouped.setdefault(row_groups[sub], []).append(sub)

    groups = [group for group in category_order if group in grouped]
    groups += sorted(set(grouped) - set(groups))

    rows = [sub for group in groups for sub in grouped[group]]
    boundaries = []
    filled = 0
    for group in groups[:-1]:
        filled += len(grouped[group])
        boundaries.append(filled - 0.5)
    return rows, boundaries


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
    fontsize = STATS_FONTSIZE * spec.font_scale
    widths = [
        max(text_width_pt(row[column], fontsize, "bold") for row in rows)
        for column in range(len(rows[0]))
    ]
    edges: list[tuple[float, float]] = []
    cursor = 0.0
    for column, width in enumerate(widths):
        if column:
            gap = STATS_GROUP_GAP_PT if column % len(STATS_ALIGN) == 0 else STATS_COLUMN_GAP_PT
            cursor += gap * spec.font_scale
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
    """Format a p-value in scientific notation, e.g. ``p=1.2E-9`` (step 06 style)."""
    mantissa, exponent = f"{pvalue:.1e}".split("e")
    return f"p={mantissa.removesuffix('.0')}E{int(exponent)}"


def log_category_stats(order: list[str], panels: list[Panel]) -> None:
    """Log per-category sample size and median for every proxy."""
    for panel in panels:
        stats = panel.frame.groupby(ROW_COL)[VALUE_COL].agg(["size", "median"])
        logger.info(f"  {panel.xlabel} — n / median per category")
        for row in order:
            if row in stats.index:
                logger.info(f"    {row:<60} n={int(stats.loc[row, 'size']):<5} med={stats.loc[row, 'median']:.3f}")


def log_pvalues(pairs: list[tuple[str, str]], panels: list[Panel]) -> None:
    """Log the two-sided Mann-Whitney U p-value of every annotated pair and panel."""
    for panel in panels:
        logger.info(f"  {panel.xlabel} — two-sided Mann-Whitney U (no correction)")
        for first, second in pairs:
            logger.info(f"    {first:<32} vs {second:<32} p={pair_pvalue(panel.frame, first, second):.3e}")


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
    font_scale: float,
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
        color = PVALUE_SIGNIFICANT_COLOR if significant else PVALUE_NEUTRAL_COLOR
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
            color=color, fontsize=PVALUE_FONTSIZE * font_scale,
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
    fontsize = STATS_FONTSIZE * spec.font_scale
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


def draw_group_legend(ax: Axes, spec: FigureSpec, colors: dict[str, str]) -> None:
    """Draw the merged-category colour key for rows that are coloured by group."""
    seen = list(dict.fromkeys(spec.row_groups[row] for row in spec.row_order))
    handles = [
        Patch(facecolor=colors[group], edgecolor="black", linewidth=0.4, label=group)
        for group in seen
    ]
    ax.legend(handles=handles, frameon=False)
    cns.take_legend_out(title="Merged category", ax=ax)
    # take_legend_out anchors the legend above the axes top; multipanel reads
    # that overhang as a top decoration of this panel and shifts its axes down,
    # which would leave the two panels' rows unaligned. Level it with the axes
    # top instead — the legend still hangs outside to the right.
    ax.get_legend().set_bbox_to_anchor((1, 1.0))


def scale_text(ax: Axes, scale: float) -> None:
    """Scale every text artist cnsplots created on one axes, keeping it in proportion."""
    if scale == 1.0:
        return
    for axis in ("x", "y"):
        size = ax.get_xticklabels()[0].get_fontsize() if axis == "x" else ax.get_yticklabels()[0].get_fontsize()
        ax.tick_params(axis=axis, which="both", labelsize=size * scale)
    for label in (ax.xaxis.label, ax.yaxis.label):
        label.set_fontsize(label.get_fontsize() * scale)
    legend = ax.get_legend()
    if legend is not None:
        for text in (*legend.get_texts(), legend.get_title()):
            text.set_fontsize(text.get_fontsize() * scale)


def scale_panel_labels(fig: Figure, labels: list[str], scale: float) -> None:
    """Scale multipanel's own A/B/C lettering, which lives outside the axes."""
    if scale == 1.0:
        return
    for text in fig.findobj(match=Text):
        if text.get_text() in labels:
            text.set_fontsize(text.get_fontsize() * scale)


def use_single_face() -> None:
    """Set both figures in ``FONT_FAMILY`` — one face for every element.

    cnsplots rewrites matplotlib's sans-serif list from its own settings each
    time it builds a figure, so the face is set there rather than in rcParams.
    At save time it also re-fonts any text holding a non-ASCII character to
    DejaVu Sans — the em dash in the axis labels and the U+2212 minus signs in
    the DR ticks are enough to trigger it — which would split a figure across
    two faces, so that fallback is pointed at the figure font as well.
    """
    cns.settings.font_sans_serif = (FONT_FAMILY, "Helvetica", "DejaVu Sans")
    cns._utils.apply_unicode_font = partial(cns.apply_unicode_font, font=FONT_FAMILY)


def draw_figure(spec: FigureSpec, colors: dict[str, str], group_order: list[str]) -> None:
    """Draw one cnsplots figure — violin panels plus the stats table — and save it."""
    plot_height = spec.row_height_pt * len(spec.row_order)
    label_width = max(len(row) for row in spec.row_order) * LABEL_CHAR_PT
    # Panel C is sized by its own table: the columns are measured in the face the
    # rows are drawn in, at this figure's text size. Leftovers here show up as a
    # gap before the legend, which hangs off this panel's right edge.
    columns, stats_width = stats_columns(spec)
    stats_width += STATS_SLACK_PT
    margin_right = PVALUE_MARGIN_PT if spec.pairs else MARGIN_RIGHT_PT
    canvas_width = (
        2 * (spec.plot_width_pt + PAD_LEFT_PT + label_width + margin_right)
        + stats_width
        + STATS_PAD_PT
        + CANVAS_SLACK_PT
    )

    mp = cns.multipanel(max_width=canvas_width)
    axes = []
    for index, panel in enumerate(spec.panels):
        ax = mp.panel(
            chr(ord("A") + index),
            width=spec.plot_width_pt,
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
            # Group membership still drives the colours via the palette below.
            hue=ROW_COL,
            hue_order=spec.row_order,
            palette=[colors[spec.row_groups[row]] for row in spec.row_order],
            legend=False,
            ax=ax,
        )
        ax.set_xlabel(panel.xlabel)
        ax.set_ylabel("")
        ax.set_xlim(*VALUE_RANGE)
        scale_text(ax, spec.font_scale)
        if index:
            ax.tick_params(labelleft=False)
        axes.append(ax)

    # Shared row labels: seaborn places the categories at 0..n-1 on the group axis.
    axes[0].set_yticks(range(len(spec.row_order)), labels=spec.row_order)

    for boundary in spec.boundaries:
        for ax in axes:
            ax.axhline(boundary, color="0.45", linewidth=0.6, linestyle=(0, (4, 2)), zorder=0)

    for index, panel in enumerate(spec.panels):
        if spec.pairs:
            annotate_pairs(
                axes[index], spec.row_order, panel.frame, spec.pairs,
                spec.plot_width_pt, spec.row_height_pt, spec.font_scale,
            )

    # Sample sizes and medians live in their own trailing panel, so the tick
    # labels stay one line and the rows can be packed more tightly. It carries a
    # label like the others: an empty one would reserve no label height and
    # multipanel would place this axes higher than the violin panels.
    stats_ax = mp.panel(
        chr(ord("A") + len(spec.panels)),
        width=stats_width,
        height=plot_height,
        pad_left=STATS_PAD_PT,
        pad_top=PAD_TOP_PT,
        margin_right=MARGIN_RIGHT_PT,
    )
    draw_stats_panel(stats_ax, spec, axes[0], columns, stats_width)

    if spec.show_legend:
        draw_group_legend(stats_ax, spec, colors)
        # After the legend exists, so its own text is zoomed with the rest.
        scale_text(stats_ax, spec.font_scale)

    scale_panel_labels(mp.fig, [chr(ord("A") + index) for index in range(len(spec.panels) + 1)], spec.font_scale)

    output_path = spec.output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf"):
        cns.savefig(output_path.with_suffix(suffix))
    logger.success(f"Saved: {output_path.with_suffix('.png')} and {output_path.with_suffix('.pdf')}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Plot DR/um distributions by phenotype category with cnsplots")
    parser.add_argument("--merged", type=Path, default=DEFAULT_MERGED, help=f"Merged categories xlsx (default: {DEFAULT_MERGED})")
    parser.add_argument("--dit-hap", type=Path, default=DEFAULT_DIT_HAP, help=f"DIT-HAP TSV (default: {DEFAULT_DIT_HAP})")
    parser.add_argument("--grna", type=Path, default=DEFAULT_GRNA, help=f"gRNA TSV (default: {DEFAULT_GRNA})")
    parser.add_argument("--revisions", type=Path, default=DEFAULT_CONFIG, help=f"Category config JSON (default: {DEFAULT_CONFIG})")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("--verbose", action="store_true", help="Enable DEBUG level logging")
    return parser.parse_args()


def build_original_spec(
    merged: pd.DataFrame,
    dit_hap: pd.DataFrame,
    grna: pd.DataFrame,
    config: CategoryConfig,
    output_dir: Path,
) -> FigureSpec:
    """Build the fine-grained Sub_category figure spec grouped by revised category."""
    working = merged.rename(columns={SUB_CATEGORY_COL: ROW_COL}).assign(
        **{GROUP_COL: lambda df: df[ROW_COL].map(lambda sub: config.revisions.get(sub, sub))}
    )
    panels = [
        Panel(
            frame=build_long_frame(working, dit_hap, "DR", row_col=ROW_COL, group_col=GROUP_COL),
            xlabel="DR — DIT-HAP",
            short_label="DR",
        ),
        Panel(
            frame=build_long_frame(working, grna, "um", row_col=ROW_COL, group_col=GROUP_COL),
            xlabel="um — gRNA",
            short_label="um",
        ),
    ]
    row_groups = working.drop_duplicates(ROW_COL).set_index(ROW_COL)[GROUP_COL].to_dict()
    row_order, boundaries = build_sub_category_rows(
        panels, config.sub_category_order, config.category_order, row_groups
    )
    return FigureSpec(
        name="Sub_category",
        row_order=row_order,
        row_groups=row_groups,
        panels=panels,
        boundaries=boundaries,
        pairs=[],
        output_path=output_dir / "DR_um_distribution_original_cnsplots.png",
        show_legend=True,
        pvalue_note=None,
        row_height_pt=ORIGINAL_ROW_HEIGHT_PT,
        plot_width_pt=ORIGINAL_PLOT_WIDTH_PT,
        font_scale=ORIGINAL_FONT_SCALE,
    )


def build_revised_spec(
    merged: pd.DataFrame,
    dit_hap: pd.DataFrame,
    grna: pd.DataFrame,
    config: CategoryConfig,
    output_dir: Path,
) -> FigureSpec:
    """Build the merged Category figure spec with adjacent-pair p-value brackets."""
    panels = [
        Panel(
            frame=build_long_frame(merged, dit_hap, "DR", row_col=CATEGORY_COL, group_col=CATEGORY_COL),
            xlabel="DR — DIT-HAP",
            short_label="DR",
        ),
        Panel(
            frame=build_long_frame(merged, grna, "um", row_col=CATEGORY_COL, group_col=CATEGORY_COL),
            xlabel="um — gRNA",
            short_label="um",
        ),
    ]
    present = set().union(*(set(panel.frame[ROW_COL]) for panel in panels))
    row_order = [cat for cat in config.category_order if cat in present]
    row_order += sorted(present - set(row_order))
    return FigureSpec(
        name="Merged Category",
        row_order=row_order,
        row_groups={cat: cat for cat in row_order},
        panels=panels,
        boundaries=[],
        pairs=adjacent_pairs(row_order, panels),
        output_path=output_dir / "DR_um_distribution_revised_cnsplots.png",
        show_legend=False,
        pvalue_note=PVALUE_NOTE,
        row_height_pt=REVISED_ROW_HEIGHT_PT,
        plot_width_pt=REVISED_PLOT_WIDTH_PT,
        font_scale=1.0,
    )


def main() -> int:
    """Main orchestrator."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    use_single_face()

    try:
        logger.info("Loading merged categories…")
        merged = load_merged(args.merged)
        logger.info(f"  {len(merged)} genes, {merged[CATEGORY_COL].nunique()} merged categories")

        logger.info("Loading DIT-HAP and gRNA…")
        dit_hap = load_dit_hap(args.dit_hap)
        grna = load_grna(args.grna)
        logger.info(f"  {len(dit_hap)} genes with DR, {len(grna)} genes with um")

        config = load_category_config(args.revisions)
        logger.info(f"  {len(config.revisions)} sub-category merges configured")

        specs = [
            build_original_spec(merged, dit_hap, grna, config, args.output_dir),
            build_revised_spec(merged, dit_hap, grna, config, args.output_dir),
        ]

        all_frames = [panel.frame for spec in specs for panel in spec.panels]
        group_order = build_group_order(all_frames, config.category_order)
        colors = build_group_colors(group_order)
        logger.info(f"  {len(group_order)} merged categories: {', '.join(group_order)}")

        for spec in specs:
            logger.info(f"=== {spec.name} figure — {len(spec.row_order)} rows ===")
            log_category_stats(spec.row_order, spec.panels)
            if spec.pairs:
                log_pvalues(spec.pairs, spec.panels)
            draw_figure(spec, colors, group_order)
    except Exception as error:
        logger.exception(f"An unexpected error occurred: {error}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
