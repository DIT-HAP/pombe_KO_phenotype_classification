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
   ``Morphological classification == "WT"``

Every adjacent pair is annotated with a two-sided Mann-Whitney U p-value
(raw, no multiple-comparison correction — the same choice as step 06).
The bracket column, label format and row layout match the revised cnsplots
figure. Colours come from the Cell journal qualitative palette: significant
p-values (p < 0.05) in the palette terracotta and the rest in the palette
grey. The five row colours are the first five Cell colours in draw order,
which keeps them maximally distinct. Both violin panels repeat the row
labels, each with that panel's own sample size in parentheses on the same
line, e.g. ``spores (n = 143)``.

Input
-----
- ``data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx``
  (sheet ``All genes``) — ``Systematic ID``, ``Sub_category``, and
  ``Morphological classification``.
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
from matplotlib.text import Text
from matplotlib.textpath import TextPath
from scipy.stats import mannwhitneyu

# 4. Local Imports
from pipeline_utils import load_dit_hap, load_grna, setup_logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================

# The face every element of the figure is drawn in.
FONT_FAMILY = "Arial"

DEFAULT_MERGED = Path("data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx")
DEFAULT_DIT_HAP = Path("data/raw/gene_level_fitting_statistics_minusDR.tsv")
DEFAULT_GRNA = Path("data/raw/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv")
DEFAULT_OUTPUT_DIR = Path("results")

MERGED_SHEET = "All genes"
ID_COL = "Systematic ID"
SUB_CATEGORY_COL = "Sub_category"
MORPHOLOGY_COL = "Morphological classification"

# Internal long-form column names handed to cnsplots.
VALUE_COL = "Value"
ROW_COL = "Row label"

# Cell journal qualitative palette. The five rows take the first five colours
# in draw order, which keeps adjacent rows maximally distinct.
PALETTE = [mcolors.to_hex(tuple(color)) for color in cns.palettes("Cell")]
PVALUE_SIGNIFICANT_COLOR = PALETTE[0]
PVALUE_NEUTRAL_COLOR = PALETTE[9]

# Panel geometry in points (1/72 inch), matching the revised figure of step 06.
ROW_HEIGHT_PT = 16.89
PLOT_WIDTH_PT = 69.25
PAD_LEFT_PT = 12
PAD_TOP_PT = 12
MARGIN_RIGHT_PT = 16
FONT_SCALE = 1.0

# Panel C is a table: one probe / count / median column per proxy, one row per
# category. Columns are sized from their widest cell at STATS_FONTSIZE.
STATS_ALIGN = ("left", "right", "right")
STATS_COLUMN_GAP_PT = 5
STATS_GROUP_GAP_PT = 13
STATS_PAD_PT = 12
STATS_SLACK_PT = 10
STATS_FONTSIZE = 7
LABEL_CHAR_PT = 4.4

# Annotated panels reserve room for the bracket column and its labels.
PVALUE_MARGIN_PT = 60
CANVAS_SLACK_PT = 190

# Matches the sign-flipped range covered by both proxies.
VALUE_RANGE = (-2.2, 0.3)

# P-value column: one column of brackets right of the axes, never staircased.
PVALUE_GAP_PT = 6
PVALUE_TICK_PT = 3
PVALUE_TEXT_PT = 4
PVALUE_INSET_PT = 3
PVALUE_FONTSIZE = 7
PVALUE_ALPHA = 0.05
PVALUE_NOTE = "p: two-sided Mann-Whitney U (uncorrected)"


# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================


@dataclass(kw_only=True, slots=True, frozen=True)
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


@dataclass(kw_only=True, slots=True, frozen=True)
class Panel:
    """One proxy panel: its long-form frame, value-axis label, and short name."""

    frame: pd.DataFrame
    xlabel: str
    short_label: str


@dataclass(kw_only=True, slots=True, frozen=True)
class FigureSpec:
    """Everything the figure needs besides the loaded data and colours."""

    name: str
    row_order: list[str]
    panels: list[Panel]
    pairs: list[tuple[str, str]]
    output_path: Path
    pvalue_note: str
    row_height_pt: float
    plot_width_pt: float
    font_scale: float


# =============================================================================
# DATA PREPARATION
# =============================================================================


def load_merged(path: Path) -> pd.DataFrame:
    """Load the merged classification with sub-category and morphology columns."""
    df = pd.read_excel(path, sheet_name=MERGED_SHEET)
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


def build_long_frame(
    selected: pd.DataFrame,
    values: pd.DataFrame,
    value_col: str,
) -> pd.DataFrame:
    """Join one proxy's values onto the selected genes and reshape to long form."""
    joined = selected.merge(values, on=ID_COL, how="inner")
    frame = pd.DataFrame(
        {
            VALUE_COL: joined[value_col].to_numpy(),
            ROW_COL: joined[ROW_COL].to_numpy(),
        }
    )
    return frame.dropna(subset=[VALUE_COL])


def build_row_colors() -> dict[str, str]:
    """Assign the first Cell-palette colours to the rows in draw order."""
    if len(PURE_ROWS) > len(PALETTE):
        logger.warning(
            f"{len(PURE_ROWS)} rows exceed the {len(PALETTE)}-colour Cell "
            "palette; colours repeat."
        )
    return {
        row.label: PALETTE[index % len(PALETTE)]
        for index, row in enumerate(PURE_ROWS)
    }


def row_tick_labels(order: list[str], frame: pd.DataFrame) -> list[str]:
    """Tick labels with this panel's own sample size, e.g. ``spores (n = 143)``."""
    counts = frame[ROW_COL].value_counts()
    return [f"{row} (n = {int(counts.get(row, 0))})" for row in order]


def build_stats_rows(spec: FigureSpec) -> list[list[str]]:
    """Panel C's cells: one list of that row's column strings per category."""
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
    """Panel C's column edges in points — ``(left, right)`` each — and its width."""
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
    """Draw the per-row sample sizes and medians as its own axes-less panel."""
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
    ax.text(
        0, len(rows), spec.pvalue_note, ha="left", va="center",
        family=FONT_FAMILY, fontsize=fontsize, fontweight="normal",
    )
    # Copy the categorical limits (including their inversion) so the text rows
    # line up with the violins; this panel carries no axes of its own.
    ax.set_xlim(0, 1)
    ax.set_ylim(*row_axis.get_ylim())
    ax.axis("off")


def scale_text(ax: Axes, scale: float) -> None:
    """Scale every text artist cnsplots created on one axes, keeping it in proportion."""
    if scale == 1.0:
        return
    for axis in ("x", "y"):
        size = (
            ax.get_xticklabels()[0].get_fontsize()
            if axis == "x"
            else ax.get_yticklabels()[0].get_fontsize()
        )
        ax.tick_params(axis=axis, which="both", labelsize=size * scale)
    for label in (ax.xaxis.label, ax.yaxis.label):
        label.set_fontsize(label.get_fontsize() * scale)


def scale_panel_labels(fig: Figure, labels: list[str], scale: float) -> None:
    """Scale multipanel's own A/B/C lettering, which lives outside the axes."""
    if scale == 1.0:
        return
    for text in fig.findobj(match=Text):
        if text.get_text() in labels:
            text.set_fontsize(text.get_fontsize() * scale)


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
    cns._utils.apply_unicode_font = partial(cns.apply_unicode_font, font=FONT_FAMILY)


def draw_figure(spec: FigureSpec, colors: dict[str, str]) -> None:
    """Draw the five-row violin figure plus its stats table and save it."""
    plot_height = spec.row_height_pt * len(spec.row_order)
    tick_labels = [
        row_tick_labels(spec.row_order, panel.frame) for panel in spec.panels
    ]
    label_width = max(len(label) for labels in tick_labels for label in labels) * LABEL_CHAR_PT
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
        scale_text(ax, spec.font_scale)
        axes.append(ax)

    for index, panel in enumerate(spec.panels):
        if spec.pairs:
            annotate_pairs(
                axes[index], spec.row_order, panel.frame, spec.pairs,
                spec.plot_width_pt, spec.row_height_pt, spec.font_scale,
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

    scale_panel_labels(
        mp.fig,
        [chr(ord("A") + index) for index in range(len(spec.panels) + 1)],
        spec.font_scale,
    )

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
    parser = argparse.ArgumentParser(
        description="Plot DR/um distributions for five pure phenotype categories with cnsplots"
    )
    parser.add_argument("--merged", type=Path, default=DEFAULT_MERGED, help=f"Merged categories xlsx (default: {DEFAULT_MERGED})")
    parser.add_argument("--dit-hap", type=Path, default=DEFAULT_DIT_HAP, help=f"DIT-HAP TSV (default: {DEFAULT_DIT_HAP})")
    parser.add_argument("--grna", type=Path, default=DEFAULT_GRNA, help=f"gRNA TSV (default: {DEFAULT_GRNA})")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})")
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
            frame=build_long_frame(selected, dit_hap, "DR"),
            xlabel="DR — DIT-HAP",
            short_label="DR",
        ),
        Panel(
            frame=build_long_frame(selected, grna, "um"),
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
        plot_width_pt=PLOT_WIDTH_PT,
        font_scale=FONT_SCALE,
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
