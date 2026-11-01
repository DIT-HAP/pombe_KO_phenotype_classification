"""
DR / um Distribution by Phenotype Category
==========================================

Reads the merged phenotype classification and joins it with two external
datasets — DIT-HAP (DR/DL) and gRNA (um/lam) — then plots the distribution
of depletion-rate proxies (DR for DIT-HAP, um for gRNA) across phenotype
categories using ``cnsplots`` violin + box plots with statistics annotations.

Both proxies are sign-flipped so that **more negative = more depleted**.
The DIT-HAP file already stores flipped ``DR`` values; gRNA ``um`` is
negated here on load for consistency. The x-axis range covers both datasets.

Each figure holds three panels that share one category axis:

- **A** — DIT-HAP depletion ``DR`` (sign-flipped: more negative = depleted).
- **B** — gRNA depletion ``um`` (negated on load, same convention).
- **C** — sample size and median per row for both proxies, one table row each.
  Keeping them off the tick labels leaves
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
  correction). Every bracket sits on one vertical
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
- ``data/raw/gene_level_fitting_statistics_minusDR.tsv``
  — ``Systematic ID`` and sign-flipped ``DR``.
- ``data/raw/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv``
  — ``Systematic ID`` and ``um`` (negated on load).
- ``data/4_categorized_genes/category_revisions.json``
  — draw orders plus the ``Sub_category -> Category`` merges.

Output
------
- ``results/DR_um_distribution_original.{png,pdf}`` — 2621 x 1825 px
- ``results/DR_um_distribution_revised.{png,pdf}`` — 2017 x 594 px

Each violin axes is 200 px wide at the 288 dpi savefig default (50 pt); height
grows with the row count. Text keeps ``cnsplots``' default sizes. Width, unlike
height, is content-tight: the export crops to the drawn extent, so the row
labels and the stats table set the final canvas width.

Every element of both figures — axis labels, tick labels and panel C's table
alike — is set in one face, Arial (``FONT_FAMILY``).

Usage
-----
    mamba run -n cnsplots python scripts/06_plot_dr_distribution.py
    mamba run -n cnsplots python scripts/06_plot_dr_distribution.py --verbose

Author:   Yusheng Yang (guidance) + Claude Code (implementation)
Date:     2026-10-13
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================

# 1. Standard Library Imports
import argparse
from pathlib import Path
import sys

# 2. Third-party Imports
import cnsplots as cns
from cnsplots._utils import apply_unicode_font as _cns_apply_unicode_font
from loguru import logger
import matplotlib as mpl
from matplotlib.axes import Axes
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

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
from dr_plotting import (  # noqa: E402
    CANVAS_SLACK_PT,
    DEFAULT_DIT_HAP,
    DEFAULT_GRNA,
    DEFAULT_MERGED,
    DEFAULT_OUTPUT_DIR,
    FONT_FAMILY,
    GROUP_COL,
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
    category_color_map,
    draw_stats_panel,
    log_category_stats,
    log_pvalues,
    stats_columns,
)
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

CATEGORY_COL = "Category"
REQUIRED_COLUMNS: tuple[str, ...] = (ID_COL, SUB_CATEGORY_COL, CATEGORY_COL)

# Panel geometry. The violin axes are 200 px wide in the 288 dpi export
# (PANEL_WIDTH_PT = 50 pt); height grows with the number of rows. Text keeps
# cnsplots' default sizes, so a row only has to clear one 7 pt tick label.
ORIGINAL_ROW_HEIGHT_PT = 14.2
REVISED_ROW_HEIGHT_PT = 16.89

# Cell journal qualitative palette. Each merged category takes the slot its
# canonical category order fixes, so a category keeps one colour in every figure.
PALETTE = [mcolors.to_hex(tuple(color)) for color in cns.palettes("Cell")]
PVALUE_SIGNIFICANT_COLOR = PALETTE[0]
PVALUE_NEUTRAL_COLOR = PALETTE[9]


# =============================================================================
# DATA PREPARATION
# =============================================================================


def load_merged(path: Path) -> pd.DataFrame:
    """Load the merged classification with systematic ID, sub-category and category."""
    df = pd.read_excel(path, sheet_name=MERGED_SHEET)
    require_columns(df, REQUIRED_COLUMNS, path)
    return df[[ID_COL, SUB_CATEGORY_COL, CATEGORY_COL]].copy()


def build_group_order(frames: list[pd.DataFrame], category_order: list[str]) -> list[str]:
    """Order merged categories as configured, appending any unlisted ones."""
    present = set().union(*(set(frame[GROUP_COL]) for frame in frames))
    ordered = [cat for cat in category_order if cat in present]
    return ordered + sorted(present - set(ordered))


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


# =============================================================================
# PLOTTING
# =============================================================================


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
    legend = ax.get_legend()
    if legend is not None:
        legend.set_bbox_to_anchor((1, 1.0))


def _apply_unicode_font(ax: Axes | None = None, font: str = FONT_FAMILY) -> None:
    """Apply the figure face as the fallback for characters the face lacks."""
    _cns_apply_unicode_font(ax, font=font)


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
    # cnsplots re-fonts non-ASCII text to its default fallback at save time;
    # rebind that module-level hook so the figure face is the fallback instead.
    vars(cns._utils)["apply_unicode_font"] = _apply_unicode_font


def draw_violin_panels(
    mp: cns.multipanel,
    spec: FigureSpec,
    colors: dict[str, str],
    plot_height: float,
    margin_right: float,
) -> list[Axes]:
    """Draw one violin axes per proxy and return them left to right."""
    axes: list[Axes] = []
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
        if index:
            ax.tick_params(labelleft=False)
        axes.append(ax)

    # Shared row labels: seaborn places the categories at 0..n-1 on the group axis.
    axes[0].set_yticks(range(len(spec.row_order)), labels=spec.row_order)
    return axes


def draw_figure(spec: FigureSpec, colors: dict[str, str]) -> None:
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
        2 * (PANEL_WIDTH_PT + PAD_LEFT_PT + label_width + margin_right)
        + stats_width
        + STATS_PAD_PT
        + CANVAS_SLACK_PT
    )

    mp = cns.multipanel(max_width=canvas_width)
    axes = draw_violin_panels(mp, spec, colors, plot_height, margin_right)

    for boundary in spec.boundaries:
        for ax in axes:
            ax.axhline(boundary, color="0.45", linewidth=0.6, linestyle=(0, (4, 2)), zorder=0)

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
        description="Plot DR/um distributions by phenotype category with cnsplots"
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
        "--revisions",
        type=Path,
        default=DEFAULT_CONFIG,
        help=f"Category config JSON (default: {DEFAULT_CONFIG})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
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
        **{GROUP_COL: lambda df: df[ROW_COL].map(config.apply)}
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
        output_path=output_dir / "DR_um_distribution_original.png",
        show_legend=True,
        pvalue_note=None,
        row_height_pt=ORIGINAL_ROW_HEIGHT_PT,
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
            frame=build_long_frame(
                merged, dit_hap, "DR", row_col=CATEGORY_COL, group_col=CATEGORY_COL
            ),
            xlabel="DR — DIT-HAP",
            short_label="DR",
        ),
        Panel(
            frame=build_long_frame(
                merged, grna, "um", row_col=CATEGORY_COL, group_col=CATEGORY_COL
            ),
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
        output_path=output_dir / "DR_um_distribution_revised.png",
        show_legend=False,
        pvalue_note=PVALUE_NOTE,
        row_height_pt=REVISED_ROW_HEIGHT_PT,
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
        colors = category_color_map(group_order, PALETTE, config.category_order)
        logger.info(f"  {len(group_order)} merged categories: {', '.join(group_order)}")

        for spec in specs:
            logger.info(f"=== {spec.name} figure — {len(spec.row_order)} rows ===")
            log_category_stats(spec.row_order, spec.panels)
            if spec.pairs:
                log_pvalues(spec.pairs, spec.panels)
            draw_figure(spec, colors)
    except Exception as error:
        logger.exception(f"An unexpected error occurred: {error}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
