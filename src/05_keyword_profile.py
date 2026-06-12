"""
Growth-tier Keyword Profile Analysis
=====================================

For each Growth_tier (1–5), extracts and categorises words from the
phenotype descriptions — growth signal words, morphological words,
modifier words — and reports their frequency.  This helps validate
whether descriptions assigned to the same tier share similar
vocabulary, and whether tiers are genuinely distinct.

The output is an inspection table, not a word cloud.  Each row is a
word; columns show its frequency per Growth_tier, allowing quick
assessment of which words drive the classification.

Input
-----
- ``results/Hayles_2013_OB_merged_categories.xlsx`` (final merged output)

Output
------
- ``data/5_keyword_profile/growth_tier_keyword_profile.xlsx``
  6 sheets:
    - ``Tier sizes`` — gene count per Growth_tier
    - ``Total vocab`` — all words, sorted by total frequency
    - ``Growth signals`` — only canonical signal words
      (spores, germinated, microcolonies, divide, etc.)
    - ``Morphology`` — shape/state words
      (long, misshapen, stubby, rounded, skittle, etc.)
    - ``Modifiers`` — frequency/degree words
      (occasionally, often, some, slightly, very, etc.)
    - ``Other words`` — remaining words not in any category above

Usage
-----
    mamba run -n bioinformatics python src/05_keyword_profile.py
    mamba run -n bioinformatics python src/05_keyword_profile.py --verbose

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-06-11
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
import re
from collections import Counter
from pathlib import Path

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_INPUT = Path("results/Hayles_2013_OB_merged_categories.xlsx")
DEFAULT_OUTPUT = Path("data/5_keyword_profile/growth_tier_keyword_profile.xlsx")

# ——— Word categories ——————————————————————————————————————————————————————————

# 1. Canonical growth-signal words (from growth_signals.py)
GROWTH_SIGNAL_WORDS = {
    "spores", "spore", "germinated", "germinate", "germination",
    "microcolonies", "microcolony",
    "divide", "division", "divided", "divides",
    "colonies", "colony",
}

# 2. Modifier / frequency words (signal that the main phenotype is partial)
MODIFIER_WORDS = {
    "occasionally", "often", "occasional", "may", "some",
    "sometimes", "mostly", "rarely", "frequently", "possible",
    "occasionally,", "often,", "some,", "possibly",
    "slightly", "very", "highly", "barely", "barely",
    "lots", "many", "more", "once", "twice", "several",
    "few", "multiple", "a",
}

# 3. Morphological words (from Hayles 2013 and the FYPO table)
MORPHOLOGY_WORDS = {
    "long", "short", "branched", "curved", "stubby", "rounded", "skittle",
    "misshapen", "swollen", "septated", "septum", "multiseptated",
    "small", "wide", "wider", "longer", "shorter",
    "thin", "narrow", "aberrantly", "abnormal",
    "tapered", "t-shaped", "dumbbell", "shaped",
    "lyses", "lysis", "dead", "die", "dying", "inviable",
    "vacuolated", "dark", "piled", "up",
    "diploids", "diploidising", "diploidises", "suppressors",
    "reverting", "revertants",
    "misplaced", "misplace", "misplaced",
    "colour", "edged", "edges", "wavy",
    "septated", "multiseptated",
}

# 4. Common stop words — removed entirely
STOP_WORDS = {
    "at", "in", "of", "to", "and", "or", "the", "a", "an",
    "cells", "cell", "are", "is", "was", "were", "be",
    "with", "for", "after", "before", "more", "most",
    "25,32", "25,32,", "32,", "25,", "'",
    "25", "32", "36", "48h", "h",
    "wt", "viable", "essential",
    "but", "not", "no", "then", "give", "gives",
    "less", "almost", "wee", "one", "two", "three",
    "than", "as", "by", "from", "has", "had",
}

# Combined ignore set for the 'all other words' table
IGNORED = GROWTH_SIGNAL_WORDS | MODIFIER_WORDS | MORPHOLOGY_WORDS | STOP_WORDS

# ——— Tokenisation ————————————————————————————————————————————————————————————

# Split on whitespace and common punctuation
_TOKEN_RE = re.compile(r"[,\s;:()]+")

# When reading from the final Excel, the description column name
DESC_COL = "Deletion mutant phenotype description"
TIER_COL = "Growth_tier"

# =============================================================================
# LOGGING SETUP
# =============================================================================


def setup_logger(log_level: str = "INFO") -> None:
    """Configure the Loguru logger."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level:<8} | {message}",
        level=log_level,
    )


setup_logger()

# =============================================================================
# CORE LOGIC
# =============================================================================


def tokenise(text: str) -> list[str]:
    """Split a description into lowercase word tokens, stripping punctuation."""
    return [t.strip(".,;:!?'") for t in _TOKEN_RE.split(text.lower()) if t.strip(".,;:!?'")]


def tier_label(tier: int) -> str:
    """Human-readable tier name."""
    return {
        1: "Spores",
        2: "Germinated",
        3: "Microcolonies",
        4: "Small colonies",
        5: "WT",
    }.get(tier, f"Tier {tier}")


@logger.catch
def build_term_table(
    term_counts: dict[int, Counter],
    term_category: str,
    selected_terms: set[str] | None = None,
    top_n: int = 50,
) -> pd.DataFrame:
    """Build a DataFrame: rows = terms, columns = tier frequencies."""
    rows: list[dict] = []
    for tier in sorted(term_counts):
        for word, count in term_counts[tier].most_common():
            if selected_terms and word not in selected_terms:
                continue
            rows.append({"term": word, "tier": tier, "count": count})

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)

    # Pivot: rows = term, columns = tier
    pivot = (
        df.groupby(["term", "tier"])["count"]
        .sum()
        .reset_index()
        .pivot(index="term", columns="tier", values="count")
        .fillna(0)
        .astype(int)
    )

    # Rename columns to human-readable
    pivot.columns = [f"Tier_{c}_{tier_label(c)}" for c in pivot.columns]

    # Add total column and sort
    pivot["Total"] = pivot.sum(axis=1)
    pivot = pivot.sort_values("Total", ascending=False).head(top_n)

    return pivot


# =============================================================================
# MAIN EXECUTION
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build keyword profile per Growth_tier.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Merged classification xlsx (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output xlsx path (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument("--verbose", action="store_true", help="DEBUG logging")
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
    logger.info(f"Loading classifications: {input_path}")
    df = pd.read_excel(input_path)
    logger.info(f"Loaded {len(df):,} genes, {df[TIER_COL].nunique()} tiers")

    # ------------------------------------------------------------------
    # 2. Tokenise per tier
    # ------------------------------------------------------------------
    tier_word_counts: dict[int, Counter] = {}
    for tier in sorted(df[TIER_COL].unique()):
        texts = df[df[TIER_COL] == tier][DESC_COL].dropna()
        counter: Counter = Counter()
        for text in texts:
            counter.update(tokenise(text))
        tier_word_counts[tier] = counter
        logger.debug(
            f"  Tier {tier} ({tier_label(tier)}): "
            f"{len(texts):,} genes, {len(counter):,} unique words",
        )

    # ------------------------------------------------------------------
    # 3. Build output tables
    # ------------------------------------------------------------------

    # 3a — All words (total vocab, top 100)
    logger.info("Building vocabulary tables …")
    total_vocab = build_term_table(tier_word_counts, "all", top_n=100)
    total_vocab.index.name = "word"

    # 3b — Growth signals only
    growth_table = build_term_table(
        tier_word_counts,
        "growth",
        selected_terms=GROWTH_SIGNAL_WORDS,
        top_n=50,
    )
    if not growth_table.empty:
        growth_table.index.name = "word"

    # 3c — Morphology words only
    morph_table = build_term_table(
        tier_word_counts,
        "morphology",
        selected_terms=MORPHOLOGY_WORDS,
        top_n=80,
    )
    if not morph_table.empty:
        morph_table.index.name = "word"

    # 3d — Modifier words only
    mod_table = build_term_table(
        tier_word_counts,
        "modifiers",
        selected_terms=MODIFIER_WORDS,
        top_n=80,
    )
    if not mod_table.empty:
        mod_table.index.name = "word"

    # ------------------------------------------------------------------
    # 4. Save
    # ------------------------------------------------------------------
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_path) as writer:
        # Summary — tier sizes
        summary = (
            df[TIER_COL]
            .value_counts()
            .sort_index()
            .rename("gene_count")
            .to_frame()
        )
        summary.index.name = "Growth_tier"
        summary.to_excel(writer, sheet_name="Tier sizes")

        # Vocab tables
        total_vocab.to_excel(writer, sheet_name="Total vocab")
        if not growth_table.empty:
            growth_table.to_excel(writer, sheet_name="Growth signals")
        if not morph_table.empty:
            morph_table.to_excel(writer, sheet_name="Morphology")
        if not mod_table.empty:
            mod_table.to_excel(writer, sheet_name="Modifiers")

        # Bonus — all OTHER words not in any category (potential discovery)
        other_terms: set[str] = set()
        for counter in tier_word_counts.values():
            other_terms.update(counter.keys())
        other_terms -= IGNORED
        other_table = build_term_table(
            tier_word_counts,
            "other",
            selected_terms=other_terms,
            top_n=50,
        )
        if not other_table.empty:
            other_table.index.name = "word"
            other_table.to_excel(writer, sheet_name="Other words")

    logger.success(f"Keyword profile saved: {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
