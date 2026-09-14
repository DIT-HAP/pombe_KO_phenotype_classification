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
- ``data/2_keyword_profile/growth_tier_keyword_profile.xlsx``
  8 sheets:
    - ``Tier sizes`` — gene count per Growth_tier
    - ``Total vocab`` — all words with Category column, sorted by
      Category then frequency
    - ``Growth signals`` — canonical growth‑stage words
      (spores, germinated, microcolonies, small colonies, etc.)
    - ``Morphology`` — shape / structural words
      (long, curved, stubby, rounded, branched, misshapen, etc.)
    - ``Viability`` — cell death / survival words
      (lysis, dead, die, inviable)
    - ``Process`` — genetic / cellular process words
      (diploidising, suppressors, revertants, stationary, etc.)
    - ``Modifiers`` — frequency, degree, and quantity words
      (occasionally, often, slightly, some, many, once, …)
    - ``Other words`` — remaining words not in any category above

Usage
-----
    mamba run -n bioinformatics python scripts/02_keyword_profile.py
    mamba run -n bioinformatics python scripts/02_keyword_profile.py --verbose

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
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# 4. Local Imports
from growth_signals import MODIFIER_DEGREE, MODIFIER_FREQ, MODIFIER_QUANT
from pipeline_utils import setup_logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_INPUT = Path("results/Hayles_2013_OB_merged_categories.xlsx")
DEFAULT_OUTPUT = Path("data/2_keyword_profile/growth_tier_keyword_profile.xlsx")

# ——— Word categories ——————————————————————————————————————————————————————————

# 1. Canonical growth-signal words (from growth_signals.py)
GROWTH_SIGNAL_WORDS = {
    "spores", "spore", "germinated", "germinate", "germination",
    "microcolonies", "microcolony",
    "divide", "division", "divided", "divides", "divisions",
    "colonies", "colony",
    "small-colonies", "small-colony",
    "very-small-colonies", "very-small-colony",
}
GROWTH_STEMS = {"spore", "germinated", "germinate", "germination",
                "microcoloni", "divide", "division", "colon"}

# 2. Morphology — shape / structural words only
MORPHOLOGY_WORDS = {
    "long", "short", "branched", "curved", "stubby", "rounded", "skittle", "skittles",
    "misshapen", "swollen", "septated", "septum", "multiseptated",
    "small", "wide", "wider", "longer", "shorter", "large", "larger",
    "thin", "narrow", "aberrantly", "abnormal",
    "tapered", "t-shaped", "dumbbell", "shaped",
    "vacuolated", "dark", "piled", "up",
    "misplaced", "misplace",
    "colour", "colored", "edged", "edges", "wavy",
}
MORPH_STEMS = {"long", "short", "branch", "curved", "curv", "stubby", "round",
               "skittle", "swollen", "sept", "misshapen",
               "small", "wide", "thin", "narrow", "aberrant", "abnormal", "larg",
               "taper", "dumbbell", "shaped", "centr",
               "vacuol", "dark",
               "misplac", "chain",
               "colour", "color", "edg", "wav"}

# 3. Viability — cell death / survival
VIABILITY_WORDS = {
    "lyses", "lysis", "lysed",
    "dead", "die", "dying", "dies",
    "inviable",
}
VIABILITY_STEMS = {"lys", "dead", "die", "dyi", "inviable"}

# 4. Process — genetic / cellular processes
PROCESS_WORDS = {
    "diploids", "diploidising", "diploidises", "diploid",
    "suppressors", "supressors", "suppressor",
    "reverting", "revertants", "revert", "reverts",
    "stationary",
    "phase",
    "sporulating", "sporulation",
}
PROCESS_STEMS = {"diploid", "suppressor", "supressor",
                 "revert",
                 "stationar", "phase",
                 "sporulat"}

# 5. Modifier sub‑categories — imported from growth_signals, as sets for
# membership/union tests.
MOD_FREQ_WORDS = set(MODIFIER_FREQ)
MOD_DEGREE_WORDS = set(MODIFIER_DEGREE)
MOD_QUANT_WORDS = set(MODIFIER_QUANT)

MOD_FREQ_STEMS = {"occasionally", "occasion", "often", "sometimes",
                  "mostly", "rare", "frequent", "possible", "may",
                  "possibl", "rapid", "initial"}

MOD_DEGREE_STEMS = {"slight", "very", "high", "barely", "weak"}

MOD_QUANT_STEMS = {"some", "many", "few", "lots", "several", "multiple",
                   "more", "multi", "many",
                   "once", "twice"}

# 6. Stop words — removed entirely
STOP_WORDS = {
    "at", "in", "of", "to", "and", "or", "the", "a", "an",
    "cells", "cell", "are", "is", "was", "were", "be",
    "with", "for", "after", "before", "more", "most",
    "25,32", "25,32,", "32,", "25,", "'",
    "25", "32", "36", "48h", "h",
    "wt", "viable", "essential", "yes",
    "but", "not", "no", "then", "give", "gives",
    "less", "almost", "wee", "one", "two", "three",
    "than", "as", "by", "from", "has", "had",
    "so", "on", "have",
}
STOP_STEMS = {"at", "in", "of", "to", "and", "or", "the", "a", "an",
               "cell", "are", "is", "was", "were", "be",
               "with", "for", "after", "before",
               "25", "32", "36", "48h", "h",
               "wt", "viable", "essential", "yes",
               "but", "not", "no", "then",
               "so", "on", "have", "has", "had",
               "than", "as", "by", "from"}

# 5. Multi‑word phrases that should be kept as single tokens.
#    Keys are the raw phrase; values are the hyphenated token.
PHRASE_PATTERNS: dict[str, str] = {
    "very small colonies": "very-small-colonies",
    "very small colony":   "very-small-colony",
    "small colonies":      "small-colonies",
    "small colony":        "small-colony",
}

# Reverse mapping: hyphenated token → original display phrase
DISPLAY_NAMES: dict[str, str] = {v: k for k, v in PHRASE_PATTERNS.items()}

# Combined ignore set for the 'all other words' table
ALL_STEMS = (
    GROWTH_STEMS | MORPH_STEMS | VIABILITY_STEMS | PROCESS_STEMS
    | MOD_FREQ_STEMS | MOD_DEGREE_STEMS | MOD_QUANT_STEMS | STOP_STEMS
)
IGNORED: set[str] = set(ALL_STEMS) | set(PHRASE_PATTERNS.values())

# ——— Tokenisation ————————————————————————————————————————————————————————————

# When reading from the final Excel, the description column name
DESC_COL = "Deletion mutant phenotype description"
TIER_COL = "Growth_tier"

# =============================================================================
# LOGGING SETUP
# =============================================================================


setup_logger()

# =============================================================================
# CORE LOGIC
# =============================================================================


def tokenise(text: str) -> list[str]:
    """Split a description into lowercase word tokens.

    Multi‑word phrases in ``PHRASE_PATTERNS`` are replaced with hyphenated
    tokens before splitting, so that ``small colonies`` becomes the single
    token ``small-colonies``.
    """
    t = text.lower()
    # Replace multi‑word phrases first (longest match wins by iteration order)
    for phrase, replacement in PHRASE_PATTERNS.items():
        t = t.replace(phrase, replacement)
    return [tok.strip(".,;:!?'") for tok in re.split(r"[,\s;:()]+", t) if tok.strip(".,;:!?'")]


def tier_label(tier: int) -> str:
    """Human-readable tier name."""
    return {
        1: "Spores",
        2: "Germinated",
        3: "Microcolonies",
        4: "Small colonies",
        5: "WT-like",
    }.get(tier, f"Tier {tier}")


def _classify_word(word: str) -> str:
    """Return the category label for *word*.

    Exact membership in the original category sets is checked first,
    then stem‑prefix matching as a fallback.
    """
    # 1. Exact match in original sets (handles hyphenated / compound tokens)
    if word in GROWTH_SIGNAL_WORDS:
        return "Growth signal"
    if word in MORPHOLOGY_WORDS:
        return "Morphology"
    if word in VIABILITY_WORDS:
        return "Viability"
    if word in PROCESS_WORDS:
        return "Process"
    if word in MOD_FREQ_WORDS:
        return "Modifier:Frequency"
    if word in MOD_DEGREE_WORDS:
        return "Modifier:Degree"
    if word in MOD_QUANT_WORDS:
        return "Modifier:Quantity"
    if word in STOP_WORDS:
        return "Stop word"

    # 2. Stem‑prefix fallback
    for stem in GROWTH_STEMS:
        if word.startswith(stem) and len(stem) >= 3:
            return "Growth signal"
    for stem in MORPH_STEMS:
        if word.startswith(stem) and len(stem) >= 3:
            return "Morphology"
    for stem in VIABILITY_STEMS:
        if word.startswith(stem) and len(stem) >= 3:
            return "Viability"
    for stem in PROCESS_STEMS:
        if word.startswith(stem) and len(stem) >= 3:
            return "Process"
    for stem in MOD_FREQ_STEMS:
        if word.startswith(stem) and len(stem) >= 3:
            return "Modifier:Frequency"
    for stem in MOD_DEGREE_STEMS:
        if word.startswith(stem) and len(stem) >= 3:
            return "Modifier:Degree"
    for stem in MOD_QUANT_STEMS:
        if word.startswith(stem) and len(stem) >= 3:
            return "Modifier:Quantity"
    for stem in STOP_STEMS:
        if word == stem or word.rstrip("s") == stem or word.rstrip("d") == stem:
            return "Stop word"

    return "Other"


@logger.catch
def build_term_table(
    term_counts: dict[int, Counter],
    term_category: str | None = None,
    selected_terms: set[str] | None = None,
    top_n: int | None = 50,
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

    # Add category column (before the numeric columns) — do this BEFORE
    # renaming the index so the classifier sees the hyphenated tokens.
    pivot.insert(0, "Category", [_classify_word(w) for w in pivot.index])

    # Restore original display names for hyphenated compound tokens
    pivot = pivot.rename(index=DISPLAY_NAMES)

    # Add total column and sort
    pivot["Total"] = pivot.iloc[:, 1:].sum(axis=1)  # skip Category col
    # Sort by Category (alphabetical), then Total descending
    pivot = pivot.sort_values(["Category", "Total"], ascending=[True, False])
    if top_n is not None:
        pivot = pivot.head(top_n)

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

    # 3a — All words (full list with category labels)
    logger.info("Building vocabulary tables …")
    total_vocab = build_term_table(tier_word_counts, "all", top_n=None)
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

    # 3d — Viability words
    viability_table = build_term_table(
        tier_word_counts, "viability", selected_terms=VIABILITY_WORDS, top_n=30,
    )
    if not viability_table.empty:
        viability_table.index.name = "word"

    # 3e — Process words
    process_table = build_term_table(
        tier_word_counts, "process", selected_terms=PROCESS_WORDS, top_n=30,
    )
    if not process_table.empty:
        process_table.index.name = "word"

    # 3f — Modifier words (all sub‑categories combined)
    mod_table = build_term_table(
        tier_word_counts,
        "modifiers",
        selected_terms=MOD_FREQ_WORDS | MOD_DEGREE_WORDS | MOD_QUANT_WORDS,
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
        if not viability_table.empty:
            viability_table.to_excel(writer, sheet_name="Viability")
        if not process_table.empty:
            process_table.to_excel(writer, sheet_name="Process")
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
