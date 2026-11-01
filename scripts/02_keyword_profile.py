"""
Growth-tier Keyword Profile Analysis
=====================================

For each Growth_tier, extracts and categorises words from the phenotype
descriptions — growth signal words, morphological words, modifier words —
and reports their frequency.  This helps validate whether descriptions
assigned to the same tier share similar vocabulary, and whether tiers are
genuinely distinct.

Tier labels are derived from the input table's ``Category`` column, so the
column names follow the current DR-ranked categories instead of a hardcoded
1-5 mapping.

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
    - ``Growth signals`` — canonical growth-stage words
      (spores, germinated, microcolonies, small colonies, etc.)
    - ``Morphology`` — shape / structural words
      (long, curved, stubby, rounded, branched, misshapen, etc.)
    - ``Viability`` — cell death / survival words
      (lysis, dead, die, inviable)
    - ``Process`` — genetic / cellular process words
      (diploidising, suppressors, revertants, stationary, etc.)
    - ``Modifiers`` — frequency, degree, and quantity words
      (occasionally, often, slightly, some, many, once, …)
    - ``Other words`` — words whose exact token or stem is not in the
      ignore sets; inflected forms of category words can still appear

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
from collections import Counter
from pathlib import Path
import re
import sys

# 2. Third-party Imports
from loguru import logger

# 3. Data Processing Imports
import pandas as pd

# 4. Project Path Setup
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# 5. Local Imports
from growth_signals import MODIFIER_DEGREE, MODIFIER_FREQ, MODIFIER_QUANT  # noqa: E402
from phenotype import PHENOTYPE_DESC_COL  # noqa: E402
from pipeline_utils import require_columns, setup_logger  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_INPUT = Path("results/Hayles_2013_OB_merged_categories.xlsx")
DEFAULT_OUTPUT = Path("data/2_keyword_profile/growth_tier_keyword_profile.xlsx")

# ——— Word categories ——————————————————————————————————————————————————————————

# 1. Growth-stage vocabulary for this profile. It is broader than
# growth_signals.py: it also counts "germination" and colony variants that
# the classifier intentionally does not treat as growth signals.
GROWTH_SIGNAL_WORDS = {
    "spores",
    "spore",
    "germinated",
    "germinate",
    "germination",
    "microcolonies",
    "microcolony",
    "divide",
    "division",
    "divided",
    "divides",
    "divisions",
    "colonies",
    "colony",
    "small-colonies",
    "small-colony",
    "very-small-colonies",
    "very-small-colony",
}
GROWTH_STEMS = {
    "spore",
    "germinated",
    "germinate",
    "germination",
    "microcoloni",
    "divide",
    "division",
    "colon",
}

# 2. Morphology — shape / structural words only
MORPHOLOGY_WORDS = {
    "long",
    "short",
    "branched",
    "curved",
    "stubby",
    "rounded",
    "skittle",
    "skittles",
    "misshapen",
    "swollen",
    "septated",
    "septum",
    "multiseptated",
    "small",
    "wide",
    "wider",
    "longer",
    "shorter",
    "large",
    "larger",
    "thin",
    "narrow",
    "aberrantly",
    "abnormal",
    "tapered",
    "t-shaped",
    "dumbbell",
    "shaped",
    "vacuolated",
    "dark",
    "piled",
    "up",
    "misplaced",
    "misplace",
    "colour",
    "colored",
    "edged",
    "edges",
    "wavy",
}
MORPH_STEMS = {
    "long",
    "short",
    "branch",
    "curved",
    "curv",
    "stubby",
    "round",
    "skittle",
    "swollen",
    "sept",
    "misshapen",
    "small",
    "wide",
    "thin",
    "narrow",
    "aberrant",
    "abnormal",
    "larg",
    "taper",
    "dumbbell",
    "shaped",
    "centr",
    "vacuol",
    "dark",
    "misplac",
    "chain",
    "colour",
    "color",
    "edg",
    "wav",
}

# 3. Viability — cell death / survival
VIABILITY_WORDS = {
    "lyses",
    "lysis",
    "lysed",
    "dead",
    "die",
    "dying",
    "dies",
    "inviable",
}
VIABILITY_STEMS = {"lys", "dead", "die", "dyi", "inviable"}

# 4. Process — genetic / cellular processes
PROCESS_WORDS = {
    "diploids",
    "diploidising",
    "diploidises",
    "diploid",
    "suppressors",
    "supressors",
    "suppressor",
    "reverting",
    "revertants",
    "revert",
    "reverts",
    "stationary",
    "phase",
    "sporulating",
    "sporulation",
}
PROCESS_STEMS = {"diploid", "suppressor", "supressor", "revert", "stationar", "phase", "sporulat"}

# 5. Modifier sub-categories — imported from growth_signals, as sets for
# membership/union tests.
MOD_FREQ_WORDS = set(MODIFIER_FREQ)
MOD_DEGREE_WORDS = set(MODIFIER_DEGREE)
MOD_QUANT_WORDS = set(MODIFIER_QUANT)

MOD_FREQ_STEMS = {
    "occasionally",
    "occasion",
    "often",
    "sometimes",
    "mostly",
    "rare",
    "frequent",
    "possible",
    "may",
    "possibl",
    "rapid",
    "initial",
}

MOD_DEGREE_STEMS = {"slight", "very", "high", "barely", "weak"}

MOD_QUANT_STEMS = {
    "some",
    "many",
    "few",
    "lots",
    "several",
    "multiple",
    "more",
    "multi",
    "once",
    "twice",
}

# 6. Stop words — removed entirely
STOP_WORDS = {
    "at",
    "in",
    "of",
    "to",
    "and",
    "or",
    "the",
    "a",
    "an",
    "cells",
    "cell",
    "are",
    "is",
    "was",
    "were",
    "be",
    "with",
    "for",
    "after",
    "before",
    "more",
    "most",
    "25,32",
    "25,32,",
    "32,",
    "25,",
    "'",
    "25",
    "32",
    "36",
    "48h",
    "h",
    "wt",
    "viable",
    "essential",
    "yes",
    "but",
    "not",
    "no",
    "then",
    "give",
    "gives",
    "less",
    "almost",
    "wee",
    "one",
    "two",
    "three",
    "than",
    "as",
    "by",
    "from",
    "has",
    "had",
    "so",
    "on",
    "have",
}
STOP_STEMS = {
    "at",
    "in",
    "of",
    "to",
    "and",
    "or",
    "the",
    "a",
    "an",
    "cell",
    "are",
    "is",
    "was",
    "were",
    "be",
    "with",
    "for",
    "after",
    "before",
    "25",
    "32",
    "36",
    "48h",
    "h",
    "wt",
    "viable",
    "essential",
    "yes",
    "but",
    "not",
    "no",
    "then",
    "so",
    "on",
    "have",
    "has",
    "had",
    "than",
    "as",
    "by",
    "from",
}

# 5. Multi-word phrases that should be kept as single tokens.
#    Keys are the raw phrase; values are the hyphenated token.
PHRASE_PATTERNS: dict[str, str] = {
    "very small colonies": "very-small-colonies",
    "very small colony": "very-small-colony",
    "small colonies": "small-colonies",
    "small colony": "small-colony",
}

# Reverse mapping: hyphenated token → original display phrase
DISPLAY_NAMES: dict[str, str] = {v: k for k, v in PHRASE_PATTERNS.items()}

# Combined ignore set for the 'all other words' table
ALL_STEMS = (
    GROWTH_STEMS
    | MORPH_STEMS
    | VIABILITY_STEMS
    | PROCESS_STEMS
    | MOD_FREQ_STEMS
    | MOD_DEGREE_STEMS
    | MOD_QUANT_STEMS
    | STOP_STEMS
)
IGNORED: set[str] = set(ALL_STEMS) | set(PHRASE_PATTERNS.values())

# Exact token -> category, checked before the stem fallback.
EXACT_WORD_CATEGORIES: tuple[tuple[str, set[str]], ...] = (
    ("Growth signal", GROWTH_SIGNAL_WORDS),
    ("Morphology", MORPHOLOGY_WORDS),
    ("Viability", VIABILITY_WORDS),
    ("Process", PROCESS_WORDS),
    ("Modifier:Frequency", MOD_FREQ_WORDS),
    ("Modifier:Degree", MOD_DEGREE_WORDS),
    ("Modifier:Quantity", MOD_QUANT_WORDS),
    ("Stop word", STOP_WORDS),
)

# Stem prefix -> category, checked in this order after the exact matches.
PREFIX_WORD_CATEGORIES: tuple[tuple[str, set[str]], ...] = (
    ("Growth signal", GROWTH_STEMS),
    ("Morphology", MORPH_STEMS),
    ("Viability", VIABILITY_STEMS),
    ("Process", PROCESS_STEMS),
    ("Modifier:Frequency", MOD_FREQ_STEMS),
    ("Modifier:Degree", MOD_DEGREE_STEMS),
    ("Modifier:Quantity", MOD_QUANT_STEMS),
)

# ——— Tokenisation ————————————————————————————————————————————————————————————

# When reading from the final Excel, the description column name
DESC_COL = PHENOTYPE_DESC_COL
TIER_COL = "Growth_tier"
CATEGORY_COL = "Category"
REQUIRED_COLUMNS: tuple[str, ...] = (TIER_COL, CATEGORY_COL, DESC_COL)

# =============================================================================
# CORE LOGIC
# =============================================================================


def tokenise(text: str) -> list[str]:
    """Split a description into lowercase word tokens.

    Multi-word phrases in ``PHRASE_PATTERNS`` are replaced with hyphenated
    tokens before splitting, so that ``small colonies`` becomes the single
    token ``small-colonies``.
    """
    t = text.lower()
    # Replace multi-word phrases first (longest match wins by iteration order)
    for phrase, replacement in PHRASE_PATTERNS.items():
        t = t.replace(phrase, replacement)
    return [tok.strip(".,;:!?'") for tok in re.split(r"[,\s;:()]+", t) if tok.strip(".,;:!?'")]


def tier_labels_from_frame(df: pd.DataFrame) -> dict[int, str]:
    """Map each Growth_tier to its most common category in the input table."""
    labels: dict[int, str] = {}
    for tier in df[TIER_COL].dropna().unique():
        categories = df.loc[df[TIER_COL] == tier, CATEGORY_COL]
        mode = categories.mode()
        labels[int(tier)] = str(mode.iloc[0]) if not mode.empty else f"Tier {tier}"
    return labels


def tier_label(tier: int, labels: dict[int, str]) -> str:
    """Human-readable tier name, derived from the input table."""
    return labels.get(int(tier), f"Tier {tier}")


def _is_stop_word(word: str) -> bool:
    """Return True for a stop word or its simple plural/past inflection."""
    return word in STOP_STEMS or word.rstrip("s") in STOP_STEMS or word.rstrip("d") in STOP_STEMS


def _classify_word(word: str) -> str:
    """Return the category label for *word*, exact match before stem fallback."""
    for label, words in EXACT_WORD_CATEGORIES:
        if word in words:
            return label
    for label, stems in PREFIX_WORD_CATEGORIES:
        if any(word.startswith(stem) and len(stem) >= 3 for stem in stems):
            return label
    if _is_stop_word(word):
        return "Stop word"
    return "Other"


def build_term_table(
    term_counts: dict[int, Counter],
    selected_terms: set[str] | None = None,
    top_n: int | None = 50,
    tier_labels: dict[int, str] | None = None,
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
        .pivot_table(
            index="term",
            columns="tier",
            values="count",
            aggfunc="sum",
            fill_value=0,
        )
        .astype(int)
    )

    # Rename columns to human-readable
    pivot.columns = [f"Tier_{c}_{tier_label(c, tier_labels or {})}" for c in pivot.columns]

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


def count_words_per_tier(
    df: pd.DataFrame,
    tier_labels: dict[int, str],
) -> dict[int, Counter]:
    """Tokenise every tier's descriptions and count each word."""
    tier_word_counts: dict[int, Counter] = {}
    for tier in sorted(df[TIER_COL].dropna().unique()):
        tier_value = int(tier)
        texts = df.loc[df[TIER_COL] == tier, DESC_COL].dropna()
        counter: Counter = Counter()
        for text in texts:
            counter.update(tokenise(text))
        tier_word_counts[tier_value] = counter
        logger.debug(
            f"  Tier {tier_value} ({tier_label(tier_value, tier_labels)}): "
            f"{len(texts):,} genes, {len(counter):,} unique words",
        )
    return tier_word_counts


def build_vocab_tables(
    tier_word_counts: dict[int, Counter],
    tier_labels: dict[int, str],
) -> dict[str, pd.DataFrame]:
    """Build the vocabulary sheets, keyed by output sheet name."""
    tables = {
        "Total vocab": build_term_table(
            tier_word_counts,
            top_n=None,
            tier_labels=tier_labels,
        ),
        "Growth signals": build_term_table(
            tier_word_counts,
            selected_terms=GROWTH_SIGNAL_WORDS,
            top_n=50,
            tier_labels=tier_labels,
        ),
        "Morphology": build_term_table(
            tier_word_counts,
            selected_terms=MORPHOLOGY_WORDS,
            top_n=80,
            tier_labels=tier_labels,
        ),
        "Viability": build_term_table(
            tier_word_counts,
            selected_terms=VIABILITY_WORDS,
            top_n=30,
            tier_labels=tier_labels,
        ),
        "Process": build_term_table(
            tier_word_counts,
            selected_terms=PROCESS_WORDS,
            top_n=30,
            tier_labels=tier_labels,
        ),
        "Modifiers": build_term_table(
            tier_word_counts,
            selected_terms=MOD_FREQ_WORDS | MOD_DEGREE_WORDS | MOD_QUANT_WORDS,
            top_n=80,
            tier_labels=tier_labels,
        ),
    }
    for table in tables.values():
        if not table.empty:
            table.index.name = "word"
    return tables


def build_other_table(
    tier_word_counts: dict[int, Counter],
    tier_labels: dict[int, str],
) -> pd.DataFrame:
    """Build the sheet of words outside every curated vocabulary set."""
    other_terms: set[str] = set()
    for counter in tier_word_counts.values():
        other_terms.update(counter.keys())
    other_terms -= IGNORED
    table = build_term_table(
        tier_word_counts,
        selected_terms=other_terms,
        top_n=50,
        tier_labels=tier_labels,
    )
    if not table.empty:
        table.index.name = "word"
    return table


def write_workbook(
    output_path: Path,
    df: pd.DataFrame,
    tables: dict[str, pd.DataFrame],
    other_table: pd.DataFrame,
) -> None:
    """Write the summary, vocabulary sheets, and other-words sheet."""
    with pd.ExcelWriter(output_path) as writer:
        summary = df[TIER_COL].value_counts().sort_index().rename("gene_count").to_frame()
        summary.index.name = TIER_COL
        summary.to_excel(writer, sheet_name="Tier sizes")

        for sheet_name, table in tables.items():
            if not table.empty:
                table.to_excel(writer, sheet_name=sheet_name)
        if not other_table.empty:
            other_table.to_excel(writer, sheet_name="Other words")


def build_profile(input_path: Path, output_path: Path) -> None:
    """Load the merged table, profile each Growth_tier, and write the workbook."""
    logger.info(f"Loading classifications: {input_path}")
    df = pd.read_excel(input_path)
    require_columns(df, REQUIRED_COLUMNS, input_path)
    logger.info(f"Loaded {len(df):,} genes, {df[TIER_COL].nunique()} tiers")
    tier_labels = tier_labels_from_frame(df)

    logger.info("Building vocabulary tables ...")
    tier_word_counts = count_words_per_tier(df, tier_labels)
    tables = build_vocab_tables(tier_word_counts, tier_labels)
    other_table = build_other_table(tier_word_counts, tier_labels)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_workbook(output_path, df, tables, other_table)
    logger.success(f"Keyword profile saved: {output_path}")


def main() -> int:
    """Run the keyword profile and propagate fatal errors as exit code 1."""
    args = parse_args()
    setup_logger(log_level="DEBUG" if args.verbose else "INFO")
    try:
        build_profile(args.input.resolve(), args.output.resolve())
    except Exception as error:
        logger.exception(f"Keyword profile build failed: {error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
