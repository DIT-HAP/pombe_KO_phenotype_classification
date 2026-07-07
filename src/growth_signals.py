"""
Growth Signal Definitions for S. pombe Phenotype Classification
===============================================================

Shared module defining growth-related keyword signals, their severity
tiers, and the ``classify_growth()`` engine that replaces the brittle
if-elif chains in the original categorize scripts.

Signals are detected independently from a phenotype description, then
resolved into a fine-grained ``Category`` label and a coarse
``Growth_tier`` (1–5) for downstream statistical analysis.

Input
-----
- A phenotype description string (e.g. from ``Deletion mutant phenotype
  description`` column in the Hayles table).

Output
------
- ``tuple[str, int]`` — (category_name, growth_tier).

Usage
-----
    from src.growth_signals import classify_growth

    cat, tier = classify_growth("VIABLE WT cells at 25,32")
    # cat == "WT-like", tier == 5

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-06-09
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import re
from dataclasses import dataclass

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================


@dataclass(frozen=True, slots=True)
class GrowthSignal:
    """One detectable keyword pattern and its corresponding growth category."""

    keyword: str  # substring to search for (case-insensitive)
    category: str  # fine-grained category label
    tier: int  # growth severity (1 = worst, 5 = best)


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

# Growth tier labels
TIER_LABELS: dict[int, str] = {
    1: "Spores (no germination)",
    2: "Germinated (limited or no division)",
    3: "Microcolonies (severely limited growth)",
    4: "Small colonies (visible but reduced)",
    5: "WT-like (normal growth)",
}

# Canonical signal table, ordered from most to least specific.
# When multiple keywords match, category names are concatenated (alphabetically
# sorted, deduplicated) and the worst tier is assigned.
GROWTH_SIGNALS: list[GrowthSignal] = [
    # tier 4 — most specific phrases first to avoid substring collision
    GrowthSignal(keyword="very small colon", category="very small colonies", tier=4),
    GrowthSignal(keyword="small colon", category="small colonies", tier=4),
    # tier 3
    GrowthSignal(keyword="microcolonies", category="microcolonies", tier=3),
    # tier 2
    GrowthSignal(keyword="germinated", category="germinated", tier=2),
    # Note: "germination" (noun) is NOT a growth signal — in the data it
    # only appears in morphological contexts ("germination long", "poor
    # germination", "from germination"), never as a standalone growth state.
    GrowthSignal(keyword="divide", category="germinated and divided", tier=2),
    GrowthSignal(keyword="division", category="germinated and divided", tier=2),
    # tier 1
    GrowthSignal(keyword="spores", category="spores", tier=1),
]

# Modifier words — when a comma‑separated segment starts with one of these,
# it is a secondary description (modifier / morphology supplement) and is
# excluded from growth‑signal detection.
MODIFIER_WORDS: tuple[str, ...] = (
    # Frequency
    "occasionally", "often", "occasional", "sometimes",
    "mostly", "rarely", "frequently", "frequency", "rare",
    "possible", "may", "possibly", "rapidly", "initially",
    # Degree
    "slightly", "very", "highly", "barely", "slight", "high", "weak",
    # Quantity
    "some", "many", "few", "lots", "several", "multiple",
    "once", "twice", "more", "multi",
)

# Growth‑signal keywords used for segment‑level analysis.
GROWTH_KEYWORDS: tuple[str, ...] = (
    "spores", "germinated", "microcolonies",
    "small colon", "very small colon",
    "divide", "division",
)

# =============================================================================
# CORE LOGIC
# =============================================================================


def _is_modifier_segment(seg: str) -> bool:
    """Return True if *seg* starts with a modifier word."""
    sl = seg.lower().strip()
    return any(sl.startswith(m) for m in MODIFIER_WORDS)


def _has_growth_signal(seg: str) -> bool:
    """Return True if *seg* contains a growth‑signal keyword."""
    sl = seg.lower()
    return any(kw in sl for kw in GROWTH_KEYWORDS)


def filter_primary_segments(text: str) -> str:
    """Remove modifier‑led comma segments, keep the rest.

    A comma‑separated segment that starts with a modifier word (e.g.
    ``some``, ``occasionally``, ``slightly``) is treated as a secondary
    description and excluded.  All other segments — whether they contain
    growth signals or pure morphology — are kept for classification.

    If every segment is modifier‑led, the first segment is returned as a
    fallback so the description is never empty.
    """
    if not isinstance(text, str) or "," not in text:
        return text
    segments = [seg.strip() for seg in text.split(",")]
    primary = [seg for seg in segments if not _is_modifier_segment(seg)]
    return ", ".join(primary) if primary else segments[0]


def count_growth_segments(text: str) -> int:
    """Count comma segments that contain a growth signal (excluding modifier‑led ones)."""
    if not isinstance(text, str) or "," not in text:
        return 1 if (isinstance(text, str) and _has_growth_signal(text)) else 0
    segments = [seg.strip() for seg in text.split(",")]
    return sum(1 for seg in segments
               if not _is_modifier_segment(seg) and _has_growth_signal(seg))


def _has_standalone_spores(text: str) -> bool:
    """Check if ``spores`` appears anywhere NOT as part of ``germinated spores``.

    A standalone occurrence (comma-separated, at the start of a sentence,
    etc.) indicates a mixed population where both dormant and germinated
    spores coexist.
    """
    for m in re.finditer(r"\bspores?\b", text):
        start = m.start()
        # Look at the 12 characters before 'spores'
        before = text[max(0, start - 12) : start].strip().rstrip(",")
        if before != "germinated":
            return True
    return False


def _detect_segment_signals(seg: str) -> list[GrowthSignal]:
    """Return all growth signals found in *seg*."""
    sl = seg.lower()
    return [s for s in GROWTH_SIGNALS if s.keyword in sl]


def _get_modifier_prefix(seg: str) -> str | None:
    """If *seg* starts with a modifier word, return that word; else None."""
    sl = seg.lower().strip()
    for m in MODIFIER_WORDS:
        if sl.startswith(m):
            return m
    return None


def _signal_to_name(signal: GrowthSignal, modifier: str | None) -> str:
    """Convert a signal + modifier into a category name.

    - No modifier → plain signal category (e.g. ``germinated``).
    - With modifier → prefixed (e.g. ``some germinated``, ``often divided``).
    """
    if modifier is None:
        return signal.category
    # "germinated and divided" → "often divided" (germinated already implied)
    if signal.category == "germinated and divided":
        return f"{modifier} divided"
    return f"{modifier} {signal.category}"


def classify_growth(description: str) -> tuple[str, int]:
    """Return (category, growth_tier) for a phenotype description string.

    The description is split into comma segments.  Each segment is scanned
    for growth signals:

    - **Primary segments** (not starting with a modifier) contribute plain
      signal names (e.g. ``germinated``, ``spores``).
    - **Modifier‑led segments** (starting with ``some``, ``often``, etc.)
      contribute prefixed signal names (e.g. ``some germinated``,
      ``often divided``) so that low‑frequency events are distinguished
      from definitive ones.

    If a signal appears in both a primary and a modifier-led segment, the
    plain (primary) version is kept — the definitive occurrence dominates.

    Tier is the best (highest) tier among all matched signals.  If no
    signals are found, returns ``("WT-like", 5)``.
    """
    if not isinstance(description, str):
        return ("WT-like", 5)

    desc_lower = description.lower()
    segments = [seg.strip() for seg in description.split(",")]

    matched: list[GrowthSignal] = []
    primary_cats: set[str] = set()    # plain signal names from primary segments
    modifier_cats: set[str] = set()   # prefixed names from modifier-led segments

    for seg in segments:
        seg_signals = _detect_segment_signals(seg)
        if not seg_signals:
            continue
        mod = _get_modifier_prefix(seg)
        for s in seg_signals:
            matched.append(s)
            name = _signal_to_name(s, mod)
            if mod is None:
                primary_cats.add(name)
            else:
                modifier_cats.add(name)

    # No growth signal → WT-like
    if not matched:
        return ("WT-like", 5)

    # Merge: primary always wins over modifier for the same base signal.
    # Track which base signals already have a plain (primary) entry.
    primary_bases: set[str] = set()
    for seg in segments:
        if _get_modifier_prefix(seg) is not None:
            continue
        for s in _detect_segment_signals(seg):
            primary_bases.add(s.category)

    all_cats: set[str] = set(primary_cats)
    for mcat in modifier_cats:
        # Extract base signal name by removing modifier prefix
        base = mcat
        for mod_word in MODIFIER_WORDS:
            if mcat.startswith(mod_word + " "):
                base = mcat[len(mod_word) + 1:]
                break
        # Map "divided" back to "germinated and divided" for comparison
        if base == "divided":
            base = "germinated and divided"
        if base not in primary_bases:
            all_cats.add(mcat)

    # Sort categories: primary (non-modifier) first by tier ascending,
    # then modifier-led categories after their primary counterpart.
    # e.g. "spores" (tier 1) before "germinated" (tier 2);
    #      "spores" before "some germinated" (modifier-led, secondary).
    # Map each category back to its base signal's tier for sorting.
    # Build a lookup: category_name → sort_key (tier, is_modifier)
    cat_sort_key: dict[str, tuple[int, int]] = {}
    for s in GROWTH_SIGNALS:
        cat_sort_key[s.category] = (s.tier, 0)
    for cat in all_cats:
        if cat not in cat_sort_key:
            # Modifier-prefixed category — find its base
            base = cat
            for mod_word in MODIFIER_WORDS:
                if cat.startswith(mod_word + " "):
                    base = cat[len(mod_word) + 1:]
                    break
            if base == "divided":
                base = "germinated and divided"
            base_tier = next(
                (s.tier for s in GROWTH_SIGNALS if s.category == base), 5
            )
            cat_sort_key[cat] = (base_tier, 1)

    unique_categories: list[str] = sorted(all_cats, key=lambda c: cat_sort_key.get(c, (5, 1)))

    # Post-process: remove redundant broader categories

    # 1. "spores" is implied by "germinated" when it only appears as part of
    #    "germinated spores" (i.e. no standalone spores).
    has_germinated = any("germinated" in cat for cat in unique_categories)
    if "spores" in unique_categories and has_germinated:
        if not _has_standalone_spores(desc_lower):
            unique_categories.remove("spores")

    # 2. "germinated" is implied by any "divided" variant
    germinated_variants = {
        c for c in unique_categories
        if "divided" in c and "germinated" in c
    }
    if germinated_variants and "germinated" in unique_categories:
        unique_categories.remove("germinated")

    # 3. "small colonies" is implied by "very small colonies"
    if "very small colonies" in unique_categories and "small colonies" in unique_categories:
        unique_categories.remove("small colonies")

    # Assign growth tier: best (highest) tier among matched signals
    growth_tier = max(s.tier for s in matched)

    composed = ", ".join(unique_categories)
    return (composed, growth_tier)


def tier_label(tier: int) -> str:
    """Return a human-readable label for a growth tier number."""
    return TIER_LABELS.get(tier, f"Unknown tier {tier}")


def classify_growth_batch(descriptions: list[str]) -> list[tuple[str, int]]:
    """Apply ``classify_growth`` to a list of descriptions."""
    return [classify_growth(d) for d in descriptions]
