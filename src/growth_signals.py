"""
Growth Signal Definitions for S. pombe Phenotype Classification
===============================================================

Shared module defining growth-related keyword signals, their severity
tiers, and the ``classify_growth()`` engine that replaces the brittle
if-elif chains in the original categorize scripts.

Signals are detected independently from a phenotype description, then
resolved into a fine-grained ``Category`` label and a **signal tier**
(1-5). The signal tier reflects biological hierarchy (1=spores, 2=germinated,
3=microcolonies, 4=small colonies, 5=WT-like) and is used internally for
category sorting. The final ``Growth_tier`` in pipeline output is re-ranked
by DIT-HAP DR median (see ``04_categorize_phenotypes.py`` and
``05_merge_categories.py``).

Input
-----
- A phenotype description string (e.g. from ``Deletion mutant phenotype
  description`` column in the Hayles table).

Output
------
- ``tuple[str, int]`` - (category_name, signal_tier).

Usage
-----
    from scripts.growth_signals import classify_growth

    cat, tier = classify_growth("VIABLE WT cells at 25,32")
    # cat == "WT-like", tier == 5

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-06-09
Version:  1.1.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from dataclasses import dataclass
import re

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================


@dataclass(frozen=True)
class GrowthSignal:
    """One detectable keyword pattern and its corresponding growth category."""

    keyword: str  # substring to search for (case-insensitive)
    category: str  # fine-grained category label
    tier: int  # growth severity (1 = worst, 5 = best)


# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

# Canonical signal table, ordered from most to least specific.
# When multiple keywords match, category names are concatenated (alphabetically
# sorted, deduplicated) and the worst tier is assigned.
GROWTH_SIGNALS: list[GrowthSignal] = [
    # tier 4 - most specific phrases first to avoid substring collision
    GrowthSignal(keyword="very small colon", category="very small colonies", tier=4),
    GrowthSignal(keyword="small colon", category="small colonies", tier=4),
    # tier 3
    GrowthSignal(keyword="microcolonies", category="microcolonies", tier=3),
    # tier 2
    GrowthSignal(keyword="germinated", category="germinated", tier=2),
    # Note: "germination" (noun) is NOT a growth signal - in the data it
    # only appears in morphological contexts ("germination long", "poor
    # germination", "from germination"), never as a standalone growth state.
    GrowthSignal(keyword="divide", category="germinated and divided", tier=2),
    GrowthSignal(keyword="division", category="germinated and divided", tier=2),
    # tier 1
    GrowthSignal(keyword="spores", category="spores", tier=1),
]

# Modifier words - when a comma-separated segment starts with one of these,
# it is a secondary description (modifier / morphology supplement) and is
# excluded from growth-signal detection. Kept grouped by type so callers
# (e.g. 02_keyword_profile) can label them; MODIFIER_WORDS is their union.
MODIFIER_FREQ: tuple[str, ...] = (
    "occasionally",
    "often",
    "occasional",
    "sometimes",
    "mostly",
    "rarely",
    "frequently",
    "frequency",
    "rare",
    "possible",
    "may",
    "possibly",
    "rapidly",
    "initially",
)
MODIFIER_DEGREE: tuple[str, ...] = (
    "slightly",
    "very",
    "highly",
    "barely",
    "slight",
    "high",
    "weak",
)
MODIFIER_QUANT: tuple[str, ...] = (
    "some",
    "many",
    "few",
    "lots",
    "several",
    "multiple",
    "once",
    "twice",
    "more",
    "multi",
)
MODIFIER_WORDS: tuple[str, ...] = MODIFIER_FREQ + MODIFIER_DEGREE + MODIFIER_QUANT

# Growth-signal keywords used for segment-level analysis.
GROWTH_KEYWORDS: tuple[str, ...] = (
    "spores",
    "germinated",
    "microcolonies",
    "small colon",
    "very small colon",
    "divide",
    "division",
)

# Within-tier order: germinated (0) before germinated and divided (1).
SUB_ORDER: dict[str, int] = {
    "germinated": 0,
    "germinated and divided": 1,
    "spores": 0,
    "microcolonies": 0,
    "small colonies": 0,
    "very small colonies": 0,
}

# Exact category -> (is_modifier, tier, within-tier order) sort key.
SIGNAL_SORT_KEYS: dict[str, tuple[int, int, int]] = {
    signal.category: (0, signal.tier, SUB_ORDER.get(signal.category, 0))
    for signal in GROWTH_SIGNALS
}

# Category -> signal tier, used to rank modifier-prefixed categories.
SIGNAL_TIERS: dict[str, int] = {signal.category: signal.tier for signal in GROWTH_SIGNALS}

WT_LIKE = "WT-like"
WT_LIKE_TIER = 5

# =============================================================================
# CORE LOGIC
# =============================================================================


def _split_segments(text: str) -> list[str]:
    """Split a description into comma segments, stripping each segment."""
    return [seg.strip() for seg in text.split(",")]


def _is_modifier_segment(seg: str) -> bool:
    """Return True if *seg* starts with a modifier word."""
    stripped = seg.lower().strip()
    return any(stripped.startswith(modifier) for modifier in MODIFIER_WORDS)


def _has_growth_signal(seg: str) -> bool:
    """Return True if *seg* contains a growth-signal keyword."""
    lowered = seg.lower()
    return any(keyword in lowered for keyword in GROWTH_KEYWORDS)


def filter_primary_segments(text: str) -> str:
    """Remove modifier-led comma segments, keep the rest.

    A comma-separated segment that starts with a modifier word (e.g.
    ``some``, ``occasionally``, ``slightly``) is treated as a secondary
    description and excluded. All other segments - whether they contain
    growth signals or pure morphology - are kept for classification.

    If every segment is modifier-led, the first segment is returned as a
    fallback so the description is never empty.
    """
    if not isinstance(text, str) or "," not in text:
        return text
    segments = _split_segments(text)
    primary = [seg for seg in segments if not _is_modifier_segment(seg)]
    return ", ".join(primary) if primary else segments[0]


def count_growth_segments(text: str) -> int:
    """Count comma segments that contain a growth signal (excluding modifier-led ones)."""
    if not isinstance(text, str):
        return 0
    if "," not in text:
        return 1 if _has_growth_signal(text) else 0
    segments = _split_segments(text)
    return sum(1 for seg in segments if not _is_modifier_segment(seg) and _has_growth_signal(seg))


def _has_standalone_spores(text: str) -> bool:
    """Check if ``spores`` appears anywhere NOT as part of ``germinated spores``.

    A standalone occurrence (comma-separated, at the start of a sentence,
    etc.) indicates a mixed population where both dormant and germinated
    spores coexist.
    """
    for match in re.finditer(r"\bspores?\b", text):
        start = match.start()
        # Look at the 12 characters before 'spores'
        before = text[max(0, start - 12) : start].strip().rstrip(",")
        if before != "germinated":
            return True
    return False


def _detect_segment_signals(seg: str) -> list[GrowthSignal]:
    """Return all growth signals found in *seg*."""
    lowered = seg.lower()
    return [signal for signal in GROWTH_SIGNALS if signal.keyword in lowered]


def _get_modifier_prefix(seg: str) -> str | None:
    """If *seg* starts with a modifier word, return that word; else None."""
    stripped = seg.lower().strip()
    for modifier in MODIFIER_WORDS:
        if stripped.startswith(modifier):
            return modifier
    return None


def _normalise_modifier(modifier: str) -> str:
    """Normalise modifier variants to a canonical form."""
    return {"occasional": "occasionally"}.get(modifier, modifier)


def _signal_to_name(signal: GrowthSignal, modifier: str | None) -> str:
    """Convert a signal and modifier into a category name.

    - No modifier -> plain signal category (e.g. ``germinated``).
    - With modifier -> prefixed (e.g. ``some germinated``, ``often divided``).
    """
    if modifier is None:
        return signal.category
    normalised = _normalise_modifier(modifier)
    # "germinated and divided" -> "often divided" (germinated already implied)
    if signal.category == "germinated and divided":
        return f"{normalised} divided"
    return f"{normalised} {signal.category}"


def _collect_categories(
    segments: list[str],
) -> tuple[list[GrowthSignal], set[str], set[str]]:
    """Collect matched signals, primary categories, and modifier-led categories."""
    matched: list[GrowthSignal] = []
    primary_categories: set[str] = set()
    modifier_categories: set[str] = set()
    for seg in segments:
        signals = _detect_segment_signals(seg)
        if not signals:
            continue
        modifier = _get_modifier_prefix(seg)
        target = primary_categories if modifier is None else modifier_categories
        for signal in signals:
            matched.append(signal)
            target.add(_signal_to_name(signal, modifier))
    return matched, primary_categories, modifier_categories


def _primary_bases(segments: list[str]) -> set[str]:
    """Return the base signal categories contributed by non-modifier segments."""
    return {
        signal.category
        for seg in segments
        if _get_modifier_prefix(seg) is None
        for signal in _detect_segment_signals(seg)
    }


def _strip_modifier_prefix(category: str) -> str:
    """Return *category* without a leading modifier word."""
    for modifier in MODIFIER_WORDS:
        if category.startswith(modifier + " "):
            return category[len(modifier) + 1 :]
    return category


def _base_category(category: str) -> str:
    """Return the underlying signal category of a possibly modifier-prefixed name."""
    base = _strip_modifier_prefix(category)
    return "germinated and divided" if base == "divided" else base


def _category_sort_key(category: str) -> tuple[int, int, int]:
    """Sort modifier-led categories after primaries, then by tier and sub-order."""
    exact = SIGNAL_SORT_KEYS.get(category)
    if exact is not None:
        return exact
    base = _base_category(category)
    return (1, SIGNAL_TIERS.get(base, WT_LIKE_TIER), SUB_ORDER.get(base, 0))


def _promote_severest_modifier(keys: dict[str, tuple[int, int, int]]) -> None:
    """Promote one modifier-led category to primary when no primary exists."""
    min_tier = min(key[1] for key in keys.values())
    for category in sorted(keys):
        if keys[category][1] == min_tier:
            keys[category] = (0, min_tier, keys[category][2])
            return


def _merge_categories(
    primary_categories: set[str],
    modifier_categories: set[str],
    primary_bases: set[str],
) -> list[str]:
    """Merge categories, letting a primary signal override its modifier variant."""
    merged = set(primary_categories)
    merged.update(
        category
        for category in modifier_categories
        if _base_category(category) not in primary_bases
    )
    keys = {category: _category_sort_key(category) for category in merged}
    if not any(key[0] == 0 for key in keys.values()):
        _promote_severest_modifier(keys)
    return sorted(keys, key=keys.__getitem__)


def _drop_implied(categories: list[str], desc_lower: str) -> list[str]:
    """Remove categories implied by a more specific matched category."""
    kept = list(categories)

    # 1. "spores" is implied by "germinated" when it only appears as part of
    #    "germinated spores" (i.e. no standalone spores).
    has_germinated = any("germinated" in category for category in kept)
    if "spores" in kept and has_germinated and not _has_standalone_spores(desc_lower):
        kept.remove("spores")

    # 2. "germinated" is implied by any "divided" variant that also names it.
    has_divided_variant = any(
        "divided" in category and "germinated" in category for category in kept
    )
    if has_divided_variant and "germinated" in kept:
        kept.remove("germinated")

    # 3. "small colonies" is implied by "very small colonies".
    if "very small colonies" in kept and "small colonies" in kept:
        kept.remove("small colonies")

    return kept


def classify_growth(description: str) -> tuple[str, int]:
    """Return (category, growth_tier) for a phenotype description string.

    The description is split into comma segments. Each segment is scanned
    for growth signals:

    - **Primary segments** (not starting with a modifier) contribute plain
      signal names (e.g. ``germinated``, ``spores``).
    - **Modifier-led segments** (starting with ``some``, ``often``, etc.)
      contribute prefixed signal names (e.g. ``some germinated``,
      ``often divided``) so that low-frequency events are distinguished
      from definitive ones.

    If a signal appears in both a primary and a modifier-led segment, the
    plain (primary) version is kept - the definitive occurrence dominates.

    Tier is the best (highest) tier among all matched signals. If no
    signals are found, returns ``("WT-like", 5)``.
    """
    if not isinstance(description, str):
        return (WT_LIKE, WT_LIKE_TIER)

    segments = _split_segments(description)
    matched, primary_categories, modifier_categories = _collect_categories(segments)
    if not matched:
        return (WT_LIKE, WT_LIKE_TIER)

    categories = _merge_categories(
        primary_categories,
        modifier_categories,
        _primary_bases(segments),
    )
    kept = _drop_implied(categories, description.lower())
    return (", ".join(kept), max(signal.tier for signal in matched))
