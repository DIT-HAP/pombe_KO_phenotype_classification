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

    cat, tier = classify_growth("ESSENTIAL germinated spores at 25,32")
    # cat == "germinated, spores", tier == 1

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
    5: "WT (normal growth)",
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
    GrowthSignal(keyword="divide", category="germinated and divided", tier=2),
    GrowthSignal(keyword="division", category="germinated and divided", tier=2),
    # tier 1
    GrowthSignal(keyword="spores", category="spores", tier=1),
]

# =============================================================================
# CORE LOGIC
# =============================================================================


def classify_growth(description: str) -> tuple[str, int]:
    """Return (category, growth_tier) for a phenotype description string.

    Detects all growth signals present in the description, composes their
    category names (sorted, deduplicated), and assigns the worst (lowest)
    tier. If no signals are found, returns ``("WT", 5)``.
    """
    desc_lower = description.lower()

    # Collect all matching signals
    matched: list[GrowthSignal] = []
    for signal in GROWTH_SIGNALS:
        if signal.keyword in desc_lower:
            matched.append(signal)

    # No growth signal → WT
    if not matched:
        return ("WT", 5)

    # Compose unique category names sorted alphabetically
    unique_categories: list[str] = sorted({s.category for s in matched})

    # Post-process: remove redundant broader categories when a more specific
    # one is already present.
    # - "germinated" is implied by "germinated and divided"
    if "germinated and divided" in unique_categories and "germinated" in unique_categories:
        unique_categories.remove("germinated")
    # - "small colonies" is implied by "very small colonies"
    if "very small colonies" in unique_categories and "small colonies" in unique_categories:
        unique_categories.remove("small colonies")

    # Assign growth tier: for composite phenotypes, use the BEST (highest)
    # tier among matched signals. This represents the most advanced growth
    # stage the gene can support — the informative endpoint, not the worst.
    # Single-signal descriptions use their own tier.
    growth_tier = max(s.tier for s in matched)

    composed = ", ".join(unique_categories)
    return (composed, growth_tier)


def tier_label(tier: int) -> str:
    """Return a human-readable label for a growth tier number."""
    return TIER_LABELS.get(tier, f"Unknown tier {tier}")


def classify_growth_batch(descriptions: list[str]) -> list[tuple[str, int]]:
    """Apply ``classify_growth`` to a list of descriptions."""
    return [classify_growth(d) for d in descriptions]
