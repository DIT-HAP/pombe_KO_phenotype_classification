"""End-to-end classification tests using real phenotype descriptions.

Each test case is a single line:
    (description, expected_category, expected_tier)

All descriptions are taken verbatim from the Hayles 2013 supplementary
table — the same text the pipeline processes in production.
"""

from __future__ import annotations

import pytest
from growth_signals import classify_growth


# =============================================================================
# (description, expected_category, expected_tier)
#
# All descriptions are verbatim from the source data.
# =============================================================================

CASES: list[tuple[str, str, int]] = [

    # ── WT-like (tier 5) — no growth signal ──────────────────────────────
    ("VIABLE WT cells at 25,32",                                          "WT-like", 5),
    ("VIABLE misshapen cells at 25,32",                                   "WT-like", 5),
    ("VIABLE slightly long cells at 25,32",                               "WT-like", 5),
    # "some germination long" — germination is NOT a signal, modifier-only segment
    ("VIABLE WT cells, some germination long at 25,32,",                  "WT-like", 5),
    ("VIABLE  WT cells, some germination long at 25,32",                  "WT-like", 5),
    # "initially branched" is modifier → no growth signal in any segment
    ("VIABLE slightly misshapen cells, initially branched septated slightly long, abnormal colony morphology at 25,32",
     "WT-like", 5),

    # ── Spores (tier 1) ──────────────────────────────────────────────────
    ("ESSENTIAL spores  at 25,32",                                        "spores", 1),
    ("ESSENTIAL spores at 25,32",                                         "spores", 1),

    # ── Germinated (tier 2) ──────────────────────────────────────────────
    ("ESSENTIAL germinated spores long at 25,32",                         "germinated", 2),
    ("ESSENTIAL germinated spores at 25,32",                              "germinated", 2),
    # "some division" is modifier-led → "some divided" (distinguished from "divided")
    ("ESSENTIAL germinated spores slightly misshapen, some division at 25,32",
     "germinated, some divided", 2),
    ("ESSENTIAL germinated spores long, some division at 25,32",          "germinated, some divided", 2),
    # "occasionally long and then divide" → "occasionally divided"
    ("ESSENTIAL germinated spores, occasionally long and then divide at 25,32",
     "germinated, occasionally divided", 2),
    # "often divide once" → "often divided"
    ("ESSENTIAL germinated spores, often divide once to give one slightly misshapen cell and one dead cell at 25,32",
     "germinated, often divided", 2),

    # ── Germinated and divided (tier 2) — definitive division ────────────
    ("ESSENTIAL germinated spores slightly misshapen and divide once at 25,32",
     "germinated and divided", 2),
    # "divide once or twice" — divide is NOT a modifier → definitive
    ("ESSENTIAL germinated spores slightly misshapen, divide once or twice at 25,32",
     "germinated and divided", 2),

    # ── Microcolonies (tier 3) ───────────────────────────────────────────
    ("ESSENTIAL microcolonies misshapen cells at 25,32",                  "microcolonies", 3),
    ("ESSENTIAL microcolonies skittle cells at 25,32",                    "microcolonies", 3),
    # "some long cells" is modifier but no growth signal → microcolonies only
    ("ESSENTIAL microcolonies slightly misshapen cells, some long cells at 25,32",
     "microcolonies", 3),

    # ── Small colonies (tier 4) ──────────────────────────────────────────
    ("VIABLE small colonies slightly misshapen cells at 25,32",           "small colonies", 4),
    ("VIABLE small colonies long cells at 25,32",                         "small colonies", 4),
    # "possibly diploidising" — no growth signal in modifier segment
    ("VIABLE small colonies long cells, possibly diploidising at 25,32",  "small colonies", 4),

    # ── Very small colonies (tier 4) ─────────────────────────────────────
    ("VIABLE very small colonies rounded cells at 25,32",                 "very small colonies", 4),

    # ── Combined: spores + germinated (mixed population, tier 2) ─────────
    ("ESSENTIAL spores, germinated spores at 25,32",                      "spores, germinated", 2),
    # "may divide once" — may is mid-segment (not segment start) → definitive divided
    ("ESSENTIAL spores, germinated spores slightly misshapen may divide once at 25,32",
     "spores, germinated and divided", 2),

    # ── Combined: spores + some germinated (low-frequency germination) ───
    # "some germinated spores" → some germinated (distinguished from definitive)
    ("ESSENTIAL spores, some germinated spores at 25,32",                 "spores, some germinated", 2),

    # ── Combined: three parallel growth signals (tier 3) ─────────────────
    ("ESSENTIAL spores, germinated spores, microcolonies misshapen cells at 25,32",
     "spores, germinated, microcolonies", 3),
    # "occasionally misshapen branched" — modifier but no growth signal → 3 signals remain
    ("ESSENTIAL spores, germinated spores, microcolonies long cells, occasionally misshapen branched at 25,32",
     "spores, germinated, microcolonies", 3),
    # "microcolonies slightly misshapen" — microcolonies starts segment → primary
    ("ESSENTIAL spores, germinated spores, microcolonies slightly misshapen cells at 25,32",
     "spores, germinated, microcolonies", 3),

    # ── Combined: germinated + microcolonies (tier 3) ────────────────────
    ("ESSENTIAL misshapen germinated spores and microcolonies misshapen cells at 32",
     "germinated, microcolonies", 3),

    # ── Combined: germinated + spores + small colonies (tier 4) ──────────
    ("ESSENTIAL spores, germinated spores, small colonies long cells at 25,32",
     "spores, germinated, small colonies", 4),

    # ── Combined: microcolonies + small colonies (tier 4) ────────────────
    ("ESSENTIAL microcolonies skittle cells, small colonies WT cells at 25,32",
     "microcolonies, small colonies", 4),

    # ── Combined: spores + microcolonies (tier 3) ────────────────────────
    ("ESSENTIAL spores, microcolonies misshapen cells at 25,32",          "spores, microcolonies", 3),

    # ── "occasional" normalised to "occasionally" ───────────────────────
    ("ESSENTIAL spores, germinated spores, occasional microcolonies of WT/ rounded cells at 25,32",
     "spores, germinated, occasionally microcolonies", 3),
    ("ESSENTIAL spores, germinated spores, occasionally microcolonies misshapen cells at 25,32",
     "spores, germinated, occasionally microcolonies", 3),

    # ── Multi-temperature: "germinated spores and colonies" → unified ───
    ("ESSENTIAL spores, germinated spores, microcolonies slightly misshapen cells at 32, spores, germinated spores at 25",
     "spores, germinated, microcolonies", 3),
]


@pytest.mark.parametrize("description, expected_cat, expected_tier", CASES)
def test_classify_growth(description: str, expected_cat: str, expected_tier: int):
    """Given a real phenotype description, the classification matches expectation."""
    cat, tier = classify_growth(description)
    assert cat == expected_cat, f"{description!r}: expected {expected_cat!r}, got {cat!r}"
    assert tier == expected_tier, f"{description!r}: expected tier {expected_tier}, got {tier}"