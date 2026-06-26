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
    # "some germination long" is a modifier segment → filtered → WT-like
    ("VIABLE WT cells, some germination long at 25,32,",                  "WT-like", 5),
    ("VIABLE  WT cells, some germination long at 25,32",                  "WT-like", 5),
    # "initially branched" is a modifier segment → filtered → WT-like
    ("VIABLE slightly misshapen cells, initially branched septated slightly long, abnormal colony morphology at 25,32",
     "WT-like", 5),

    # ── Spores (tier 1) ──────────────────────────────────────────────────
    ("ESSENTIAL spores  at 25,32",                                        "spores", 1),
    ("ESSENTIAL spores at 25,32",                                         "spores", 1),
    # "some germinated spores" is modifier-led → Single → spores only
    ("ESSENTIAL spores, some germinated spores at 25,32",                 "spores", 1),

    # ── Germinated (tier 2) ──────────────────────────────────────────────
    ("ESSENTIAL germinated spores long at 25,32",                         "germinated", 2),
    ("ESSENTIAL germinated spores at 25,32",                              "germinated", 2),
    # "some division" is modifier-led → filtered → germinated only
    ("ESSENTIAL germinated spores slightly misshapen, some division at 25,32",
     "germinated", 2),
    ("ESSENTIAL germinated spores long, some division at 25,32",          "germinated", 2),
    # "occasionally long and then divide" is modifier-led → filtered → germinated
    ("ESSENTIAL germinated spores, occasionally long and then divide at 25,32",
     "germinated", 2),
    # "often divide once" is modifier-led → filtered → germinated
    ("ESSENTIAL germinated spores, often divide once to give one slightly misshapen cell and one dead cell at 25,32",
     "germinated", 2),

    # ── Germinated and divided (tier 2) ──────────────────────────────────
    ("ESSENTIAL germinated spores slightly misshapen and divide once at 25,32",
     "germinated and divided", 2),
    # "divide once or twice" — divide is NOT a modifier → kept → divided
    ("ESSENTIAL germinated spores slightly misshapen, divide once or twice at 25,32",
     "germinated and divided", 2),

    # ── Microcolonies (tier 3) ───────────────────────────────────────────
    ("ESSENTIAL microcolonies misshapen cells at 25,32",                  "microcolonies", 3),
    ("ESSENTIAL microcolonies skittle cells at 25,32",                    "microcolonies", 3),
    # "some long cells" is modifier-led → filtered → microcolonies
    ("ESSENTIAL microcolonies slightly misshapen cells, some long cells at 25,32",
     "microcolonies", 3),

    # ── Small colonies (tier 4) ──────────────────────────────────────────
    ("VIABLE small colonies slightly misshapen cells at 25,32",           "small colonies", 4),
    ("VIABLE small colonies long cells at 25,32",                         "small colonies", 4),
    # "possibly diploidising" is modifier-led → filtered → small colonies
    ("VIABLE small colonies long cells, possibly diploidising at 25,32",  "small colonies", 4),

    # ── Very small colonies (tier 4) ─────────────────────────────────────
    ("VIABLE very small colonies rounded cells at 25,32",                 "very small colonies", 4),

    # ── Combined: spores + germinated (mixed population, tier 2) ─────────
    ("ESSENTIAL spores, germinated spores at 25,32",                      "germinated, spores", 2),
    # "may divide once" — may is a modifier → that segment is modifier-led?
    # Actually "germinated spores slightly misshapen may divide once" is one
    # segment (no comma before "may"), so may does not start a segment.
    # Result: spores + germinated + divided
    ("ESSENTIAL spores, germinated spores slightly misshapen may divide once at 25,32",
     "germinated and divided, spores", 2),

    # ── Combined: spores + germinated + microcolonies (tier 3) ───────────
    ("ESSENTIAL spores, germinated spores, microcolonies misshapen cells at 25,32",
     "germinated, microcolonies, spores", 3),
    # "occasionally misshapen branched" is modifier-led → filtered, 3 signals remain
    ("ESSENTIAL spores, germinated spores, microcolonies long cells, occasionally misshapen branched at 25,32",
     "germinated, microcolonies, spores", 3),
    # "microcolonies slightly misshapen" — microcolonies starts segment → kept
    ("ESSENTIAL spores, germinated spores, microcolonies slightly misshapen cells at 25,32",
     "germinated, microcolonies, spores", 3),

    # ── Combined: germinated + microcolonies (tier 3) ────────────────────
    ("ESSENTIAL misshapen germinated spores and microcolonies misshapen cells at 32",
     "germinated, microcolonies", 3),

    # ── Combined: germinated + spores + small colonies (tier 4) ──────────
    ("ESSENTIAL spores, germinated spores, small colonies long cells at 25,32",
     "germinated, small colonies, spores", 4),

    # ── Combined: microcolonies + small colonies (tier 4) ────────────────
    ("ESSENTIAL microcolonies skittle cells, small colonies WT cells at 25,32",
     "microcolonies, small colonies", 4),

    # ── Combined: spores + microcolonies (tier 3) ────────────────────────
    ("ESSENTIAL spores, microcolonies misshapen cells at 25,32",          "microcolonies, spores", 3),
]


@pytest.mark.parametrize("description, expected_cat, expected_tier", CASES)
def test_classify_growth(description: str, expected_cat: str, expected_tier: int):
    """Given a real phenotype description, the classification matches expectation."""
    cat, tier = classify_growth(description)
    assert cat == expected_cat, f"{description!r}: expected {expected_cat!r}, got {cat!r}"
    assert tier == expected_tier, f"{description!r}: expected tier {expected_tier}, got {tier}"