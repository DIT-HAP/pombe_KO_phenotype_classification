"""Unit tests for the growth-signal detection engine.

Tests :func:`growth_signals.classify_growth` exhaustively:
  - each canonical signal in isolation
  - combined / mixed signals (spores + germinated, etc.)
  - multi-word phrases (small colonies, very small colonies)
  - germinated-spores relationship (standalone vs implied)
  - default WT-like fallback
  - case insensitivity
"""

from __future__ import annotations

import pytest
from growth_signals import classify_growth


# ── Helpers ─────────────────────────────────────────────────────────────────────


def assert_result(desc: str, *, cat: str, tier: int):
    """Assert that *classify_growth* returns the expected (cat, tier)."""
    got_cat, got_tier = classify_growth(desc)
    assert got_cat == cat, f"{desc!r}: expected cat={cat!r}, got {got_cat!r}"
    assert got_tier == tier, f"{desc!r}: expected tier={tier}, got {got_tier}"


# ── WT-like (no growth signal) ──────────────────────────────────────────────────


class TestNoSignal:
    def test_wt_cells(self):
        assert_result("VIABLE WT cells", cat="WT-like", tier=5)

    def test_empty_string(self):
        assert_result("", cat="WT-like", tier=5)

    def test_only_morphology(self):
        assert_result("VIABLE misshapen cells at 25,32", cat="WT-like", tier=5)

    def test_only_modifier(self):
        assert_result("VIABLE slightly long cells", cat="WT-like", tier=5)

    def test_viable_essential_no_signal(self):
        assert_result("VIABLE cells at 25,32", cat="WT-like", tier=5)


# ── Single signals ──────────────────────────────────────────────────────────────


class TestSpores:
    def test_basic(self):
        assert_result("ESSENTIAL spores", cat="spores", tier=1)

    def test_with_modifiers(self):
        assert_result("ESSENTIAL spores misshapen cells", cat="spores", tier=1)

    def test_at_temperature(self):
        assert_result("ESSENTIAL spores at 25,32", cat="spores", tier=1)


class TestGerminated:
    def test_basic(self):
        assert_result("ESSENTIAL germinated spores", cat="germinated", tier=2)

    def test_with_morphology(self):
        assert_result("ESSENTIAL germinated spores long", cat="germinated", tier=2)

    def test_long_form(self):
        """``germination`` (noun form) should also be detected."""
        assert_result("VIABLE germination at 25", cat="germinated", tier=2)

    def test_some_germination_filtered_to_wt_like(self):
        """``classify_growth`` filters modifier segments internally.

        ``some germination long`` starts with ``some`` → dropped, leaving
        ``VIABLE WT cells`` → no signal → WT-like.
        """
        assert_result(
            "VIABLE WT cells, some germination long",
            cat="WT-like", tier=5,
        )

    def test_germinated_spores_no_separate_spores(self):
        """``germinated spores`` alone — spores is implied, not separate."""
        assert_result("ESSENTIAL germinated spores", cat="germinated", tier=2)


class TestGerminatedAndDivided:
    def test_divide_word(self):
        assert_result("ESSENTIAL germinated spores divide",
                       cat="germinated and divided", tier=2)

    def test_division_word(self):
        assert_result("ESSENTIAL germinated spores division",
                       cat="germinated and divided", tier=2)

    def test_divides_plural(self):
        assert_result("germinated spores divides once",
                       cat="germinated and divided", tier=2)

    def test_combined_division_instance(self):
        assert_result("germinated spores, division",
                       cat="germinated and divided", tier=2)


class TestMicrocolonies:
    def test_basic(self):
        assert_result("ESSENTIAL microcolonies", cat="microcolonies", tier=3)

    def test_with_morphology(self):
        assert_result("ESSENTIAL microcolonies misshapen cells",
                       cat="microcolonies", tier=3)


class TestSmallColonies:
    def test_phrase(self):
        assert_result("VIABLE small colonies", cat="small colonies", tier=4)

    def test_phrase_with_morphology(self):
        assert_result("VIABLE small colonies long cells",
                       cat="small colonies", tier=4)


class TestVerySmallColonies:
    def test_phrase(self):
        assert_result("VIABLE very small colonies", cat="very small colonies", tier=4)

    def test_phrase_with_morphology(self):
        assert_result("VIABLE very small colonies misshapen cells",
                       cat="very small colonies", tier=4)


# ── Multi-signal / Combined ─────────────────────────────────────────────────────


class TestCombinedSignals:
    def test_spores_and_germinated(self):
        """Standalone ``spores`` + ``germinated`` → mixed population."""
        assert_result("spores, germinated spores",
                       cat="germinated, spores", tier=2)

    def test_spores_germinated_microcolonies(self):
        desc = "spores, germinated spores, microcolonies"
        assert_result(desc, cat="germinated, microcolonies, spores", tier=3)

    def test_spores_germinated_microcolonies_modifier(self):
        """``occasionally misshapen branched`` is a modifier segment → ignored."""
        desc = "ESSENTIAL spores, germinated spores, microcolonies long cells, occasionally misshapen branched"
        # After filter_primary_segments: "spores, germinated spores, microcolonies long cells"
        assert_result(desc, cat="germinated, microcolonies, spores", tier=3)

    def test_spores_germinated_mixed(self):
        assert_result("ESSENTIAL spores, germinated spores",
                       cat="germinated, spores", tier=2)

    def test_germinated_microcolonies(self):
        assert_result("germinated spores, microcolonies",
                       cat="germinated, microcolonies", tier=3)


# ── Standalone Spores (the _has_standalone_spores logic) ────────────────────────


class TestStandaloneSpores:
    def test_germinated_spores_implies_no_standalone(self):
        """``germinated spores`` without a separate ``spores`` term → only germinated."""
        assert_result("germinated spores long",
                       cat="germinated", tier=2)

    def test_spores_and_germinated_spores_standalone(self):
        """``spores, germinated spores`` → mixed population."""
        assert_result("spores, germinated spores",
                       cat="germinated, spores", tier=2)


# ── Single / Multiple equivalence (classify_growth is independent of PC) ─────────


class TestCrossCutting:
    def test_single_growth_signal_whole_descriptions(self):
        """A single continuous description with one signal → correct category."""
        assert_result("VIABLE small colonies long cells",
                       cat="small colonies", tier=4)

    def test_microcolonies_with_morphology(self):
        """Morphology words (misshapen, etc.) do NOT trigger growth signals."""
        assert_result("ESSENTIAL microcolonies misshapen cells",
                       cat="microcolonies", tier=3)


# ── String matching / case handling ─────────────────────────────────────────────


class TestRobustness:
    def test_case_insensitivity(self):
        assert_result("Essential Spores", cat="spores", tier=1)

    def test_leading_trailing_whitespace(self):
        assert_result("  ESSENTIAL spores  ", cat="spores", tier=1)

    @pytest.mark.parametrize(
        "desc, expected_cat, expected_tier",
        [
            ("germinated spores", "germinated", 2),
            ("germinated, microcolonies", "germinated, microcolonies", 3),
            ("spores only here", "spores", 1),
            ("small colonies", "small colonies", 4),
            ("WT-like", "WT-like", 5),
        ],
    )
    def test_parametrized(self, desc, expected_cat, expected_tier):
        assert_result(desc, cat=expected_cat, tier=expected_tier)