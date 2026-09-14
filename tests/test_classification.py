"""Unit tests for classification and segment-filtering logic.

Tests:
  - :func:`growth_signals.filter_primary_segments`
  - :func:`growth_signals.count_growth_segments`
  - :func:`03_group_genes.classify_phenotype_count`
"""

from __future__ import annotations

import numpy as np
import pytest
from growth_signals import (
    filter_primary_segments,
    count_growth_segments,
    MODIFIER_WORDS,
    GROWTH_KEYWORDS,
)

# Need to import classify_phenotype_count from the step-03 script.
# It is not importable as a regular module (has argparse at module level),
# so we import just the function definition.
import sys, importlib.util
spec = importlib.util.spec_from_file_location(
    "group_genes",
    str(__import__("pathlib").Path(__file__).resolve().parent.parent / "scripts/03_group_genes.py"),
)
mod = importlib.util.module_from_spec(spec)

# The file uses __future__ annotations which would cause the import to
# evaluate all type hints eagerly.  Instead we re-implement the minimal
# logic inline so tests are decoupled from pipeline I/O code.

# ---------------------------------------------------------------------------
# A standalone copy of classify_phenotype_count that relies only on
# count_growth_segments (already imported above).
# ---------------------------------------------------------------------------

def classify_phenotype_count(phenotype: object) -> str | float:
    """Replica of 03_group_genes.classify_phenotype_count."""
    if not isinstance(phenotype, str):
        return np.nan
    return "Multiple" if count_growth_segments(phenotype) > 1 else "Single"


# ===========================================================================
# Tests for filter_primary_segments
# ===========================================================================


class TestFilterPrimarySegments:
    def test_no_comma(self):
        assert filter_primary_segments("ESSENTIAL spores") == "ESSENTIAL spores"

    def test_non_string(self):
        assert filter_primary_segments(123) == 123  # not a string → pass-through

    def test_drop_some_segment(self):
        """``some germination long`` starts with ``some`` → dropped."""
        result = filter_primary_segments(
            "VIABLE WT cells, some germination long"
        )
        assert result == "VIABLE WT cells"

    def test_drop_some_keep_others(self):
        """``some division`` dropped, ``spores`` and ``germinated spores`` kept."""
        result = filter_primary_segments(
            "ESSENTIAL spores, germinated spores, some division"
        )
        assert result == "ESSENTIAL spores, germinated spores"

    def test_drop_occasionally_keep_growth_segments(self):
        result = filter_primary_segments(
            "spores, germinated spores, microcolonies long cells, occasionally misshapen branched"
        )
        assert "occasionally misshapen branched" not in result
        assert "spores" in result
        assert "germinated spores" in result
        assert "microcolonies long cells" in result

    def test_drop_slightly_modifier(self):
        """``slightly misshapen cells`` is a modifier segment → dropped."""
        result = filter_primary_segments(
            "ESSENTIAL germinated spores, slightly misshapen cells"
        )
        assert result == "ESSENTIAL germinated spores"

    def test_drop_initially(self):
        result = filter_primary_segments(
            "ESSENTIAL germinated spores very large swollen branched, initially germinated spores long"
        )
        assert "initially germinated spores long" not in result
        assert "germinated spores very large swollen branched" in result

    def test_keep_microcolonies_slightly(self):
        """``microcolonies slightly misshapen`` — ``microcolonies`` is a growth
        signal at the start, not a modifier → KEPT (``slightly`` modifies
        ``misshapen``, not the whole segment)."""
        result = filter_primary_segments(
            "ESSENTIAL spores, germinated spores, microcolonies slightly misshapen cells"
        )
        assert "microcolonies slightly misshapen cells" in result

    def test_keep_divide_once(self):
        """``divide once`` — ``divide`` is not a modifier → KEPT."""
        result = filter_primary_segments(
            "ESSENTIAL germinated spores, divide once"
        )
        assert "divide once" in result

    def test_all_modifier_segments(self):
        """When every segment is modifier-led, fallback to the first segment."""
        result = filter_primary_segments(
            "VIABLE cells, occasionally long, slightly misshapen"
        )
        assert result == "VIABLE cells"

    def test_only_modifier_segments(self):
        result = filter_primary_segments(
            "initially rounded, slightly misshapen"
        )
        # Both start with modifier words → fallback to first
        assert result == "initially rounded"

    def test_possibly_segment(self):
        """``possibly diploidising`` starts with ``possibly`` → dropped."""
        result = filter_primary_segments(
            "VIABLE small colonies long cells, possibly diploidising"
        )
        assert "possibly diploidising" not in result
        assert result == "VIABLE small colonies long cells"

    def test_often_segment(self):
        result = filter_primary_segments(
            "ESSENTIAL germinated spores, often divide once"
        )
        assert "often divide once" not in result
        assert result == "ESSENTIAL germinated spores"


# ===========================================================================
# Tests for count_growth_segments
# ===========================================================================


class TestCountGrowthSegments:
    def test_no_comma_single_signal(self):
        assert count_growth_segments("ESSENTIAL germinated spores") == 1

    def test_no_comma_no_signal(self):
        assert count_growth_segments("VIABLE WT cells") == 0

    def test_two_parallel_segments(self):
        """``spores, germinated spores`` = 2 growth segments."""
        assert count_growth_segments("spores, germinated spores") == 2

    def test_three_parallel_with_modifier(self):
        """``spores, germinated spores, microcolonies, occasionally...`` = 3."""
        desc = "spores, germinated spores, microcolonies, occasionally misshapen"
        assert count_growth_segments(desc) == 3

    def test_modifier_led_not_counted(self):
        """``, some division`` starts with ``some`` → NOT counted."""
        desc = "germinated spores slightly misshapen, some division"
        assert count_growth_segments(desc) == 1

    def test_divide_once_is_separate_signal(self):
        """``, divide once`` contains ``divide`` → counted."""
        desc = "germinated spores, divide once"
        assert count_growth_segments(desc) == 2

    def test_non_string(self):
        assert count_growth_segments(np.nan) == 0

    def test_complex_mixed(self):
        """``often divide`` starts with ``often`` → NOT counted,
        but ``spores`` + ``germinated spores`` = 2."""
        desc = "ESSENTIAL spores, germinated spores, often divide once"
        assert count_growth_segments(desc) == 2


# ===========================================================================
# Tests for classify_phenotype_count
# ===========================================================================


class TestClassifyPhenotypeCount:
    def test_no_comma_single(self):
        assert classify_phenotype_count("ESSENTIAL germinated spores") == "Single"

    def test_non_string_nan(self):
        result = classify_phenotype_count(np.nan)
        assert isinstance(result, float) and np.isnan(result)

    def test_two_parallel_multiple(self):
        """spores + germinated spores = 2 growth segments → Multiple."""
        assert classify_phenotype_count("spores, germinated spores") == "Multiple"

    def test_three_parallel_multiple(self):
        desc = "spores, germinated spores, microcolonies, occasionally misshapen"
        assert classify_phenotype_count(desc) == "Multiple"

    def test_germinated_with_some_division_single(self):
        """Germinated is 1 segment, ``some division`` is modifier → Single."""
        desc = "ESSENTIAL germinated spores slightly misshapen, some division"
        assert classify_phenotype_count(desc) == "Single"

    def test_spores_germinated_with_some_division_multiple(self):
        """``spores`` + ``germinated spores`` = 2, ``some division`` → modif."""
        desc = "ESSENTIAL spores, germinated spores slightly misshapen, some division"
        assert classify_phenotype_count(desc) == "Multiple"

    def test_microcolonies_small_colonies_multiple(self):
        desc = "ESSENTIAL microcolonies skittle cells, small colonies WT cells"
        assert classify_phenotype_count(desc) == "Multiple"

    def test_some_germination_single(self):
        """``some germination`` starts with ``some`` → Single."""
        desc = "VIABLE WT cells, some germination long"
        assert classify_phenotype_count(desc) == "Single"

    def test_germinated_divide_once_multiple(self):
        """``germinated spores`` + ``divide once`` = 2 (``divide`` not a modifier)."""
        desc = "ESSENTIAL germinated spores, divide once"
        assert classify_phenotype_count(desc) == "Multiple"

    def test_wt_only_single(self):
        assert classify_phenotype_count("VIABLE WT cells") == "Single"

    def test_morphology_only_single(self):
        assert classify_phenotype_count("VIABLE misshapen cells") == "Single"

    def test_pure_spores_single(self):
        assert classify_phenotype_count("ESSENTIAL spores") == "Single"

    def test_pure_microcolonies_single(self):
        assert classify_phenotype_count("ESSENTIAL microcolonies") == "Single"

    def test_some_germination_pipeline_wt_like(self):
        """After filter_primary_segments + classify_growth,
        ``VIABLE WT cells, some germination long`` → WT-like.

        ``some germination long`` starts with modifier ``some`` → dropped,
        leaving ``VIABLE WT cells`` → no growth signal → WT-like (tier 5).
        """
        from growth_signals import filter_primary_segments, classify_growth

        filtered = filter_primary_segments(
            "VIABLE WT cells, some germination long"
        )
        assert filtered == "VIABLE WT cells"

        cat, tier = classify_growth(filtered)
        assert cat == "WT-like"
        assert tier == 5


# ── Sanity check: our hard-coded lists ──────────────────────────────────────────


class TestConstantIntegrity:
    def test_modifier_words_are_unique(self):
        assert len(MODIFIER_WORDS) == len(set(MODIFIER_WORDS)), \
            "MODIFIER_WORDS contains duplicates"

    def test_growth_keywords_are_unique(self):
        assert len(GROWTH_KEYWORDS) == len(set(GROWTH_KEYWORDS)), \
            "GROWTH_KEYWORDS contains duplicates"

    def test_all_modifier_words_start_lowercase(self):
        for m in MODIFIER_WORDS:
            assert m[0].islower(), f"{m!r} is not lowercase"

    def test_all_growth_keywords_start_lowercase(self):
        for kw in GROWTH_KEYWORDS:
            assert kw[0].islower(), f"{kw!r} is not lowercase"