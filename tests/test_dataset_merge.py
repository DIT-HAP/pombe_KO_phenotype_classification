"""Tests for the Kim 2010 + Hayles 2013 union table."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

from dataset_merge import (
    CHANGE_CLASS_COL,
    CONSISTENCY_COL,
    DEFAULT_HAYLES2013,
    DEFAULT_HAYLES_DOCX,
    DEFAULT_KIM2010,
    DEFAULT_KIM_TABLE_S3,
    DISCREPANCY_COL,
    HAYLES_STUDY_COL,
    PRESENCE_COL,
    build_merged_dataset,
    load_hayles2013,
    load_kim2010,
    load_kim_table_s3,
    normalize_call,
    normalize_id,
    parse_table_s3,
    read_docx_tables,
    reconstruct_gene_name_from_date,
    validate_kim_table_s3,
)

ROOT = Path(__file__).resolve().parent.parent


def test_normalize_id_strips_all_whitespace():
    assert normalize_id("  SPAC1F3.02c ") == "SPAC1F3.02c"
    assert normalize_id("SPAC 1002.01") == "SPAC1002.01"
    assert normalize_id("") is None
    assert normalize_id(pd.NA) is None


def test_normalize_call_maps_words_and_stray_spaces():
    assert normalize_call("viable") == "V"
    assert normalize_call("essential") == "E"
    assert normalize_call(" E ") == "E"
    assert normalize_call("V ") == "V"
    assert normalize_call("unknown") is None


@pytest.mark.parametrize(
    ("date", "expected"),
    [
        (datetime(2001, 10, 1), "Oct-1"),
        (datetime(2001, 9, 1), "Sep-1"),
        (datetime(2011, 9, 1), "Sep-11"),
        (datetime(2001, 10, 5), "Oct-5"),
    ],
)
def test_reconstruct_gene_name_from_date(date: datetime, expected: str):
    assert reconstruct_gene_name_from_date(date) == expected


def test_parse_table_s3_repairs_truncated_id():
    tables = [
        [["Cosmid ID Newly constructed deletions", "Gene name"]],
        [
            [
                "Cosmid ID",
                "Gene name",
                "Kim et al 2010 dispensability deletion dataset",
                "Kim et al 2010 dispensability analysis dataset",
                "DispensabilityThis study",
                "Comments",
            ],
            ["SPAC13D6.0", "", "V", "E", "E", "Re-analysed -same as pub"],
        ],
    ]
    table = parse_table_s3(tables)
    assert table["id"].tolist() == ["SPAC13D6.05"]
    assert table["hayles2013_study_s3"].tolist() == ["E"]
    assert table["strain_status"].tolist() == ["Re-analysed"]
    assert table["dispensability_change"].tolist() == ["same (E)"]


def _kim_frame(rows: list[dict]) -> pd.DataFrame:
    """Build a canonical Kim 2010 frame with the loader's columns."""
    defaults = {
        "gene_name_2010": None,
        "pct_ko_2010": 90.0,
        "kim2010_deletion": None,
        "kim2010_analysis": None,
        "published_phenotype": None,
        "constructed": True,
        "construction_method": None,
        "gene_product_2010": None,
        "paralogues_2010": None,
        "sc_orthologue_2010": None,
        "ts_remake": None,
        "source_format": "full_table_s1",
        "kim2010_call": None,
    }
    return pd.DataFrame([{**defaults, **row} for row in rows])


def _hayles_frame(rows: list[dict]) -> pd.DataFrame:
    """Build a canonical Hayles 2013 frame with the loader's columns."""
    defaults = {
        "gene_name_2013": None,
        "working_code": None,
        "pct_ko_2013": 90.0,
        "hayles2013_study": None,
        "phenotype_description": None,
        "phenotype_classification": "WT",
        "gene_product_2013": None,
        "paralogues_2013": None,
        "cell_length": None,
    }
    return pd.DataFrame([{**defaults, **row} for row in rows])


def _merged_fixture() -> pd.DataFrame:
    kim = _kim_frame(
        [
            # both, consistent, 2010 columns agree
            {"id": "A", "gene_name_2010": "a10", "kim2010_deletion": "V",
             "kim2010_analysis": "V", "kim2010_call": "V"},
            # both, inconsistent, documented in Table S3
            {"id": "B", "gene_name_2010": "b10", "kim2010_deletion": "V",
             "kim2010_analysis": "V", "kim2010_call": "V"},
            # both, 2013 matches the 2010 deletion call, not the analysis call
            {"id": "C", "gene_name_2010": "c10", "kim2010_deletion": "V",
             "kim2010_analysis": "E", "kim2010_call": "E"},
            # both, 2013 agrees with the analysis call; the two 2010 columns disagree
            {"id": "G", "gene_name_2010": "g10", "kim2010_deletion": "E",
             "kim2010_analysis": "V", "kim2010_call": "V"},
            # 2013-only: Kim candidate was never constructed
            {"id": "D", "gene_name_2010": "d10", "constructed": False},
            # 2010-only: constructed but absent from 2013
            {"id": "E", "gene_name_2010": "e10", "kim2010_deletion": "V",
             "kim2010_analysis": "V", "kim2010_call": "V"},
            # neither: candidate not constructed, no 2013 row
            {"id": "F", "gene_name_2010": "f10", "constructed": False},
        ]
    )
    hayles = _hayles_frame(
        [
            {"id": "A", "gene_name_2013": "a13", "hayles2013_study": "V",
             "hayles2013_kim2010_deletion": "V",
             "hayles2013_kim2010_analysis": "V"},
            {"id": "B", "gene_name_2013": "b13", "hayles2013_study": "E",
             "hayles2013_kim2010_deletion": "V",
             "hayles2013_kim2010_analysis": "V"},
            {"id": "C", "gene_name_2013": "c13", "hayles2013_study": "V",
             "hayles2013_kim2010_deletion": "V",
             "hayles2013_kim2010_analysis": "E"},
            {"id": "G", "gene_name_2013": "g13", "hayles2013_study": "V",
             "hayles2013_kim2010_deletion": "E",
             "hayles2013_kim2010_analysis": "V"},
            {"id": "D", "gene_name_2013": "d13", "hayles2013_study": "E"},
        ]
    )
    table_s3 = pd.DataFrame({"id": ["B"], "s3_comment": ["Re-analysed"]})
    return build_merged_dataset(kim, hayles, table_s3)


def test_merged_presence_column():
    merged = _merged_fixture().set_index("Systematic ID")

    assert merged.loc["A", PRESENCE_COL] == "both"
    assert merged.loc["B", PRESENCE_COL] == "both"
    assert merged.loc["C", PRESENCE_COL] == "both"
    assert merged.loc["D", PRESENCE_COL] == "2013_only"
    assert merged.loc["E", PRESENCE_COL] == "2010_only"
    assert merged.loc["F", PRESENCE_COL] == "neither"
    assert len(merged) == 7


def test_merged_consistency_column():
    merged = _merged_fixture().set_index("Systematic ID")

    assert merged.loc["A", CONSISTENCY_COL] == "consistent"
    assert merged.loc["B", CONSISTENCY_COL] == "inconsistent"
    assert merged.loc["C", CONSISTENCY_COL] == "inconsistent"
    assert merged.loc["G", CONSISTENCY_COL] == "consistent"
    assert merged.loc["D", CONSISTENCY_COL] == "not_comparable"
    assert merged.loc["E", CONSISTENCY_COL] == "not_comparable"
    assert merged.loc["F", CONSISTENCY_COL] == "not_comparable"


def test_merged_change_classes_and_discrepancy():
    merged = _merged_fixture().set_index("Systematic ID")

    assert merged.loc["A", CHANGE_CLASS_COL] == "no_change"
    assert merged.loc["B", CHANGE_CLASS_COL] == "documented_change"
    assert merged.loc["C", CHANGE_CLASS_COL] == (
        "matches_kim2010_deletion_not_analysis"
    )
    assert merged.loc["G", CHANGE_CLASS_COL] == (
        "kim2010_deletion_analysis_discrepancy"
    )
    assert merged.loc["D", CHANGE_CLASS_COL] == "hayles2013_addition"
    assert merged.loc["E", CHANGE_CLASS_COL] == "no_hayles2013_row"
    assert merged.loc["F", CHANGE_CLASS_COL] == "no_data_either"

    assert merged.loc["A", DISCREPANCY_COL] == "same (V)"
    assert "Kim2010_analysis=V vs Hayles2013_study=E" in (
        merged.loc["B", DISCREPANCY_COL]
    )
    assert "Table S3: Re-analysed" in merged.loc["B", DISCREPANCY_COL]
    assert merged.loc["C", DISCREPANCY_COL].startswith(
        "Kim2010_analysis=E vs Hayles2013_study=V; "
        "matches_kim2010_deletion_not_analysis"
    )
    assert merged.loc["G", DISCREPANCY_COL] == (
        "same (V); Kim2010 deletion=E vs analysis=V"
    )
    assert merged.loc["E", DISCREPANCY_COL] == "no 2013 screen row"
    assert merged.loc["F", DISCREPANCY_COL] == (
        "no 2010 analysis call and no 2013 row"
    )


def test_merged_gene_name_prefers_2013_then_falls_back_to_2010():
    merged = _merged_fixture().set_index("Systematic ID")
    assert merged.loc["A", "Gene name"] == "a13"
    assert merged.loc["E", "Gene name"] == "e10"
    assert merged.loc["E", "Kim2010_Gene_name"] == "e10"


def test_merged_keeps_source_columns():
    merged = _merged_fixture()
    for col in (
        "Systematic ID",
        "Gene name",
        PRESENCE_COL,
        CONSISTENCY_COL,
        DISCREPANCY_COL,
        CHANGE_CLASS_COL,
        "Kim2010_Deletion_dataset",
        "Kim2010_Analysis_dataset",
        HAYLES_STUDY_COL,
        "Phenotypic classification used for analysis",
    ):
        assert col in merged.columns

    # The 2013 copies of Kim's two dispensability columns are not carried.
    assert "Kim et al 2010 Deletion dataset" not in merged.columns
    assert "Kim et al 2010 Analysis dataset" not in merged.columns


def test_merged_without_table_s3_marks_changed_rows_accordingly():
    kim = _kim_frame(
        [
            {"id": "A", "kim2010_deletion": "V", "kim2010_analysis": "V",
             "kim2010_call": "V"},
        ]
    )
    hayles = _hayles_frame(
        [{"id": "A", "hayles2013_study": "E",
          "hayles2013_kim2010_deletion": "V",
          "hayles2013_kim2010_analysis": "V"}]
    )
    merged = build_merged_dataset(kim, hayles, None).set_index("Systematic ID")

    assert merged.loc["A", CONSISTENCY_COL] == "inconsistent"
    assert merged.loc["A", CHANGE_CLASS_COL] == "changed_table_s3_not_available"


def test_real_inputs_build_expected_union():
    kim = load_kim2010(ROOT / DEFAULT_KIM2010)
    hayles = load_hayles2013(ROOT / DEFAULT_HAYLES2013)

    assert len(kim) == 4914
    assert int(kim["constructed"].sum()) == 4836
    assert len(hayles) == 4843

    # Excel converted these three gene names into dates; the loader repairs them.
    names = hayles.set_index("id")["gene_name_2013"]
    assert names["SPAC1F3.10c"] == "Oct-1"
    assert names["SPAC5D6.05"] == "Sep-11"
    assert names["SPBC4C3.12"] == "Sep-1"

    docx = ROOT / DEFAULT_HAYLES_DOCX
    table_s3 = (
        parse_table_s3(read_docx_tables(docx)) if docx.exists() else None
    )
    merged = build_merged_dataset(kim, hayles, table_s3)

    assert len(merged) == 4914
    assert merged["Systematic ID"].is_unique
    assert merged[PRESENCE_COL].value_counts().to_dict() == {
        "both": 4825,
        "neither": 60,
        "2013_only": 18,
        "2010_only": 11,
    }
    assert merged[CONSISTENCY_COL].value_counts().to_dict() == {
        "consistent": 4796,
        "not_comparable": 89,
        "inconsistent": 29,
    }

    classes = merged[CHANGE_CLASS_COL].value_counts().to_dict()
    assert classes["no_change"] == 4789
    assert classes["no_data_either"] == 60
    assert classes["hayles2013_addition"] == 18
    assert classes["no_hayles2013_row"] == 11

    if table_s3 is not None:
        assert len(table_s3) == 21
        assert classes["documented_change"] == 14
        assert classes["documented_discrepancy_in_kim2010_columns"] == 7
        assert classes["matches_kim2010_deletion_not_analysis"] == 11
        assert classes["undocumented_change"] == 4

        assert merged["Hayles2013_Strain_status"].value_counts().to_dict() == {
            "Re-analysed": 12,
            "Deletion remade": 9,
        }
        indexed = merged.set_index("Systematic ID")
        assert indexed.loc["SPAC144.08", "Hayles2013_Strain_status"] == (
            "Deletion remade"
        )
        assert indexed.loc[
            "SPAC144.08", "Hayles2013_2010_dispensability_change"
        ] == "V -> E"
        assert indexed.loc["SPAC13D6.05", "Hayles2013_TableS3_comment"] == (
            "Re-analysed - same as pub"
        )


def test_kim_table_s3_matches_table_s1():
    kim = load_kim2010(ROOT / DEFAULT_KIM2010)
    s3 = load_kim_table_s3(ROOT / DEFAULT_KIM_TABLE_S3)

    assert len(s3) == 40
    validate_kim_table_s3(s3, kim)


def test_kim_table_s3_validation_rejects_transcription_errors(tmp_path):
    kim = load_kim2010(ROOT / DEFAULT_KIM2010)
    source = pd.read_csv(
        ROOT / DEFAULT_KIM_TABLE_S3, sep="\t", dtype=str
    ).fillna("")

    def drop_last_row(frame: pd.DataFrame) -> None:
        frame.drop(frame.index[-1], inplace=True)

    def use_non_disagreement_id(frame: pd.DataFrame) -> None:
        frame.loc[0, "systematic_id"] = "SPAC1002.03c"
        frame.loc[0, ["published_data", "our_data", "final_judgment"]] = [
            "V", "V", "V",
        ]

    def change_criteria(frame: pd.DataFrame) -> None:
        frame.loc[0, "criteria_for_judgment"] = "Disagreement"

    def clear_footnote(frame: pd.DataFrame) -> None:
        frame.loc[3, "criteria_detail"] = ""

    cases = [
        ("dropped_row", drop_last_row, "numbering"),
        ("non_disagreement", use_non_disagreement_id, "disagreement set"),
        ("criteria_drift", change_criteria, "criteria distribution"),
        ("footnote_drift", clear_footnote, "footnote rows"),
    ]

    for label, mutate, message in cases:
        mutated = source.copy()
        mutate(mutated)
        path = tmp_path / f"{label}.tsv"
        mutated.to_csv(path, sep="\t", index=False)
        with pytest.raises(ValueError, match=message):
            validate_kim_table_s3(load_kim_table_s3(path), kim)
