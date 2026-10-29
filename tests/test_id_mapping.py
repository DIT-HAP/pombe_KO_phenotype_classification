"""Tests for the shared PomBase ID-mapping module."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from id_mapping import (
    PomBaseCol,
    build_paralog_columns,
    build_update_maps,
    collect_mapping_ids,
    combine_mapping_notes,
    normalize_sysid,
    read_pombase_annotation,
    update_sysids,
)


@pytest.fixture()
def annotation() -> pd.DataFrame:
    """Small PomBase-like annotation with every lookup layer."""
    return pd.DataFrame(
        {
            PomBaseCol.GENE_SYSID: [
                "SPAC1002.01c",
                "SPAC823.17",
                "SPNCRNA.9001",
                "SPCC622.17",
            ],
            PomBaseCol.GENE_NAME: ["mrx11", "tom6", "prt2", "fake"],
            PomBaseCol.GENE_PRODUCT: ["product A", "product B", "lncRNA", "product D"],
            PomBaseCol.GENE_TYPE: [
                "protein coding gene",
                "protein coding gene",
                "lncRNA gene",
                "pseudogene",
            ],
            PomBaseCol.SYNONYMS: [
                "SPAC1002.01, alias1",
                "SPAC823.02",
                "SPBC8E4.02c",
                None,
            ],
        }
    )


def test_normalize_sysid_case_and_suffix():
    assert normalize_sysid("spac1002.01C") == "SPAC1002.01c"
    assert normalize_sysid("SPAC1556.06.1") == "SPAC1556.06"
    assert normalize_sysid(None) is None


def test_build_update_maps_layers(annotation):
    sysids, names, synonyms = build_update_maps(annotation)

    assert "SPAC1002.01c" in sysids
    assert "SPNCRNA.9001" not in sysids
    assert names.loc["mrx11"] == "SPAC1002.01c"
    assert synonyms.loc["SPAC823.02"] == "SPAC823.17"


def test_update_sysids_all_layers(annotation):
    ids = [
        "SPAC1002.01c",
        "mrx11",
        "SPAC823.02",
        "SPBC8E4.02c",
        "SPAC0000.00",
    ]
    updated, notes = update_sysids(ids, annotation)

    assert updated == [
        "SPAC1002.01c",
        "SPAC1002.01c",
        "SPAC823.17",
        "SPBC8E4.02c",
        "SPAC0000.00",
    ]
    assert notes[0] == "NO_CHANGE"
    assert notes[1].startswith("VIA_NAME")
    assert notes[2].startswith("VIA_SYNONYM")
    assert notes[3].startswith("RECLASSIFIED")
    assert "lncRNA" in notes[3]
    assert notes[4] == "NOT_FOUND"


def test_collect_mapping_ids_deduplicates_main_and_paralog():
    kim = pd.DataFrame(
        {
            "id": ["SPAC1002.01c", "SPAC1002.01c", "SPAC823.17"],
            "paralogues_2010": ["SPAC1111.01|SPAC2222.02", "none", "being curated"],
        }
    )
    hayles = pd.DataFrame({"id": ["SPAC1002.01c", "SPCC622.17"]})

    ids = collect_mapping_ids(kim, hayles)

    assert ids == [
        "SPAC1002.01c",
        "SPAC823.17",
        "SPCC622.17",
        "SPAC1111.01",
        "SPAC2222.02",
    ]


def test_build_paralog_columns_updates_and_drops_self():
    raw = pd.Series(["SPAC823.02|SPAC1002.01C", "being curated", None])
    main_ids = pd.Series(["SPAC1002.01c", "SPAC823.17", "SPCC622.17"])
    id_map = {"SPAC823.02": "SPAC823.17", "SPAC1002.01C": "SPAC1002.01c"}
    name_map = {"SPAC823.17": "tom6", "SPAC1002.01c": "mrx11"}

    ids, names = build_paralog_columns(raw, main_ids, id_map, name_map)

    assert ids.iloc[0] == "SPAC823.17"
    assert ids.iloc[1] == "being curated"
    assert pd.isna(ids.iloc[2])
    assert names.iloc[0] == "tom6"
    assert pd.isna(names.iloc[1])
    assert pd.isna(names.iloc[2])


def test_combine_mapping_notes_single_and_collapsed():
    notes = {"SPAC27D7.09c": "NO_CHANGE", "SPAC27D7.10c": "VIA_SYNONYM → SPAC27D7.09c"}

    assert combine_mapping_notes("SPAC27D7.09c", notes) == "NO_CHANGE"
    assert combine_mapping_notes(
        "SPAC27D7.09c; SPAC27D7.10c", notes
    ).startswith("COLLAPSED: SPAC27D7.09c: NO_CHANGE")


def test_read_pombase_annotation_strips_whitespace(tmp_path: Path):
    path = tmp_path / "annotation.tsv"
    path.write_text(
        "gene_systematic_id\tgene_name\tgene_product\tgene_type\tsynonyms\n"
        " SPAC1002.01c \t mrx11 \t product A \t protein coding gene \t alias1\n"
    )

    df = read_pombase_annotation(path)

    assert df.loc[0, PomBaseCol.GENE_SYSID] == "SPAC1002.01c"
    assert df.loc[0, PomBaseCol.GENE_NAME] == "mrx11"
    assert df.loc[0, PomBaseCol.SYNONYMS] == "alias1"
