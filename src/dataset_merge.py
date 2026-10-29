"""
Shared Kim 2010 / Hayles 2013 dataset merge
===========================================

Library module used by ``01_format_hayles_phenotypes.py`` for the source
loaders and duplicate-ID collapse, and by the archived
``01_format_and_update_ids.py`` for the union builder and Table S3 parsing.
It loads the two source tables, canonicalises their column names, and builds
a union table with:

- one row per current systematic ID (the caller updates IDs before calling
  this module and uses ``collapse_duplicate_ids()`` for any collisions);
- a ``Presence`` column: ``both``, ``2010_only``, ``2013_only`` or
  ``neither`` (data availability, not just ID membership);
- a ``Dispensability_consistency`` column for genes present in both:
  ``consistent``, ``inconsistent`` or ``not_comparable``;
- a human-readable ``Dispensability_discrepancy`` column describing every
  inconsistent / non-comparable case;
- a machine-readable ``Dispensability_change_class`` column using the
  taxonomy from the 2013 Table S3 reconciliation;
- Table S3 per-gene annotations for the 21 documented changes:
  ``Hayles2013_Strain_status`` (deletion remade / re-analysed),
  ``Hayles2013_2010_dispensability_change`` (Kim 2010 analysis → 2013 study),
  and the verbatim ``Hayles2013_TableS3_comment``.

Presence is defined by **data availability**:

- ``both``: Kim 2010 has an ``Analysis dataset`` V/E call *and* the gene is
  in the 2013 screen (4,825 genes).
- ``2010_only``: Kim 2010 has a call but the gene is absent from 2013 (11
  constructed strains).
- ``2013_only``: the gene is in the 2013 screen but Kim 2010 has no
  ``Analysis dataset`` call (18 Hayles additions).
- ``neither``: Kim lists the gene as a candidate but never constructed it,
  and it is absent from 2013 (60 genes, no dispensability data in either).

The union therefore has 4,914 rows: 4,825 + 11 + 18 + 60.

Input
-----
- Kim 2010 full Supplementary Table 1 ``.xls``.
- Kim 2010 Supplementary Table 3 ``.tsv`` (hand-transcribed: the 40
  published-vs-own disagreements and the paper's judgment criteria).
- Hayles 2013 main ``.xlsx``.
- Optional Hayles Supplementary Material 1 ``.docx`` for Table S3.

Output
------
- Canonical DataFrames for the two sources, plus a merged union table with
  the columns described above.

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-10-28
Version: 1.3.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import re
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime
from pathlib import Path

# 2. Data Processing Imports
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

DEFAULT_KIM2010 = Path("data/references/kim2010/kim2010_supp_table1.xls")
DEFAULT_KIM_TABLE_S3 = Path("data/references/kim2010/kim2010_supp_table3.tsv")
DEFAULT_HAYLES2013 = Path(
    "data/references/hayles2013/hayles2013_supp_table3.xlsx"
)
DEFAULT_HAYLES_DOCX = Path("data/references/hayles2013/hayles2013_supp.docx")

KIM_ID = "Systemic ID"  # sic: the published table spells it "Systemic"
KIM_DELETION = "Deletion dataset"
KIM_ANALYSIS = "Analysis dataset"
HAYLES_ID = "Systematic ID"
HAYLES_STUDY = "Gene dispensability. This study"
HAYLES_STUDY_COL = "Hayles2013_Gene_dispensability"
HAYLES_PHENOTYPE = "Phenotypic classification used for analysis"
S3_STRAIN_STATUS_COL = "Hayles2013_Strain_status"
S3_CHANGE_COL = "Hayles2013_2010_dispensability_change"
S3_COMMENT_COL = "Hayles2013_TableS3_comment"

# WordprocessingML namespace for the supplementary .docx tables.
WORD_NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}

# Table S3 of the 2013 supplement contains one truncated locus ID.
TABLE_S3_ID_REPAIRS = {"SPAC13D6.0": "SPAC13D6.05"}

# Canonical column names in the merged table.
PRESENCE_COL = "Presence"
CONSISTENCY_COL = "Dispensability_consistency"
DISCREPANCY_COL = "Dispensability_discrepancy"
CHANGE_CLASS_COL = "Dispensability_change_class"

# Kim 2010 Supplementary Table 3 (published-vs-own disagreements).
KIM_S3_COLUMNS = (
    "table_s3_no",
    "systematic_id",
    "gene_name",
    "published_data",
    "our_data",
    "criteria_for_judgment",
    "criteria_detail",
    "final_judgment",
)
KIM_S3_NO_COL = "Kim2010_TableS3_no"
KIM_S3_CRITERIA_COL = "Kim2010_Judgment_criteria"
KIM_S3_DETAIL_COL = "Kim2010_Judgment_criteria_detail"

# Expected shape of the hand-transcribed Kim Table S3, used for validation.
KIM_S3_EXPECTED_ROWS = 40
KIM_S3_CRITERIA_COUNTS = {
    "Disagreement": 19,
    "More of CDS deleted in this study": 13,
    "Essential in our assay conditions": 4,
    "Different interpretation due to slow growth or dying cells": 2,
    "Published data not confirmed": 2,
}
KIM_S3_DETAIL_IDS = {
    "SPAC15A10.01",  # atm1: KO%-dependent
    "SPAC1B3.14",    # vma3: media pH
    "SPAC343.05",    # vma1: media pH
    "SPAC56F8.10",   # met9: methionine
    "SPCC1753.04",   # tol1: sulfite
}

# Source-column renames shared by the union and 01's per-source outputs.
KIM_COLUMN_RENAMES = {
    "gene_name_2010": "Kim2010_Gene_name",
    "pct_ko_2010": "Kim2010_KO_pct",
    "kim2010_deletion": "Kim2010_Deletion_dataset",
    "kim2010_analysis": "Kim2010_Analysis_dataset",
    "published_phenotype": "Kim2010_Published_phenotype",
    "constructed": "Kim2010_Constructed",
    "construction_method": "Kim2010_Construction_method",
    "gene_product_2010": "Kim2010_Gene_function",
    "paralogues_2010": "Kim2010_Paralog",
    "sc_orthologue_2010": "Kim2010_Sc_ortholog",
    "ts_remake": "Kim2010_ts_remake",
}

HAYLES_COLUMN_RENAMES = {
    "gene_name_2013": "Hayles2013_Gene_name",
    "working_code": "Working code",
    "pct_ko_2013": "% KO",
    "hayles2013_study": HAYLES_STUDY_COL,
    "phenotype_description": "Deletion mutant phenotype description",
    "phenotype_classification": HAYLES_PHENOTYPE,
    "gene_product_2013": "Gene product description",
    "paralogues_2013": "Paralogues",
    "cell_length": "Cell length >15.6 (st dev) for 34 novel long genes",
}


# =============================================================================
# VALUE NORMALIZATION
# =============================================================================


def normalize_id(value: object) -> str | None:
    """Strip all whitespace from a systematic ID and drop missing values."""
    if value is None or value is pd.NA:
        return None
    text = re.sub(r"\s+", "", str(value))
    if not text or text.lower() == "nan":
        return None
    return text


def normalize_call(value: object) -> str | None:
    """Map viable/essential (and stray-spaced V/E) to ``V`` or ``E``."""
    if value is None or value is pd.NA:
        return None
    text = str(value).strip().upper()
    if not text or text == "NAN":
        return None
    if text in {"V", "VIABLE"}:
        return "V"
    if text in {"E", "ESSENTIAL"}:
        return "E"
    logger.warning(f"Unrecognised dispensability call {value!r}; treating as missing")
    return None


def normalize_published_call(value: object) -> str | None:
    """Map Kim's published phenotype text to ``V`` / ``E``, silently ignoring ``nodata``."""
    if value is None or value is pd.NA:
        return None
    text = str(value).strip().upper()
    if text in {"V", "VIABLE"}:
        return "V"
    if text in {"E", "ESSENTIAL"}:
        return "E"
    return None


def clean_text_column(series: pd.Series) -> pd.Series:
    """Strip a text column, turning empty strings into missing values."""
    return series.astype("string").str.strip().replace("", pd.NA)


def optional_column(df: pd.DataFrame, name: str) -> pd.Series:
    """Return a column, or an all-missing column of the same length."""
    if name in df.columns:
        return df[name]
    return pd.Series(pd.NA, index=df.index, dtype="string")


def reconstruct_gene_name_from_date(value: datetime) -> str:
    """Rebuild an Excel auto-converted gene name such as ``Oct-1`` from its date.

    Excel reads ``Oct-1`` as 1 Oct of the current year, but reads ``Sep-11``
    as 1 Sep 2011 (the 11 becomes a year). When the day is 1 and the year is
    after 2001, the intended day is therefore ``year - 2000``.
    """
    day = value.year - 2000 if value.day == 1 and value.year > 2001 else value.day
    return f"{value.strftime('%b')}-{day}"


# =============================================================================
# LOADERS
# =============================================================================


@logger.catch(reraise=True)
def load_kim2010(path: Path) -> pd.DataFrame:
    """Load Kim 2010 full Supplementary Table 1 into canonical columns.

    Every constructed gene has a ``kim2010_analysis`` call; the 78
    ``Constructed? = no`` rows mostly have no call (two have a deletion-only
    partial call).
    """
    df = pd.read_excel(path, sheet_name=0)
    df.columns = [str(c).strip() for c in df.columns]

    if "Constructed?" not in df.columns:
        raise ValueError(
            f"{path}: this is not the full Kim 2010 Supplementary Table 1 "
            "(missing 'Constructed?'). The 1:1-ortholog Table 11 subset is "
            "not sufficient for the complete reconciliation."
        )
    return _load_kim_full(df, path)


def _load_kim_full(df: pd.DataFrame, path: Path) -> pd.DataFrame:
    """Canonicalize Kim 2010 full Supplementary Table 1."""
    required = [KIM_ID, "Constructed?", KIM_DELETION, KIM_ANALYSIS]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    out = pd.DataFrame(
        {
            "id": df[KIM_ID].map(normalize_id),
            "gene_name_2010": clean_text_column(df["Gene name"]),
            "pct_ko_2010": pd.to_numeric(df["KO%"], errors="coerce"),
            "kim2010_deletion": df[KIM_DELETION].map(normalize_call),
            "kim2010_analysis": df[KIM_ANALYSIS].map(normalize_call),
            "published_phenotype": clean_text_column(df["Published phenotype"]),
            "constructed": clean_text_column(df["Constructed?"])
            .str.lower()
            .eq("yes"),
            "construction_method": clean_text_column(
                optional_column(df, "Construction method")
            ),
            "gene_product_2010": clean_text_column(
                optional_column(df, "Gene function")
            ),
            "paralogues_2010": clean_text_column(optional_column(df, "Paralog")),
            "sc_orthologue_2010": clean_text_column(
                optional_column(df, "Sc ortholog")
            ),
            "ts_remake": clean_text_column(
                optional_column(
                    df,
                    "Heterozygous diploid containing ts mutation (being remade)",
                )
            ),
            "source_format": "full_table_s1",
        }
    )
    if out["id"].isna().any():
        raise ValueError(f"{path}: {out['id'].isna().sum()} rows have no systematic ID")
    if out["id"].duplicated().any():
        dup = out.loc[out["id"].duplicated(), "id"].tolist()
        raise ValueError(f"{path}: duplicate systematic IDs {dup[:10]}")

    out["kim2010_call"] = out["kim2010_analysis"].where(
        out["kim2010_analysis"].notna(), out["kim2010_deletion"]
    )
    logger.info(
        f"Kim 2010 full Table S1: {len(out):,} rows; "
        f"{int(out['constructed'].sum()):,} constructed"
    )
    return out


@logger.catch(reraise=True)
def load_hayles2013(path: Path) -> pd.DataFrame:
    """Load the Hayles 2013 main supplementary table."""
    df = pd.read_excel(path, sheet_name=0)
    df.columns = [str(c).strip() for c in df.columns]
    required = [
        HAYLES_ID,
        HAYLES_STUDY,
        HAYLES_PHENOTYPE,
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    # Excel auto-converted three gene names into dates (Oct-1, Sep-11, Sep-1);
    # reconstruct them before the column is cast to strings.
    gene_names = df["Gene name"].map(
        lambda value: reconstruct_gene_name_from_date(value)
        if isinstance(value, datetime)
        else value
    )

    out = pd.DataFrame(
        {
            "id": df[HAYLES_ID].map(normalize_id),
            "gene_name_2013": clean_text_column(gene_names),
            "working_code": clean_text_column(df["Working code"]),
            "pct_ko_2013": pd.to_numeric(df["% KO"], errors="coerce"),
            "hayles2013_study": df[HAYLES_STUDY].map(normalize_call),
            "phenotype_description": clean_text_column(
                df["Deletion mutant phenotype description"]
            ),
            "phenotype_classification": clean_text_column(df[HAYLES_PHENOTYPE]),
            "gene_product_2013": clean_text_column(df["Gene product description"]),
            "paralogues_2013": clean_text_column(df["Paralogues"]),
            "cell_length": clean_text_column(
                optional_column(
                    df, "Cell length >15.6 (st dev) for 34 novel long genes"
                )
            ),
        }
    )
    if out["id"].isna().any():
        raise ValueError(f"{path}: {out['id'].isna().sum()} rows have no systematic ID")
    if out["id"].duplicated().any():
        dup = out.loc[out["id"].duplicated(), "id"].tolist()
        raise ValueError(f"{path}: duplicate systematic IDs {dup[:10]}")
    logger.info(f"Hayles 2013: {len(out):,} genes")
    return out


# =============================================================================
# SUPPLEMENTARY DOCX TABLES (S3)
# =============================================================================


@logger.catch(reraise=True)
def read_docx_tables(path: Path) -> list[list[list[str]]]:
    """Extract every Word table as a list of rows of cell strings."""
    with zipfile.ZipFile(path) as archive:
        document = ET.fromstring(archive.read("word/document.xml"))
    tables: list[list[list[str]]] = []
    for table in document.findall(".//w:tbl", WORD_NS):
        rows: list[list[str]] = []
        for row in table.findall("./w:tr", WORD_NS):
            cells = []
            for cell in row.findall("./w:tc", WORD_NS):
                text = "".join(
                    node.text or "" for node in cell.findall(".//w:t", WORD_NS)
                )
                cells.append(text.strip())
            rows.append(cells)
        tables.append(rows)
    return tables


def _cell(row: list[str], index: int) -> str:
    """Return a row cell or an empty string when the row is short."""
    return row[index] if index < len(row) else ""


def _strain_status(comment: str | None) -> str | None:
    """Classify a Table S3 comment as a remade or re-analysed strain."""
    if not comment:
        return None
    if comment.startswith("Deletion remade"):
        return "Deletion remade"
    if comment.startswith("Re-analysed"):
        return "Re-analysed"
    return None


def _dispensability_change(before: str | None, after: str | None) -> str | None:
    """Summarise the Table S3 change from Kim 2010 to the 2013 study."""
    if before is None or after is None:
        return None
    if before == after:
        return f"same ({before})"
    return f"{before} -> {after}"


@logger.catch(reraise=True)
def parse_table_s3(tables: list[list[list[str]]]) -> pd.DataFrame:
    """Parse Hayles Table S3 (documented dispensability changes)."""
    table = next(
        (
            t
            for t in tables
            if t
            and "Comments" in _cell(t[0], -1)
            and any("Dispensability" in cell for cell in t[0])
        ),
        None,
    )
    if table is None:
        raise ValueError("Table S3 not found in the supplementary docx")

    rows = []
    for cells in table[1:]:
        gene_id = normalize_id(_cell(cells, 0))
        if gene_id is None:
            continue
        if gene_id in TABLE_S3_ID_REPAIRS:
            repaired = TABLE_S3_ID_REPAIRS[gene_id]
            logger.warning(f"Table S3: repairing truncated ID {gene_id} -> {repaired}")
            gene_id = repaired
        comment = _cell(cells, 5) or None
        analysis = normalize_call(_cell(cells, 3))
        study = normalize_call(_cell(cells, 4))
        rows.append(
            {
                "id": gene_id,
                "gene_name_s3": _cell(cells, 1) or None,
                "kim2010_deletion_s3": normalize_call(_cell(cells, 2)),
                "kim2010_analysis_s3": analysis,
                "hayles2013_study_s3": study,
                "s3_comment": comment,
                "strain_status": _strain_status(comment),
                "dispensability_change": _dispensability_change(analysis, study),
            }
        )
    logger.info(f"Table S3: {len(rows)} documented changes parsed")
    return pd.DataFrame(rows)


# =============================================================================
# KIM 2010 TABLE S3 (LITERATURE DISAGREEMENTS)
# =============================================================================


@logger.catch(reraise=True)
def load_kim_table_s3(path: Path = DEFAULT_KIM_TABLE_S3) -> pd.DataFrame:
    """Load the hand-transcribed Kim 2010 Supplementary Table 3."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    missing = [col for col in KIM_S3_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    out = df[list(KIM_S3_COLUMNS)].copy()
    for col in out.columns:
        out[col] = out[col].astype("string").str.strip().replace("", pd.NA)
    out["table_s3_no"] = pd.to_numeric(
        out["table_s3_no"], errors="raise"
    ).astype(int)
    out["systematic_id"] = out["systematic_id"].map(normalize_id)
    if out["systematic_id"].isna().any():
        raise ValueError(f"{path}: rows without a systematic ID")
    if out["systematic_id"].duplicated().any():
        duplicates = out.loc[
            out["systematic_id"].duplicated(), "systematic_id"
        ].tolist()
        raise ValueError(f"{path}: duplicate systematic IDs {duplicates}")

    for col in ("published_data", "our_data", "final_judgment"):
        out[col] = out[col].map(normalize_call)

    logger.info(f"Kim Table S3: {len(out)} published-vs-own disagreements")
    return out


@logger.catch(reraise=True)
def validate_kim_table_s3(s3: pd.DataFrame, kim2010: pd.DataFrame) -> None:
    """Check Table S3's row set, calls, numbering, criteria and details."""
    table = kim2010.set_index("id")

    expected_rows = list(range(1, KIM_S3_EXPECTED_ROWS + 1))
    if s3["table_s3_no"].tolist() != expected_rows:
        raise ValueError(
            f"Kim Table S3 numbering is not 1..{KIM_S3_EXPECTED_ROWS}: "
            f"{s3['table_s3_no'].tolist()}"
        )

    missing = [
        gene_id for gene_id in s3["systematic_id"] if gene_id not in table.index
    ]
    if missing:
        raise ValueError(f"Kim Table S3 IDs absent from Table S1: {missing}")

    expected = {
        "published_data": s3["systematic_id"].map(
            lambda gene_id: normalize_published_call(
                table.at[gene_id, "published_phenotype"]
            )
        ),
        "our_data": s3["systematic_id"].map(
            lambda gene_id: table.at[gene_id, "kim2010_deletion"]
        ),
        "final_judgment": s3["systematic_id"].map(
            lambda gene_id: table.at[gene_id, "kim2010_analysis"]
        ),
    }
    for column, values in expected.items():
        mismatches = s3.loc[s3[column] != values, ["systematic_id", column]]
        if not mismatches.empty:
            detail = mismatches.assign(expected=values.loc[mismatches.index])
            raise ValueError(
                f"Kim Table S3 {column} disagrees with Table S1:\n"
                f"{detail.to_string(index=False)}"
            )

    # The S3 rows must be exactly the Table S1 published-vs-own disagreements.
    published = table["published_phenotype"].map(normalize_published_call)
    own = table["kim2010_deletion"]
    disagreements = set(
        table.index[published.notna() & own.notna() & (published != own)]
    )
    listed = set(s3["systematic_id"])
    if disagreements != listed:
        raise ValueError(
            "Kim Table S3 does not match the Table S1 disagreement set; "
            f"missing={sorted(disagreements - listed)}, "
            f"extra={sorted(listed - disagreements)}"
        )

    criteria_counts = s3["criteria_for_judgment"].value_counts().to_dict()
    if criteria_counts != KIM_S3_CRITERIA_COUNTS:
        raise ValueError(
            "Kim Table S3 criteria distribution changed: "
            f"{criteria_counts} != {KIM_S3_CRITERIA_COUNTS}"
        )

    detail_ids = set(s3.loc[s3["criteria_detail"].notna(), "systematic_id"])
    if detail_ids != KIM_S3_DETAIL_IDS:
        raise ValueError(
            "Kim Table S3 footnote rows changed: "
            f"missing={sorted(KIM_S3_DETAIL_IDS - detail_ids)}, "
            f"extra={sorted(detail_ids - KIM_S3_DETAIL_IDS)}"
        )

    logger.info("Kim Table S3 row set, calls, criteria and details agree with Table S1")


# =============================================================================
# MERGE
# =============================================================================


def _merge_annotations(
    kim2010: pd.DataFrame,
    hayles2013: pd.DataFrame,
    table_s3: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build the presence/consistency/discrepancy annotations for the union."""
    k = kim2010.set_index("id")
    h = hayles2013.set_index("id")
    s3 = table_s3.set_index("id").to_dict("index") if table_s3 is not None else {}
    s3_available = table_s3 is not None

    # 2013 order first, then Kim-only rows in Kim's table order.
    h_ids = list(hayles2013["id"])
    h_set = set(h_ids)
    extra = [gene_id for gene_id in kim2010["id"] if gene_id not in h_set]

    rows = []
    for gene_id in h_ids + extra:
        in_2010 = gene_id in k.index
        in_2013 = gene_id in h.index
        deletion = k.at[gene_id, "kim2010_deletion"] if in_2010 else None
        analysis = k.at[gene_id, "kim2010_analysis"] if in_2010 else None
        constructed = k.at[gene_id, "constructed"] if in_2010 else None
        study = h.at[gene_id, "hayles2013_study"] if in_2013 else None
        s3_entry = s3.get(gene_id)
        has_2010 = pd.notna(analysis)

        if has_2010 and in_2013:
            presence = "both"
            if analysis == study:
                consistency = "consistent"
                if pd.notna(deletion) and deletion != analysis:
                    detail = (
                        f"same ({analysis}); Kim2010 deletion={deletion} "
                        f"vs analysis={analysis}"
                    )
                    if s3_entry is not None:
                        change_class = "documented_discrepancy_in_kim2010_columns"
                    else:
                        change_class = "kim2010_deletion_analysis_discrepancy"
                else:
                    detail = f"same ({analysis})"
                    change_class = "no_change"
            else:
                consistency = "inconsistent"
                if s3_entry is not None:
                    change_class = "documented_change"
                elif study == deletion:
                    change_class = "matches_kim2010_deletion_not_analysis"
                elif not s3_available:
                    change_class = "changed_table_s3_not_available"
                else:
                    change_class = "undocumented_change"
                detail = (
                    f"Kim2010_analysis={analysis} vs "
                    f"Hayles2013_study={study}; {change_class}"
                )
                if s3_entry and s3_entry.get("s3_comment"):
                    detail += f"; Table S3: {s3_entry['s3_comment']}"
        elif has_2010 and not in_2013:
            presence = "2010_only"
            consistency = "not_comparable"
            change_class = "no_hayles2013_row"
            detail = "no 2013 screen row"
        elif in_2013 and not has_2010:
            presence = "2013_only"
            consistency = "not_comparable"
            change_class = "hayles2013_addition"
            if pd.notna(deletion):
                detail = f"no 2010 analysis call (Kim2010 deletion-only={deletion})"
            elif constructed is False:
                detail = "no 2010 analysis call (Kim candidate not constructed)"
            else:
                detail = "no 2010 analysis call"
        else:
            presence = "neither"
            consistency = "not_comparable"
            change_class = "no_data_either"
            detail = "no 2010 analysis call and no 2013 row"

        rows.append(
            {
                "id": gene_id,
                PRESENCE_COL: presence,
                CONSISTENCY_COL: consistency,
                DISCREPANCY_COL: detail,
                CHANGE_CLASS_COL: change_class,
                S3_STRAIN_STATUS_COL: (
                    s3_entry.get("strain_status") if s3_entry else None
                ),
                S3_CHANGE_COL: (
                    s3_entry.get("dispensability_change") if s3_entry else None
                ),
                S3_COMMENT_COL: s3_entry.get("s3_comment") if s3_entry else None,
            }
        )
    return pd.DataFrame(rows)


def collapse_duplicate_ids(
    df: pd.DataFrame,
    id_col: str = "id",
    original_id_col: str = "original_id",
) -> pd.DataFrame:
    """Collapse rows sharing a current ID, combining their original IDs."""
    if original_id_col not in df.columns:
        df = df.assign(**{original_id_col: df[id_col]})
    if df[id_col].is_unique:
        return df

    rows = []
    for gene_id, group in df.groupby(id_col, sort=False):
        if len(group) == 1:
            rows.append(group.iloc[0])
            continue

        # Prefer the row whose original ID is already the current ID.
        primary_mask = (
            group[original_id_col].astype("string").str.upper()
            == str(gene_id).upper()
        )
        primary = (
            group.loc[primary_mask].iloc[0] if primary_mask.any() else group.iloc[0]
        )
        row = primary.copy()

        orig_ids = list(dict.fromkeys(str(value) for value in group[original_id_col]))
        row[original_id_col] = "; ".join(orig_ids)

        if "constructed" in group.columns:
            row["constructed"] = bool(group["constructed"].fillna(False).any())

        # Take a non-null gene name / ts-remake flag if the primary row lacks it.
        for col in ("gene_name_2010", "updated_gene_name_2010", "ts_remake"):
            if col in group.columns and pd.isna(row[col]):
                non_null = group[col].dropna()
                if not non_null.empty:
                    row[col] = non_null.iloc[0]

        for col in ("kim2010_deletion", "kim2010_analysis"):
            if col in group.columns:
                values = group[col].dropna().unique().tolist()
                if len(values) > 1:
                    logger.warning(f"Conflicting {col} for {gene_id}: {values}")

        logger.warning(f"Collapsed {len(group)} rows onto {gene_id}: {orig_ids}")
        rows.append(row)

    return pd.DataFrame(rows).reset_index(drop=True)


def build_merged_dataset(
    kim2010: pd.DataFrame,
    hayles2013: pd.DataFrame,
    table_s3: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build the union table used by 01 (2013 columns + prefixed 2010 columns)."""
    annotations = _merge_annotations(kim2010, hayles2013, table_s3)

    kim = kim2010.rename(columns=KIM_COLUMN_RENAMES).drop(
        columns=["kim2010_call", "source_format"],
        errors="ignore",
    )

    hayles = hayles2013.rename(columns=HAYLES_COLUMN_RENAMES)

    merged = annotations.merge(kim, on="id", how="left").merge(
        hayles, on="id", how="left"
    )
    merged = merged.rename(columns={"id": "Systematic ID"})

    def _first_available(df: pd.DataFrame, columns: list[str]) -> pd.Series:
        """First non-null value across the named columns."""
        result = pd.Series(pd.NA, index=df.index, dtype="string")
        for col in columns:
            if col in df.columns:
                result = result.where(result.notna(), df[col].astype("string"))
        return result

    # Current names are supplied pre-merge by 01; raw names are the fallback
    # when build_merged_dataset() is used directly on loader output.
    merged["Gene name"] = _first_available(
        merged,
        [
            "updated_gene_name_2013",
            "updated_gene_name_2010",
            "Hayles2013_Gene_name",
            "Kim2010_Gene_name",
        ],
    ).fillna(merged["Systematic ID"])

    # Original identifiers are the pre-mapping values carried by 01.
    merged["Original Gene name"] = _first_available(
        merged, ["Hayles2013_Gene_name", "Kim2010_Gene_name"]
    )

    def _combine_original_ids(row: pd.Series) -> str:
        """Combine the 2010/2013 pre-mapping IDs, dropping duplicates."""
        parts: list[str] = []
        for col in ("original_id_2010", "original_id_2013"):
            value = row.get(col)
            if pd.notna(value):
                for part in str(value).split(";"):
                    part = part.strip()
                    if part and part not in parts:
                        parts.append(part)
        return "; ".join(parts) if parts else str(row["Systematic ID"])

    merged["Original Systematic ID"] = merged.apply(_combine_original_ids, axis=1)

    merged = merged.drop(
        columns=[
            "original_id_2010",
            "original_id_2013",
            "updated_gene_name_2010",
            "updated_gene_name_2013",
        ],
        errors="ignore",
    )

    ordered = [
        "Systematic ID",
        "Gene name",
        PRESENCE_COL,
        CONSISTENCY_COL,
        DISCREPANCY_COL,
        CHANGE_CLASS_COL,
        S3_STRAIN_STATUS_COL,
        S3_CHANGE_COL,
        S3_COMMENT_COL,
        "Working code",
        "% KO",
        "Deletion mutant phenotype description",
        HAYLES_PHENOTYPE,
        HAYLES_STUDY_COL,
        "Gene product description",
        "Paralogues",
        "Cell length >15.6 (st dev) for 34 novel long genes",
        "Hayles2013_Gene_name",
        "Kim2010_Gene_name",
        "Kim2010_KO_pct",
        "Kim2010_Deletion_dataset",
        "Kim2010_Analysis_dataset",
        "Kim2010_Published_phenotype",
        "Kim2010_Constructed",
        "Kim2010_Construction_method",
        "Kim2010_Gene_function",
        "Kim2010_Paralog",
        "Kim2010_Sc_ortholog",
        "Kim2010_ts_remake",
        "Original Systematic ID",
        "Original Gene name",
    ]
    merged = merged[ordered]
    logger.info(
        f"Merged union: {len(merged):,} rows; "
        + ", ".join(
            f"{k}={v}" for k, v in merged[PRESENCE_COL].value_counts().items()
        )
    )
    return merged
