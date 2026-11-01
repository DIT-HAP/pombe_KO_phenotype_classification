"""
Hayles 2013 Phenotype Table
===========================

Builds the gene-level phenotype table from Hayles 2013. The Kim 2010 table
contributes only its updated systematic IDs and processed Paralog list.
Kim phenotype data are not merged. The Hayles table contributes its
phenotype description, morphological classification, and dispensability
call.

Steps
-----
1. Load Kim's Systemic ID and Paralog columns.
2. Update every main and Paralog ID through the current PomBase annotation.
   Convert each Paralog ID to its current gene name; keep the ID when the
   annotation has no gene name.
3. Write the Kim paralog pre-merge table with mapping notes.
4. Load Hayles' four columns and rename them to clear names.
5. Clean the Hayles phenotype text (whitespace, trailing period, temperature
   notation).
6. Update the Hayles IDs and write the Hayles pre-merge table with notes.
7. Left-join the Kim paralog pre-merge onto the Hayles pre-merge by the
   updated Systematic ID.
8. Add current PomBase Gene name, Gene product, and Gene type. Fill an
   empty Gene name with the updated Systematic ID.

Input
-----
- Kim 2010 full Supplementary Table 1 (.xls): Systemic ID, Paralog.
- Hayles 2013 main table (.xlsx): Systematic ID, Deletion mutant phenotype
  description, Phenotypic classification used for analysis,
  Gene dispensability. This study.
- PomBase annotation TSV.

Output
------
- data/1_formatted/Kim2010_paralog_premerge.xlsx:
  Systematic ID, Kim_ID, Kim_Note, Kim_Paralog, Kim_Paralog_gene_name
  (4,913 rows).
- data/1_formatted/Hayles2013_phenotypes_premerge.xlsx:
  Systematic ID, Hayles_ID, Hayles_Note, Phenotype_description,
  Morphological_classification, Gene_dispensability (4,843 rows).
- data/1_formatted/Hayles2013_phenotypes_annotated.xlsx:
  Hayles pre-merge plus Kim paralog columns and Gene name, Gene product,
  Gene type (4,843 rows).

Usage
-----
    cd <project root>
    mamba run -n bioinformatics python scripts/01_format_hayles_phenotypes.py
    mamba run -n bioinformatics python scripts/01_format_hayles_phenotypes.py --verbose

Author:   Yusheng Yang (guidance) + Claude Code (implementation)
Date:     2026-10-29
Version:  2.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
import sys

# 2. Third-party Imports
from loguru import logger

# 3. Data Processing Imports
import pandas as pd

# 4. Project Path Setup
SCRIPT_DIR = Path(__file__).resolve().parent
SRC_DIR = SCRIPT_DIR.parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# 5. Project Imports
from dataset_merge import (  # noqa: E402
    collapse_duplicate_ids,
    load_hayles2013,
    load_kim2010,
)
from id_mapping import (  # noqa: E402
    PomBaseCol,
    build_gene_name_map,
    build_paralog_columns,
    collect_mapping_ids,
    combine_mapping_notes,
    read_pombase_annotation,
    update_sysids,
)
from phenotype import (  # noqa: E402
    DISPENSABILITY_COL,
    MORPHOLOGY_COL,
    PHENOTYPE_DESC_COL,
    RAW_DISPENSABILITY_COL,
    RAW_MORPHOLOGY_COL,
    RAW_PHENOTYPE_DESC_COL,
    format_hayles_phenotypes,
)
from pipeline_utils import setup_logger  # noqa: E402

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================


class PhenotypeCol(StrEnum):
    """Column names in the pre-merge and final tables."""

    SYSTEMATIC_ID = "Systematic ID"
    GENE_NAME = "Gene name"
    GENE_PRODUCT = "Gene product"
    GENE_TYPE = "Gene type"
    KIM_ID = "Kim_ID"
    HAYLES_ID = "Hayles_ID"
    KIM_NOTE = "Kim_Note"
    HAYLES_NOTE = "Hayles_Note"
    KIM_PARALOG = "Kim_Paralog"
    KIM_PARALOG_GENE_NAME = "Kim_Paralog_gene_name"
    PHENOTYPE_DESCRIPTION = PHENOTYPE_DESC_COL
    MORPHOLOGICAL_CLASSIFICATION = MORPHOLOGY_COL
    GENE_DISPENSABILITY = DISPENSABILITY_COL


PROJECT_ROOT = SCRIPT_DIR.parent

DEFAULT_KIM_PATH = PROJECT_ROOT / "data/references/kim2010/kim2010_supp_table1.xls"
DEFAULT_HAYLES_PATH = PROJECT_ROOT / "data/references/hayles2013/hayles2013_supp_table3.xlsx"
DEFAULT_ANNOTATION_PATH = (
    PROJECT_ROOT / "data/references/pombase-2026-06-01_gene_IDs_names_products.tsv"
)

DEFAULT_KIM_PREMERGE_OUTPUT = PROJECT_ROOT / "data/1_formatted/Kim2010_paralog_premerge.xlsx"
DEFAULT_HAYLES_PREMERGE_OUTPUT = (
    PROJECT_ROOT / "data/1_formatted/Hayles2013_phenotypes_premerge.xlsx"
)
DEFAULT_OUTPUT = PROJECT_ROOT / "data/1_formatted/Hayles2013_phenotypes_annotated.xlsx"

KIM_PREMERGE_COLUMNS: tuple[PhenotypeCol, ...] = (
    PhenotypeCol.SYSTEMATIC_ID,
    PhenotypeCol.KIM_ID,
    PhenotypeCol.KIM_NOTE,
    PhenotypeCol.KIM_PARALOG,
    PhenotypeCol.KIM_PARALOG_GENE_NAME,
)

HAYLES_PREMERGE_COLUMNS: tuple[PhenotypeCol, ...] = (
    PhenotypeCol.SYSTEMATIC_ID,
    PhenotypeCol.HAYLES_ID,
    PhenotypeCol.HAYLES_NOTE,
    PhenotypeCol.PHENOTYPE_DESCRIPTION,
    PhenotypeCol.MORPHOLOGICAL_CLASSIFICATION,
    PhenotypeCol.GENE_DISPENSABILITY,
)

FINAL_COLUMNS: tuple[PhenotypeCol, ...] = (
    PhenotypeCol.SYSTEMATIC_ID,
    PhenotypeCol.GENE_NAME,
    PhenotypeCol.GENE_PRODUCT,
    PhenotypeCol.GENE_TYPE,
    PhenotypeCol.HAYLES_ID,
    PhenotypeCol.KIM_ID,
    PhenotypeCol.HAYLES_NOTE,
    PhenotypeCol.KIM_NOTE,
    PhenotypeCol.PHENOTYPE_DESCRIPTION,
    PhenotypeCol.MORPHOLOGICAL_CLASSIFICATION,
    PhenotypeCol.GENE_DISPENSABILITY,
    PhenotypeCol.KIM_PARALOG,
    PhenotypeCol.KIM_PARALOG_GENE_NAME,
)

# =============================================================================
# CONFIGURATION & DATACLASSES
# =============================================================================


@dataclass(kw_only=True, frozen=True)
class MergeConfig:
    """Runtime paths for the Hayles 2013 phenotype table."""

    kim_path: Path
    hayles_path: Path
    annotation_path: Path
    kim_premerge_path: Path
    hayles_premerge_path: Path
    output_path: Path


@dataclass(kw_only=True, frozen=True)
class MergeResult:
    """The two pre-merge tables and the merged output."""

    kim_premerge: pd.DataFrame
    hayles_premerge: pd.DataFrame
    merged: pd.DataFrame


# =============================================================================
# CORE LOGIC
# =============================================================================


def load_sources(
    config: MergeConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load the two source tables and the PomBase annotation."""
    kim = load_kim2010(config.kim_path)
    hayles = format_hayles_phenotypes(load_hayles2013(config.hayles_path))
    annotation = read_pombase_annotation(config.annotation_path)
    return kim, hayles, annotation


def map_systematic_ids(
    kim: pd.DataFrame,
    hayles: pd.DataFrame,
    annotation: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str], dict[str, str]]:
    """Update both ID columns and keep the pre-update IDs and notes."""
    mapping_ids = collect_mapping_ids(kim, hayles)
    mapped_ids, notes = update_sysids(mapping_ids, annotation)
    id_map = dict(zip(mapping_ids, mapped_ids, strict=True))
    note_map = dict(zip(mapping_ids, notes, strict=True))

    kim = kim.copy()
    hayles = hayles.copy()
    kim[PhenotypeCol.KIM_ID] = kim["id"]
    hayles[PhenotypeCol.HAYLES_ID] = hayles["id"]
    kim["id"] = kim["id"].map(id_map)
    hayles["id"] = hayles["id"].map(id_map)
    return kim, hayles, id_map, note_map


def prepare_kim(
    kim: pd.DataFrame,
    id_map: dict[str, str],
    note_map: dict[str, str],
    annotation: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse duplicate Kim IDs, attach notes, and process Paralog lists."""
    kim = collapse_duplicate_ids(
        kim,
        id_col="id",
        original_id_col=PhenotypeCol.KIM_ID,
    )
    kim[PhenotypeCol.KIM_NOTE] = kim[PhenotypeCol.KIM_ID].map(
        lambda original_ids: combine_mapping_notes(original_ids, note_map)
    )
    paralog_ids, paralog_names = build_paralog_columns(
        kim["paralogues_2010"],
        kim["id"],
        id_map,
        name_map=build_gene_name_map(annotation),
    )
    kim[PhenotypeCol.KIM_PARALOG] = paralog_ids
    kim[PhenotypeCol.KIM_PARALOG_GENE_NAME] = paralog_names
    return kim


def prepare_hayles(
    hayles: pd.DataFrame,
    note_map: dict[str, str],
) -> pd.DataFrame:
    """Collapse duplicate Hayles IDs and attach ID mapping notes."""
    hayles = collapse_duplicate_ids(
        hayles,
        id_col="id",
        original_id_col=PhenotypeCol.HAYLES_ID,
    )
    hayles[PhenotypeCol.HAYLES_NOTE] = hayles[PhenotypeCol.HAYLES_ID].map(
        lambda original_ids: combine_mapping_notes(original_ids, note_map)
    )
    return hayles


def build_kim_premerge(kim: pd.DataFrame) -> pd.DataFrame:
    """Kim's updated ID, original ID, note, and processed Paralog columns."""
    return kim[
        [
            "id",
            PhenotypeCol.KIM_ID,
            PhenotypeCol.KIM_NOTE,
            PhenotypeCol.KIM_PARALOG,
            PhenotypeCol.KIM_PARALOG_GENE_NAME,
        ]
    ].rename(columns={"id": PhenotypeCol.SYSTEMATIC_ID})[list(KIM_PREMERGE_COLUMNS)]


def build_hayles_premerge(hayles: pd.DataFrame) -> pd.DataFrame:
    """Hayles' updated ID, original ID, note, and four renamed columns."""
    return hayles[
        [
            "id",
            PhenotypeCol.HAYLES_ID,
            PhenotypeCol.HAYLES_NOTE,
            RAW_PHENOTYPE_DESC_COL,
            RAW_MORPHOLOGY_COL,
            RAW_DISPENSABILITY_COL,
        ]
    ].rename(
        columns={
            "id": PhenotypeCol.SYSTEMATIC_ID,
            RAW_PHENOTYPE_DESC_COL: PhenotypeCol.PHENOTYPE_DESCRIPTION,
            RAW_MORPHOLOGY_COL: PhenotypeCol.MORPHOLOGICAL_CLASSIFICATION,
            RAW_DISPENSABILITY_COL: PhenotypeCol.GENE_DISPENSABILITY,
        }
    )[list(HAYLES_PREMERGE_COLUMNS)]


def merge_tables(
    hayles_premerge: pd.DataFrame,
    kim_premerge: pd.DataFrame,
) -> pd.DataFrame:
    """Left-join the Kim paralog pre-merge onto the Hayles pre-merge by updated ID."""
    merged = hayles_premerge.merge(
        kim_premerge,
        on=PhenotypeCol.SYSTEMATIC_ID,
        how="left",
        sort=False,
    )
    missing = int(merged[PhenotypeCol.KIM_ID].isna().sum())
    if missing:
        logger.warning(f"{missing:,} Hayles rows have no Kim ID match")
    return merged.sort_values(PhenotypeCol.SYSTEMATIC_ID).reset_index(drop=True)


def add_annotation(merged: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    """Add current PomBase gene name, product, and type by updated systematic ID.

    A direct systematic ID wins; a Hayles ID that is only a current PomBase
    synonym (for example ``SPBC8E4.02c``) resolves through the synonym table.
    """
    by_id = annotation.drop_duplicates(PomBaseCol.GENE_SYSID).set_index(PomBaseCol.GENE_SYSID)
    synonyms = (
        annotation[[PomBaseCol.GENE_SYSID, PomBaseCol.SYNONYMS]]
        .dropna(subset=[PomBaseCol.SYNONYMS])
        .assign(synonym=lambda df: df[PomBaseCol.SYNONYMS].str.split(","))
        .explode("synonym")
    )
    synonym_to_id = (
        synonyms.assign(synonym=lambda df: df["synonym"].str.strip())
        .loc[lambda df: df["synonym"].ne("")]
        .drop_duplicates("synonym")
        .set_index("synonym")[PomBaseCol.GENE_SYSID]
    )

    ids = merged[PhenotypeCol.SYSTEMATIC_ID]
    canonical = ids.where(ids.isin(by_id.index), ids.map(synonym_to_id))
    gene_name = canonical.map(by_id[PomBaseCol.GENE_NAME]).astype("string").str.strip()
    gene_name = gene_name.mask(gene_name.isna() | gene_name.eq(""), ids)

    annotated = merged.copy()
    annotated[PhenotypeCol.GENE_NAME] = gene_name
    annotated[PhenotypeCol.GENE_PRODUCT] = canonical.map(by_id[PomBaseCol.GENE_PRODUCT])
    annotated[PhenotypeCol.GENE_TYPE] = canonical.map(by_id[PomBaseCol.GENE_TYPE])
    return annotated


def build_phenotype_table(config: MergeConfig) -> MergeResult:
    """Load, map, prepare both sources, merge, and annotate."""
    kim, hayles, annotation = load_sources(config)
    kim, hayles, id_map, note_map = map_systematic_ids(kim, hayles, annotation)
    kim = prepare_kim(kim, id_map, note_map, annotation)
    hayles = prepare_hayles(hayles, note_map)

    kim_premerge = build_kim_premerge(kim)
    hayles_premerge = build_hayles_premerge(hayles)
    merged = merge_tables(hayles_premerge, kim_premerge)
    merged = add_annotation(merged, annotation)
    merged = merged[list(FINAL_COLUMNS)]

    logger.info(f"Kim paralog pre-merge: {len(kim_premerge):,} rows")
    logger.info(f"Hayles pre-merge: {len(hayles_premerge):,} rows")
    logger.info(f"Hayles phenotype table: {len(merged):,} rows")
    logger.info(f"Hayles rows without Kim ID: {merged[PhenotypeCol.KIM_ID].isna().sum():,}")
    return MergeResult(
        kim_premerge=kim_premerge,
        hayles_premerge=hayles_premerge,
        merged=merged,
    )


# =============================================================================
# MAIN EXECUTION
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build the Hayles 2013 phenotype table with Kim Paralog.",
    )
    parser.add_argument(
        "--kim2010",
        type=Path,
        default=DEFAULT_KIM_PATH,
        help=f"Kim 2010 full Table S1 (default: {DEFAULT_KIM_PATH})",
    )
    parser.add_argument(
        "--hayles2013",
        type=Path,
        default=DEFAULT_HAYLES_PATH,
        help=f"Hayles 2013 main table (default: {DEFAULT_HAYLES_PATH})",
    )
    parser.add_argument(
        "--annotation",
        type=Path,
        default=DEFAULT_ANNOTATION_PATH,
        help=(f"PomBase annotation TSV (default: {DEFAULT_ANNOTATION_PATH})"),
    )
    parser.add_argument(
        "--kim-premerge-output",
        type=Path,
        default=DEFAULT_KIM_PREMERGE_OUTPUT,
        help=(f"Kim paralog pre-merge workbook (default: {DEFAULT_KIM_PREMERGE_OUTPUT})"),
    )
    parser.add_argument(
        "--hayles-premerge-output",
        type=Path,
        default=DEFAULT_HAYLES_PREMERGE_OUTPUT,
        help=(f"Hayles pre-merge workbook (default: {DEFAULT_HAYLES_PREMERGE_OUTPUT})"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Merged output workbook (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    return parser.parse_args()


def main() -> int:
    """Write the two pre-merge tables and the annotated phenotype table."""
    args = parse_args()
    setup_logger("DEBUG" if args.verbose else "INFO")
    config = MergeConfig(
        kim_path=args.kim2010,
        hayles_path=args.hayles2013,
        annotation_path=args.annotation,
        kim_premerge_path=args.kim_premerge_output,
        hayles_premerge_path=args.hayles_premerge_output,
        output_path=args.output,
    )

    try:
        result = build_phenotype_table(config)
        outputs = (
            (config.kim_premerge_path, result.kim_premerge),
            (config.hayles_premerge_path, result.hayles_premerge),
            (config.output_path, result.merged),
        )
        for path, frame in outputs:
            path.parent.mkdir(parents=True, exist_ok=True)
            frame.to_excel(path, index=False)
            logger.success(f"Saved: {path} ({len(frame):,} rows)")
    except Exception as exc:
        logger.exception(f"Phenotype table build failed: {exc}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
