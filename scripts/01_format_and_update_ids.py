#!/usr/bin/env python3
"""
Format Hayles 2013 Phenotype Data and Update Gene Systematic IDs
=================================================================

Two-phase pipeline that (1) cleans and normalises the raw Hayles 2013
supplementary table, and (2) maps all gene systematic IDs against a
current PomBase annotation to correct identifiers changed by genome
annotation updates since 2013.

Phase 1 — Data Formatting
  - Strips whitespace from Systematic ID, phenotype description,
    classification, and dispensability columns
  - Normalises temperature notation: coerces variants such as
    "25, 32", "25 32", "32, 25" to the canonical "25,32"
  - Drops spurious columns (e.g. 'Unnamed: 12')
  - Prints summary statistics of phenotypic classification and
    gene dispensability categories

Phase 2 — Gene ID Mapping
  - Loads the PomBase gene annotation (gene_IDs_names_products.tsv)
  - Builds lookup tables from current systematic IDs, gene names,
    and synonyms (restricted to protein-coding genes by default)
  - Maps each input identifier through three layers:
      1. Direct systematic ID match
      2. Gene name → systematic ID lookup
      3. Synonym → systematic ID lookup
  - Identifies genes reclassified as pseudogene or ncRNA
  - Looks up current gene names for all entries

This script supersedes the original ``format_phenotype_input.py``
(which only performed Phase 1).

Input
-----
- ``data/raw/rsob130053supp2.xlsx``
  Hayles 2013 supplementary table (12 columns, 4,843 rows).
  Required columns: 'Systematic ID', 'Gene name'.

- PomBase annotation TSV (downloaded by 00_download_pombase_annotation.py)
  File: ``data/references/pombase-{release}_gene_IDs_names_products.tsv``

Output
------
- ``data/1_formatted/Hayles_2013_OB_formatted_phenotypes.xlsx``
  Cleaned phenotype data with updated Systematic ID and Gene name
  columns (12 columns, 4,843 rows).

- ``data/1_formatted/gene_id_mapping_changelog.xlsx``
  Change log with columns: original_systematic_id, original_gene_name,
  updated_systematic_id, updated_gene_name, note, update_type.

- The main output xlsx now also preserves the original identifiers in
  columns ``Original Systematic ID``, ``Original Gene name``, and the
  per‑gene mapping outcome in column ``note``.

Usage
-----
    mamba run -n bioinformatics python src/01_format_and_update_ids.py
    mamba run -n bioinformatics python src/01_format_and_update_ids.py --verbose
    mamba run -n bioinformatics python src/01_format_and_update_ids.py \\
        --gene-filter "gene_type != 'pseudogene'"

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-06-08
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
from enum import StrEnum
import re
import sys
from pathlib import Path
from typing import Any

# 2. Data Processing Imports
import numpy as np
import pandas as pd

# 3. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS & ENUMS
# =============================================================================


class PomBaseCol(StrEnum):
    """Column names in the PomBase gene annotation TSV."""

    GENE_SYSID = "gene_systematic_id"
    GENE_NAME = "gene_name"
    GENE_TYPE = "gene_type"
    SYNONYMS = "synonyms"


# PomBase monthly release version — single source of truth for this pipeline
DEFAULT_RELEASE = "2026-06-01"

# Paths relative to project root
DEFAULT_RAW_DATA = Path("data/raw/rsob130053supp2.xlsx")
DEFAULT_ANNOTATION = (
    Path("data/references")
    / f"pombase-{DEFAULT_RELEASE}_gene_IDs_names_products.tsv"
)
DEFAULT_OUTPUT = Path("data/1_formatted/Hayles_2013_OB_formatted_phenotypes.xlsx")
DEFAULT_CHANGELOG = Path("data/1_formatted/gene_id_mapping_changelog.xlsx")

# =============================================================================
# LOGGING SETUP
# =============================================================================


def setup_logger(log_level: str = "INFO") -> None:
    """Configure the Loguru logger."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level:<8} | {message}",
        level=log_level,
    )


setup_logger()

# =============================================================================
# CORE LOGIC
# =============================================================================


@logger.catch
def format_phenotype_data(df: pd.DataFrame) -> pd.DataFrame:
    """Clean whitespace and normalise temperature notation in phenotype columns."""
    # Strip whitespace from key columns
    df["Systematic ID"] = df["Systematic ID"].str.strip()
    df["Deletion mutant phenotype description"] = (
        df["Deletion mutant phenotype description"]
        .str.strip()
        .str.strip(".")
        .str.strip()
    )
    df["Phenotypic classification used for analysis"] = (
        df["Phenotypic classification used for analysis"].str.strip()
    )
    df["Gene dispensability. This study"] = (
        df["Gene dispensability. This study"].str.strip()
    )

    # Normalise temperature notation: coerce all variants of "25, 32" to "25,32"
    df.replace(
        to_replace="25, 32", value="25,32", regex=True, inplace=True,
    )
    df.replace(to_replace="25 32", value="25,32", regex=True, inplace=True)
    df.replace(
        to_replace="32, 25", value="25,32", regex=True, inplace=True,
    )

    # Summary statistics
    logger.info(
        "Classification categories:\n{}",
        df["Phenotypic classification used for analysis"].value_counts().to_string(),
    )
    logger.info(
        "Dispensability categories:\n{}",
        df["Gene dispensability. This study"].value_counts().to_string(),
    )

    return df


@logger.catch
def read_pombase_annotation(path: Path) -> pd.DataFrame:
    """Load and preprocess the PomBase annotation TSV."""
    logger.info(f"Loading PomBase annotation: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str)

    # Strip whitespace from every string column
    for col in df.columns:
        if df[col].dtype == "object":
            df[col] = df[col].str.strip()

    logger.info(
        f"Loaded {len(df):,} rows, {df[PomBaseCol.GENE_TYPE].nunique()} feature types",
    )
    for t, c in df[PomBaseCol.GENE_TYPE].value_counts().items():
        logger.debug(f"  {t}: {c:,}")

    return df


@logger.catch
def build_update_maps(
    annotation: pd.DataFrame,
    gene_filter: str = f"{PomBaseCol.GENE_TYPE} == 'protein coding gene'",
) -> tuple[set[str], pd.Series, pd.Series]:
    """Build the three lookup structures (sysID set, name→sysID, synonym→sysID)."""
    filtered = annotation.query(gene_filter).copy()
    logger.info(
        f"Filtered to {len(filtered):,} rows matching: {gene_filter}",
    )

    # Fallback: fill NaN gene names with systematic ID
    filtered[PomBaseCol.GENE_NAME] = filtered[PomBaseCol.GENE_NAME].fillna(filtered[PomBaseCol.GENE_SYSID])

    # 1. Current systematic IDs (unique set)
    sysids_now: set[str] = set(filtered[PomBaseCol.GENE_SYSID].unique())
    logger.debug(f"Unique systematic IDs in filtered set: {len(sysids_now):,}")

    # 2. Gene name → systematic ID
    names2id = (
        filtered.set_index(PomBaseCol.GENE_NAME)[PomBaseCol.GENE_SYSID]
        .drop_duplicates()
    )
    logger.debug(f"Gene name → SysID mappings: {len(names2id):,}")

    # 3. Synonyms → systematic ID (explode comma-separated list)
    syn_rows = filtered[[PomBaseCol.GENE_SYSID, PomBaseCol.SYNONYMS]].dropna(subset=[PomBaseCol.SYNONYMS])
    exploded = (
        syn_rows.assign(
            **{PomBaseCol.SYNONYMS: syn_rows[PomBaseCol.SYNONYMS].str.split(",")}
        )
        .explode(PomBaseCol.SYNONYMS)
    )
    exploded[PomBaseCol.SYNONYMS] = exploded[PomBaseCol.SYNONYMS].str.strip()
    exploded = exploded[exploded[PomBaseCol.SYNONYMS] != ""]

    synonyms2id = (
        exploded.set_index(PomBaseCol.SYNONYMS)[PomBaseCol.GENE_SYSID]
    )
    logger.debug(f"Synonym → SysID mappings: {len(synonyms2id):,}")

    return sysids_now, names2id, synonyms2id


@logger.catch
def normalize_sysid(gene: Any) -> str | Any:
    """Normalise a systematic ID to consistent case (e.g. SPAC1002.01c)."""
    if not isinstance(gene, str):
        return gene
    gene = gene.strip()
    if "." in gene:
        parts = gene.split(".")
        return f"{parts[0].upper()}.{parts[1].lower()}"
    return gene


@logger.catch
def update_sysids(
    genes: list[Any],
    annotation: pd.DataFrame,
    gene_filter: str = f"{PomBaseCol.GENE_TYPE} == 'protein coding gene'",
) -> tuple[list[Any], list[str]]:
    """Map gene identifiers to current systematic IDs via three-layer lookup (SysID→Name→Synonym)."""
    sysids_now, names2id, synonyms2id = build_update_maps(annotation, gene_filter)

    updated_ids: list[Any] = []
    notes: list[str] = []

    for i, gene in enumerate(genes):
        if i < 5 or (i + 1) % 500 == 0:
            logger.debug(f"Processing gene {i + 1}/{len(genes)}: {gene!r}")

        # Handle NaN / missing
        if isinstance(gene, float) and np.isnan(gene):
            updated_ids.append(gene)
            notes.append("NA_INPUT")
            continue

        # Normalize string
        gene_str: str = normalize_sysid(gene)

        # Layer 1: direct systematic ID match
        if gene_str in sysids_now:
            updated_ids.append(gene_str)
            notes.append("NO_CHANGE")
            continue

        # Layer 2: gene name → systematic ID
        if gene_str in names2id.index:
            result = names2id.loc[gene_str]
            if isinstance(result, pd.Series):
                updated_ids.append(np.nan)
                notes.append(f"MULTIPLE_NAMES: {result.tolist()}")
                logger.warning(
                    f"Gene name {gene_str!r} maps to multiple SysIDs: "
                    f"{result.tolist()}",
                )
            else:
                updated_ids.append(result)
                notes.append(f"VIA_NAME → {result}")
                logger.info(f"{gene_str} updated via name → {result}")
            continue

        # Layer 3: synonym → systematic ID
        if gene_str in synonyms2id.index:
            result = synonyms2id.loc[gene_str]
            if isinstance(result, pd.Series):
                updated_ids.append(np.nan)
                notes.append(f"MULTIPLE_SYNONYMS: {result.tolist()}")
                logger.warning(
                    f"Synonym {gene_str!r} maps to multiple SysIDs: "
                    f"{result.tolist()}",
                )
            else:
                updated_ids.append(result)
                notes.append(f"VIA_SYNONYM → {result}")
                logger.info(f"{gene_str} updated via synonym → {result}")
            continue

        # Not found in filtered set — check full annotation for diagnosis
        full_mask = annotation[PomBaseCol.GENE_SYSID] == gene_str
        if full_mask.any():
            match_idx = full_mask.idxmax()
            actual_type = annotation.loc[match_idx, PomBaseCol.GENE_TYPE]
            updated_ids.append(gene_str)
            notes.append(f"RECLASSIFIED → {actual_type}")
            logger.info(
                f"{gene_str} is now annotated as '{actual_type}' "
                f"(was protein-coding in Hayles 2013)",
            )
        else:
            # Also check if this ID appears as a synonym (ID was retired,
            # target gene reclassified)
            syn_mask = annotation[PomBaseCol.SYNONYMS].str.contains(
                rf"(?:^|,)\s*{re.escape(gene_str)}\s*(?:,|$)",
                na=False,
                regex=True,
            )
            if syn_mask.any():
                match_idx = syn_mask.idxmax()
                actual_type = annotation.loc[match_idx, PomBaseCol.GENE_TYPE]
                current_id = annotation.loc[match_idx, PomBaseCol.GENE_SYSID]
                updated_ids.append(gene_str)
                notes.append(f"RECLASSIFIED → {actual_type} (synonym of {current_id})")
                logger.info(
                    f"{gene_str} is now a synonym of {current_id} "
                    f"(annotated as '{actual_type}')",
                )
            else:
                updated_ids.append(gene_str)
                notes.append("NOT_FOUND")
                logger.warning(
                    f"{gene_str} not found in any annotation layer",
                )

    return updated_ids, notes


@logger.catch
def lookup_gene_names(
    sysids: list[Any],
    annotation: pd.DataFrame,
) -> pd.Series:
    """Look up current gene names for a list of systematic IDs."""
    # Build name lookup from full annotation (not filtered — we want names
    # even for non-protein-coding genes if they ended up there)
    name_map = (
        annotation
        .dropna(subset=[PomBaseCol.GENE_SYSID])
        .set_index(PomBaseCol.GENE_SYSID)[PomBaseCol.GENE_NAME]
    )

    names = []
    for sid in sysids:
        if isinstance(sid, float) and np.isnan(sid):
            names.append(np.nan)
        elif sid in name_map.index:
            val = name_map.loc[sid]
            names.append(val if pd.notna(val) else np.nan)
        else:
            logger.warning(f"Systematic ID {sid!r} not found in annotation for name lookup")
            names.append(np.nan)

    return pd.Series(names, dtype="object")


@logger.catch
def build_changelog(
    original_ids: list[Any],
    original_names: list[Any],
    updated_ids: list[Any],
    updated_names: list[Any],
    notes: list[str],
) -> pd.DataFrame:
    """Build a changelog DataFrame summarising every mapping result."""
    df = pd.DataFrame({
        "original_systematic_id": original_ids,
        "original_gene_name": original_names,
        "updated_systematic_id": updated_ids,
        "updated_gene_name": updated_names,
        "note": notes,
    })

    # Categorise the note into a simpler update_type
    def categorize(note: str) -> str:
        if note == "NA_INPUT":
            return "NA"
        if note == "NO_CHANGE":
            return "NO_CHANGE"
        if note.startswith("VIA_NAME"):
            return "VIA_NAME"
        if note.startswith("VIA_SYNONYM"):
            return "VIA_SYNONYM"
        if note.startswith("MULTIPLE"):
            return "CONFLICT"
        if note.startswith("RECLASSIFIED"):
            return "RECLASSIFIED"
        if note == "NOT_FOUND":
            return "NOT_FOUND"
        return "OTHER"

    df["update_type"] = df["note"].apply(categorize)
    return df


# =============================================================================
# MAIN EXECUTION
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Format raw Hayles 2013 phenotypes and update systematic IDs "
        "against current PomBase annotation.",
    )
    parser.add_argument(
        "--raw",
        type=Path,
        default=DEFAULT_RAW_DATA,
        help=f"Hayles raw input (default: {DEFAULT_RAW_DATA})",
    )
    parser.add_argument(
        "--annotation",
        type=Path,
        default=DEFAULT_ANNOTATION,
        help=f"PomBase annotation TSV (default: {DEFAULT_ANNOTATION})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output formatted phenotypes xlsx (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--changelog",
        type=Path,
        default=DEFAULT_CHANGELOG,
        help=f"Output changelog xlsx (default: {DEFAULT_CHANGELOG})",
    )
    parser.add_argument(
        "--gene-filter",
        type=str,
        default=f"{PomBaseCol.GENE_TYPE} == 'protein coding gene'",
        help="Pandas query to filter annotation rows for mapping",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable DEBUG level logging",
    )
    return parser.parse_args()


def main() -> int:
    """Main orchestrator."""
    args = parse_args()

    if args.verbose:
        setup_logger(log_level="DEBUG")

    # Resolve paths
    raw_path = args.raw.resolve()
    annot_path = args.annotation.resolve()
    output_path = args.output.resolve()
    changelog_path = args.changelog.resolve()

    # ------------------------------------------------------------------
    # 1. Load raw Hayles data
    # ------------------------------------------------------------------
    logger.info(f"Loading raw phenotypes: {raw_path}")
    raw = pd.read_excel(raw_path)

    # Drop the spurious 'Unnamed: 12' column if present
    unnamed = [c for c in raw.columns if c.startswith("Unnamed:")]
    if unnamed:
        raw.drop(columns=unnamed, inplace=True)
        logger.debug(f"Dropped spurious column(s): {unnamed}")

    # ------------------------------------------------------------------
    # 2. Format phenotype data (strip whitespace, normalise temperature)
    # ------------------------------------------------------------------
    logger.info("Formatting phenotype data…")
    raw = format_phenotype_data(raw)

    # ------------------------------------------------------------------
    # 3. Extract original identifiers before mapping
    # ------------------------------------------------------------------
    original_ids: list[Any] = raw["Systematic ID"].tolist()
    original_names: list[Any] = raw["Gene name"].tolist()
    logger.info(f"Loaded {len(raw):,} genes, {raw['Gene name'].isna().sum():,} missing gene names")

    # ------------------------------------------------------------------
    # 4. Load PomBase annotation
    # ------------------------------------------------------------------
    if not annot_path.exists():
        logger.error(
            f"Annotation file not found: {annot_path}\n"
            "Run 00_download_pombase_annotation.py first.",
        )
        return 1

    annotation = read_pombase_annotation(annot_path)

    # ------------------------------------------------------------------
    # 5. Update systematic IDs
    # ------------------------------------------------------------------
    logger.info("Mapping systematic IDs through PomBase annotation…")
    updated_ids, notes = update_sysids(
        original_ids,
        annotation,
        gene_filter=args.gene_filter,
    )

    # Count changes
    n_changed = sum(1 for n in notes if n.startswith("VIA"))
    n_not_found = sum(1 for n in notes if n == "NOT_FOUND")
    n_conflict = sum(1 for n in notes if n.startswith("MULTIPLE"))
    logger.info(
        f"Mapping complete: {n_changed} updated, "
        f"{n_not_found} not found, {n_conflict} conflicts",
    )

    # ------------------------------------------------------------------
    # 6. Update gene names
    # ------------------------------------------------------------------
    logger.info("Looking up current gene names…")
    updated_names_series = lookup_gene_names(updated_ids, annotation)
    updated_names = updated_names_series.tolist()

    # ------------------------------------------------------------------
    # 7. Build changelog
    # ------------------------------------------------------------------
    changelog_df = build_changelog(
        original_ids,
        original_names,
        updated_ids,
        updated_names,
        notes,
    )
    logger.info(f"Changelog summary:\n{changelog_df['update_type'].value_counts().to_string()}")

    # ------------------------------------------------------------------
    # 8. Update the raw DataFrame and save
    # ------------------------------------------------------------------
    raw["Original Systematic ID"] = original_ids
    raw["Original Gene name"] = original_names
    raw["Systematic ID"] = updated_ids
    raw["Gene name"] = updated_names
    raw["note"] = notes

    # Fill remaining NaN gene names with the corresponding Systematic ID
    raw["Gene name"] = raw["Gene name"].fillna(raw["Systematic ID"])

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    changelog_path.parent.mkdir(parents=True, exist_ok=True)

    raw.to_excel(output_path, index=False)
    logger.success(f"Updated phenotypes saved: {output_path} ({len(raw):,} rows)")

    changelog_df.to_excel(changelog_path, index=False)
    logger.success(f"Changelog saved: {changelog_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
