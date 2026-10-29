"""
Systematic ID Mapping
=====================

Shared ID-update engine for the Kim 2010 and Hayles 2013 tables. It maps
pre-update identifiers to current PomBase systematic IDs through three
lookup layers:

1. current systematic ID match;
2. gene name → systematic ID;
3. synonym → systematic ID.

IDs that do not match any layer are diagnosed against the full annotation:
a row can be `RECLASSIFIED` (for example, a former protein-coding gene now
annotated as a pseudogene) or `NOT_FOUND`.

The module also collects the IDs that need mapping, processes Kim's Paralog
lists with the same maps, and combines per-ID notes for collapsed rows.

Usage
-----
    from id_mapping import (
        collect_mapping_ids,
        read_pombase_annotation,
        update_sysids,
    )

    annotation = read_pombase_annotation(annotation_path)
    ids = collect_mapping_ids(kim, hayles)
    mapped_ids, notes = update_sysids(ids, annotation)

Author:   Yusheng Yang (guidance) + Claude Code (implementation)
Date:     2026-10-29
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
from enum import StrEnum
from pathlib import Path
import re
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
    GENE_PRODUCT = "gene_product"
    GENE_TYPE = "gene_type"
    SYNONYMS = "synonyms"


DEFAULT_GENE_FILTER = f"{PomBaseCol.GENE_TYPE} == 'protein coding gene'"

# =============================================================================
# POMBASE ANNOTATION
# =============================================================================


@logger.catch
def read_pombase_annotation(path: Path) -> pd.DataFrame:
    """Load and preprocess the PomBase annotation TSV."""
    logger.info(f"Loading PomBase annotation: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str)

    # Strip whitespace from every string column
    for col in df.columns:
        if pd.api.types.is_string_dtype(df[col]):
            df[col] = df[col].str.strip()

    logger.info(
        f"Loaded {len(df):,} rows, {df[PomBaseCol.GENE_TYPE].nunique()} feature types",
    )
    for gene_type, count in df[PomBaseCol.GENE_TYPE].value_counts().items():
        logger.debug(f"  {gene_type}: {count:,}")

    return df


@logger.catch
def build_update_maps(
    annotation: pd.DataFrame,
    gene_filter: str = DEFAULT_GENE_FILTER,
) -> tuple[set[str], pd.Series, pd.Series]:
    """Build the three lookup structures (sysID set, name→sysID, synonym→sysID)."""
    filtered = annotation.query(gene_filter).copy()
    logger.info(
        f"Filtered to {len(filtered):,} rows matching: {gene_filter}",
    )

    # Fallback: fill NaN gene names with systematic ID
    filtered[PomBaseCol.GENE_NAME] = filtered[PomBaseCol.GENE_NAME].fillna(
        filtered[PomBaseCol.GENE_SYSID]
    )

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
    syn_rows = filtered[
        [PomBaseCol.GENE_SYSID, PomBaseCol.SYNONYMS]
    ].dropna(subset=[PomBaseCol.SYNONYMS])
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


# =============================================================================
# ID UPDATE
# =============================================================================


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
    gene_filter: str = DEFAULT_GENE_FILTER,
) -> tuple[list[Any], list[str]]:
    """Map identifiers to current systematic IDs via the three lookup layers."""
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


# =============================================================================
# ID COLLECTION & PARALOGS
# =============================================================================


@logger.catch
def collect_mapping_ids(
    kim2010: pd.DataFrame,
    hayles2013: pd.DataFrame,
) -> list[str]:
    """Collect unique main-table and Kim-paralog IDs before the ID update."""
    seen: set[str] = set()
    ids: list[str] = []

    def add(value: object) -> None:
        if pd.isna(value):
            return
        text = str(value).strip()
        if not text or text in seen:
            return
        seen.add(text)
        ids.append(text)

    for value in kim2010["id"]:
        add(value)
    for value in hayles2013["id"]:
        add(value)
    for value in kim2010["paralogues_2010"].dropna():
        text = str(value).strip()
        if text.lower() in {"none", "being curated"}:
            continue
        for token in text.split("|"):
            add(token)

    logger.info(f"Collected {len(ids):,} unique IDs (main + paralog) for mapping")
    return ids


@logger.catch
def build_paralog_columns(
    raw_paralogs: pd.Series,
    main_ids: pd.Series,
    id_map: dict[str, str],
    name_map: dict[str, Any],
) -> tuple[pd.Series, pd.Series]:
    """Map Kim paralog IDs/names to current values and drop self-references."""
    paralog_ids: list[str | None] = []
    paralog_names: list[str | None] = []

    for raw, main_id in zip(raw_paralogs, main_ids):
        if pd.isna(raw):
            paralog_ids.append(None)
            paralog_names.append(None)
            continue

        text = str(raw).strip()
        if not text or text.lower() == "none":
            paralog_ids.append(None)
            paralog_names.append(None)
            continue
        if text.lower() == "being curated":
            paralog_ids.append("being curated")
            paralog_names.append(None)
            continue

        ids: list[str] = []
        names: list[str] = []
        for token in text.split("|"):
            token = token.strip()
            if not token:
                continue
            new_id = id_map.get(token, token)
            if pd.isna(new_id):
                new_id = token
            if str(new_id).upper() == str(main_id).upper():
                continue
            if new_id in ids:
                continue
            ids.append(new_id)
            name = name_map.get(new_id)
            names.append(name if pd.notna(name) and name else new_id)

        paralog_ids.append("|".join(ids) if ids else None)
        paralog_names.append("|".join(names) if names else None)

    return (
        pd.Series(paralog_ids, dtype="string"),
        pd.Series(paralog_names, dtype="string"),
    )


def combine_mapping_notes(original_ids: object, note_map: dict[str, str]) -> str:
    """Combine per-original-ID mapping notes for one row."""
    parts = [
        part.strip()
        for part in str(original_ids).split(";")
        if part.strip()
    ]
    if len(parts) > 1:
        detail = "; ".join(
            f"{gene_id}: {note_map.get(gene_id, 'NOT_FOUND')}"
            for gene_id in parts
        )
        return f"COLLAPSED: {detail}"
    gene_id = parts[0] if parts else str(original_ids)
    return note_map.get(gene_id, "NOT_FOUND")
