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
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
import re
from typing import Any

# 2. Third-party Imports
from loguru import logger

# 3. Data Processing Imports
import numpy as np
import pandas as pd

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


def read_pombase_annotation(path: Path) -> pd.DataFrame:
    """Load and preprocess the PomBase annotation TSV."""
    logger.info(f"Loading PomBase annotation: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str)

    missing = set(PomBaseCol) - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns {sorted(missing)}")

    df = df[list(PomBaseCol)].copy()
    for col in PomBaseCol:
        df[col] = df[col].str.strip()

    logger.info(
        f"Loaded {len(df):,} rows, {df[PomBaseCol.GENE_TYPE].nunique()} feature types",
    )
    for gene_type, count in df[PomBaseCol.GENE_TYPE].value_counts().items():
        logger.debug(f"  {gene_type}: {count:,}")

    return df


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
    names2id = filtered.set_index(PomBaseCol.GENE_NAME)[PomBaseCol.GENE_SYSID].drop_duplicates()
    logger.debug(f"Gene name → SysID mappings: {len(names2id):,}")

    # 3. Synonyms → systematic ID (explode comma-separated list)
    syn_rows = filtered[[PomBaseCol.GENE_SYSID, PomBaseCol.SYNONYMS]].dropna(
        subset=[PomBaseCol.SYNONYMS]
    )
    exploded = syn_rows.assign(
        **{PomBaseCol.SYNONYMS: syn_rows[PomBaseCol.SYNONYMS].str.split(",")}
    ).explode(PomBaseCol.SYNONYMS)
    exploded[PomBaseCol.SYNONYMS] = exploded[PomBaseCol.SYNONYMS].str.strip()
    exploded = exploded[exploded[PomBaseCol.SYNONYMS] != ""]

    synonyms2id = exploded.set_index(PomBaseCol.SYNONYMS)[PomBaseCol.GENE_SYSID]
    logger.debug(f"Synonym → SysID mappings: {len(synonyms2id):,}")

    return sysids_now, names2id, synonyms2id


# =============================================================================
# ID UPDATE
# =============================================================================


def normalize_sysid(gene: Any) -> Any:
    """Normalise a systematic ID to consistent case (e.g. SPAC1002.01c)."""
    if not isinstance(gene, str):
        return gene
    gene = gene.strip()
    if "." in gene:
        parts = gene.split(".")
        return f"{parts[0].upper()}.{parts[1].lower()}"
    return gene


@dataclass(frozen=True)
class _LookupMaps:
    """The three lookup layers built from one PomBase annotation."""

    sysids: set[str]
    names: pd.Series
    synonyms: pd.Series


@dataclass(frozen=True)
class _IdOutcome:
    """One identifier's updated value together with its mapping note."""

    updated: Any
    note: str


def _is_missing_gene(gene: Any) -> bool:
    """Return True for a float NaN input."""
    return isinstance(gene, float) and np.isnan(gene)


def _lookup_in_mapping(
    gene_str: str,
    mapping: pd.Series,
    kind: str,
) -> _IdOutcome | None:
    """Resolve *gene_str* through one name or synonym layer, if present."""
    if gene_str not in mapping.index:
        return None
    result = mapping.loc[gene_str]
    if isinstance(result, pd.Series):
        # Keep the original ID rather than dropping the row.
        return _IdOutcome(gene_str, f"MULTIPLE_{kind}S: {result.tolist()}")
    return _IdOutcome(result, f"VIA_{kind} → {result}")


def _diagnose_unmapped(gene_str: str, annotation: pd.DataFrame) -> _IdOutcome:
    """Explain why *gene_str* is absent from the filtered annotation."""
    by_sysid = annotation[PomBaseCol.GENE_SYSID] == gene_str
    if by_sysid.any():
        actual_type = annotation.loc[by_sysid.idxmax(), PomBaseCol.GENE_TYPE]
        return _IdOutcome(gene_str, f"RECLASSIFIED → {actual_type}")

    synonym_pattern = rf"(?:^|,)\s*{re.escape(gene_str)}\s*(?:,|$)"
    by_synonym = annotation[PomBaseCol.SYNONYMS].str.contains(
        synonym_pattern,
        na=False,
        regex=True,
    )
    if not by_synonym.any():
        return _IdOutcome(gene_str, "NOT_FOUND")

    row = annotation.loc[by_synonym.idxmax()]
    note = f"RECLASSIFIED → {row[PomBaseCol.GENE_TYPE]} (synonym of {row[PomBaseCol.GENE_SYSID]})"
    return _IdOutcome(gene_str, note)


def _resolve_gene(
    gene: Any,
    maps: _LookupMaps,
    annotation: pd.DataFrame,
) -> _IdOutcome:
    """Resolve one identifier through the three lookup layers and diagnosis."""
    gene_str = normalize_sysid(gene)
    if gene_str in maps.sysids:
        return _IdOutcome(gene_str, "NO_CHANGE")
    for kind, mapping in (("NAME", maps.names), ("SYNONYM", maps.synonyms)):
        outcome = _lookup_in_mapping(gene_str, mapping, kind)
        if outcome is not None:
            return outcome
    return _diagnose_unmapped(gene_str, annotation)


def _log_outcome(gene: Any, outcome: _IdOutcome) -> None:
    """Log one resolution with the severity its note implies."""
    if outcome.note.startswith("MULTIPLE"):
        logger.warning(f"{gene!r} maps to multiple SysIDs: {outcome.note}")
    elif outcome.note.startswith("VIA"):
        logger.info(f"{gene} updated: {outcome.note}")
    elif outcome.note.startswith("RECLASSIFIED"):
        logger.info(f"{gene} reclassified: {outcome.note}")
    elif outcome.note == "NOT_FOUND":
        logger.warning(f"{gene} not found in any annotation layer")


def update_sysids(
    genes: list[Any],
    annotation: pd.DataFrame,
    gene_filter: str = DEFAULT_GENE_FILTER,
) -> tuple[list[Any], list[str]]:
    """Map identifiers to current systematic IDs via the three lookup layers."""
    sysids_now, names2id, synonyms2id = build_update_maps(annotation, gene_filter)
    maps = _LookupMaps(sysids=sysids_now, names=names2id, synonyms=synonyms2id)

    updated_ids: list[Any] = []
    notes: list[str] = []
    for index, gene in enumerate(genes):
        if index < 5 or (index + 1) % 500 == 0:
            logger.debug(f"Processing gene {index + 1}/{len(genes)}: {gene!r}")

        if _is_missing_gene(gene):
            updated_ids.append(gene)
            notes.append("NA_INPUT")
            continue

        outcome = _resolve_gene(gene, maps, annotation)
        _log_outcome(gene, outcome)
        updated_ids.append(outcome.updated)
        notes.append(outcome.note)

    return updated_ids, notes


# =============================================================================
# ID COLLECTION & PARALOGS
# =============================================================================


def _split_paralog_tokens(value: object) -> list[str]:
    """Split Kim's pipe-separated Paralog field, dropping field sentinels."""
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text or text.lower() in {"none", "being curated"}:
        return []
    return [token.strip() for token in text.split("|") if token.strip()]


def collect_mapping_ids(
    kim2010: pd.DataFrame,
    hayles2013: pd.DataFrame,
) -> list[str]:
    """Collect unique main-table and Kim-paralog IDs before the ID update."""
    main_ids = pd.concat([kim2010["id"], hayles2013["id"]], ignore_index=True)
    paralog_ids = kim2010["paralogues_2010"].dropna().map(_split_paralog_tokens).explode()
    ids = pd.concat([main_ids, paralog_ids], ignore_index=True).dropna().astype(str).str.strip()
    ids = ids[ids.ne("")].drop_duplicates().tolist()

    logger.info(f"Collected {len(ids):,} unique IDs (main + paralog) for mapping")
    return ids


def build_gene_name_map(annotation: pd.DataFrame) -> dict[str, str]:
    """Map current systematic IDs to gene names, falling back to the ID."""
    by_id = (
        annotation.dropna(subset=[PomBaseCol.GENE_SYSID])
        .drop_duplicates(PomBaseCol.GENE_SYSID)
        .set_index(PomBaseCol.GENE_SYSID)[PomBaseCol.GENE_NAME]
    )

    cleaned = by_id.astype("string").str.strip()
    fallback = pd.Series(by_id.index.astype(str), index=by_id.index)
    resolved = cleaned.where(cleaned.notna() & cleaned.ne(""), fallback)
    return {str(key): str(value) for key, value in resolved.items()}


def _map_paralog_tokens(
    tokens: list[str],
    main_id: Any,
    id_map: dict[str, str],
    name_map: dict[str, Any],
) -> tuple[list[str], list[str]]:
    """Map one row's paralog tokens, dropping self-references and duplicates."""
    ids: list[str] = []
    names: list[str] = []
    for token in tokens:
        new_id = id_map.get(token, token)
        if str(new_id).upper() == str(main_id).upper() or new_id in ids:
            continue
        ids.append(new_id)
        names.append(name_map.get(new_id, new_id))
    return ids, names


def build_paralog_columns(
    raw_paralogs: pd.Series,
    main_ids: pd.Series,
    id_map: dict[str, str],
    name_map: dict[str, Any],
) -> tuple[pd.Series, pd.Series]:
    """Map Kim paralog IDs/names to current values and drop self-references."""
    paralog_ids: list[str | None] = []
    paralog_names: list[str | None] = []

    for raw, main_id in zip(raw_paralogs, main_ids, strict=True):
        text = "" if pd.isna(raw) else str(raw).strip()
        if text.lower() == "being curated":
            paralog_ids.append("being curated")
            paralog_names.append("being curated")
            continue

        ids, names = _map_paralog_tokens(
            _split_paralog_tokens(text),
            main_id,
            id_map,
            name_map,
        )
        paralog_ids.append("|".join(ids) if ids else None)
        paralog_names.append("|".join(names) if names else None)

    return (
        pd.Series(paralog_ids, dtype="string"),
        pd.Series(paralog_names, dtype="string"),
    )


def combine_mapping_notes(original_ids: object, note_map: dict[str, str]) -> str:
    """Combine per-original-ID mapping notes for one row."""
    parts = [part.strip() for part in str(original_ids).split(";") if part.strip()]
    if len(parts) > 1:
        detail = "; ".join(f"{gene_id}: {note_map.get(gene_id, 'NOT_FOUND')}" for gene_id in parts)
        return f"COLLAPSED: {detail}"
    gene_id = parts[0] if parts else str(original_ids)
    return note_map.get(gene_id, "NOT_FOUND")
