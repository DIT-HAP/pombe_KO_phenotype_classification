#!/usr/bin/env python3
"""
Download PomBase Gene Annotation File
======================================

Download a specific monthly release of the PomBase gene annotation file
(gene_IDs_names_products.tsv) and cache it locally. This file serves as
the reference for updating gene systematic identifiers in downstream
pipeline steps.

The PomBase annotation covers all annotated S. pombe feature types
(protein-coding genes, tRNAs, rRNAs, snoRNAs, snRNAs, lncRNAs,
pseudogenes, etc.) in 8 columns — systematic IDs, gene symbols,
synonyms, feature types, and cross-references.

Input
-----
- URL to PomBase monthly release (constructed from ``--release``).
  Example: https://www.pombase.org/monthly_releases/2026/
           pombase-2026-06-01/gene_names_and_identifiers/
           gene_IDs_names_products.tsv

Output
------
- ``data/references/pombase-{release}_gene_IDs_names_products.tsv``
  Tab-separated file; the 4 columns used by downstream scripts are:
    - gene_systematic_id (str): PomBase systematic identifier
    - gene_name (str): Standard gene symbol (may be empty)
    - gene_type (str): Feature type (e.g. 'protein coding gene')
    - synonyms (str): Comma-separated alternative names

Usage
-----
    mamba run -n bioinformatics python scripts/00_download_pombase_annotation.py
    mamba run -n bioinformatics python scripts/00_download_pombase_annotation.py --release 2026-05-01
    mamba run -n bioinformatics python scripts/00_download_pombase_annotation.py --force --verbose

Author:   Yusheng Yang (guidance) + Hermes (implementation)
Date:     2026-06-08
Version:  1.0.0
"""

# =============================================================================
# IMPORTS
# =============================================================================
# 1. Standard Library Imports
import argparse
import sys
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen

# 2. Third-party Imports
from loguru import logger

# =============================================================================
# GLOBAL CONSTANTS
# =============================================================================

# Default release version (YYYY-MM-DD format, matching PomBase monthly release)
DEFAULT_RELEASE = "2026-06-01"

# Default output directory (relative to project root)
DEFAULT_OUTDIR = Path("data/references")

# Base URL pattern for PomBase monthly release gene annotation files.
# The release directory uses a ``pombase-`` prefix in the path.
POMBASE_URL_TEMPLATE = (
    "https://www.pombase.org/monthly_releases/{year}/pombase-{release}/"
    "gene_names_and_identifiers/gene_IDs_names_products.tsv"
)

# =============================================================================
# LOGGING SETUP
# =============================================================================


def setup_logger(log_level: str = "INFO") -> None:
    """Configure the Loguru logger with stdout output."""
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
def build_url(release: str) -> str:
    """Build the full PomBase download URL for the given release version."""
    year = release[:4]
    return POMBASE_URL_TEMPLATE.format(year=year, release=release)


@logger.catch
def build_output_path(outdir: Path, release: str) -> Path:
    """Build the local cache file path from release version."""
    filename = f"pombase-{release}_gene_IDs_names_products.tsv"
    return outdir / filename


@logger.catch
def download_file(url: str, dst: Path) -> int:
    """Download a file from *url* to *dst*; returns bytes downloaded."""
    logger.info(f"Downloading from: {url}")

    try:
        with urlopen(url, timeout=120) as response:
            status = response.status
            if status != 200:
                msg = f"HTTP {status} — {response.reason}"
                raise RuntimeError(msg)

            content = response.read()
    except URLError as e:
        raise RuntimeError(f"Network error while downloading: {e}") from e

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(content)

    return len(content)


@logger.catch
def inspect_file(path: Path) -> dict:
    """Print file statistics; returns metadata dict for downstream use."""
    if not path.exists():
        logger.warning(f"File not found: {path}")
        return {}

    file_size = path.stat().st_size
    logger.info(f"File size: {file_size:,} bytes ({file_size / 1024 / 1024:.2f} MB)")

    # Peek at the file to find the header line
    with path.open("r", encoding="utf-8") as f:
        lines = [f.readline() for _ in range(5)]

    header_line = None
    comment_lines = 0
    for line in lines:
        if line.startswith("#"):
            comment_lines += 1
            continue
        header_line = line
        break

    columns: list[str] = []
    if header_line:
        columns = header_line.strip().split("\t")
        logger.info(f"Columns ({len(columns)}): {', '.join(columns)}")
    else:
        logger.warning("Could not find header line in the file")

    logger.info(f"Comment lines at top: {comment_lines}")

    # Count total lines efficiently (avoids loading the whole file)
    with path.open("r", encoding="utf-8") as f:
        total_lines = sum(1 for _ in f)
    data_lines = total_lines - comment_lines - 1  # subtract comments + header
    logger.info(f"Total lines: {total_lines:,} (≈{data_lines:,} data rows)")

    return {
        "path": str(path),
        "file_size": file_size,
        "total_lines": total_lines,
        "data_lines": data_lines,
        "columns": columns,
        "release": path.stem.replace("pombase-", "").replace("_gene_IDs_names_products", ""),
    }


# =============================================================================
# MAIN EXECUTION
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Download PomBase gene annotation (gene_IDs_names_products.tsv) "
        "for a specific monthly release.",
    )
    parser.add_argument(
        "--release",
        type=str,
        default=DEFAULT_RELEASE,
        help=f"PomBase monthly release version (default: {DEFAULT_RELEASE})",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help=f"Output directory for cached files (default: {DEFAULT_OUTDIR})",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download even if the file already exists locally",
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

    outdir = args.outdir.resolve()
    url = build_url(args.release)
    dst = build_output_path(outdir, args.release)

    logger.info(f"Release: {args.release}")
    logger.info(f"Output:  {dst}")

    # Cache hit — skip unless --force
    if dst.exists() and not args.force:
        logger.info("File already cached locally. Use --force to re-download.")
        inspect_file(dst)
        return 0

    # Download
    try:
        bytes_downloaded = download_file(url, dst)
        logger.info(f"Downloaded: {bytes_downloaded:,} bytes ({bytes_downloaded / 1024 / 1024:.2f} MB)")

        metadata = inspect_file(dst)
        logger.info(
            f"Annotation file saved to: {dst} "
            f"({len(metadata.get('columns', []))} columns, {metadata.get('data_lines', 0)} rows)",
        )
    except (RuntimeError, URLError, OSError) as e:
        logger.error(f"Download failed: {e}")
        if dst.exists():
            dst.unlink()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
