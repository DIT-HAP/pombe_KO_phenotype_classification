"""Tests for the shell pipeline runner (run_pipeline.sh).

These exercise only the dry-run/selection path, so they need neither mamba
nor the pipeline data files. They pin the execution order (which is not the
numeric order) and the explicit-path wiring between steps.
"""

from __future__ import annotations

from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "run_pipeline.sh"


def run_pipeline(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["bash", str(SCRIPT), *args],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def selected_steps(stdout: str) -> list[str]:
    return re.findall(r"^\[(?:STEP|DRY-RUN) \d+/\d+\] (\d{2})", stdout, re.MULTILINE)


def test_dry_run_uses_dependency_order():
    result = run_pipeline("--dry-run")
    assert result.returncode == 0
    assert selected_steps(result.stdout) == [
        "00",
        "01",
        "03",
        "04",
        "05",
        "02",
        "06",
        "07",
    ]


def test_only_is_sorted_in_execution_order():
    result = run_pipeline("--dry-run", "--only", "06,00")
    assert result.returncode == 0
    assert selected_steps(result.stdout) == ["00", "06"]


def test_from_and_to_select_a_contiguous_execution_window():
    result = run_pipeline("--dry-run", "--from", "04", "--to", "02")
    assert result.returncode == 0
    assert selected_steps(result.stdout) == ["04", "05", "02"]


def test_skip_removes_steps():
    result = run_pipeline("--dry-run", "--skip", "06,07")
    assert result.returncode == 0
    assert selected_steps(result.stdout) == ["00", "01", "03", "04", "05", "02"]


def test_only_rejects_combined_range_options():
    result = run_pipeline("--dry-run", "--only", "05", "--from", "04")
    assert result.returncode == 2
    assert "--only cannot be combined" in result.stderr


def test_unknown_step_is_rejected():
    result = run_pipeline("--dry-run", "--only", "99")
    assert result.returncode == 2
    assert "unknown step" in result.stderr


def test_annotation_path_follows_release():
    result = run_pipeline("--dry-run", "--only", "00,01", "--release", "2026-05-01")
    assert result.returncode == 0
    assert "pombase-2026-05-01_gene_IDs_names_products.tsv" in result.stdout
    assert "--release 2026-05-01" in result.stdout


def test_list_shows_all_steps_in_execution_order():
    result = run_pipeline("--list")
    assert result.returncode == 0
    listed = re.findall(r"^  (\d{2})\s+\S", result.stdout, re.MULTILINE)
    assert listed == ["00", "01", "03", "04", "05", "02", "06", "07"]
