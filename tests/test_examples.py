"""Smoke tests for examples."""

from __future__ import annotations

import runpy
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).resolve().parents[1] / "examples"
BASIC_DIR = EXAMPLES_DIR / "basic"
DATABASE_DIR = EXAMPLES_DIR / "database"
VALIDATION_DIR = EXAMPLES_DIR / "validation"
CLI_DIR = EXAMPLES_DIR / "cli"


def _list_scripts(directory: Path) -> list[Path]:
    return sorted(path for path in directory.glob("*.py") if path.name != "__init__.py")


CORE_SCRIPTS = _list_scripts(BASIC_DIR)
DATABASE_SCRIPTS = _list_scripts(DATABASE_DIR)
VALIDATION_SCRIPTS = _list_scripts(VALIDATION_DIR)
CLI_SCRIPTS = sorted(CLI_DIR.glob("*.sh"))


@pytest.mark.parametrize("script_path", CLI_SCRIPTS, ids=lambda p: p.name)
def test_cli_example_runs_successfully(script_path: Path) -> None:
    proc = subprocess.run(
        ["bash", str(script_path)],
        cwd=str(EXAMPLES_DIR.parent),
        env={"PYTHON": sys.executable, "PATH": "/usr/bin:/bin"},
        check=False,
        text=True,
        capture_output=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


@pytest.mark.parametrize("script_path", CORE_SCRIPTS, ids=lambda p: p.name)
def test_basic_example_runs_successfully(script_path: Path) -> None:
    runpy.run_path(str(script_path), run_name="__main__")


@pytest.mark.parametrize("script_path", DATABASE_SCRIPTS, ids=lambda p: p.name)
def test_database_example_runs_successfully(script_path: Path) -> None:
    runpy.run_path(str(script_path), run_name="__main__")


@pytest.mark.parametrize("script_path", VALIDATION_SCRIPTS, ids=lambda p: p.name)
def test_validation_example_runs_successfully(
    script_path: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("thermo")
    monkeypatch.setenv("CHEMTHERMO_OUTDIR", str(tmp_path))
    runpy.run_path(str(script_path), run_name="__main__")
