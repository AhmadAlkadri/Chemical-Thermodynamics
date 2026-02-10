"""Smoke tests for examples."""

from __future__ import annotations

import runpy
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).resolve().parents[1] / "examples"
BASIC_DIR = EXAMPLES_DIR / "basic"
VALIDATION_DIR = EXAMPLES_DIR / "validation"


def _list_scripts(directory: Path) -> list[Path]:
    return sorted(path for path in directory.glob("*.py") if path.name != "__init__.py")


CORE_SCRIPTS = _list_scripts(BASIC_DIR)
VALIDATION_SCRIPTS = _list_scripts(VALIDATION_DIR)


@pytest.mark.parametrize("script_path", CORE_SCRIPTS, ids=lambda p: p.name)
def test_basic_example_runs_successfully(script_path: Path) -> None:
    runpy.run_path(str(script_path), run_name="__main__")


@pytest.mark.parametrize("script_path", VALIDATION_SCRIPTS, ids=lambda p: p.name)
def test_validation_example_runs_successfully(
    script_path: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytest.importorskip("thermo")
    monkeypatch.setenv("CHEMTHERMO_OUTDIR", str(tmp_path))
    runpy.run_path(str(script_path), run_name="__main__")
