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
def test_validation_example_runs_successfully(script_path: Path, tmp_path: Path) -> None:
    pytest.importorskip("thermo")

    script_copy = tmp_path / script_path.name
    script_copy.write_text(script_path.read_text(encoding="utf-8"), encoding="utf-8")

    runpy.run_path(str(script_copy), run_name="__main__")
