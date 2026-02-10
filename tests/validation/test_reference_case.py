from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "examples" / "validation" / "00_reference_case.py"


def test_validation_reference_case_script_passes() -> None:
    pytest.importorskip("thermo")

    proc = subprocess.run(
        [sys.executable, str(SCRIPT)],
        cwd=str(REPO_ROOT),
        check=False,
        text=True,
        capture_output=True,
    )

    assert proc.returncode == 0, proc.stderr
    assert "PASS" in proc.stdout
