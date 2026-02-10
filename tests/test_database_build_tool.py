from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from chemthermo.schemas import Database

REPO_ROOT = Path(__file__).resolve().parents[1]
BUILD_SCRIPT = REPO_ROOT / "tools" / "build_database.py"


def _run_build(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(BUILD_SCRIPT), *args],
        cwd=str(cwd or REPO_ROOT),
        check=False,
        text=True,
        capture_output=True,
    )


def test_build_database_writes_schema_valid_payload(tmp_path: Path) -> None:
    output_path = tmp_path / "components.generated.json"

    proc = _run_build("--output", str(output_path))
    assert proc.returncode == 0, proc.stderr
    assert output_path.exists()

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    db = Database(**payload)
    assert db.schema_version == 1
    assert len(db.components) > 0


def test_build_database_check_passes_against_canonical_runtime_db() -> None:
    proc = _run_build("--check")
    assert proc.returncode == 0, proc.stderr
    assert "Canonical runtime database is up to date" in proc.stdout


def test_build_database_check_detects_drift() -> None:
    # Restricting to one component intentionally differs from canonical payload.
    proc = _run_build("--check", "--names", "Methane")
    assert proc.returncode == 1
    assert "out of sync" in proc.stderr


def test_build_database_writes_non_runtime_mirror(tmp_path: Path) -> None:
    output_path = tmp_path / "components.generated.json"
    mirror_path = tmp_path / "components.mirror.json"

    proc = _run_build(
        "--output",
        str(output_path),
        "--write-mirror",
        "--mirror-output",
        str(mirror_path),
    )
    assert proc.returncode == 0, proc.stderr
    assert output_path.exists()
    assert mirror_path.exists()
