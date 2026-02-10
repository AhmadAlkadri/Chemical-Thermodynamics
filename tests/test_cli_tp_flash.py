from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PATH = REPO_ROOT / "tests" / "fixtures" / "cli" / "tp_flash_v1.json"


def _run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "chemthermo", *args],
        cwd=str(REPO_ROOT),
        check=False,
        text=True,
        capture_output=True,
    )


def test_cli_tp_flash_json_output_matches_fixture_contract() -> None:
    proc = _run_cli(
        "tp-flash",
        "--components",
        "Methane,Ethane,Propane",
        "--z",
        "0.5,0.3,0.2",
        "--temperature-k",
        "240",
        "--pressure-pa",
        "3000000",
        "--format",
        "json",
    )
    assert proc.returncode == 0, proc.stderr

    payload = json.loads(proc.stdout)
    fixture = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))

    assert set(payload) == set(fixture)
    assert payload["cli_schema_version"] == fixture["cli_schema_version"]
    assert payload["command"] == fixture["command"]
    assert payload["solver"] == fixture["solver"]

    assert payload["inputs"]["components"] == fixture["inputs"]["components"]
    assert payload["inputs"]["z_mole"] == fixture["inputs"]["z_mole"]
    assert payload["result"]["component_order"] == fixture["result"]["component_order"]
    assert payload["result"]["phase_names"] == fixture["result"]["phase_names"]

    assert payload["result"]["vapor_fraction"] == pytest.approx(
        fixture["result"]["vapor_fraction"], rel=1e-9, abs=1e-12
    )

    assert np.allclose(
        payload["result"]["phases"]["liquid"]["fractions"],
        fixture["result"]["phases"]["liquid"]["fractions"],
        rtol=1e-9,
        atol=1e-12,
    )
    assert np.allclose(
        payload["result"]["phases"]["vapor"]["fractions"],
        fixture["result"]["phases"]["vapor"]["fractions"],
        rtol=1e-9,
        atol=1e-12,
    )

    for key, expected in fixture["diagnostics"].items():
        actual = payload["diagnostics"][key]
        if isinstance(expected, float):
            assert actual == pytest.approx(expected, rel=1e-9, abs=1e-12)
        else:
            assert actual == expected


def test_cli_tp_flash_text_output_contains_core_fields() -> None:
    proc = _run_cli(
        "tp-flash",
        "--components",
        "Methane,Ethane,Propane",
        "--z",
        "0.5,0.3,0.2",
        "--temperature-k",
        "240",
        "--pressure-pa",
        "3000000",
    )
    assert proc.returncode == 0, proc.stderr
    assert "TP flash (chemthermo CLI)" in proc.stdout
    assert "Vapor fraction beta" in proc.stdout
    assert "Diagnostics:" in proc.stdout


def test_cli_tp_flash_validation_error_returns_exit_code_1() -> None:
    proc = _run_cli(
        "tp-flash",
        "--components",
        "Methane,Ethane,Propane",
        "--z",
        "0.5,0.5",
        "--temperature-k",
        "240",
        "--pressure-pa",
        "3000000",
    )
    assert proc.returncode == 1
    assert "same number" in proc.stderr


def test_cli_tp_flash_usage_error_returns_exit_code_2() -> None:
    proc = _run_cli("tp-flash")
    assert proc.returncode == 2
    assert "usage:" in proc.stderr


def test_cli_tp_flash_nonconvergence_returns_exit_code_3() -> None:
    proc = _run_cli(
        "tp-flash",
        "--components",
        "Methane,Ethane",
        "--z",
        "0.5,0.5",
        "--temperature-k",
        "240",
        "--pressure-pa",
        "3000000",
        "--max-iter",
        "1",
        "--tol",
        "1e-12",
    )
    assert proc.returncode == 3
    assert "did not converge" in proc.stderr
