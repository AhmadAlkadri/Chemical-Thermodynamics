from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_PATH = REPO_ROOT / "tests" / "fixtures" / "cli" / "tp_flash_v1.json"
GAMMA_FIXTURE_PATH = REPO_ROOT / "tests" / "fixtures" / "cli" / "tp_flash_gamma_phi_v1.json"


def _run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "chemthermo", *args],
        cwd=str(REPO_ROOT),
        check=False,
        text=True,
        capture_output=True,
    )


def test_cli_tp_flash_json_output_matches_fixture_contract() -> None:
    """Phi-phi CLI contract against the pre-`flash-auto-phase-detection` fixture.

    The fixture's `result` block is deliberately still the pre-slice one. Its
    numbers were produced by the Wilson-seeded iteration; the tangent-plane
    path reaches the *same* equilibrium from a different starting K, so the two
    differ only by the width of the K-update tolerance (`tol = 1e-8`). Achieved
    agreement at this state: 1.11e-9 relative on the vapor fraction and 1.8e-9
    relative on the compositions, hence `rel=1e-7` here instead of the previous
    `rel=1e-9`. The `diagnostics` block and `solver.algorithm` were regenerated
    (iteration path and new keys); see ADR-0008 and validation Case F-1.
    """
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
        fixture["result"]["vapor_fraction"], rel=1e-7, abs=1e-12
    )

    assert np.allclose(
        payload["result"]["phases"]["liquid"]["fractions"],
        fixture["result"]["phases"]["liquid"]["fractions"],
        rtol=1e-7,
        atol=1e-12,
    )
    assert np.allclose(
        payload["result"]["phases"]["vapor"]["fractions"],
        fixture["result"]["phases"]["vapor"]["fractions"],
        rtol=1e-7,
        atol=1e-12,
    )

    for key, expected in fixture["diagnostics"].items():
        actual = payload["diagnostics"][key]
        if isinstance(expected, float):
            assert actual == pytest.approx(expected, rel=1e-9, abs=1e-12)
        else:
            assert actual == expected


def test_cli_tp_flash_json_diagnostics_carry_the_phase_detection_keys() -> None:
    """New diagnostics keys serialize through the CLI without a schema bump.

    ADR-0003/ADR-0004 require `cli_schema_version` to be bumped for a
    *structural* break. `diagnostics` is a free-form mapping whose keys are
    documented as implementation details, so adding keys inside it removes
    nothing and changes no type: the version stays 1 (ADR-0008).
    """
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

    assert payload["cli_schema_version"] == 1
    diagnostics = payload["diagnostics"]
    for key in (
        "phase_detection",
        "stability_status",
        "tpd_min",
        "feed_branch",
        "k_seed",
        "mass_balance_residual",
        "fugacity_residual",
        "delta_g_split_rt",
        "stability_trials",
    ):
        assert key in diagnostics, key
    assert diagnostics["phase_detection"] == "tangent-plane"
    assert diagnostics["stability_status"] == "unstable"
    assert diagnostics["k_seed"] == "stability"
    # Round-trips through json.dumps already (the CLI printed it), so every
    # value is a JSON scalar.
    assert json.loads(json.dumps(diagnostics)) == diagnostics


def test_cli_tp_flash_json_diagnostics_carry_the_post_split_keys() -> None:
    """The post-split stability keys serialize through the CLI (ADR-0009).

    `cli_schema_version` stays 1 for the same reason as in ADR-0008: these keys
    live inside the free-form `diagnostics` mapping, which removes nothing and
    changes no type. The numeric per-phase tangent-plane distances are asserted
    by magnitude rather than pinned in the fixture: they are near-cancellation
    quantities of order 1e-10 whose last digits are not a contract.
    """
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

    assert payload["cli_schema_version"] == 1
    diagnostics = payload["diagnostics"]
    assert diagnostics["post_split_checked"] is True
    assert diagnostics["post_split_stable"] is True
    assert diagnostics["post_split_status"] == "stable"
    assert diagnostics["phase_stability_liquid"] == "stable"
    assert diagnostics["phase_stability_vapor"] == "stable"
    for key in (
        "post_split_tpd_min",
        "phase_stability_tpd_min_liquid",
        "phase_stability_tpd_min_vapor",
    ):
        assert abs(float(diagnostics[key])) < 1e-8, key
    assert json.loads(json.dumps(diagnostics)) == diagnostics


def test_cli_tp_flash_gamma_phi_json_output_matches_fixture_contract() -> None:
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
        "--flash-mode",
        "gamma-phi",
        "--format",
        "json",
    )
    assert proc.returncode == 0, proc.stderr

    payload = json.loads(proc.stdout)
    fixture = json.loads(GAMMA_FIXTURE_PATH.read_text(encoding="utf-8"))

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


def test_cli_tp_flash_gamma_phi_text_output_contains_core_fields() -> None:
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
        "--flash-mode",
        "gamma-phi",
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


def test_cli_tp_flash_gamma_phi_missing_pair_returns_exit_code_1() -> None:
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
        "--flash-mode",
        "gamma-phi",
    )
    assert proc.returncode == 1
    assert "Missing NRTL parameters" in proc.stderr


def test_cli_tp_flash_usage_error_returns_exit_code_2() -> None:
    proc = _run_cli("tp-flash")
    assert proc.returncode == 2
    assert "usage:" in proc.stderr


def test_cli_tp_flash_invalid_mode_returns_exit_code_2() -> None:
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
        "--flash-mode",
        "bad-mode",
    )
    assert proc.returncode == 2
    assert "invalid choice" in proc.stderr


def test_cli_tp_flash_nonconvergence_returns_exit_code_3() -> None:
    """A flash that cannot converge exits 3 and says so on stderr.

    The starved state is run in ``gamma-phi`` mode because that mode goes
    through the legacy Wilson-heuristic path, which ADR-0016 deliberately left
    alone. The same starved ``phi-phi`` state is now *rescued* by the
    second-order stage (see the companion test below), so it is no longer a
    non-convergence example; the exit-code contract itself is unchanged.
    """
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
        "--flash-mode",
        "gamma-phi",
        "--max-iter",
        "1",
        "--tol",
        "1e-12",
    )
    assert proc.returncode == 3
    assert "did not converge" in proc.stderr


def test_cli_tp_flash_phi_phi_survives_a_starved_iteration_budget() -> None:
    """ADR-0016: the phi-phi second-order stage finishes a starved split."""
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
        "--format",
        "json",
    )
    assert proc.returncode == 0, proc.stderr
    payload = json.loads(proc.stdout)
    assert payload["diagnostics"]["converged_stage"] == "second-order"
    assert payload["diagnostics"]["ssi_iterations"] == 1
