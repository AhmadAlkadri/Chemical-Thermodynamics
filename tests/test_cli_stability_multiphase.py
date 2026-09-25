"""CLI contract of ADR-0033: `stability-tp`, `tp-flash --eos`, `--max-phases`.

Golden fixtures in `tests/fixtures/cli/` were written by the CLI itself at the
commit that introduced them. They are compared the way ADR-0032 compares any
capture from one machine: every discrete field exactly (keys, names,
statuses, branches, counts, settings) and floats to a stated bound -
`1e-10` absolute plus relative here, two orders inside the flash tolerance
(`tol = 1e-8`) these answers converged to, because this is a contract test,
not a bit-identity guard. Diagnostics keys that name *which* of several tied
trials was reported (`minimizing_trial*`, ADR-0032) are excluded from the
exact comparison; the rest of the diagnostics mapping is compared.

The ADR-0003/0004 fixtures and their tests (`tests/test_cli_tp_flash.py`) are
untouched: existing invocations produce byte-identical output (ADR-0033
decision 1, checked when this contract was introduced).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
from _capture_identity import capture_deviations

import chemthermo as ct
from chemthermo import cli

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = REPO_ROOT / "tests" / "fixtures" / "cli"

#: Diagnostics that name a tied trial or count its iterations (ADR-0032).
TIE_SENSITIVE_PREFIX = "minimizing_trial"

GOLDEN = [
    (
        "stability_tp_v1_pr_stable.json",
        ["stability-tp", "--components", "Methane,Ethane", "--z", "0.5,0.5",
         "--temperature-k", "300", "--pressure-pa", "100000"],
    ),
    (
        "stability_tp_v1_pr_unstable.json",
        ["stability-tp", "--components", "Methane,n-Hexane", "--z", "0.5,0.5",
         "--temperature-k", "300", "--pressure-pa", "2000000"],
    ),
    (
        "stability_tp_v1_pcsaft_unstable.json",
        ["stability-tp", "--components", "Methane,n-Hexane", "--z", "0.5,0.5",
         "--temperature-k", "300", "--pressure-pa", "2000000", "--eos", "pc-saft"],
    ),
    (
        "tp_flash_v1_pcsaft_single_phase.json",
        ["tp-flash", "--components", "Methane,Ethane", "--z", "0.5,0.5",
         "--temperature-k", "300", "--pressure-pa", "100000", "--eos", "pc-saft"],
    ),
    (
        "tp_flash_v1_pcsaft_vle.json",
        ["tp-flash", "--components", "Methane,n-Hexane", "--z", "0.5,0.5",
         "--temperature-k", "300", "--pressure-pa", "2000000", "--eos", "pc-saft"],
    ),
    (
        "tp_flash_v1_pr_three_liquid.json",
        ["tp-flash", "--components", "Water,Ethanol,n-Hexane", "--z", "0.2,0.4,0.4",
         "--temperature-k", "280", "--pressure-pa", "101325", "--max-phases", "3"],
    ),
]  # fmt: skip


def _run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "chemthermo", *args],
        cwd=str(REPO_ROOT),
        check=False,
        text=True,
        capture_output=True,
    )


def _json(*args: str) -> dict[str, Any]:
    proc = _run_cli(*args, "--format", "json")
    assert proc.returncode == 0, proc.stderr
    payload: dict[str, Any] = json.loads(proc.stdout)
    return payload


def _without_tie_keys(diagnostics: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in diagnostics.items() if not k.startswith(TIE_SENSITIVE_PREFIX)}


@pytest.mark.parametrize(("fixture_name", "args"), GOLDEN, ids=[name for name, _ in GOLDEN])
def test_golden_output(fixture_name: str, args: list[str]) -> None:
    payload = _json(*args)
    fixture = json.loads((FIXTURES / fixture_name).read_text(encoding="utf-8"))

    assert set(payload) == set(fixture)
    for key in ("cli_schema_version", "command", "inputs", "solver"):
        assert payload[key] == fixture[key], key
    violations, _ = capture_deviations(
        fixture["result"], payload["result"], rtol=1e-10, atol=1e-10, path="result"
    )
    assert not violations, violations
    assert set(payload["diagnostics"]) == set(fixture["diagnostics"])
    violations, _ = capture_deviations(
        _without_tie_keys(fixture["diagnostics"]),
        _without_tie_keys(payload["diagnostics"]),
        rtol=1e-10,
        atol=1e-10,
        path="diagnostics",
    )
    assert not violations, violations


# ---------------------------------------------------------------------------
# stability-tp
# ---------------------------------------------------------------------------


def test_stability_scope_is_stated_in_every_payload_and_in_help() -> None:
    payload = _json(*GOLDEN[0][1])
    assert payload["result"]["status"] == "stable"
    assert payload["result"]["stable"] is True
    assert payload["result"]["stability_scope"] == "bounded-trial-set"
    assert payload["result"]["trial_composition"] is None

    help_text = _run_cli("stability-tp", "--help").stdout
    assert "not a global proof" in " ".join(help_text.split())

    text = _run_cli(*GOLDEN[0][1]).stdout
    assert "Status: stable" in text
    assert "not a global proof" in text


def test_an_unstable_feed_reports_its_incipient_phase() -> None:
    payload = _json(*GOLDEN[1][1])
    result = payload["result"]
    assert result["status"] == "unstable"
    assert result["stable"] is False
    assert result["tpd_min"] < -1e-8
    assert result["phase_branch"] == "vapor"
    assert result["feed_branch"] == "liquid"
    w = result["trial_composition"]
    assert sum(w) == pytest.approx(1.0, abs=1e-12)
    # K = w / z, feed -> incipient (StabilityResult.k_values).
    for k_i, w_i, z_i in zip(result["k_values"], w, payload["inputs"]["z_mole"]):
        assert k_i == pytest.approx(w_i / z_i, rel=1e-12)
    assert result["non_trivial_trial_count"] >= 1


def test_pc_saft_and_peng_robinson_agree_on_the_verdict_but_not_the_number() -> None:
    pr = _json(*GOLDEN[1][1])["result"]
    saft = _json(*GOLDEN[2][1])["result"]
    assert (pr["status"], saft["status"]) == ("unstable", "unstable")
    assert abs(pr["tpd_min"] - saft["tpd_min"]) > 1e-3


def test_the_cli_and_the_library_report_the_same_numbers() -> None:
    payload = _json(*GOLDEN[1][1])
    direct = ct.stability_tp(
        ct.Mixture.from_database(["Methane", "n-Hexane"], [0.5, 0.5]),
        temperature_K=300.0,
        pressure_Pa=2.0e6,
        eos=ct.PengRobinsonEOS(),
    )
    assert payload["result"]["tpd_min"] == direct.tpd_min
    assert payload["result"]["trial_composition"] == list(direct.trial_composition or ())


def test_an_inconclusive_verdict_exits_3_and_still_prints_the_payload(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    def no_verdict(mixture: ct.Mixture, **kwargs: Any) -> ct.StabilityResult:
        return ct.StabilityResult(
            temperature_K=kwargs["temperature_K"],
            pressure_Pa=kwargs["pressure_Pa"],
            feed_composition=tuple(mixture.composition.fractions),
            stable=False,
            status="inconclusive",
            tpd_min=0.0,
            diagnostics={"tm_at_stationary_point": float("-inf")},
        )

    monkeypatch.setattr(cli, "stability_tp", no_verdict)
    code = cli.main([*GOLDEN[0][1], "--format", "json"])
    assert code == 3
    payload = json.loads(capsys.readouterr().out)  # valid JSON despite the -inf
    assert payload["result"]["status"] == "inconclusive"
    assert payload["diagnostics"]["tm_at_stationary_point"] == "-inf"


def test_stability_errors_keep_the_exit_code_contract() -> None:
    base = ["stability-tp", "--z", "0.5,0.5", "--temperature-k", "300", "--pressure-pa", "1e5"]
    assert _run_cli(*base, "--components", "Methane,Unobtainium").returncode == 1
    assert _run_cli(*base, "--components", "Methane").returncode == 1
    assert _run_cli(*base, "--components", "Methane,Ethane", "--eos", "srk").returncode == 2
    assert _run_cli("stability-tp", "--components", "Methane").returncode == 2


# ---------------------------------------------------------------------------
# tp-flash --eos / --max-phases
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fixture_name", [name for name, _ in GOLDEN if "tp_flash" in name])
def test_flash_answers_conserve_mass_from_the_json_alone(fixture_name: str) -> None:
    args = dict(GOLDEN)[fixture_name]
    payload = _json(*args)
    result = payload["result"]
    z = payload["inputs"]["z_mole"]
    fractions = result["phase_fractions"]
    assert list(fractions) == result["phase_names"] == list(result["phases"])
    assert sum(fractions.values()) == pytest.approx(1.0, abs=1e-12)
    assert all(0.0 < beta <= 1.0 for beta in fractions.values())
    for i, z_i in enumerate(z):
        recombined = sum(
            fractions[name] * result["phases"][name]["fractions"][i] for name in fractions
        )
        assert recombined == pytest.approx(z_i, abs=1e-10)
    for name in fractions:
        assert sum(result["phases"][name]["fractions"]) == pytest.approx(1.0, abs=1e-12)
    if "vapor" in fractions:
        assert result["vapor_fraction"] == fractions["vapor"]
    else:
        assert result["vapor_fraction"] is None


def test_three_liquids_come_back_under_the_v1_layout() -> None:
    payload = _json(*dict(GOLDEN)["tp_flash_v1_pr_three_liquid.json"])
    assert payload["cli_schema_version"] == 1
    assert payload["result"]["phase_names"] == ["liquid1", "liquid2", "liquid3"]
    assert payload["solver"]["settings"]["max_phases"] == 3
    assert payload["diagnostics"]["post_split_status"] == "stable"


def test_max_phases_is_written_only_when_given() -> None:
    args = ["tp-flash", "--components", "Methane,Ethane,Propane", "--z", "0.5,0.3,0.2",
            "--temperature-k", "240", "--pressure-pa", "3000000"]  # fmt: skip
    assert "max_phases" not in _json(*args)["solver"]["settings"]
    assert _json(*args, "--max-phases", "2")["solver"]["settings"]["max_phases"] == 2


def test_a_phase_budget_below_the_answer_is_a_solver_refusal() -> None:
    args = dict(GOLDEN)["tp_flash_v1_pr_three_liquid.json"][:-2]
    proc = _run_cli(*args, "--max-phases", "2")
    assert proc.returncode == 3
    assert "max_phases" in proc.stderr


def test_pc_saft_is_labelled_and_restricted_to_phi_phi() -> None:
    payload = _json(*dict(GOLDEN)["tp_flash_v1_pcsaft_vle.json"])
    assert payload["solver"]["eos"] == "pc_saft"
    assert payload["solver"]["method"] == "phi-phi"
    assert payload["result"]["phase_names"] == ["liquid", "vapor"]

    base = ["tp-flash", "--components", "Methane,Ethane", "--z", "0.5,0.5",
            "--temperature-k", "240", "--pressure-pa", "3000000"]  # fmt: skip
    assert _run_cli(*base, "--eos", "pc-saft", "--flash-mode", "gamma-phi").returncode == 2
    assert _run_cli(*base, "--max-phases", "0").returncode == 2
    assert _run_cli(*base, "--max-phases", "two").returncode == 2
