"""The benchmark harness runs, and its comparison refuses a changed answer (ADR-0023).

Deliberately small. The harness is a maintainer's instrument, not a library
feature, so what is worth pinning is that it *works* and that its one hard rule
holds: a comparison that finds two records disagreeing on a case's
``result_hash`` reports a failure rather than a speedup. The committed records
under ``benchmarks/`` are the measurement; this file is the guard on the
instrument.

The timing case chosen here is the cheapest in the workload
(``pr-stability-ternary``, a few milliseconds) and it is run with one repeat,
so this module costs a fraction of a second.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from chemthermo.bench import CASES, CASES_BY_ID, compare, main, run_case
from chemthermo.bench._record import SCHEMA, StateOutcome, result_hash

CHEAP_CASE = "pr-stability-ternary"


def test_every_case_id_is_unique_and_addressable() -> None:
    ids = [case.id for case in CASES]
    assert len(ids) == len(set(ids))
    assert set(CASES_BY_ID) == set(ids)


def test_one_case_runs_and_records_the_campaign_fields() -> None:
    entry = run_case(CASES_BY_ID[CHEAP_CASE], repeats=1)

    assert entry["status"] == "ok"
    assert entry["id"] == CHEAP_CASE
    assert entry["derivative_mode"] == "analytic"
    assert entry["model"] == "Peng-Robinson"
    assert entry["components"] == ["Methane", "Ethane", "Propane"]
    # Convergence criteria, initialization and iteration counts are all
    # present: a wall time without them is not a reproducible record.
    assert entry["settings"]["stability_tol"] == pytest.approx(1e-10)
    assert entry["wall_time_s"]["median"] > 0.0
    assert entry["wall_time_s"]["min"] <= entry["wall_time_s"]["median"]
    assert entry["wall_time_s"]["max"] >= entry["wall_time_s"]["median"]
    assert entry["peak_memory_bytes"] > 0
    assert entry["result_hash"].startswith("sha256:")

    outcome = entry["outcomes"][0]
    assert outcome["temperature_K"] == 240.0
    assert outcome["pressure_Pa"] == 3.0e6
    assert outcome["status"] == "ok"
    assert outcome["iterations"]["stability_trials"] >= 1
    assert outcome["initialization"]["phase_detection"] == "tangent-plane"

    # The record is JSON, all the way down.
    json.dumps(entry)


def test_the_result_hash_moves_with_the_answer_and_not_with_the_clock() -> None:
    base = StateOutcome(
        temperature_K=300.0,
        pressure_Pa=1.0e5,
        composition=(0.4, 0.6),
        phases=("liquid", "vapor"),
        phase_compositions=((0.7, 0.3), (0.1, 0.9)),
        phase_fractions=(0.25, 0.75),
    )
    same_answer_other_timing = StateOutcome(
        temperature_K=300.0,
        pressure_Pa=1.0e5,
        composition=(0.4, 0.6),
        phases=("liquid", "vapor"),
        phase_compositions=((0.7, 0.3), (0.1, 0.9)),
        phase_fractions=(0.25, 0.75),
        iterations={"iterations": 99},
    )
    moved = StateOutcome(
        temperature_K=300.0,
        pressure_Pa=1.0e5,
        composition=(0.4, 0.6),
        phases=("liquid", "vapor"),
        phase_compositions=((0.7, 0.3), (0.1, 0.9)),
        phase_fractions=(0.25 + 1e-9, 0.75),
    )

    assert result_hash([base]) == result_hash([same_answer_other_timing])
    assert result_hash([base]) != result_hash([moved])


def _record(case_id: str, *, median: float, digest: str) -> dict:
    return {
        "schema": SCHEMA,
        "cases": [
            {
                "id": case_id,
                "status": "ok",
                "wall_time_s": {"median": median, "min": median, "max": median},
                "result_hash": digest,
            }
        ],
    }


def test_compare_reports_a_speedup_when_the_answer_held() -> None:
    before = _record("x", median=2.0, digest="sha256:abc")
    after = _record("x", median=1.0, digest="sha256:abc")

    lines, identical = compare(before, after)

    assert identical is True
    assert any("2.00x" in line for line in lines)
    assert any("identical" in line for line in lines)


def test_compare_fails_when_a_result_hash_changed() -> None:
    before = _record("x", median=2.0, digest="sha256:abc")
    after = _record("x", median=1.0, digest="sha256:def")

    lines, identical = compare(before, after)

    assert identical is False
    assert any("CHANGED" in line for line in lines)
    assert any("not an optimization" in line for line in lines)


def test_the_cli_compare_exit_status_follows_the_hashes(tmp_path: Path) -> None:
    before_path = tmp_path / "before.json"
    after_path = tmp_path / "after.json"
    before = _record("x", median=2.0, digest="sha256:abc")
    before_path.write_text(json.dumps(before), encoding="utf-8")

    after_path.write_text(json.dumps(_record("x", median=1.0, digest="sha256:abc")), "utf-8")
    assert main(["--compare", str(before_path), str(after_path)]) == 0

    changed = copy.deepcopy(before)
    changed["cases"][0]["result_hash"] = "sha256:def"
    after_path.write_text(json.dumps(changed), encoding="utf-8")
    assert main(["--compare", str(before_path), str(after_path)]) == 1


def test_the_cli_writes_a_record_with_an_environment_and_a_commit(tmp_path: Path) -> None:
    out = tmp_path / "record.json"

    assert main(["--repeats", "1", "--case", CHEAP_CASE, "--out", str(out)]) == 0

    record = json.loads(out.read_text(encoding="utf-8"))
    assert record["schema"] == SCHEMA
    assert record["environment"]["python"]
    assert record["environment"]["numpy"]
    assert "sha" in record["git"]
    assert [entry["id"] for entry in record["cases"]] == [CHEAP_CASE]
