"""Reproducible benchmark records for the reference flash paths (ADR-0023).

**Internal package.** Nothing here is re-exported from ``chemthermo`` and
nothing here is part of the public API (ADR-0001): it is a maintainer's
instrument, not a library feature.

What it is for
--------------
An optimization in this repository is accepted on two pieces of evidence: the
answer did not move, and the clock did. This package produces both in one
artefact. ``python -m chemthermo.bench --out record.json`` runs the fixed
workload of :mod:`chemthermo.bench._cases` and writes, per case:

- the model, components, state ``(T, P)`` and overall composition;
- the phase count obtained, and the phase names, compositions and fractions;
- the convergence criteria in force (every ``FlashSettings`` /
  ``StabilitySettings`` field that can change an answer);
- the initialization (``k_seed``, ``phase_detection``) the solver reported;
- the iteration counts the solver reported;
- the derivative mode (``"analytic"`` throughout - chemthermo differentiates no
  model by finite difference);
- the **median** wall time over ``--repeats`` timed runs after a warm-up run,
  with the minimum and maximum alongside;
- the peak allocation of one separate ``tracemalloc``-instrumented run (never
  the timed one: the tracer roughly triples the wall time it observes);
- the hardware and interpreter the numbers were taken on, and the git commit;
- a ``result_hash`` over the accepted thermodynamic answer, and the failure
  status when a state was refused.

Comparing two records
---------------------
``python -m chemthermo.bench --compare before.json after.json`` prints the
per-case median before, after and ratio, and **fails** (exit status 1) if any
case's ``result_hash`` differs. That is the whole acceptance rule for a
performance change here: same hash, better clock. See ``benchmarks/README.md``.

Why the median, and why a warm-up
---------------------------------
The first call through any of these paths pays for databank reads, parameter
resolution and numpy's own first-touch costs, none of which is what an
optimization moves; the workload is therefore *prepared* outside the timed
region and run once untimed before the repeats begin. The median rather than
the minimum because a macOS laptop's scheduler produces occasional slow
repeats and no fast ones, so the minimum flatters and the mean is dragged;
min and max are recorded so that spread stays visible rather than hidden
behind one number.
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
import tracemalloc
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from ._cases import CASES, CASES_BY_ID, DERIVATIVE_MODE, BenchCase, CaseSkipped
from ._record import SCHEMA, environment, git_state, result_hash

__all__ = [
    "CASES",
    "CASES_BY_ID",
    "BenchCase",
    "compare",
    "main",
    "run_all",
    "run_case",
]

#: Timed repeats per case when ``--repeats`` is not given.
DEFAULT_REPEATS = 5


def run_case(case: BenchCase, *, repeats: int = DEFAULT_REPEATS) -> dict[str, Any]:
    """Run one case and return its record entry.

    Args:
        case: The workload to run.
        repeats: Timed repeats after the warm-up. Must be at least one.

    Returns:
        A JSON-ready mapping; see the module docstring for the fields.
    """
    if repeats < 1:
        raise ValueError("repeats must be at least 1.")

    entry: dict[str, Any] = {
        "id": case.id,
        "description": case.description,
        "model": case.model,
        "route": case.route,
        "components": list(case.components),
        "settings": dict(case.settings),
        "derivative_mode": DERIVATIVE_MODE,
        "repeats": repeats,
    }

    try:
        payload = case.prepare()
    except CaseSkipped as exc:
        entry.update({"status": "skipped", "reason": str(exc)})
        return entry

    # Warm-up: not timed, and its outcome is what the record describes only
    # after the timed repeats have confirmed it does not change.
    outcomes = case.invoke(payload)

    timings: list[float] = []
    for _ in range(repeats):
        started = time.perf_counter()
        repeat_outcomes = case.invoke(payload)
        timings.append(time.perf_counter() - started)
        if result_hash(repeat_outcomes) != result_hash(outcomes):
            raise RuntimeError(
                f"Case {case.id!r} is not deterministic: a repeat produced a different "
                "result hash than the warm-up run."
            )

    tracemalloc.start()
    case.invoke(payload)
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    entry.update(
        {
            "status": "ok",
            "states": len(outcomes),
            "refusals": sum(1 for outcome in outcomes if outcome.status != "ok"),
            "phase_counts": [len(outcome.phases) for outcome in outcomes],
            "wall_time_s": {
                "median": statistics.median(timings),
                "min": min(timings),
                "max": max(timings),
                "samples": list(timings),
            },
            "peak_memory_bytes": int(peak),
            "result_hash": result_hash(outcomes),
            "outcomes": [outcome.as_json() for outcome in outcomes],
        }
    )
    return entry


def run_all(*, repeats: int = DEFAULT_REPEATS, only: Sequence[str] | None = None) -> dict[str, Any]:
    """Run the workload and return a complete record.

    Args:
        repeats: Timed repeats per case.
        only: Case ids to run; ``None`` runs every case in :data:`CASES`.

    Returns:
        A JSON-ready record: schema tag, timestamp, git state, environment and
        one entry per case.

    Raises:
        KeyError: If ``only`` names a case id that does not exist.
    """
    selected = CASES if only is None else tuple(CASES_BY_ID[name] for name in only)
    return {
        "schema": SCHEMA,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git": git_state(),
        "environment": environment(),
        "repeats": repeats,
        "cases": [run_case(case, repeats=repeats) for case in selected],
    }


def compare(before: Mapping[str, Any], after: Mapping[str, Any]) -> tuple[list[str], bool]:
    """Compare two records; return the report lines and whether results held.

    The second element is ``True`` only when every case present in both records
    has the same ``result_hash``. A case that is missing from one side, or that
    was skipped on one side, is reported and does not on its own make the
    comparison fail - a changed *answer* does.
    """
    before_cases = {entry["id"]: entry for entry in before.get("cases", [])}
    after_cases = {entry["id"]: entry for entry in after.get("cases", [])}
    lines: list[str] = []
    identical = True

    lines.append(f"{'case':<28} {'before / s':>12} {'after / s':>12} {'speedup':>9}  result")
    lines.append("-" * 78)
    for case_id in sorted(set(before_cases) | set(after_cases)):
        left = before_cases.get(case_id)
        right = after_cases.get(case_id)
        if left is None or right is None:
            side = "after only" if left is None else "before only"
            lines.append(f"{case_id:<28} {'-':>12} {'-':>12} {'-':>9}  {side}")
            continue
        if left.get("status") != "ok" or right.get("status") != "ok":
            state = f"{left.get('status')} -> {right.get('status')}"
            lines.append(f"{case_id:<28} {'-':>12} {'-':>12} {'-':>9}  {state}")
            continue
        left_time = float(left["wall_time_s"]["median"])
        right_time = float(right["wall_time_s"]["median"])
        ratio = left_time / right_time if right_time > 0.0 else float("inf")
        same = left["result_hash"] == right["result_hash"]
        identical = identical and same
        verdict = "identical" if same else "CHANGED"
        lines.append(f"{case_id:<28} {left_time:12.6f} {right_time:12.6f} {ratio:8.2f}x  {verdict}")

    lines.append("")
    lines.append(
        "result hashes identical across every shared case"
        if identical
        else "RESULT HASHES DIFFER: this is not an optimization"
    )
    return lines, identical


def _write(path: Path, record: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=False)
        handle.write("\n")


def _summary(record: Mapping[str, Any]) -> list[str]:
    lines = [f"{'case':<28} {'median / s':>12} {'peak / MiB':>11}  states  result"]
    lines.append("-" * 78)
    for entry in record["cases"]:
        if entry.get("status") != "ok":
            lines.append(f"{entry['id']:<28} {'-':>12} {'-':>11}  {'-':>6}  {entry['status']}")
            continue
        lines.append(
            f"{entry['id']:<28} {entry['wall_time_s']['median']:12.6f} "
            f"{entry['peak_memory_bytes'] / 1048576.0:11.2f}  {entry['states']:6d}  "
            f"{entry['result_hash'][7:19]}"
        )
    return lines


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point; see ``python -m chemthermo.bench --help``."""
    parser = argparse.ArgumentParser(
        prog="python -m chemthermo.bench",
        description="Run the chemthermo benchmark workload, or compare two records.",
    )
    parser.add_argument("--out", type=Path, default=None, help="write the record to this path")
    parser.add_argument(
        "--repeats",
        type=int,
        default=DEFAULT_REPEATS,
        help=f"timed repeats per case after a warm-up run (default {DEFAULT_REPEATS})",
    )
    parser.add_argument(
        "--case",
        action="append",
        dest="cases",
        default=None,
        help="run only this case id; repeatable",
    )
    parser.add_argument("--list", action="store_true", help="list the case ids and exit")
    parser.add_argument(
        "--compare",
        nargs=2,
        metavar=("BEFORE", "AFTER"),
        type=Path,
        default=None,
        help="compare two records instead of running; exits 1 if any result hash differs",
    )
    args = parser.parse_args(argv)

    if args.list:
        for case in CASES:
            print(f"{case.id:<28} {case.description}")
        return 0

    if args.compare is not None:
        before_path, after_path = args.compare
        with before_path.open("r", encoding="utf-8") as handle:
            before = json.load(handle)
        with after_path.open("r", encoding="utf-8") as handle:
            after = json.load(handle)
        lines, identical = compare(before, after)
        print("\n".join(lines))
        return 0 if identical else 1

    unknown = [name for name in (args.cases or ()) if name not in CASES_BY_ID]
    if unknown:
        parser.error(f"unknown case id(s): {', '.join(unknown)}")

    record = run_all(repeats=args.repeats, only=args.cases)
    print("\n".join(_summary(record)))
    if args.out is not None:
        _write(args.out, record)
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":  # pragma: no cover - exercised through __main__.py
    sys.exit(main())
