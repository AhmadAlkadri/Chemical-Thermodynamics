"""The off-platform comparison of ADR-0032 is not vacuous.

`tests/_capture_identity.py` compares bit for bit on the capture platform and
with a stated bound elsewhere. These checks run its bounded path directly, on
every platform, so what it lets through and what it catches is pinned wherever
the suite runs.
"""

from __future__ import annotations

import math

from _capture_identity import capture_deviations

PINNED = {
    "phase_names": ["liquid", "vapor"],
    "phases": {"liquid": [0.25, 0.75], "vapor": [0.9, 0.1]},
    "vapor_fraction": 0.4,
    "diagnostics": {"iterations": 12, "converged": True, "tpd_min": -0.5, "residual": 1e-16},
}


def _moved(**changes: object) -> dict[str, object]:
    moved = {key: value for key, value in PINNED.items()}
    moved["diagnostics"] = dict(PINNED["diagnostics"])  # type: ignore[arg-type]
    for key, value in changes.items():
        if key in moved:
            moved[key] = value
        else:
            moved["diagnostics"][key] = value  # type: ignore[index]
    return moved


def test_identical_values_pass_with_zero_deviation() -> None:
    violations, worst = capture_deviations(PINNED, _moved(), rtol=1e-12, atol=1e-12)
    assert violations == []
    assert worst == 0.0


def test_a_last_bit_move_passes() -> None:
    violations, worst = capture_deviations(
        PINNED, _moved(tpd_min=math.nextafter(-0.5, 0.0), residual=-3e-16), rtol=1e-12, atol=1e-12
    )
    assert violations == []
    assert 0.0 < worst < 1.0


def test_a_move_beyond_the_bound_fails() -> None:
    violations, _ = capture_deviations(PINNED, _moved(tpd_min=-0.5 + 1e-11), rtol=1e-12, atol=1e-12)
    assert len(violations) == 1
    assert "tpd_min" in violations[0]


def test_every_discrete_field_is_exact() -> None:
    for changes in (
        {"phase_names": ["vapor", "liquid"]},
        {"iterations": 13},
        {"converged": False},
        {"converged": 1},
        {"vapor_fraction": None},
        {"extra_key": 0.0},
        {"phases": {"liquid": [0.25, 0.75]}},
        {"phases": {"liquid": [0.25, 0.75, 0.0], "vapor": [0.9, 0.1]}},
    ):
        violations, _ = capture_deviations(PINNED, _moved(**changes), rtol=1.0, atol=1.0)
        assert violations, changes


def test_non_finite_floats_are_exact() -> None:
    violations, _ = capture_deviations({"x": math.inf}, {"x": 1e308}, rtol=1.0, atol=1.0)
    assert violations
    assert capture_deviations({"x": math.nan}, {"x": math.nan}, rtol=0.0, atol=0.0)[0] == []


def test_key_order_does_not_matter_but_key_sets_do() -> None:
    reordered = {"b": 2.0, "a": 1.0}
    assert capture_deviations({"a": 1.0, "b": 2.0}, reordered, rtol=0.0, atol=0.0)[0] == []
