"""The runtime gate and the bounded comparison of ADR-0032 are not vacuous.

`tests/_capture_identity.py` compares bit for bit on the capture runtime and
with a stated bound elsewhere. These checks pin, on every runtime, which
runtimes get the exact comparison and what the bounded path lets through and
what it catches.
"""

from __future__ import annotations

import functools
import math
import operator
import platform
import sys

import pytest
from _capture_identity import (
    CAPTURE_RUNTIME,
    Runtime,
    capture_deviations,
    current_runtime,
    is_capture_runtime,
    on_capture_runtime,
)

# --------------------------------------------------------------------------
# Which runtime is the capture runtime
# --------------------------------------------------------------------------


def test_the_capture_runtime_is_macos_arm64_cpython_3_11() -> None:
    assert CAPTURE_RUNTIME == Runtime("Darwin", "arm64", "CPython", (3, 11))
    assert is_capture_runtime(Runtime("Darwin", "arm64", "CPython", (3, 11)))


@pytest.mark.parametrize(
    "runtime",
    [
        # Same machine, newer interpreter: CPython 3.12 changed float sum().
        Runtime("Darwin", "arm64", "CPython", (3, 12)),
        Runtime("Darwin", "arm64", "CPython", (3, 13)),
        Runtime("Darwin", "arm64", "CPython", (3, 10)),
        # Same interpreter version, another implementation.
        Runtime("Darwin", "arm64", "PyPy", (3, 11)),
        # Same interpreter, another CPU or OS.
        Runtime("Darwin", "x86_64", "CPython", (3, 11)),
        Runtime("Linux", "x86_64", "CPython", (3, 11)),
        Runtime("Linux", "aarch64", "CPython", (3, 11)),
        Runtime("Windows", "AMD64", "CPython", (3, 11)),
    ],
    ids=str,
)
def test_every_other_runtime_gets_the_bounded_comparison(runtime: Runtime) -> None:
    assert not is_capture_runtime(runtime)


def test_the_current_runtime_is_read_from_the_interpreter() -> None:
    runtime = current_runtime()
    assert runtime.system == platform.system()
    assert runtime.machine == platform.machine()
    assert runtime.implementation == platform.python_implementation()
    assert runtime.python == sys.version_info[:2]
    assert on_capture_runtime() is is_capture_runtime(runtime)


def test_the_interpreter_change_the_gate_tracks() -> None:
    """Why the minor version is in the gate: from CPython 3.12 the built-in
    ``sum()`` of floats is compensated, so it no longer equals a left fold.
    Captures made under 3.11 used the left fold."""
    values = [1.0, 1e100, 1.0, -1e100]
    left_fold = functools.reduce(operator.add, values, 0)
    assert left_fold == 0.0
    compensated = platform.python_implementation() == "CPython" and sys.version_info >= (3, 12)
    assert sum(values) == (2.0 if compensated else 0.0)


# --------------------------------------------------------------------------
# The bounded comparison used off the capture runtime
# --------------------------------------------------------------------------

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
