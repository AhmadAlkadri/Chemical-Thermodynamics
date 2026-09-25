"""Platform-aware comparison against floats captured on one machine (ADR-0032).

Several guards pin values that were *captured* on the development machine
(macOS arm64) and stored as literals or JSON. Their job is to prove that a
refactor or a speed-up changed no arithmetic. On the capture platform they
still do exactly that: every float is compared with ``==``.

On any other platform the same inputs give floats that differ in the last few
bits - ``exp``/``log``/``pow`` and numpy's SIMD kernels are not bit-reproducible
across CPUs or C libraries, and two Linux x86_64 hosts have been measured to
disagree with each other as well as with macOS (ledger Case P-11,
"cross-platform"). There a float passes when

    |actual - expected| <= atol + rtol * max(|actual|, |expected|)

with the bound stated by the caller and justified in ADR-0032, and **every
discrete field is still exact**: dict keys, list lengths, strings, ints, bools,
``None``, and non-finite floats. A phase name, a verdict, a status, an
iteration count or a diagnostics key that moves fails on every platform.

In-process A/B identity guards (the same state computed twice, two ways, in
one run) do not use this module: they are platform-independent and stay exact
everywhere.
"""

from __future__ import annotations

import math
import numbers
import platform
from typing import Any

#: Where the pinned fixtures were captured. Exact comparison applies here.
#: Reproduced exactly on 2026-09-25 with CPython 3.11.6 and numpy 2.4.6
#: (``.agents/handoffs/cloud-continuation.md`` section 7); the numpy version is
#: informational, not part of the gate.
CAPTURE_SYSTEM = "Darwin"
CAPTURE_MACHINE = "arm64"


def on_capture_platform() -> bool:
    """True where the pinned floats must be reproduced bit for bit."""
    return platform.system() == CAPTURE_SYSTEM and platform.machine() == CAPTURE_MACHINE


def _kind(value: Any) -> str:
    """What ``==`` would compare: a bool is not an int, a list is a tuple."""
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, numbers.Integral):
        return "int"
    if isinstance(value, float):
        return "float"
    if isinstance(value, dict):
        return "dict"
    if isinstance(value, (list, tuple)):
        return "sequence"
    return type(value).__name__


def capture_deviations(
    expected: Any, actual: Any, *, rtol: float, atol: float, path: str = ""
) -> tuple[list[str], float]:
    """Walk two JSON-shaped values; return (violations, worst float deviation).

    A violation is any discrete mismatch or any float outside the bound. The
    worst deviation is the largest ``|actual - expected| / (atol + rtol * scale)``
    seen over all finite float pairs, so ``<= 1`` means inside the bound; it is
    reported so a failure (or a curious reader) sees how much headroom is left.
    """
    violations: list[str] = []
    worst = 0.0

    def walk(e: Any, a: Any, where: str) -> None:
        nonlocal worst
        if isinstance(e, float) and isinstance(a, float):
            if not (math.isfinite(e) and math.isfinite(a)):
                if not (e == a or (math.isnan(e) and math.isnan(a))):
                    violations.append(f"{where}: non-finite {e!r} != {a!r}")
                return
            allowed = atol + rtol * max(abs(e), abs(a))
            deviation = abs(a - e)
            if allowed > 0.0:
                worst = max(worst, deviation / allowed)
            if deviation > allowed:
                violations.append(f"{where}: {a!r} vs pinned {e!r} (|diff| {deviation:.3g})")
            return
        if _kind(e) != _kind(a):
            violations.append(f"{where}: type {type(a).__name__} vs pinned {type(e).__name__}")
            return
        if isinstance(e, dict):
            if set(e) != set(a):
                violations.append(f"{where}: keys {sorted(a)} vs pinned {sorted(e)}")
                return
            for key in e:
                walk(e[key], a[key], f"{where}/{key}")
            return
        if isinstance(e, (list, tuple)):
            if len(e) != len(a):
                violations.append(f"{where}: length {len(a)} vs pinned {len(e)}")
                return
            for index, (ei, ai) in enumerate(zip(e, a)):
                walk(ei, ai, f"{where}/{index}")
            return
        if e != a:
            violations.append(f"{where}: {a!r} vs pinned {e!r}")

    walk(expected, actual, path)
    return violations, worst


def assert_matches_capture(
    expected: Any, actual: Any, *, rtol: float, atol: float, label: str
) -> None:
    """Exact on the capture platform; discrete-exact and float-bounded elsewhere."""
    if on_capture_platform():
        assert actual == expected, label
        return
    violations, worst = capture_deviations(expected, actual, rtol=rtol, atol=atol, path=label)
    assert not violations, (
        f"{len(violations)} deviation(s) from the {CAPTURE_SYSTEM} {CAPTURE_MACHINE} capture "
        f"(off-platform bound rtol={rtol:g}, atol={atol:g}; worst float at "
        f"{worst:.2f} of the bound): " + "; ".join(violations[:10])
    )
