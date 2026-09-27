"""Runtime-aware comparison against floats captured on one machine (ADR-0032).

Several guards pin values that were *captured* on the development machine and
stored as literals or JSON. Their job is to prove that a refactor or a
speed-up changed no arithmetic. On the capture runtime they still do exactly
that: every float is compared with ``==``.

The capture runtime is three things, not one (ADR-0032, amended 2026-09-26):

- **operating system and architecture**: macOS (``Darwin``) on ``arm64``;
- **interpreter**: CPython, **3.11** (major.minor; the patch level is not part
  of the gate). CPython 3.12 made the built-in ``sum()`` of floats
  compensated (Neumaier summation), so the same code on the same machine
  rounds differently from 3.12 on: under 3.12 and 3.13 the flash fixture
  moves in the last bits of four near-zero ``tpd`` diagnostics and is
  reproduced bit for bit again once ``sum`` is replaced by a plain left fold;
- **numpy**: *not* part of the gate. The capture was reproduced exactly with
  numpy 2.4.2 and 2.4.6 under CPython 3.11, and with numpy 2.5.3 under 3.12 /
  3.13 once ``sum`` is a left fold.

Off the capture runtime - another OS or CPU, another CPython minor version,
another interpreter - the same inputs can give floats that differ in the last
few bits: ``exp``/``log``/``pow`` and numpy's SIMD kernels are not
bit-reproducible across CPUs or C libraries, two Linux x86_64 hosts have been
measured to disagree with each other as well as with macOS (ledger Case P-11,
"cross-platform"), and the interpreter's own float arithmetic can change
between versions. There a float passes when

    |actual - expected| <= atol + rtol * max(|actual|, |expected|)

with the bound stated by the caller and justified in ADR-0032, and **every
discrete field is still exact**: dict keys, list lengths, strings, ints, bools,
``None``, and non-finite floats. A phase name, a verdict, a status, an
iteration count or a diagnostics key that moves fails on every runtime.

In-process A/B identity guards (the same state computed twice, two ways, in
one run) do not use this module: they are runtime-independent and stay exact
everywhere.
"""

from __future__ import annotations

import math
import numbers
import platform
import sys
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class Runtime:
    """What decides whether captured floats can be reproduced bit for bit."""

    system: str  #: ``platform.system()``, e.g. ``"Darwin"``, ``"Linux"``
    machine: str  #: ``platform.machine()``, e.g. ``"arm64"``, ``"x86_64"``
    implementation: str  #: ``platform.python_implementation()``, e.g. ``"CPython"``
    python: tuple[int, int]  #: interpreter (major, minor)

    def __str__(self) -> str:
        major, minor = self.python
        return f"{self.system} {self.machine} {self.implementation} {major}.{minor}"


#: Where the pinned fixtures were captured; exact comparison applies here only.
#: Reproduced exactly on 2026-09-25 with CPython 3.11.6 and numpy 2.4.6
#: (``.agents/handoffs/cloud-continuation.md`` section 7) and on 2026-09-26
#: with CPython 3.11.6 and numpy 2.4.2; the numpy version is informational.
CAPTURE_RUNTIME = Runtime(
    system="Darwin", machine="arm64", implementation="CPython", python=(3, 11)
)


def current_runtime() -> Runtime:
    """The runtime this process is on, in the terms of :data:`CAPTURE_RUNTIME`."""
    return Runtime(
        system=platform.system(),
        machine=platform.machine(),
        implementation=platform.python_implementation(),
        python=(sys.version_info.major, sys.version_info.minor),
    )


def is_capture_runtime(runtime: Runtime) -> bool:
    """True when ``runtime`` is the one the pinned floats were captured on."""
    return runtime == CAPTURE_RUNTIME


def on_capture_runtime() -> bool:
    """True where the pinned floats must be reproduced bit for bit."""
    return is_capture_runtime(current_runtime())


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
    """Exact on the capture runtime; discrete-exact and float-bounded elsewhere."""
    if on_capture_runtime():
        assert actual == expected, label
        return
    violations, worst = capture_deviations(expected, actual, rtol=rtol, atol=atol, path=label)
    assert not violations, (
        f"{len(violations)} deviation(s) from the {CAPTURE_RUNTIME} capture on "
        f"{current_runtime()} (bound rtol={rtol:g}, atol={atol:g}; worst float at "
        f"{worst:.2f} of the bound): " + "; ".join(violations[:10])
    )
