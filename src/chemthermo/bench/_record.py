"""Record shapes, result hashing and environment capture for the benchmark harness.

Internal to :mod:`chemthermo.bench` (ADR-0001, ADR-0023): nothing here is
exported from ``chemthermo``.

The one idea worth stating: a benchmark record is only useful next to the
*answer* the run produced. A faster run that moved a mole fraction is not an
optimization, so every case carries a :func:`result_hash` over the accepted
thermodynamic result - phase names, compositions and phase fractions, each
formatted to twelve significant digits - and
:func:`chemthermo.bench.compare` refuses to report a speedup when two records
disagree on it.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

#: Digits kept when a float enters the result hash. Twelve is well inside the
#: solvers' own convergence tolerances (``1e-08`` on the flash, ``1e-10`` on the
#: stability test) and well outside the last-bit noise a pure code motion can
#: produce, so the hash changes when an *answer* changes and not when a sum is
#: reassociated. A bit-identity claim is proved by
#: ``tests/test_flash_refactor_bit_identity.py`` and by
#: ``tests/test_eos_branch_reuse.py``, not by this hash.
RESULT_DIGITS = 12

#: Schema tag written into every record file, so a future field change can be
#: detected rather than silently mis-read.
SCHEMA = "chemthermo-bench/1"


def _format(value: float) -> str:
    return f"{float(value):.{RESULT_DIGITS - 1}e}"


@dataclass(frozen=True)
class StateOutcome:
    """The accepted answer at one ``(T, P, z)`` state of a case.

    Attributes:
        temperature_K: State temperature in K.
        pressure_Pa: State pressure in Pa.
        composition: Overall (feed) mole fractions.
        phases: Phase names in the order the solver reported them.
        phase_compositions: Mole fractions of each phase, in ``phases`` order.
        phase_fractions: Molar phase fractions, in ``phases`` order.
        iterations: Iteration counts read out of the solver's own diagnostics.
        initialization: What seeded the solve (``k_seed``, ``phase_detection``).
        status: ``"ok"``, or the exception class name when the solver refused
            the state. A refusal is recorded, never swallowed: a "speedup" that
            comes from a state no longer converging has to be visible.
        error: The refusal message, or None.
    """

    temperature_K: float
    pressure_Pa: float
    composition: tuple[float, ...]
    phases: tuple[str, ...] = ()
    phase_compositions: tuple[tuple[float, ...], ...] = ()
    phase_fractions: tuple[float, ...] = ()
    iterations: Mapping[str, int] = field(default_factory=dict)
    initialization: Mapping[str, str] = field(default_factory=dict)
    status: str = "ok"
    error: str | None = None

    def as_json(self) -> dict[str, Any]:
        return {
            "temperature_K": self.temperature_K,
            "pressure_Pa": self.pressure_Pa,
            "composition": list(self.composition),
            "phases": list(self.phases),
            "phase_compositions": [list(values) for values in self.phase_compositions],
            "phase_fractions": list(self.phase_fractions),
            "iterations": dict(self.iterations),
            "initialization": dict(self.initialization),
            "status": self.status,
            "error": self.error,
        }

    def hash_payload(self) -> list[Any]:
        """The part of this outcome the result hash is taken over."""
        return [
            _format(self.temperature_K),
            _format(self.pressure_Pa),
            [_format(value) for value in self.composition],
            list(self.phases),
            [[_format(value) for value in row] for row in self.phase_compositions],
            [_format(value) for value in self.phase_fractions],
            self.status,
        ]


def result_hash(outcomes: Sequence[StateOutcome]) -> str:
    """Return ``sha256:<hex>`` over the accepted results of a whole case.

    Deterministic across machines to the extent the solvers are: the payload is
    decimal text at :data:`RESULT_DIGITS` significant digits, serialized with
    sorted keys and no whitespace, so two runs of the same code on the same
    inputs hash the same and a changed answer does not.
    """
    payload = json.dumps(
        [outcome.hash_payload() for outcome in outcomes],
        sort_keys=True,
        separators=(",", ":"),
    )
    return "sha256:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _git(*args: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *args], capture_output=True, text=True, timeout=10, check=False
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def git_state() -> dict[str, Any]:
    """Return the commit the measurement was taken at, and whether it was dirty.

    ``None`` for both when git is not available or the working directory is not
    a repository; a record produced outside the repository is still a valid
    record, it just cannot name a commit.
    """
    sha = _git("rev-parse", "HEAD")
    status = _git("status", "--porcelain")
    return {
        "sha": sha,
        "dirty": None if status is None else bool(status),
    }


def environment() -> dict[str, Any]:
    """Return the hardware and interpreter the numbers were measured on.

    Wall times only mean something next to this. ``cpu`` is the most specific
    processor string the platform offers (``machdep.cpu.brand_string`` on
    macOS, ``platform.processor()`` elsewhere), because ``platform.processor()``
    is often just the architecture on Darwin.
    """
    import numpy as np

    from .. import __version__

    cpu = platform.processor() or None
    if platform.system() == "Darwin":
        brand = _sysctl("machdep.cpu.brand_string")
        if brand:
            cpu = brand
    return {
        "platform": platform.platform(),
        "system": platform.system(),
        "machine": platform.machine(),
        "cpu": cpu,
        "cpu_count": _cpu_count(),
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "numpy": np.__version__,
        "chemthermo": __version__,
    }


def _sysctl(name: str) -> str | None:
    try:
        completed = subprocess.run(
            ["sysctl", "-n", name], capture_output=True, text=True, timeout=5, check=False
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout.strip() or None


def _cpu_count() -> int | None:
    import os

    return os.cpu_count()
