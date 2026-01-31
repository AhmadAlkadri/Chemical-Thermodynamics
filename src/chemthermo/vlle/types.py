"""Public data types for the VLLE API boundary."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping, Sequence

Scalar = float | int | str | bool
Options = Mapping[str, Scalar]


@dataclass(frozen=True)
class VLLEInputs:
    """Inputs for a VLLE solve request."""

    temperature_K: float
    pressure_Pa: float
    z: Sequence[float]
    components: Sequence[str] | None = None
    model: str | None = None
    options: Options = field(default_factory=dict)


@dataclass(frozen=True)
class VLLEPhase:
    """Phase information returned by a VLLE solver."""

    name: str
    fraction: float
    composition: Sequence[float]


@dataclass(frozen=True)
class VLLEResult:
    """VLLE solve results."""

    phases: Sequence[VLLEPhase]
    diagnostics: Mapping[str, Scalar] = field(default_factory=dict)
