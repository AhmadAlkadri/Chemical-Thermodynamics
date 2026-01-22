"""Result containers for flash calculations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

from ..core import Composition
from ..exceptions import InputRangeError
from ..validation import validate_pressure, validate_temperature


@dataclass(frozen=True)
class PhaseResult:
    """Per-phase composition and properties."""

    name: str
    composition: Composition
    properties: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("Phase name must be non-empty.")


@dataclass(frozen=True)
class FlashResult:
    """Outcome of a flash calculation."""

    temperature_K: float
    pressure_Pa: float
    phases: Mapping[str, PhaseResult]
    vapor_fraction: float | None = None
    diagnostics: Mapping[str, float | str | int | bool] = field(default_factory=dict)

    def __post_init__(self) -> None:
        validate_temperature(self.temperature_K)
        validate_pressure(self.pressure_Pa)

        if not self.phases:
            raise ValueError("FlashResult must include at least one phase.")

        for key, phase in self.phases.items():
            if key != phase.name:
                raise ValueError(
                    f"Phase mapping key '{key}' does not match phase name '{phase.name}'."
                )

        if self.vapor_fraction is not None:
            if not (0.0 <= self.vapor_fraction <= 1.0):
                raise InputRangeError("vapor_fraction must be between 0 and 1.")

    def phase_names(self) -> list[str]:
        return list(self.phases.keys())
