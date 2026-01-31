"""Result containers for flash calculations."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

from ..core import Composition
from ..exceptions import InputRangeError
from ..validation import COMPOSITION_SUM_TOL, validate_pressure, validate_temperature


@dataclass(frozen=True)
class PhaseResult:
    """Per-phase composition and properties.

    Attributes:
        name: Phase name (non-empty).
        composition: Phase composition as a Composition object (fractions sum
            to 1 within COMPOSITION_SUM_TOL).
        properties: Optional scalar properties for the phase (user-defined).
    """

    name: str
    composition: Composition
    properties: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("Phase name must be non-empty.")


@dataclass(frozen=True)
class FlashResult:
    """Outcome of a flash calculation.

    Attributes:
        temperature_K: Temperature in K.
        pressure_Pa: Pressure in Pa.
        phases: Mapping from phase name to PhaseResult.
        vapor_fraction: Vapor fraction (0-1), if applicable.
        phase_fractions: Mapping from phase name to phase fraction.
        diagnostics: Diagnostic metadata from the solver.

    Invariants:
        - temperature_K > 0 and pressure_Pa > 0.
        - phases is non-empty, and mapping keys match phase.name.
        - phase_fractions values are in [0, 1] and sum to 1 within
          COMPOSITION_SUM_TOL.
        - If vapor_fraction is provided and phase_fractions includes "vapor",
          the two agree within COMPOSITION_SUM_TOL.

    Phase naming conventions:
        - VLE: "liquid", "vapor"
        - Single-phase: "liquid" or "vapor"
    """

    temperature_K: float
    pressure_Pa: float
    phases: Mapping[str, PhaseResult]
    vapor_fraction: float | None = None
    phase_fractions: Mapping[str, float] = field(default_factory=dict)
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
            if "vapor" in self.phase_fractions:
                if abs(self.phase_fractions["vapor"] - self.vapor_fraction) > COMPOSITION_SUM_TOL:
                    raise InputRangeError(
                        "vapor_fraction must match phase_fractions['vapor'] when provided."
                    )

        if self.phase_fractions:
            total = 0.0
            for name, fraction in self.phase_fractions.items():
                if name not in self.phases:
                    raise ValueError(f"Phase fraction provided for unknown phase '{name}'.")
                if not (0.0 <= fraction <= 1.0):
                    raise InputRangeError("phase_fractions must be between 0 and 1.")
                total += float(fraction)
            if abs(total - 1.0) > COMPOSITION_SUM_TOL:
                raise InputRangeError("phase_fractions must sum to 1 within tolerance.")

    def phase_names(self) -> list[str]:
        return list(self.phases.keys())
