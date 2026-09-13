"""Assemble a single- or two-phase `FlashResult` from converged compositions."""

from __future__ import annotations

import numpy as np

from ..core import Composition, Mixture
from ..validation import COMPOSITION_SUM_TOL
from .results import FlashResult, PhaseResult


def _single_phase_result(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    *,
    phase_name: str,
    vapor_fraction: float | None,
    diagnostics: dict[str, float | int | str | bool],
) -> FlashResult:
    composition = Composition(
        fractions=mixture.fractions, basis=mixture.basis, normalize=False, tol=COMPOSITION_SUM_TOL
    )
    phase = PhaseResult(name=phase_name, composition=composition)
    phase_fractions = {phase_name: 1.0}
    return FlashResult(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        phases={phase_name: phase},
        vapor_fraction=vapor_fraction,
        phase_fractions=phase_fractions,
        diagnostics=diagnostics,
    )


def _two_phase_result(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    x: np.ndarray,
    y: np.ndarray,
    beta: float,
    *,
    names: tuple[str, str],
    vapor_fraction: float | None,
    diagnostics: dict[str, float | int | str | bool],
) -> FlashResult:
    first_name, second_name = names
    first = PhaseResult(
        name=first_name,
        composition=Composition(
            fractions=tuple(x.tolist()),
            basis=mixture.basis,
            normalize=False,
            tol=COMPOSITION_SUM_TOL,
        ),
    )
    second = PhaseResult(
        name=second_name,
        composition=Composition(
            fractions=tuple(y.tolist()),
            basis=mixture.basis,
            normalize=False,
            tol=COMPOSITION_SUM_TOL,
        ),
    )
    phase_fractions = {first_name: 1.0 - float(beta), second_name: float(beta)}
    return FlashResult(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        phases={first_name: first, second_name: second},
        vapor_fraction=None if vapor_fraction is None else float(vapor_fraction),
        phase_fractions=phase_fractions,
        diagnostics=diagnostics,
    )
