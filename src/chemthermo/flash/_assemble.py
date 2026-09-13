"""Assemble a `FlashResult` with any number of phases from converged compositions."""

from __future__ import annotations

from typing import Sequence

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


def _multi_phase_result(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    *,
    compositions: Sequence[np.ndarray],
    fractions: Sequence[float] | np.ndarray,
    names: Sequence[str],
    diagnostics: dict[str, float | int | str | bool],
) -> FlashResult:
    """Assemble a `FlashResult` for any number of phases (ADR-0011).

    ``vapor_fraction`` is the fraction of the phase named ``"vapor"`` when the
    phase set contains one, and None otherwise: reporting a vapor fraction for
    a set with no vapor in it would be fiction, and that is the same rule the
    two-phase liquid-liquid paths already follow.

    Phase fractions are renormalized to sum to exactly one in floating point,
    because `FlashResult` enforces that invariant and a converged multiphase
    Rachford-Rice solution satisfies it only to round-off.
    """
    total = float(sum(float(value) for value in fractions))
    normalized = [float(value) / total for value in fractions]
    # Absorb the remaining ULPs into the largest phase, which is the least
    # sensitive to them, so the sum is exactly 1.0.
    largest = max(range(len(normalized)), key=lambda index: normalized[index])
    normalized[largest] += 1.0 - float(sum(normalized))

    phases = {
        name: PhaseResult(
            name=name,
            composition=Composition(
                fractions=tuple(np.asarray(composition, dtype=float).tolist()),
                basis=mixture.basis,
                normalize=True,
                tol=COMPOSITION_SUM_TOL,
            ),
        )
        for name, composition in zip(names, compositions)
    }
    phase_fractions = {name: value for name, value in zip(names, normalized)}
    vapor_fraction = phase_fractions.get("vapor")
    return FlashResult(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        phases=phases,
        vapor_fraction=None if vapor_fraction is None else float(vapor_fraction),
        phase_fractions=phase_fractions,
        diagnostics=diagnostics,
    )
