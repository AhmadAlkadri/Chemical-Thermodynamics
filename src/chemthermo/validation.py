"""Input validation helpers for chemthermo."""

from __future__ import annotations

import math
from typing import Iterable, Sequence

from .exceptions import CompositionError, InputRangeError

COMPOSITION_SUM_TOL = 1e-8


def validate_temperature(temperature: float) -> float:
    """Validate temperature in Kelvin (must be finite and > 0)."""

    value = float(temperature)
    if not math.isfinite(value) or value <= 0.0:
        raise InputRangeError(f"Temperature must be finite and > 0 K; got {temperature!r}.")
    return value


def validate_pressure(pressure: float) -> float:
    """Validate pressure in Pascals (must be finite and > 0)."""

    value = float(pressure)
    if not math.isfinite(value) or value <= 0.0:
        raise InputRangeError(f"Pressure must be finite and > 0 Pa; got {pressure!r}.")
    return value


def validate_fractions(
    fractions: Sequence[float] | Iterable[float],
    *,
    normalize: bool = False,
    tol: float = COMPOSITION_SUM_TOL,
) -> list[float]:
    """Validate composition fractions.

    Parameters
    ----------
    fractions:
        Sequence of mole or mass fractions.
    normalize:
        When True, normalize the fractions to sum to 1. When False (default),
        the sum must already be 1 within tolerance.
    tol:
        Absolute tolerance for the sum-to-one check.
    """

    values = [float(value) for value in fractions]
    if not values:
        raise CompositionError("Composition must contain at least one fraction.")

    for value in values:
        if not math.isfinite(value):
            raise CompositionError("Composition contains a non-finite value.")
        if value < 0.0:
            raise CompositionError("Composition fractions must be non-negative.")

    total = sum(values)
    if total <= 0.0:
        raise CompositionError("Composition fractions must sum to a positive value.")

    if normalize:
        return [value / total for value in values]

    if abs(total - 1.0) > tol:
        raise CompositionError(
            f"Composition fractions must sum to 1 within {tol}; got {total:.8f}."
        )

    return values
