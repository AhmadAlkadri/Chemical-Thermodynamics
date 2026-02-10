"""Shared helpers for VLE/flash equilibrium calculations."""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from ..core import Mixture
from ..exceptions import CompositionError, ModelError


def as_float_array(values: Sequence[float]) -> np.ndarray:
    """Convert values to a 1D float array."""
    return np.array(list(values), dtype=float)


def wilson_k(mixture: Mixture, temperature_K: float, pressure_Pa: float) -> np.ndarray:
    """Wilson K-value estimate for each component (dimensionless)."""
    values = []
    for component in mixture.components:
        ln_k = math.log(component.pc_pa / pressure_Pa) + 5.373 * (1.0 + component.omega) * (
            1.0 - component.tc_k / temperature_K
        )
        values.append(math.exp(ln_k))

    K = np.array(values, dtype=float)
    if np.any(~np.isfinite(K)) or np.any(K <= 0.0):
        raise ModelError("Wilson K-value initialization failed.")
    return K


def normalize_composition(
    values: np.ndarray,
    *,
    label: str,
    error_cls: type[Exception] = CompositionError,
) -> np.ndarray:
    """Normalize a non-negative composition vector to sum to one."""
    if np.any(~np.isfinite(values)):
        raise error_cls(f"Non-finite {label} composition encountered.")
    if np.any(values < 0.0):
        raise error_cls(f"Negative {label} composition encountered.")

    total = float(np.sum(values))
    if total <= 0.0:
        raise error_cls(f"{label.capitalize()} composition has non-positive total.")

    normalized = values / total
    if np.any(normalized < 0.0):
        raise error_cls(f"Negative normalized {label} composition encountered.")

    return normalized
