"""Shared helpers for VLE/flash equilibrium calculations."""

from __future__ import annotations

import math
from typing import NamedTuple, Sequence

import numpy as np

from ..core import Mixture
from ..exceptions import CompositionError, ModelError
from ..models import EquationOfState


def as_float_array(values: Sequence[float]) -> np.ndarray:
    """Convert values to a 1D float array."""
    return np.array(list(values), dtype=float)


#: Wilson K-value used for a component declared **non-volatile** (ADR-0022).
#:
#: The Wilson correlation is written in ``Tc``, ``Pc`` and ``omega``, and a
#: polymer has none of the three. ``1e-10`` says "essentially absent from the
#: vapour-like trial": ``z K`` puts the polymer at a mole fraction of order
#: ``1e-10`` in the ``wilson-vapor`` start and ``z / K`` makes it dominant in
#: the ``wilson-liquid`` start, which is exactly the pair of estimates a
#: polymer/solvent system needs. It is an initial **estimate** - the trial then
#: iterates on the model - so no verdict, composition or fugacity is a function
#: of this number; what it buys is that the deterministic trial set stays
#: deterministic and complete for a component with no critical constants.
NON_VOLATILE_WILSON_K = 1e-10


def wilson_k(mixture: Mixture, temperature_K: float, pressure_Pa: float) -> np.ndarray:
    """Wilson K-value estimate for each component (dimensionless).

    A component whose ``volatile`` flag is ``False`` (only reachable through
    :meth:`chemthermo.Component.custom`, ADR-0022) gets
    :data:`NON_VOLATILE_WILSON_K` instead: it has no critical constants to
    evaluate the correlation with. Every databank component is volatile, so
    this function is unchanged for every mixture built before ADR-0022.
    """
    values = []
    for component in mixture.components:
        if not component.volatile:
            values.append(NON_VOLATILE_WILSON_K)
            continue
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


class EosBranchTerms(NamedTuple):
    """One density/compressibility branch of an equation of state at one composition.

    Attributes:
        phi: The fugacity coefficients the model returned, or ``None`` when
            they are not representable as doubles (see
            :func:`eos_branch_terms`).
        ln_phi: ``ln phi_i``, always usable. It is ``np.log(phi)`` whenever
            ``phi`` is not ``None`` - the same double the pre-ADR-0022 code
            produced - and the model's own ``ln phi`` otherwise.
    """

    phi: np.ndarray | None
    ln_phi: np.ndarray


def eos_branch_terms(
    eos: EquationOfState,
    *,
    mixture: Mixture,
    temperature: float,
    pressure: float,
    composition: np.ndarray,
    phase: str,
) -> EosBranchTerms:
    """Fugacity coefficients and their logarithm on one branch (ADR-0022 guard).

    The ordinary path is unchanged and is taken first: ask the model for
    ``phi``, and when every entry is finite and positive return it together
    with ``np.log(phi)``. Those are the doubles every caller used before
    ADR-0022, so nothing that converged then moves now.

    The guard is for the case the polymer slice measured: a chain molecule can
    have ``|ln phi|`` far past 709, where ``exp`` is ``0.0`` or ``inf`` and the
    model's perfectly finite answer is destroyed by the last step of
    ``fugacity_coefficients``. A 53 000 g/mol polyethylene (``m = 1393.9``) in
    n-pentane at 453 K and 10 MPa has ``ln phi_PE = -1690.6``; before this
    guard the stability test refused the state with "non-finite or non-positive
    fugacity coefficients" and the flash could not run at all. Only then is
    :meth:`~chemthermo.models.EquationOfState.log_fugacity_coefficients` asked,
    and a model that does not implement it (the default ``None``) produces the
    same refusal as before.

    Raises:
        ModelError: If the branch is unusable - the model raised, returned the
            wrong shape, or returned a ``phi`` that is not representable while
            offering no logarithmic route.
    """
    try:
        values = as_float_array(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature,
                pressure_Pa=pressure,
                composition=composition.tolist(),
                phase=phase,
            )
        )
    except Exception as exc:  # noqa: BLE001 - a branch may be absent here
        raise ModelError(str(exc)) from exc

    if values.shape != composition.shape:
        raise ModelError("inconsistent fugacity coefficient shape")
    if np.all(np.isfinite(values)) and np.all(values > 0.0):
        return EosBranchTerms(phi=values, ln_phi=np.log(values))

    try:
        logarithms = eos.log_fugacity_coefficients(
            mixture=mixture,
            temperature_K=temperature,
            pressure_Pa=pressure,
            composition=composition.tolist(),
            phase=phase,
        )
    except Exception as exc:  # noqa: BLE001 - a branch may be absent here
        raise ModelError(str(exc)) from exc
    if logarithms is None:
        raise ModelError("non-finite or non-positive fugacity coefficients")

    ln_phi = as_float_array(logarithms)
    if ln_phi.shape != composition.shape:
        raise ModelError("inconsistent fugacity coefficient shape")
    if not np.all(np.isfinite(ln_phi)):
        raise ModelError("non-finite log fugacity coefficients")
    return EosBranchTerms(phi=None, ln_phi=ln_phi)
