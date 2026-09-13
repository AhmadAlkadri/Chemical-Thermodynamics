"""NRTL activity coefficient model.

Index convention (fixed throughout this module and
`chemthermo.parameters.nrtl`)
------------------------------------------------------------------
``tau[i, j]`` is ``tau_ij`` (the ordered, generally asymmetric binary
interaction parameter of component ``j`` acting on component ``i``),
``alpha[i, j]`` is the non-randomness factor ``alpha_ij = alpha_ji``, and the
Boltzmann-type factor is

    G_ij = exp(-alpha_ij * tau_ij),      tau_ii = 0,  G_ii = 1.

A pair entry ``(A, B, tau_12, tau_21, alpha_12, alpha_21)`` therefore stores
``tau_AB`` in ``tau_12``; `NRTLParameters.for_components` places it at
``tau[index(A), index(B)]``.

Working equation
----------------
The Renon-Prausnitz local-composition model (Renon & Prausnitz, AIChE J. 14
(1968) 135) writes the reduced molar excess Gibbs energy as a sum over
*columns* of the parameter matrices,

    g^E(x) = G^E / (R T) = sum_j x_j * C_j(x) / S_j(x)                      (1)

with the two mole-fraction weighted averages

    S_j(x) = sum_k G_kj x_k                                                 (2)
    C_j(x) = sum_k tau_kj G_kj x_k                                          (3)

Both sums run over the *first* (row) index, i.e. down column ``j``.
Differentiating ``n g^E`` with respect to ``n_i`` at fixed T, P and the other
mole numbers gives the activity coefficient

    ln gamma_i = C_i / S_i
                 + sum_j [ x_j G_ij / S_j ] * ( tau_ij - C_j / S_j )        (4)

Two properties of (4) are easy to get wrong and are both covered by tests:

* the first term has a **single** denominator ``S_i``; it is not a sum of
  per-term quotients, and
* every ``S`` and ``C`` is a **column** sum (equations 2-3), not a row sum.
  Using row sums silently produces a function that is not the derivative of
  any ``g^E``, so it violates the Gibbs-Duhem relation
  ``sum_i x_i d ln gamma_i = 0`` at constant T, P (see
  `tests/test_activity_nrtl.py`). The two forms coincide only when ``G`` is
  symmetric, which is why a symmetric-tau binary test cannot detect the
  difference.

Equation (4) is the same expression used by Tessier, Brennecke & Stadtherr,
Chem. Eng. Sci. 55 (2000) 1785, section 2.2, whose published stationary points
are reproduced in `tests/validation/test_nrtl_tessier2000.py` and
`examples/validation/07_nrtl_tessier_stationary_points.py`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from ..core import Mixture
from ..exceptions import CompositionError, ModelError
from ..parameters import NRTLParameters
from ..validation import COMPOSITION_SUM_TOL, validate_fractions, validate_temperature
from .base import ActivityModel


def _default_parameters() -> NRTLParameters:
    return NRTLParameters.load()


@dataclass(frozen=True)
class NRTL(ActivityModel):
    """Non-random two-liquid (NRTL) activity coefficient model.

    Parameters are supplied by NRTLParameters and must cover all ordered
    component pairs in the mixture. Returns positive, finite activity
    coefficients for valid inputs.

    The implemented equation, its index convention and the Gibbs-Duhem
    property it satisfies are documented in the module docstring.
    """

    parameters: NRTLParameters = field(default_factory=_default_parameters)
    name: str = "NRTL"

    def activity_coefficients(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        composition: Sequence[float],
    ) -> Sequence[float]:
        """Return activity coefficients for a liquid composition.

        Evaluates equation (4) of the module docstring,

            ln gamma_i = C_i / S_i
                         + sum_j x_j G_ij / S_j * (tau_ij - C_j / S_j)

        with the column sums S_j = sum_k G_kj x_k and
        C_j = sum_k tau_kj G_kj x_k.

        Args:
            mixture: Mixture providing component identities.
            temperature_K: Temperature in K.
            composition: Mole fractions, sum to 1 within COMPOSITION_SUM_TOL.

        Returns:
            Activity coefficients (dimensionless), one per component.

        Raises:
            CompositionError: If composition length mismatches or is invalid.
            ModelError: If parameters are missing or results are non-physical.
        """
        validate_temperature(temperature_K)

        if len(composition) != len(mixture.components):
            raise CompositionError("Composition length must match number of mixture components.")

        x = np.array(
            validate_fractions(composition, normalize=False, tol=COMPOSITION_SUM_TOL), dtype=float
        )
        count = x.size
        if count == 1:
            return [1.0]

        tau, alpha = self.parameters.for_mixture(mixture)
        if tau.shape != (count, count) or alpha.shape != (count, count):
            raise ModelError("NRTL parameters returned inconsistent array shapes.")

        G = np.exp(-alpha * tau)

        # Column sums: S[j] = sum_k G[k, j] x[k], C[j] = sum_k tau[k, j] G[k, j] x[k].
        S = G.T @ x
        if np.any(S <= 0.0):
            raise ModelError("Invalid NRTL summations encountered.")
        C = (tau * G).T @ x
        A = C / S  # A[j] = C_j / S_j

        ln_gamma = A + np.sum((x * G * (tau - A)) / S, axis=1)

        gamma = np.exp(ln_gamma)
        if np.any(~np.isfinite(gamma)) or np.any(gamma <= 0.0):
            raise ModelError("Non-positive activity coefficients from NRTL.")

        return gamma.tolist()
