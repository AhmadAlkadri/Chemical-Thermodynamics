"""NRTL activity coefficient model."""

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
        S = G @ x
        if np.any(S <= 0.0):
            raise ModelError("Invalid NRTL summations encountered.")

        A = np.zeros(count, dtype=float)
        for j in range(count):
            A[j] = float(np.sum(x * tau[:, j] * G[:, j]) / S[j])

        ln_gamma = np.zeros(count, dtype=float)
        for i in range(count):
            term1 = float(np.sum(x * tau[:, i] * G[:, i] / S))
            term2 = float(np.sum(x * G[i, :] / S[i] * (tau[i, :] - A)))
            ln_gamma[i] = term1 + term2

        gamma = np.exp(ln_gamma)
        if np.any(~np.isfinite(gamma)) or np.any(gamma <= 0.0):
            raise ModelError("Non-positive activity coefficients from NRTL.")

        return gamma.tolist()
