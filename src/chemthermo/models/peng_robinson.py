"""Peng-Robinson equation of state in SI units."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from ..core import Mixture
from ..exceptions import CompositionError, ModelError
from ..validation import (
    COMPOSITION_SUM_TOL,
    validate_fractions,
    validate_pressure,
    validate_temperature,
)
from .base import EquationOfState

R_J_PER_MOL_K = 8.314462618


@dataclass(frozen=True)
class PengRobinsonEOS(EquationOfState):
    """Peng-Robinson EOS in SI units with simple quadratic mixing.

    Supports "vapor" and "liquid" phase labels and returns fugacity
    coefficients (dimensionless).
    """

    kij: float = 0.0
    name: str = "Peng-Robinson"

    def fugacity_coefficients(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> Sequence[float]:
        """Return fugacity coefficients for a phase composition.

        Args:
            mixture: Mixture providing component properties in SI units.
            temperature_K: Temperature in K.
            pressure_Pa: Pressure in Pa.
            composition: Mole fractions, sum to 1 within COMPOSITION_SUM_TOL.
            phase: "vapor" selects the largest real root; "liquid" the smallest.

        Returns:
            Fugacity coefficients (dimensionless), one per component.

        Raises:
            CompositionError: If composition length mismatches or is invalid.
            ModelError: If EOS parameters or roots are invalid.
            ValueError: If phase is not "vapor" or "liquid".
        """

        temperature = validate_temperature(temperature_K)
        pressure = validate_pressure(pressure_Pa)

        if len(composition) != len(mixture.components):
            raise CompositionError("Composition length must match number of mixture components.")

        fractions = validate_fractions(composition, normalize=False, tol=COMPOSITION_SUM_TOL)
        y = np.array(fractions, dtype=float)

        a_i, b_i = self._component_parameters(mixture, temperature)
        aij = np.sqrt(np.outer(a_i, a_i)) * (1.0 - self.kij)

        a_mix = float(np.sum(y[:, None] * y[None, :] * aij))
        b_mix = float(np.sum(y * b_i))

        if a_mix <= 0.0 or b_mix <= 0.0:
            raise ModelError("Invalid mixture parameters for Peng-Robinson EOS.")

        A = a_mix * pressure / (R_J_PER_MOL_K**2 * temperature**2)
        B = b_mix * pressure / (R_J_PER_MOL_K * temperature)

        roots = self._compressibility_roots(A, B)
        if not roots:
            raise ModelError("No real compressibility roots found for Peng-Robinson EOS.")

        if phase == "vapor":
            Z = max(roots)
        elif phase == "liquid":
            Z = min(roots)
        else:
            raise ValueError("phase must be 'vapor' or 'liquid'.")

        if Z <= B:
            raise ModelError("Invalid state: Z <= B for Peng-Robinson EOS.")

        log_term = self._log_term(Z, B)
        sqrt2 = math.sqrt(2.0)

        sum_y_aij = np.dot(aij, y)
        log_phi = []
        for i in range(len(y)):
            term1 = b_i[i] / b_mix * (Z - 1.0) - math.log(Z - B)
            term2 = (
                (A / (2.0 * sqrt2 * B)) * (2.0 * sum_y_aij[i] / a_mix - b_i[i] / b_mix) * log_term
            )
            log_phi.append(term1 - term2)

        phi = np.exp(np.array(log_phi, dtype=float))
        return phi.tolist()

    def compressibility_factor(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> float:
        """Return the compressibility factor Z for a phase composition."""
        temperature = validate_temperature(temperature_K)
        pressure = validate_pressure(pressure_Pa)

        if len(composition) != len(mixture.components):
            raise CompositionError("Composition length must match number of mixture components.")

        fractions = validate_fractions(composition, normalize=False, tol=COMPOSITION_SUM_TOL)
        y = np.array(fractions, dtype=float)

        a_i, b_i = self._component_parameters(mixture, temperature)
        aij = np.sqrt(np.outer(a_i, a_i)) * (1.0 - self.kij)

        a_mix = float(np.sum(y[:, None] * y[None, :] * aij))
        b_mix = float(np.sum(y * b_i))

        if a_mix <= 0.0 or b_mix <= 0.0:
            raise ModelError("Invalid mixture parameters for Peng-Robinson EOS.")

        A = a_mix * pressure / (R_J_PER_MOL_K**2 * temperature**2)
        B = b_mix * pressure / (R_J_PER_MOL_K * temperature)

        roots = self._compressibility_roots(A, B)
        if not roots:
            raise ModelError("No real compressibility roots found for Peng-Robinson EOS.")

        if phase == "vapor":
            return max(roots)
        elif phase == "liquid":
            return min(roots)
        else:
            raise ValueError("phase must be 'vapor' or 'liquid'.")

    @staticmethod
    def _component_parameters(
        mixture: Mixture, temperature: float
    ) -> tuple[np.ndarray, np.ndarray]:
        a_values = []
        b_values = []
        for component in mixture.components:
            tc = component.tc_k
            pc = component.pc_pa
            omega = component.omega

            tr = temperature / tc
            sqrt_tr = math.sqrt(tr)
            kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
            alpha = (1.0 + kappa * (1.0 - sqrt_tr)) ** 2

            a_i = 0.45724 * (R_J_PER_MOL_K**2) * (tc**2) / pc * alpha
            b_i = 0.07780 * R_J_PER_MOL_K * tc / pc

            a_values.append(a_i)
            b_values.append(b_i)

        return np.array(a_values, dtype=float), np.array(b_values, dtype=float)

    @staticmethod
    def _compressibility_roots(A: float, B: float) -> list[float]:
        coeffs = [
            1.0,
            -(1.0 - B),
            A - 3.0 * B**2 - 2.0 * B,
            -(A * B - B**2 - B**3),
        ]
        roots = np.roots(coeffs)
        real_roots = [float(root.real) for root in roots if abs(root.imag) < 1e-8]
        return sorted(root for root in real_roots if root > 0.0)

    @staticmethod
    def _log_term(Z: float, B: float) -> float:
        sqrt2 = math.sqrt(2.0)
        num = Z + (1.0 + sqrt2) * B
        den = Z + (1.0 - sqrt2) * B
        if num <= 0.0 or den <= 0.0:
            raise ModelError("Invalid log term for Peng-Robinson EOS.")
        return math.log(num / den)
