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
from ._kij import KijInput, KijPairs, canonicalize_kij, kij_matrix
from .base import EquationOfState

R_J_PER_MOL_K = 8.314462618

# ``KijInput`` / ``KijPairs`` / ``canonicalize_kij`` / ``kij_matrix`` moved to
# the internal ``chemthermo.models._kij`` in the ``pcsaft-residual-helmholtz``
# slice so PC-SAFT can reuse the same per-pair contract (ADR-0006, ADR-0014).
# Nothing about their behavior changed; Peng-Robinson stays bit-identical.


@dataclass(frozen=True)
class PengRobinsonEOS(EquationOfState):
    """Peng-Robinson EOS in SI units with simple quadratic mixing.

    Supports "vapor" and "liquid" phase labels and returns fugacity
    coefficients (dimensionless).

    ``kij`` is either a scalar applied to every off-diagonal pair (never to
    the diagonal, so pure-component behavior never changes with ``kij``), or a
    mapping from an unordered pair of component names to a per-pair value,
    e.g. ``{("Methane", "n-Decane"): 0.0411}``. Names are matched via
    ``chemthermo.data.normalize_name`` and pairs missing from the mapping
    default to ``0.0``. See ADR-0006 for the rationale.

    After construction ``kij`` holds the *canonical* form: unchanged if given
    as a scalar, or normalized to a sorted ``KijPairs`` tuple if given as a
    mapping (see :func:`chemthermo.models._kij.canonicalize_kij`); the type
    annotation below includes that canonical tuple shape for that reason.
    """

    kij: KijInput | KijPairs = 0.0
    name: str = "Peng-Robinson"

    def __post_init__(self) -> None:
        object.__setattr__(self, "kij", canonicalize_kij(self.kij, model="PengRobinsonEOS"))

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

        a_i, b_i, aij, a_mix, b_mix = self._mixture_parameters(mixture, temperature, y)

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

        _a_i, _b_i, _aij, a_mix, b_mix = self._mixture_parameters(mixture, temperature, y)

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

    def _mixture_parameters(
        self, mixture: Mixture, temperature: float, y: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
        """Return ``(a_i, b_i, aij, a_mix, b_mix)`` for a phase composition ``y``.

        Single source of truth for the quadratic mixing rule shared by
        :meth:`fugacity_coefficients` and :meth:`compressibility_factor`.
        ``aij``'s diagonal always uses a ``(1 - kij) = 1`` factor (pure a_ii is
        never corrupted by a nonzero kij); off-diagonal entries use the
        resolved per-pair (or scalar) kij for this mixture's component order.
        """
        a_i, b_i = self._component_parameters(mixture, temperature)
        kij_mat = self._kij_matrix(mixture)
        aij = np.sqrt(np.outer(a_i, a_i)) * (1.0 - kij_mat)

        a_mix = float(np.sum(y[:, None] * y[None, :] * aij))
        b_mix = float(np.sum(y * b_i))
        return a_i, b_i, aij, a_mix, b_mix

    def _kij_matrix(self, mixture: Mixture) -> np.ndarray:
        """Build the dense n x n kij matrix for ``mixture``'s component order.

        Thin wrapper over the shared helper in ``chemthermo.models._kij``; see
        that module for the (unchanged) rules.
        """
        return kij_matrix(self.kij, mixture.component_names)

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
