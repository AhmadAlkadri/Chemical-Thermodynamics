"""Peng-Robinson equation of state in SI units."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from ..core import Mixture
from ..data import normalize_name
from ..exceptions import CompositionError, ModelError
from ..validation import (
    COMPOSITION_SUM_TOL,
    validate_fractions,
    validate_pressure,
    validate_temperature,
)
from .base import EquationOfState

R_J_PER_MOL_K = 8.314462618

#: Canonical stored form of a per-pair kij matrix: a sorted tuple of
#: ``((name_a, name_b), value)`` entries with ``name_a < name_b`` (both
#: normalized via ``chemthermo.data.normalize_name``). Kept immutable so the
#: frozen dataclass stays hashable/comparable and its repr is deterministic.
KijPairs = tuple[tuple[tuple[str, str], float], ...]
KijInput = float | Mapping[tuple[str, str], float]


def _canonicalize_kij(kij: KijInput | KijPairs) -> float | KijPairs:
    """Normalize constructor input for ``PengRobinsonEOS.kij``.

    A scalar is returned unchanged (as a ``float``); it is applied to every
    off-diagonal pair and never to the diagonal. A mapping from unordered
    component-name pairs to values is normalized (name case/whitespace via
    ``normalize_name``, pair order) into a sorted tuple. Both orders of a pair
    must agree if both are given; a pair naming the same component twice is
    rejected. Values are looked up per-mixture later, so pair names unknown to
    any particular mixture are simply never used.

    Idempotent: re-running this on an already-canonical tuple (as happens on
    ``dataclasses.replace``) returns it unchanged rather than re-validating,
    since it was already validated when first constructed.
    """
    if isinstance(kij, tuple):
        return kij
    if isinstance(kij, Mapping):
        canonical: dict[tuple[str, str], float] = {}
        for (name_a, name_b), value in kij.items():
            key_a = normalize_name(name_a)
            key_b = normalize_name(name_b)
            if not key_a or not key_b:
                raise ModelError("PengRobinsonEOS kij component names must be non-empty.")
            if key_a == key_b:
                raise ModelError(
                    "PengRobinsonEOS kij pair components must be distinct "
                    f"(got {name_a!r} paired with itself)."
                )
            pair_key = (key_a, key_b) if key_a < key_b else (key_b, key_a)
            value_f = float(value)
            if pair_key in canonical and canonical[pair_key] != value_f:
                raise ModelError(
                    f"Conflicting kij values given for pair {pair_key!r}: "
                    f"{canonical[pair_key]!r} vs {value_f!r}."
                )
            canonical[pair_key] = value_f
        return tuple(sorted(canonical.items()))
    return float(kij)


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
    mapping (see :func:`_canonicalize_kij`); the type annotation below
    includes that canonical tuple shape for that reason.
    """

    kij: KijInput | KijPairs = 0.0
    name: str = "Peng-Robinson"

    def __post_init__(self) -> None:
        object.__setattr__(self, "kij", _canonicalize_kij(self.kij))

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
        kij_matrix = self._kij_matrix(mixture)
        aij = np.sqrt(np.outer(a_i, a_i)) * (1.0 - kij_matrix)

        a_mix = float(np.sum(y[:, None] * y[None, :] * aij))
        b_mix = float(np.sum(y * b_i))
        return a_i, b_i, aij, a_mix, b_mix

    def _kij_matrix(self, mixture: Mixture) -> np.ndarray:
        """Build the dense n x n kij matrix for ``mixture``'s component order.

        The diagonal is always zero. A scalar ``kij`` fills every off-diagonal
        entry; a per-pair mapping fills only the pairs it names (by
        normalized component name) and defaults missing pairs to zero.
        """
        n = len(mixture.components)
        matrix = np.zeros((n, n), dtype=float)

        if isinstance(self.kij, (int, float)):
            if self.kij != 0.0:
                matrix[:, :] = self.kij
                np.fill_diagonal(matrix, 0.0)
            return matrix

        # self.kij is always the canonical KijPairs tuple here (never a raw
        # Mapping at runtime; __post_init__ already converted it). Both
        # branches call the same dict(...) constructor -- the isinstance
        # split exists only so pyright resolves a single dict() overload per
        # branch instead of over-widening a `Mapping | KijPairs` union.
        pairs: dict[tuple[str, str], float]
        if isinstance(self.kij, Mapping):
            pairs = dict(self.kij)
        else:
            pairs = dict(self.kij)
        names = [normalize_name(name) for name in mixture.component_names]
        for i in range(n):
            for j in range(i + 1, n):
                key = (names[i], names[j]) if names[i] < names[j] else (names[j], names[i])
                value = pairs.get(key, 0.0)
                matrix[i, j] = value
                matrix[j, i] = value
        return matrix

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
