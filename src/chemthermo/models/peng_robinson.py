"""Peng-Robinson equation of state in SI units."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .._eos_memo import MISS, active_memo
from ..core import Mixture
from ..exceptions import CompositionError, ModelError
from ..validation import (
    COMPOSITION_SUM_TOL,
    validate_fractions,
    validate_pressure,
    validate_temperature,
)
from ._kij import KijInput, KijPairs, canonicalize_kij, kij_matrix
from .base import KAPPA_LIQUID_THRESHOLD, EquationOfState

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

        log_phi = self._branch_log_phi(Z, A=A, B=B, y=y, b_i=b_i, aij=aij, a_mix=a_mix, b_mix=b_mix)

        phi = np.exp(np.array(log_phi, dtype=float))
        return phi.tolist()

    def ln_fugacity_branches(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
    ) -> dict[str, list[float]]:
        """Return ``ln phi_i`` on both compressibility roots from **one** cubic solve.

        The ``EquationOfState`` capability of ADR-0023. Validation, the mixing
        rule and ``np.roots`` are each done once here where two
        :meth:`fugacity_coefficients` calls at the same ``(T, P, x)`` did them
        twice, and the two branches then differ only in which root of the
        already-solved cubic they read - ``max`` for the vapour, ``min`` for the
        liquid, exactly as that method selects them.

        The values are ``ln phi`` *before* the exponential
        :meth:`fugacity_coefficients` finishes with, so ``exp`` of what comes
        back here is the same double that method returns.

        A branch whose root is inadmissible (``Z <= B``) or whose logarithmic
        term does not exist is **omitted** rather than reported, and the caller
        falls back to the per-branch route for that composition, where the
        refusal carries its own message.

        Raises:
            CompositionError: If the composition length mismatches or is invalid.
            ModelError: If the mixture parameters are invalid or the cubic has
                no positive real root.
            InputRangeError: If the temperature or pressure is out of range.
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

        branches: dict[str, list[float]] = {}
        for phase, Z in (("vapor", max(roots)), ("liquid", min(roots))):
            if Z <= B:
                continue
            try:
                branches[phase] = self._branch_log_phi(
                    Z, A=A, B=B, y=y, b_i=b_i, aij=aij, a_mix=a_mix, b_mix=b_mix
                )
            except ModelError:
                continue
        return branches

    @staticmethod
    def _branch_log_phi(
        Z: float,
        *,
        A: float,
        B: float,
        y: np.ndarray,
        b_i: np.ndarray,
        aij: np.ndarray,
        a_mix: float,
        b_mix: float,
    ) -> list[float]:
        """``ln phi_i`` on one compressibility root.

        The single source of the log-fugacity expression, shared by
        :meth:`fugacity_coefficients` and :meth:`ln_fugacity_branches`. The
        per-component loop is kept rather than vectorized: at the two- and
        three-component sizes this library's reference paths use, three numpy
        operations on a length-3 array cost more than the loop they replace
        (measured at 3.0 us against 1.4 us for the ternary of ADR-0023's
        ``pr-flash-ternary`` case), and a vectorized form would also have to
        reproduce this expression's exact operation order to stay
        bit-identical.
        """
        log_term = PengRobinsonEOS._log_term(Z, B)
        sqrt2 = math.sqrt(2.0)

        sum_y_aij = np.dot(aij, y)
        log_phi = []
        for i in range(len(y)):
            term1 = b_i[i] / b_mix * (Z - 1.0) - math.log(Z - B)
            term2 = (
                (A / (2.0 * sqrt2 * B)) * (2.0 * sum_y_aij[i] / a_mix - b_i[i] / b_mix) * log_term
            )
            log_phi.append(term1 - term2)
        return log_phi

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

    def phase_identity(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> str:
        """Return "liquid" or "vapor" from the compressibility criterion (ADR-0017).

        ``kappa = -P / (V * dP/dV)`` is the dimensionless isothermal
        compressibility times pressure, evaluated at the root ``phase``
        selects. ``dP/dV`` is the exact analytic derivative of the
        Peng-Robinson pressure equation
        ``P = RT/(V - b) - a / (V^2 + 2 b V - b^2)``, whose root in ``Z`` is
        algebraically the same cubic :meth:`compressibility_factor` already
        solves - so this needs no new root-finding and no finite difference.
        See ``chemthermo.models.base.KAPPA_LIQUID_THRESHOLD`` for the
        threshold and ADR-0017 for the measured separation between liquid and
        vapor roots.

        When the cubic has a single real root at this composition, both
        ``phase="liquid"`` and ``phase="vapor"`` name that same root and this
        method returns the same identity either way - which is exactly the
        case the historical vapor-first Gibbs tie-break (ADR-0008 decision 3,
        superseded by ADR-0017) could not tell apart.
        """
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
            Z = max(roots)
        elif phase == "liquid":
            Z = min(roots)
        else:
            raise ValueError("phase must be 'vapor' or 'liquid'.")

        if Z <= B:
            raise ModelError("Invalid state: Z <= B for Peng-Robinson EOS.")

        molar_volume = Z * R_J_PER_MOL_K * temperature / pressure
        denominator = molar_volume**2 + 2.0 * b_mix * molar_volume - b_mix**2
        dP_dV = (
            -R_J_PER_MOL_K * temperature / (molar_volume - b_mix) ** 2
            + a_mix * (2.0 * molar_volume + 2.0 * b_mix) / denominator**2
        )
        if not math.isfinite(dP_dV) or dP_dV >= 0.0:
            raise ModelError(
                "Peng-Robinson root is mechanically unstable (dP/dV >= 0); phase identity "
                "is undefined there."
            )
        kappa = -pressure / (molar_volume * dP_dV)
        return "liquid" if kappa < KAPPA_LIQUID_THRESHOLD else "vapor"

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
        """The positive real roots of the Peng-Robinson cubic in ``Z``, ascending.

        Inside a ``flash_tp`` / ``stability_tp`` call this is served from the
        call-local memo of :mod:`chemthermo._eos_memo` (ADR-0030) when the same
        ``(A, B)`` has already been solved in that call. ``A`` and ``B`` are the
        *only* inputs - the cubic's coefficients are functions of them alone -
        so the key needs neither the mixture nor the model instance, and a hit
        returns the roots the first solve produced rather than a re-solve of
        them. ``fugacity_coefficients``, ``ln_fugacity_branches``,
        ``compressibility_factor`` and ``phase_identity`` all reach this with
        the same ``(A, B)`` at one state, which is where the repeats are.
        Outside such a call there is no memo and this runs exactly as before.

        ``numpy.roots`` is what this has always used and what it still means:
        for a monic polynomial with a non-zero constant term that function is
        exactly "build the companion matrix, take its eigenvalues", and the
        companion matrix is built here instead so the eigenvalue call is
        reached without ``numpy.roots``' own polynomial bookkeeping - the
        ``atleast_1d``, the non-zero scan, the trimming, the dtype check, the
        division by a leading coefficient that is exactly ``1.0``, and the
        ``hstack`` of the trailing zeros. Measured at 8.6 us against 14.7 us on
        this cubic, roughly a sixth of one ``fugacity_coefficients`` call, and
        the eigenvalues are the **same array**: nothing about the numerics
        changed, only the route to LAPACK.

        The one case ``numpy.roots`` treats differently is a constant term of
        exactly ``0.0``, where it deflates to a quadratic and appends a zero
        root - a different eigenproblem, and so possibly different last bits.
        That case is handed back to ``numpy.roots`` rather than reproduced, and
        it is not reachable in practice anyway: the constant term is
        ``-(A B - B^2 - B^3)`` and ``B > 0``.

        No closed-form (Cardano) route is used. It would be far faster still
        and it does **not** reproduce these doubles, and bit-identity is the
        gate this slice is held to; see ADR-0023.
        """
        memo = active_memo()
        key: tuple[object, ...] | None = None
        if memo is not None:
            key = ("pr-compressibility-roots", A, B)
            cached = memo.lookup(key)
            if cached is not MISS:
                # A fresh list of the stored doubles: the caller owns its list,
                # and the memo's copy cannot be mutated underneath a later hit.
                return list(cached)

        c1 = -(1.0 - B)
        c2 = A - 3.0 * B**2 - 2.0 * B
        c3 = -(A * B - B**2 - B**3)

        if c3 == 0.0:
            roots = np.roots([1.0, c1, c2, c3])
        else:
            companion = np.zeros((3, 3), dtype=float)
            companion[0, 0] = -c1
            companion[0, 1] = -c2
            companion[0, 2] = -c3
            companion[1, 0] = 1.0
            companion[2, 1] = 1.0
            roots = np.linalg.eigvals(companion)

        real_roots = [float(root.real) for root in roots if abs(root.imag) < 1e-8]
        positive = sorted(root for root in real_roots if root > 0.0)
        if memo is not None and key is not None:
            memo.store(key, tuple(positive))
        return positive

    @staticmethod
    def _log_term(Z: float, B: float) -> float:
        sqrt2 = math.sqrt(2.0)
        num = Z + (1.0 + sqrt2) * B
        den = Z + (1.0 - sqrt2) * B
        if num <= 0.0 or den <= 0.0:
            raise ModelError("Invalid log term for Peng-Robinson EOS.")
        return math.log(num / den)
