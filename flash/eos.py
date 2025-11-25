"""Equations of state and helper utilities.

The flash package previously shipped only interface-level placeholders for
equations of state. This module adds a concrete Peng–Robinson implementation
that can compute compressibility factors, molar volumes, and densities for
pure components and mixtures using standard mixing rules.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from .system import Component, EquationOfState, Mixture


R_BAR_L_PER_MOL_K = 0.08314472  # L·bar·mol⁻¹·K⁻¹


@dataclass
class EosResult:
    """Collection of volumetric properties returned by an EOS calculation."""

    z_factors: List[float]
    molar_volumes_m3_per_mol: List[float]
    densities_kg_per_m3: List[float]
    details: Dict[str, float]

    def preferred_phase_volume(self, phase: str = "vapor") -> float:
        """Return a representative molar volume for the requested phase.

        Parameters
        ----------
        phase:
            Either ``"vapor"`` (largest real root) or ``"liquid"`` (smallest real
            root).
        """

        if not self.z_factors:
            raise ValueError("No compressibility factors available for selection.")

        if phase == "vapor":
            return self.molar_volumes_m3_per_mol[self.z_factors.index(max(self.z_factors))]
        if phase == "liquid":
            return self.molar_volumes_m3_per_mol[self.z_factors.index(min(self.z_factors))]
        raise ValueError(f"Unsupported phase label: {phase}")


class PengRobinsonEOS(EquationOfState):
    """Peng–Robinson cubic equation of state.

    The implementation follows the conventional form with the original Soave
    alpha correlation and simple quadratic mixing rules. Binary interaction
    parameters are assumed to be zero for now but can be incorporated later via
    the ``kij`` argument.
    """

    name = "Peng-Robinson"

    def __init__(self, kij: float | None = None) -> None:
        self.kij = kij or 0.0

    def _component_parameters(
        self, components: Sequence[Component], temperature: float
    ) -> List[Dict[str, float]]:
        params: List[Dict[str, float]] = []
        for comp in components:
            try:
                tc = float(comp.metadata["Tc[K]"])
                pc = float(comp.metadata["Pc[bar]"])
                omega = float(comp.metadata.get("omega", 0.0))
                mw = float(comp.metadata.get("MW[g/mol]", 0.0))
            except KeyError as exc:
                raise KeyError(
                    f"Missing critical property {exc!s} for component '{comp.name}'."
                ) from exc

            tr = temperature / tc
            kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
            sqrt_tr = np.sqrt(tr)
            alpha = (1 + kappa * (1 - sqrt_tr)) ** 2
            base = 1 + kappa * (1 - sqrt_tr)
            dalpha_dT = -kappa * base / (tc * sqrt_tr)

            a_i = 0.45724 * (R_BAR_L_PER_MOL_K**2) * (tc**2) / pc
            b_i = 0.07780 * R_BAR_L_PER_MOL_K * tc / pc
            da_i_dT = a_i * dalpha_dT

            params.append({"a": a_i * alpha, "b": b_i, "mw": mw, "da_dT": da_i_dT})
        return params

    def _apply_mixing_rules(
        self, fractions: Iterable[float], params: Sequence[Dict[str, float]]
    ) -> Dict[str, float]:
        y = np.array(list(fractions), dtype=float)
        y /= y.sum()

        a_values = np.array([p["a"] for p in params])
        b_values = np.array([p["b"] for p in params])
        da_values = np.array([p["da_dT"] for p in params])

        a_matrix = np.sqrt(np.outer(a_values, a_values)) * (1 - self.kij)
        a_mix = float(np.sum(y[:, None] * y[None, :] * a_matrix))
        b_mix = float(np.sum(y * b_values))
        mw_mix = float(np.sum(y * np.array([p["mw"] for p in params])))

        # Temperature derivative of the mixed a parameter
        da_mix = 0.0
        for i, y_i in enumerate(y):
            for j, y_j in enumerate(y):
                a_ij = (1 - self.kij) * math.sqrt(a_values[i] * a_values[j])
                term = 0.5 * a_ij * (
                    da_values[i] / a_values[i] + da_values[j] / a_values[j]
                )
                da_mix += y_i * y_j * term

        return {"a": a_mix, "b": b_mix, "mw": mw_mix, "da_dT": da_mix}

    def _solve_cubic(self, temperature: float, pressure: float, a: float, b: float) -> List[float]:
        A = a * pressure / (R_BAR_L_PER_MOL_K**2 * temperature**2)
        B = b * pressure / (R_BAR_L_PER_MOL_K * temperature)

        coeffs = [
            1,
            -(1 - B),
            A - 3 * B**2 - 2 * B,
            -(A * B - B**2 - B**3),
        ]

        roots = np.roots(coeffs)
        real_roots = sorted([float(r.real) for r in roots if abs(r.imag) < 1e-8])
        return real_roots

    def _calc_aij_matrix(self, params: Sequence[Dict[str, float]]) -> np.ndarray:
        a_values = np.array([p["a"] for p in params])
        return np.sqrt(np.outer(a_values, a_values)) * (1 - self.kij)

    def _a_mix(self, y: np.ndarray, aij: np.ndarray) -> float:
        return float(np.sum(y[:, None] * y[None, :] * aij))

    def _b_mix(self, y: np.ndarray, params: Sequence[Dict[str, float]]) -> float:
        return float(np.sum(y * np.array([p["b"] for p in params])))

    def _dA_over_dT(self, A: float, a_mix: float, da_dT: float, temperature: float) -> float:
        """Return dA/dT using A = aP/(R^2 T^2)."""

        if a_mix == 0:
            return 0.0
        return A * (da_dT / a_mix - 2.0 / temperature)

    def _t_dA_dT(self, A: float, a_mix: float, da_dT: float, temperature: float) -> float:
        """Return T * dA/dT for residual enthalpy expression."""

        if a_mix == 0:
            return 0.0
        return A * (temperature * da_dT / a_mix - 2.0)

    def _fugacity_log_coefficients(
        self,
        y: np.ndarray,
        params: Sequence[Dict[str, float]],
        Z: float,
        A: float,
        B: float,
        a_mix: float,
        b_mix: float,
        aij: np.ndarray,
    ) -> np.ndarray:
        bi = np.array([p["b"] for p in params])
        # For each component, sum_j y_j * a_ij
        sum_y_aij = np.dot(aij, y)

        log_phi = []
        log_term = math.log((Z + (1 + math.sqrt(2)) * B) / (Z + (1 - math.sqrt(2)) * B))
        for i in range(len(y)):
            first = bi[i] / b_mix * (Z - 1.0) - math.log(Z - B)
            second = (A / (2.0 * math.sqrt(2) * B)) * (
                2.0 * sum_y_aij[i] / a_mix - bi[i] / b_mix
            ) * log_term
            log_phi.append(first - second)
        return np.array(log_phi)

    def _residual_properties(
        self, Z: float, A: float, B: float, a_mix: float, da_dT: float, temperature: float
    ) -> Tuple[float, float]:
        """Return (HR_over_RT, SR_over_R)."""

        if B <= 0 or Z <= 0:
            return 0.0, 0.0
        log_term = math.log((Z + (1 + math.sqrt(2)) * B) / (Z + (1 - math.sqrt(2)) * B))
        t_dA_dT = self._t_dA_dT(A, a_mix, da_dT, temperature)
        dA_dT = self._dA_over_dT(A, a_mix, da_dT, temperature)

        HR_over_RT = Z - 1.0 + (t_dA_dT - A) / (2 * math.sqrt(2) * B) * log_term
        SR_over_R = math.log(Z - B) - dA_dT / (2 * math.sqrt(2) * B) * log_term
        return HR_over_RT, SR_over_R

    def compute_volumes(
        self, mixture: Mixture, temperature: float, pressure: float
    ) -> EosResult:
        """Compute compressibility factors and volumes for a mixture.

        Parameters
        ----------
        mixture:
            Mixture definition including component metadata.
        temperature:
            Temperature in Kelvin.
        pressure:
            Pressure in bar.
        """

        params = self._component_parameters(mixture.components, temperature)
        mixed = self._apply_mixing_rules(mixture.mole_fractions, params)
        z_factors = self._solve_cubic(temperature, pressure, mixed["a"], mixed["b"])

        molar_volumes = [
            z * R_BAR_L_PER_MOL_K * temperature / pressure / 1000 for z in z_factors
        ]  # convert L/mol to m³/mol
        densities = [
            (mixed["mw"] / 1000.0) / v if v != 0 else float("nan") for v in molar_volumes
        ]

        return EosResult(
            z_factors=z_factors,
            molar_volumes_m3_per_mol=molar_volumes,
            densities_kg_per_m3=densities,
            details={"A": mixed["a"], "B": mixed["b"], "MW_mix[g/mol]": mixed["mw"], "da_dT": mixed["da_dT"]},
        )

    def fugacity_coefficients(
        self, mixture: Mixture, composition: Sequence[float], temperature: float, pressure: float, phase_root: str = "vapor"
    ) -> Tuple[np.ndarray, float, float, float, float, float]:
        """Return fugacity coefficients and key state properties for a phase.

        Parameters
        ----------
        mixture:
            Mixture definition.
        composition:
            Phase mole fractions.
        temperature, pressure:
            State conditions (K, bar).
        phase_root:
            ``"vapor"`` selects the largest real root; ``"liquid"`` selects the smallest.

        Returns
        -------
        (phi, Z, Vm, rho, HR, SR)
            Fugacity coefficients, compressibility factor, molar volume (m3/mol),
            mass density (kg/m3), residual enthalpy (J/mol), and residual entropy (J/mol/K).
        """

        y = np.array(composition, dtype=float)
        if y.sum() <= 0:
            raise ValueError("Composition must contain at least one positive fraction.")
        y /= y.sum()

        params = self._component_parameters(mixture.components, temperature)
        aij = self._calc_aij_matrix(params)
        a_mix = self._a_mix(y, aij)
        b_mix = self._b_mix(y, params)
        da_mix = self._apply_mixing_rules(y, params)["da_dT"]  # reuse helper for derivative

        Z_roots = self._solve_cubic(temperature, pressure, a_mix, b_mix)
        if not Z_roots:
            raise ValueError("No real compressibility factors found.")
        Z = max(Z_roots) if phase_root == "vapor" else min(Z_roots)

        A = a_mix * pressure / (R_BAR_L_PER_MOL_K**2 * temperature**2)
        B = b_mix * pressure / (R_BAR_L_PER_MOL_K * temperature)
        log_phi = self._fugacity_log_coefficients(y, params, Z, A, B, a_mix, b_mix, aij)

        phi = np.exp(log_phi)
        Vm = Z * R_BAR_L_PER_MOL_K * temperature / pressure / 1000  # m3/mol
        MW_mix = float(np.sum(y * np.array([p["mw"] for p in params])))
        rho = (MW_mix / 1000.0) / Vm if Vm != 0 else float("nan")

        HR_over_RT, SR_over_R = self._residual_properties(Z, A, B, a_mix, da_mix, temperature)
        HR = HR_over_RT * 8.314462618 * temperature
        SR = SR_over_R * 8.314462618

        return phi, Z, Vm, rho, HR, SR
