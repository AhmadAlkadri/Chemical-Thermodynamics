"""Equations of state and helper utilities.

The flash package previously shipped only interface-level placeholders for
equations of state. This module adds a concrete Peng–Robinson implementation
that can compute compressibility factors, molar volumes, and densities for
pure components and mixtures using standard mixing rules.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence

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
            alpha = (1 + kappa * (1 - np.sqrt(tr))) ** 2

            a_i = 0.45724 * (R_BAR_L_PER_MOL_K**2) * (tc**2) / pc
            b_i = 0.07780 * R_BAR_L_PER_MOL_K * tc / pc

            params.append({"a": a_i * alpha, "b": b_i, "mw": mw})
        return params

    def _apply_mixing_rules(
        self, fractions: Iterable[float], params: Sequence[Dict[str, float]]
    ) -> Dict[str, float]:
        y = np.array(list(fractions), dtype=float)
        y /= y.sum()

        a_values = np.array([p["a"] for p in params])
        b_values = np.array([p["b"] for p in params])

        a_matrix = np.sqrt(np.outer(a_values, a_values)) * (1 - self.kij)
        a_mix = float(np.sum(y[:, None] * y[None, :] * a_matrix))
        b_mix = float(np.sum(y * b_values))
        mw_mix = float(np.sum(y * np.array([p["mw"] for p in params])))

        return {"a": a_mix, "b": b_mix, "mw": mw_mix}

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
            details={"A": mixed["a"], "B": mixed["b"], "MW_mix[g/mol]": mixed["mw"]},
        )

