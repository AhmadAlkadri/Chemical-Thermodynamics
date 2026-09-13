"""Interfaces for thermodynamic models."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Sequence

from ..core import Mixture

#: Compressibility-based liquid/vapor threshold (ADR-0017). An implementation
#: of :meth:`EquationOfState.phase_identity` reports ``"liquid"`` when the
#: dimensionless isothermal-compressibility ratio
#: ``kappa = (1 / rho) * (d rho / dP)_T * P`` (equivalently
#: ``kappa = -P / (V * (dP/dV)_T)``) is below this value at the evaluated
#: root, and ``"vapor"`` otherwise. ``kappa`` is exactly 1 for an ideal gas
#: (``P = rho R T``) and well below 1 for a liquid; the value below was chosen
#: so that it separates every liquid root from every vapor root, with a wide
#: margin, over both the Peng-Robinson and the PC-SAFT validation grids - see
#: ``.agents/brain/adr/0017-phase-identity-by-compressibility.md`` for the
#: measured kappa tables.
KAPPA_LIQUID_THRESHOLD: float = 0.5


class EquationOfState(ABC):
    """Interface for equation-of-state models.

    Implementations should accept SI units (K, Pa) and return dimensionless
    fugacity coefficients for a specified phase.
    """

    name: str = "generic-eos"

    @abstractmethod
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
            mixture: Mixture providing component properties.
            temperature_K: Temperature in K.
            pressure_Pa: Pressure in Pa.
            composition: Mole fractions, sum to 1 within COMPOSITION_SUM_TOL.
            phase: Phase label (e.g., "vapor" or "liquid").

        Returns:
            Fugacity coefficients (dimensionless), one per component.
        """

    def phase_identity(
        self,
        *,
        mixture: Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: Sequence[float],
        phase: str,
    ) -> str | None:
        """Return ``"liquid"`` or ``"vapor"`` from a compressibility criterion (ADR-0017).

        This identifies the *root* ``phase`` selects (the same root
        :meth:`fugacity_coefficients` would use) from the model's own
        pressure-volume behavior, rather than from which of the two labels
        happens to have lower Gibbs energy - the two agree whenever the model
        has two distinct roots, and only the compressibility criterion is
        meaningful when it has one (see
        ``chemthermo.models.base.KAPPA_LIQUID_THRESHOLD``).

        Args:
            mixture: Mixture providing component properties.
            temperature_K: Temperature in K.
            pressure_Pa: Pressure in Pa.
            composition: Mole fractions, sum to 1 within COMPOSITION_SUM_TOL.
            phase: Phase label (``"vapor"`` or ``"liquid"``) selecting which
                root to identify, with the same semantics
                :meth:`fugacity_coefficients` uses.

        Returns:
            ``"liquid"`` or ``"vapor"``, or ``None`` when this implementation
            does not support compressibility-based identification. **The
            default implementation always returns ``None``**, which is the
            documented fallback: a caller that receives ``None`` must keep its
            pre-ADR-0017 behavior (the min-Gibbs tie-break for a single-phase
            result, the Wilson volatility ranking for a phi-phi split)
            unchanged.
        """
        return None


class ActivityModel(ABC):
    """Interface for activity coefficient models.

    Implementations should accept SI units (K) and return dimensionless
    activity coefficients for the liquid phase.
    """

    name: str = "generic-activity"

    @abstractmethod
    def activity_coefficients(
        self, *, mixture: Mixture, temperature_K: float, composition: Sequence[float]
    ) -> Sequence[float]:
        """Return activity coefficients for a phase composition.

        Args:
            mixture: Mixture providing component identities.
            temperature_K: Temperature in K.
            composition: Mole fractions, sum to 1 within COMPOSITION_SUM_TOL.

        Returns:
            Activity coefficients (dimensionless), one per component.
        """
