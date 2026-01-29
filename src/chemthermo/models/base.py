"""Interfaces for thermodynamic models."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Sequence

from ..core import Mixture


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
