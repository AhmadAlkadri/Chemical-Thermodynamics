"""Interfaces for thermodynamic models."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Sequence

from ..core import Mixture


class EquationOfState(ABC):
    """Base class for equation-of-state models."""

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
        """Return fugacity coefficients for a phase composition."""


class ActivityModel(ABC):
    """Base class for activity coefficient models."""

    name: str = "generic-activity"

    @abstractmethod
    def activity_coefficients(
        self, *, mixture: Mixture, temperature_K: float, composition: Sequence[float]
    ) -> Sequence[float]:
        """Return activity coefficients for a phase composition."""
