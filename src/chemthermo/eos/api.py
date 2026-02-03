"""Protocols for residual-Helmholtz equation-of-state models."""

from __future__ import annotations

from typing import Protocol, Sequence


class EOSProtocol(Protocol):
    """Protocol for residual Helmholtz EOS implementations.

    Implementations should return the reduced residual Helmholtz energy
    (A^res / RT), given temperature, molar volume, and composition.
    """

    name: str

    def num_components(self) -> int:
        """Return the number of components the EOS instance was configured for."""
        ...

    def residual_helmholtz(
        self,
        *,
        temperature_K: float,
        volume_m3: float,
        composition: Sequence[float],
    ) -> float:
        """Return reduced residual Helmholtz energy for the given state."""
        ...
