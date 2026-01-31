"""Public interfaces for VLLE engines."""

from __future__ import annotations

from typing import Mapping, Protocol, Sequence

from .types import Scalar, VLLEResult


class VLLEEngine(Protocol):
    """Protocol for VLLE solver plugins."""

    def solve(
        self,
        *,
        temperature_K: float,
        pressure_Pa: float,
        z: Sequence[float],
        components: Sequence[str] | None = None,
        model: str | None = None,
        options: Mapping[str, Scalar] | None = None,
    ) -> VLLEResult:
        """Solve for VLLE phase split at the given state."""
        ...
