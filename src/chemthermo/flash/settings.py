"""Flash calculation settings."""

from __future__ import annotations

from dataclasses import dataclass

from ..exceptions import InputRangeError


@dataclass(frozen=True)
class FlashSettings:
    """Solver settings for flash calculations."""

    max_iter: int = 100
    tol: float = 1e-8
    damping: float | None = None

    def __post_init__(self) -> None:
        if self.max_iter <= 0:
            raise InputRangeError("max_iter must be positive.")
        if self.tol <= 0.0:
            raise InputRangeError("tol must be positive.")
        if self.damping is not None:
            if not (0.0 < self.damping <= 1.0):
                raise InputRangeError("damping must be in (0, 1] when specified.")
