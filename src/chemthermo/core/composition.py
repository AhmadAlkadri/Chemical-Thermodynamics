"""Composition representation with explicit normalization control."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

from ..exceptions import CompositionError
from ..validation import COMPOSITION_SUM_TOL, validate_fractions

CompositionBasis = Literal["mole", "mass"]


@dataclass(frozen=True)
class Composition:
    """Component fractions on a mole or mass basis."""

    fractions: Sequence[float]
    basis: CompositionBasis = "mole"
    normalize: bool = False
    tol: float = COMPOSITION_SUM_TOL

    def __post_init__(self) -> None:
        if self.basis not in ("mole", "mass"):
            raise CompositionError(f"Unsupported composition basis: {self.basis!r}.")

        values = validate_fractions(self.fractions, normalize=self.normalize, tol=self.tol)
        object.__setattr__(self, "fractions", tuple(values))

    def __len__(self) -> int:
        return len(self.fractions)

    def to_list(self) -> list[float]:
        return list(self.fractions)

    def sum(self) -> float:
        return sum(self.fractions)
