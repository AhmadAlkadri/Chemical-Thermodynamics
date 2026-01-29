"""Mixtures of components with validated compositions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from ..exceptions import CompositionError
from .component import Component
from .composition import Composition, CompositionBasis


@dataclass(frozen=True)
class Mixture:
    """A mixture of components with an associated composition.

    The number of components must match the composition length, and the
    composition must already satisfy its own validation contract.
    """

    components: tuple[Component, ...]
    composition: Composition

    def __post_init__(self) -> None:
        if not self.components:
            raise CompositionError("Mixture must contain at least one component.")
        if len(self.components) != len(self.composition):
            raise CompositionError(
                "Mixture components and composition fractions must have the same length."
            )

    @classmethod
    def from_components(
        cls,
        components: Sequence[Component] | Iterable[Component],
        fractions: Sequence[float] | Iterable[float],
        *,
        basis: CompositionBasis = "mole",
        normalize: bool = False,
    ) -> "Mixture":
        composition = Composition(fractions=tuple(fractions), basis=basis, normalize=normalize)
        return cls(components=tuple(components), composition=composition)

    @classmethod
    def from_database(
        cls,
        names: Sequence[str] | Iterable[str],
        fractions: Sequence[float] | Iterable[float],
        *,
        basis: CompositionBasis = "mole",
        normalize: bool = False,
    ) -> "Mixture":
        components = tuple(Component.from_database(name) for name in names)
        return cls.from_components(components, fractions, basis=basis, normalize=normalize)

    @property
    def fractions(self) -> tuple[float, ...]:
        return tuple(self.composition.fractions)

    @property
    def basis(self) -> CompositionBasis:
        return self.composition.basis

    @property
    def component_names(self) -> list[str]:
        return [component.name for component in self.components]
