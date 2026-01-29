"""Component metadata and property access."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from ..data import get_component_record
from ..exceptions import PropertyNotFoundError


@dataclass(frozen=True)
class Component:
    """Chemical component metadata in SI units.

    Attributes:
        name: Component common name.
        formula: Chemical formula string.
        properties: Mapping of property name to value (SI units).
        antoine: Optional Antoine coefficients mapping.

    Notes:
        Use Component.from_database to load a record from the built-in
        databank. Missing properties raise PropertyNotFoundError via
        require_property().
    """

    name: str
    formula: str
    properties: Mapping[str, float]
    antoine: Mapping[str, float] | None = None

    @classmethod
    def from_database(cls, name: str) -> "Component":
        record = get_component_record(name)
        antoine = record.get("antoine")
        return cls(
            name=str(record["name"]),
            formula=str(record["formula"]),
            properties=dict(record["properties"]),
            antoine=dict(antoine) if isinstance(antoine, dict) else None,
        )

    def require_property(self, key: str) -> float:
        try:
            value = self.properties[key]
        except KeyError as exc:
            raise PropertyNotFoundError(
                f"Property '{key}' not available for component '{self.name}'."
            ) from exc
        return float(value)

    @property
    def mw_kg_per_mol(self) -> float:
        return self.require_property("MW_kg_per_mol")

    @property
    def tc_k(self) -> float:
        return self.require_property("Tc_K")

    @property
    def pc_pa(self) -> float:
        return self.require_property("Pc_Pa")

    @property
    def omega(self) -> float:
        return self.require_property("omega")
