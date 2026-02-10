"""Component metadata and property access."""

from __future__ import annotations

from ..citations import BibliographicDatabase
from ..data import get_component_record
from ..schemas import AntoineCoefficients, ComponentData, Parameter


class Component:
    """Chemical component metadata in SI units."""

    def __init__(self, data: ComponentData):
        self._data = data
        self._bib_db: BibliographicDatabase | None = None

    @property
    def name(self) -> str:
        return self._data.name

    @property
    def formula(self) -> str:
        return self._data.formula

    @classmethod
    def from_database(cls, name: str) -> "Component":
        record = get_component_record(name)
        model = ComponentData(**record)
        return cls(model)

    def _bibliography(self) -> BibliographicDatabase:
        if self._bib_db is None:
            self._bib_db = BibliographicDatabase()
        return self._bib_db

    def get_citation(self, property_name: str) -> str:
        """Return the citation text for a given property."""
        prop = getattr(self._data, property_name, None)

        key: str | None = None
        if isinstance(prop, Parameter):
            key = prop.source_key
        elif isinstance(prop, AntoineCoefficients) and property_name == "antoine":
            key = prop.source_key

        if not key:
            return "No citation available."

        return self._bibliography().get_citation_text(key)

    @property
    def mw_kg_per_mol(self) -> float:
        return self._data.MW.value

    @property
    def tc_k(self) -> float:
        return self._data.Tc.value

    @property
    def pc_pa(self) -> float:
        return self._data.Pc.value

    @property
    def omega(self) -> float:
        return self._data.omega.value

    @property
    def antoine(self) -> AntoineCoefficients | None:
        return self._data.antoine
