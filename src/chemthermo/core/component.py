"""Component metadata and property access."""

from __future__ import annotations

from ..citations import BibliographicDatabase
from ..data import get_component_record
from ..schemas import AntoineCoefficients, ComponentData, Parameter


class Component:
    """Chemical component metadata in SI units.

    Wraps the Pydantic ComponentData model to provide property access and citation queries.
    """

    def __init__(self, data: ComponentData):
        self._data = data
        self._bib_db = BibliographicDatabase()

    @property
    def name(self) -> str:
        return self._data.name

    @property
    def formula(self) -> str:
        return self._data.formula

    @classmethod
    def from_database(cls, name: str) -> "Component":
        # Temporary: Access the raw dict from legacy get_component_record and wrap it
        # TODO: Update get_component_record to return ComponentData directly in future slices
        record = get_component_record(name)

        # We need to adapt the legacy record layout to the new Schema if we are still loading old data
        # Check if the record matches the new schema (has MW as object, etc.)
        # If it's the old dict format, we might fail or need an adapter.
        # For now, let's assume we are transitioning and might need to handle both
        # or that we assume the DB has been migrated (Slice C).
        # Actually, since Slice C (Migration) hasn't run yet, get_component_record returns old dicts.
        # BUT, the goal of this slice is to "Update Component to use new Parameter".
        # So I should probably modify get_component_record OR make this adapter here.

        if "schema_version" in record:
            # It's already new format? No, individual records don't have schema_version, the DB does.
            pass

        # Since we haven't migrated data yet, we can't fully instantiate ComponentData from legacy dicts
        # without some conversion.
        # However, for the scope of Slice B (Provenance Plumbing), we might just want to support the *capability*.
        # Let's try to construct ComponentData from the record if possible, or fail if not.

        # Adaptation for legacy dict (if needed for testing before Slice C)
        # Old dict: {"name": "Methane", "properties": {"Tc_K": ...}, ...}
        # New model: name, MW: Parameter(...), Tc: Parameter(...)

        # For now, let's just assume the input 'record' IS compliant or we adapt it.
        # Since we are "Implementing Provenance Plumbing", let's build the plumbing to work with the NEW schema.
        try:
            model = ComponentData(**record)
            return cls(model)
        except Exception:
            # Fallback/Adapter logic would go here if we wanted hybrid support.
            # but to keep it simple, let's assume valid data is passed or allow failure until migration.
            # Actually, creating a dummy adapter to allow tests to run might be good.
            raise NotImplementedError(
                "Legacy data adapter not yet implemented. Run migration (Slice C)."
            )

    def get_citation(self, property_name: str) -> str:
        """Return the citation text for a given property."""
        prop = getattr(self._data, property_name, None)

        if isinstance(prop, Parameter):
            key = prop.source_key
            if key:
                return self._bib_db.get_citation_text(key)

        if isinstance(prop, AntoineCoefficients) and property_name == "antoine":
            key = prop.source_key
            if key:
                return self._bib_db.get_citation_text(key)

        return "No citation available."

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
