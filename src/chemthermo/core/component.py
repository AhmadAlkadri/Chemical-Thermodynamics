"""Component metadata and property access."""

from __future__ import annotations

import math

from ..citations import BibliographicDatabase
from ..data import get_component_record
from ..exceptions import InputRangeError, PropertyNotFoundError
from ..schemas import (
    AntoineCoefficients,
    ComponentCore,
    ComponentData,
    CustomComponentData,
    Parameter,
)


class Component:
    """Chemical component metadata in SI units."""

    def __init__(self, data: ComponentCore):
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

    @classmethod
    def custom(
        cls,
        name: str,
        *,
        mw_kg_per_mol: float,
        formula: str = "",
        tc_k: float | None = None,
        pc_pa: float | None = None,
        omega: float | None = None,
        volatile: bool = True,
        antoine: AntoineCoefficients | None = None,
        source: str = "",
    ) -> "Component":
        """Build a component from caller-supplied data, with no databank lookup.

        This is the entry point for a compound the packaged databank does not
        carry - a polymer, in the slice that added it (ADR-0022). The molar
        mass is required because it is what converts a mass fraction to a mole
        fraction and what a segments-per-mass PC-SAFT record is multiplied by;
        the **critical constants are optional**, because a polymer has none and
        inventing one to satisfy a schema would put a fabricated number where a
        model might read it. :attr:`tc_k`, :attr:`pc_pa` and :attr:`omega`
        raise :class:`~chemthermo.exceptions.PropertyNotFoundError` when they
        were not given.

        Args:
            name: Component name. It is what a PC-SAFT parameter record is
                matched against (after ``chemthermo.data.normalize_name``), so
                it must agree with the record's ``name``.
            mw_kg_per_mol: Molar mass in **kg/mol** (SI, as everywhere else in
                this package): 16.4 for a polyethylene of ``Mw = 16400 g/mol``.
            formula: Optional chemical formula, documentation only.
            tc_k: Critical temperature in K, or None.
            pc_pa: Critical pressure in Pa, or None.
            omega: Acentric factor, or None.
            volatile: Whether a vapour-like stability trial should place this
                component in the vapour at all. ``False`` - the polymer case -
                replaces the component's Wilson K-value estimate by the fixed
                ``chemthermo.flash._common.NON_VOLATILE_WILSON_K = 1e-10``,
                which is "essentially absent from the vapour-like trial". It is
                an initial **estimate** and nothing else reads it: no verdict,
                composition or fugacity depends on it.
            antoine: Optional Antoine coefficients, for the modified-Raoult
                path. A component without them is refused there, as before.
            source: Free-text provenance for the numbers above. Recorded on the
                molar-mass parameter so ``get_citation("MW")`` can find it.

        Returns:
            A :class:`Component` usable in :class:`~chemthermo.core.Mixture`.

        Raises:
            InputRangeError: If ``mw_kg_per_mol`` is not finite and positive,
                or if a critical constant that *was* given is not.
        """
        mw = _positive(mw_kg_per_mol, "mw_kg_per_mol", name)
        data = CustomComponentData(
            name=name,
            formula=formula,
            MW=Parameter(value=mw, units="kg/mol", source_key=source or None),
            Tc=(
                None if tc_k is None else Parameter(value=_positive(tc_k, "tc_k", name), units="K")
            ),
            Pc=(
                None
                if pc_pa is None
                else Parameter(value=_positive(pc_pa, "pc_pa", name), units="Pa")
            ),
            omega=(
                None if omega is None else Parameter(value=_finite(omega, "omega", name), units="-")
            ),
            volatile=bool(volatile),
            antoine=antoine,
        )
        return cls(data)

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
        """Critical temperature in K.

        Raises:
            PropertyNotFoundError: For a :meth:`custom` component built without
                one (a polymer has no critical point).
        """
        return self._required("Tc", "critical temperature").value

    @property
    def pc_pa(self) -> float:
        """Critical pressure in Pa; see :attr:`tc_k` for the missing-value rule."""
        return self._required("Pc", "critical pressure").value

    @property
    def omega(self) -> float:
        """Acentric factor; see :attr:`tc_k` for the missing-value rule."""
        return self._required("omega", "acentric factor").value

    @property
    def volatile(self) -> bool:
        """Whether a vapour-like trial estimate should place this component in the vapour.

        Always ``True`` for a databank component. ``False`` only for a
        :meth:`custom` component declared non-volatile (ADR-0022), and read by
        exactly one place: the Wilson K-value **estimate** in
        :func:`chemthermo.flash._common.wilson_k`.
        """
        return bool(getattr(self._data, "volatile", True))

    @property
    def antoine(self) -> AntoineCoefficients | None:
        return self._data.antoine

    def _required(self, field: str, label: str) -> Parameter:
        value = getattr(self._data, field, None)
        if value is None:
            raise PropertyNotFoundError(
                f"Component {self.name!r} has no {label}: it was built with "
                f"Component.custom(...) and no {field!r} was supplied. Supply one, or use a "
                "model that does not need it."
            )
        return value


def _finite(value: float, label: str, owner: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise InputRangeError(f"Component {owner!r}: {label} must be finite (got {value!r}).")
    return number


def _positive(value: float, label: str, owner: str) -> float:
    number = _finite(value, label, owner)
    if number <= 0.0:
        raise InputRangeError(f"Component {owner!r}: {label} must be positive (got {value!r}).")
    return number
