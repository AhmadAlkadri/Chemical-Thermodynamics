"""Pydantic models for the component database schema."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field


class Parameter(BaseModel):
    """A physical property value with units and provenance."""

    value: float
    units: str
    source_key: str | None = Field(default=None, description="BibTeX key or legacy identifier.")
    uncertainty: float | None = None
    method: str | None = Field(default=None, description="Measurement or estimation method.")


class AntoineCoefficients(BaseModel):
    """Antoine equation coefficients for vapor pressure."""

    A: float
    B: float
    C: float
    Tmin_K: float
    Tmax_K: float
    units: Literal["bar", "Pa", "kPa"] = "bar"
    source_key: str | None = None


class ComponentCore(BaseModel):
    """Fields every component record carries, whatever its provenance.

    Split out of :class:`ComponentData` so that a user-supplied component
    (:class:`CustomComponentData`) shares the identity and molar-mass fields
    without relaxing anything the packaged databank promises. ``ComponentData``
    keeps ``Tc`` / ``Pc`` / ``omega`` **required**, so ``Database`` validates
    exactly as it did before (ADR-0022).
    """

    name: str = Field(..., description="Common name (e.g., 'Methane').")
    formula: str = Field(..., description="Chemical formula (e.g., 'CH4').")
    CAS: str | None = Field(default=None, description="CAS Registry Number.")

    MW: Parameter = Field(..., description="Molecular/Atomic Weight.")

    # Temperature dependent models
    antoine: AntoineCoefficients | None = None


class ComponentData(ComponentCore):
    """Schema for a single chemical component record in the packaged databank."""

    # Critical properties
    Tc: Parameter = Field(..., description="Critical Temperature.")
    Pc: Parameter = Field(..., description="Critical Pressure.")
    omega: Parameter = Field(..., description="Acentric Factor.")


class CustomComponentData(ComponentCore):
    """A component described by the caller rather than looked up (ADR-0022).

    The critical constants are **optional** here and only here: a polymer has
    no measurable critical point, and requiring one would mean inventing a
    number to satisfy a schema. ``Component.tc_k`` / ``pc_pa`` / ``omega``
    raise :class:`~chemthermo.exceptions.PropertyNotFoundError` when they are
    absent, so a model that genuinely needs them says so instead of reading a
    placeholder.

    ``volatile`` is the one piece of *modelling* information the flash and the
    stability test need from a component with no critical constants: a
    ``volatile=False`` component is given the fixed Wilson estimate
    ``K = 1e-10`` (``chemthermo.flash._common.NON_VOLATILE_WILSON_K``) instead
    of the Wilson correlation, which puts it at essentially zero in a
    vapour-like trial composition. It is an *initial estimate*, never a model
    statement: nothing downstream reads ``volatile``.

    This model is deliberately **not** a ``ComponentData``, so it cannot enter
    :class:`Database` and the databank's ``schema_version`` is untouched.
    """

    Tc: Parameter | None = None
    Pc: Parameter | None = None
    omega: Parameter | None = None
    volatile: bool = True


class Database(BaseModel):
    """Root container for the component database."""

    schema_version: int = 1
    components: list[ComponentData]
