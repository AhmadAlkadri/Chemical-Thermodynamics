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


class ComponentData(BaseModel):
    """Schema for a single chemical component record."""

    name: str = Field(..., description="Common name (e.g., 'Methane').")
    formula: str = Field(..., description="Chemical formula (e.g., 'CH4').")
    CAS: str | None = Field(default=None, description="CAS Registry Number.")

    # Critical properties
    MW: Parameter = Field(..., description="Molecular/Atomic Weight.")
    Tc: Parameter = Field(..., description="Critical Temperature.")
    Pc: Parameter = Field(..., description="Critical Pressure.")
    omega: Parameter = Field(..., description="Acentric Factor.")

    # Temperature dependent models
    antoine: AntoineCoefficients | None = None


class Database(BaseModel):
    """Root container for the component database."""

    schema_version: int = 1
    components: list[ComponentData]
