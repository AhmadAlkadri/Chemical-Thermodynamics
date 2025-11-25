"""Abstractions for equilibrium calculations.

This module defines high-level classes that model mixtures, components,
thermodynamic models, and the system-level helper that orchestrates flash
calculations. The goal is to provide clear interfaces that later
implementations can build upon.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from numbers import Number
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd


# Default location for the packaged component database.
DEFAULT_DATABASE_PATH = Path(__file__).resolve().parent.parent / "database" / "database.h5"
R_J_PER_MOL_K = 8.314462618  # J·mol⁻¹·K⁻¹


@dataclass
class Component:
    """Represents a single chemical component.

    Only a handful of fields are defined now; additional properties (critical
    temperatures, Antoine parameters, etc.) can be added as the database schema
    is formalized.
    """

    name: str
    metadata: Dict[str, float] = field(default_factory=dict)

    @classmethod
    def from_database(
        cls, name: str, database_path: Path | None = None
    ) -> "Component":
        """Load component metadata from the HDF database if available.

        The current implementation is intentionally lightweight and forgiving:
        missing data results in an empty metadata dictionary rather than an
        exception. Downstream code can make stricter assertions once the schema
        stabilizes. When no ``database_path`` is provided, the bundled
        :data:`DEFAULT_DATABASE_PATH` is used.
        """

        store_path = (database_path or DEFAULT_DATABASE_PATH).expanduser().resolve()
        if not store_path.exists():
            return cls(name=name)

        metadata: Dict[str, float] = {}
        with pd.HDFStore(store_path) as store:
            if "all" in store:
                table = store["all"]
                match = table[table["Name"].str.match(name, case=False)]
                if not match.empty:
                    row = match.iloc[0]
                    metadata = {}
                    for key, value in row.items():
                        if isinstance(value, (int, float)):
                            metadata[key] = float(value)
                        elif (
                            key == "Antoine Constants"
                            and isinstance(value, (list, tuple))
                            and len(value) == 5
                        ):
                            A, B, C, Tmin, Tmax = value
                            metadata.update(
                                {
                                    "Antoine_A": float(A),
                                    "Antoine_B": float(B),
                                    "Antoine_C": float(C),
                                    "Antoine_Tmin[K]": float(Tmin),
                                    "Antoine_Tmax[K]": float(Tmax),
                                }
                            )
            if "/Cp/gases" in store:
                cp_table = store["/Cp/gases"]
                cp_match = cp_table[cp_table["Name"].str.match(name, case=False)]
                if not cp_match.empty:
                    cp_row = cp_match.iloc[0]
                    metadata.update(
                        {
                            "CpA": float(cp_row["A"]),
                            "CpB": float(cp_row["B*(10**3)"]) * 1e-3,
                            "CpC": float(cp_row["C*(10**6)"]) * 1e-6,
                            "CpD": float(cp_row["D*(10**(-5))"]) * 1e-5,
                            "CpE": float(cp_row["E*(10**9)"]) * 1e-9,
                            "Cp_Tmin[K]": float(cp_row["Tmin"]),
                            "Cp_Tmax[K]": float(cp_row["Tmax"]),
                        }
                    )
        return cls(name=name, metadata=metadata)

    def saturation_pressure(self, temperature: float, units: str = "bar") -> float:
        """Return vapor pressure from the stored Antoine parameters.

        The database provides Antoine coefficients as :math:`\\ln P_{\\text{bar}} = A - B / (T + C)`,
        valid over ``Antoine_Tmin[K]``–``Antoine_Tmax[K]``. Results are returned in
        ``units`` (``"kPa"`` or ``"bar"``). A ``ValueError`` is raised when the required
        coefficients are missing.
        """

        units_normalized = units.lower()
        try:
            A = self.metadata["Antoine_A"]
            B = self.metadata["Antoine_B"]
            C = self.metadata["Antoine_C"]
        except KeyError as exc:
            raise ValueError(
                f"Antoine parameters unavailable for component '{self.name}'."
            ) from exc

        tmin = self.metadata.get("Antoine_Tmin[K]")
        tmax = self.metadata.get("Antoine_Tmax[K]")
        if tmin is not None and tmax is not None:
            if not (tmin <= temperature <= tmax):
                raise ValueError(
                    f"Temperature {temperature} K outside Antoine valid range "
                    f"{tmin}–{tmax} K for '{self.name}'."
                )

        pressure_bar = float(math.exp(A - B / (float(temperature) + C)))

        if units_normalized == "kpa":
            return pressure_bar * 100.0
        if units_normalized == "bar":
            return pressure_bar
        raise ValueError(f"Unsupported pressure units '{units}'. Expected 'kPa' or 'bar'.")

    def heat_capacity(self, temperature: float, units: str = "J/mol/K") -> float:
        """Return ideal-gas heat capacity using the tabulated correlation.

        The correlation follows ``Cp/R = A + B*T + C*T^2 + D*T^-2 + E*T^3`` with
        ``T`` in Kelvin. Coefficients are read from the ``/Cp/gases`` table and
        automatically rescaled from the tabulated units (see ``database/Cp_gas.txt``).
        Supported output units are ``J/mol/K`` (default), ``kJ/mol/K``, and ``R``.
        Raises ``ValueError`` if coefficients are missing or the temperature is
        outside the stated validity range.
        """

        units_normalized = units.lower()
        try:
            A = self.metadata["CpA"]
            B = self.metadata["CpB"]
            C = self.metadata["CpC"]
            D = self.metadata["CpD"]
            E = self.metadata["CpE"]
        except KeyError as exc:
            raise ValueError(
                f"Heat capacity coefficients unavailable for component '{self.name}'."
            ) from exc

        tmin = self.metadata.get("Cp_Tmin[K]")
        tmax = self.metadata.get("Cp_Tmax[K]")
        if tmin is not None and tmax is not None:
            if not (tmin <= temperature <= tmax):
                raise ValueError(
                    f"Temperature {temperature} K outside Cp correlation range "
                    f"{tmin}–{tmax} K for '{self.name}'."
                )

        T = float(temperature)
        cp_over_r = A + B * T + C * T**2 + D * T**-2 + E * T**3

        if units_normalized in {"r"}:
            return cp_over_r
        if units_normalized in {"j/mol/k", "j"}:
            return cp_over_r * R_J_PER_MOL_K
        if units_normalized in {"kj/mol/k", "kj"}:
            return cp_over_r * R_J_PER_MOL_K / 1000.0

        raise ValueError(
            f"Unsupported heat capacity units '{units}'. Expected 'J/mol/K', 'kJ/mol/K', or 'R'."
        )


@dataclass
class Mixture:
    """A collection of components with associated mole fractions."""

    components: List[Component]
    mole_fractions: List[float]

    @classmethod
    def from_names(
        cls,
        names: Iterable[str],
        mole_fractions: Iterable[float],
        database_path: Path | None = None,
    ) -> "Mixture":
        """Build a mixture by looking up each component in the database.

        If ``database_path`` is omitted, the function falls back to
        :data:`DEFAULT_DATABASE_PATH`, allowing quickstart examples to avoid
        hard-coding file locations.
        """

        comps = [Component.from_database(name, database_path) for name in names]
        fractions = list(mole_fractions)
        total = sum(fractions)
        normalized = [x / total for x in fractions] if total else fractions
        return cls(components=comps, mole_fractions=normalized)

    def as_dataframe(self) -> pd.DataFrame:
        """Return component data as a tidy DataFrame for quick inspection."""

        rows = []
        for comp, frac in zip(self.components, self.mole_fractions):
            row = {"name": comp.name, "z": frac}
            row.update(comp.metadata)
            rows.append(row)
        return pd.DataFrame(rows)


class EquationOfState:
    """Base class for cubic or other equation-of-state models."""

    name: str = "generic-EOS"

    def compute_fugacity_coefficients(self, mixture: Mixture, temperature: float, pressure: float):
        """Return placeholder fugacity coefficients.

        Concrete subclasses should provide real implementations.
        """

        raise NotImplementedError


class ActivityModel:
    """Base class for excess-Gibbs energy/activity models."""

    name: str = "generic-activity"

    def compute_activity_coefficients(self, mixture: Mixture, temperature: float):
        raise NotImplementedError


@dataclass
class FlashResult:
    """Container for the outcome of a flash calculation."""

    temperature: float
    pressure: float
    phases: Dict[str, pd.DataFrame] = field(default_factory=dict)
    model_details: Dict[str, str] = field(default_factory=dict)

    def available_specific_properties(self) -> List[str]:
        """List the specific properties computed for this flash result.

        The list combines system-level properties (temperature and pressure)
        with any per-phase columns that convey thermophysical quantities such as
        compressibility, molar volume, density, or composition. Non-physical
        identifiers like component names and phase labels are omitted.
        """

        properties = {"temperature", "pressure"}
        for phase_df in self.phases.values():
            properties.update(
                col
                for col in phase_df.columns
                if col not in {"name", "phase"} and not phase_df[col].isna().all()
            )
        return sorted(properties)

    def specific_properties(self) -> pd.DataFrame:
        """Return the specific property values in a tidy table.

        System-level properties (``temperature`` and ``pressure``) are included
        alongside phase-level quantities. If a property varies by component
        (e.g., composition ``z``), the per-component mapping is returned in the
        ``value`` field to keep the output compact. A ``units`` column is
        included when units can be inferred from the column name (e.g.
        ``Vm[m3/mol]``) or from a small set of known unitless quantities.
        """

        def _split_units(prop_name: str) -> tuple[str, Optional[str]]:
            if "[" in prop_name and prop_name.endswith("]"):
                start = prop_name.rfind("[")
                return prop_name[:start].strip(), prop_name[start + 1 : -1]
            return prop_name, None

        def _units_for(base_prop: str, parsed_unit: Optional[str]) -> Optional[str]:
            if parsed_unit:
                return parsed_unit
            return {
                "temperature": "K",
                "pressure": "bar",
                "z": "mole fraction",
                "Z": "dimensionless",
                "omega": "dimensionless",
            }.get(base_prop)

        rows = [
            {
                "phase": "system",
                "property": "temperature",
                "units": _units_for("temperature", None),
                "value": self.temperature,
            },
            {
                "phase": "system",
                "property": "pressure",
                "units": _units_for("pressure", None),
                "value": self.pressure,
            },
        ]

        for phase, phase_df in self.phases.items():
            for col in phase_df.columns:
                if col in {"name", "phase"}:
                    continue

                base_prop, parsed_unit = _split_units(col)
                units = _units_for(base_prop, parsed_unit)
                series = phase_df[col].dropna()
                if series.empty:
                    continue

                unique_values = series.unique()
                if len(unique_values) == 1:
                    value = unique_values[0]
                    if isinstance(value, Number):
                        value = float(value)
                    rows.append(
                        {"phase": phase, "property": col, "units": units, "value": value}
                    )
                    continue

                component_values = (
                    phase_df[["name", col]]
                    .dropna()
                    .set_index("name")[col]
                    .apply(lambda v: float(v) if isinstance(v, Number) else v)
                    .to_dict()
                )
                rows.append(
                    {
                        "phase": phase,
                        "property": f"{col} (by component)",
                        "units": units,
                        "value": component_values,
                    }
                )

        return pd.DataFrame(rows)


class FlashCalculator:
    """Strategy interface for performing flash calculations."""

    def flash(
        self, mixture: Mixture, temperature: float, pressure: float
    ) -> FlashResult:
        raise NotImplementedError


@dataclass
class SimpleFlashCalculator(FlashCalculator):
    """Very light-weight placeholder flash implementation.

    The purpose of this stub is to demonstrate the expected data flow without
    committing to a particular thermodynamic model. It currently assumes a
    single vapor phase and echoes back the feed composition.
    """

    eos_model: Optional[EquationOfState] = None
    activity_model: Optional[ActivityModel] = None

    def flash(
        self, mixture: Mixture, temperature: float, pressure: float
    ) -> FlashResult:
        # The eventual implementation should detect multi-phase regions and
        # allocate compositions accordingly. For now, keep a single phase.
        vapor_df = mixture.as_dataframe().copy()
        vapor_df["phase"] = "vapor"

        if self.eos_model is not None:
            eos_result = self.eos_model.compute_volumes(
                mixture=mixture, temperature=temperature, pressure=pressure
            )
            if eos_result.z_factors:
                vapor_z = max(eos_result.z_factors)
                vapor_df["Z"] = vapor_z
                vapor_df["Vm[m3/mol]"] = eos_result.preferred_phase_volume("vapor")
                density = eos_result.densities_kg_per_m3[
                    eos_result.z_factors.index(vapor_z)
                ]
                vapor_df["rho[kg/m3]"] = density
                vapor_df["v_specific[m3/kg]"] = float("nan") if density == 0 else 1.0 / density

        model_details: Dict[str, str] = {}
        if self.eos_model:
            model_details["eos"] = self.eos_model.name
        if self.activity_model:
            model_details["activity_model"] = self.activity_model.name

        return FlashResult(
            temperature=temperature,
            pressure=pressure,
            phases={"vapor": vapor_df},
            model_details=model_details,
        )


@dataclass
class NPTEquilibriumSystem:
    """Represents a system held at constant pressure and temperature.

    This class bundles together the state variables, the mixture description,
    and the flash calculation strategy. Real calculations can swap in more
    sophisticated flash calculators or EOS/activity models without altering
    client code.
    """

    temperature: float
    pressure: float
    mixture: Mixture
    flash_calculator: FlashCalculator = field(default_factory=SimpleFlashCalculator)

    def set_state(self, temperature: Optional[float] = None, pressure: Optional[float] = None) -> None:
        """Update temperature and/or pressure in place."""

        if temperature is not None:
            self.temperature = float(temperature)
        if pressure is not None:
            self.pressure = float(pressure)

    def flash(self) -> FlashResult:
        """Perform an isothermal, isobaric flash calculation.

        The default implementation delegates to the configured
        :class:`FlashCalculator`. Advanced calculators can compute phase splits,
        enforce material balances, and iterate until fugacity/activity criteria
        are satisfied.
        """

        return self.flash_calculator.flash(
            mixture=self.mixture, temperature=self.temperature, pressure=self.pressure
        )
