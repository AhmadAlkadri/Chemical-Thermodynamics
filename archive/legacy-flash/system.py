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
import numpy as np


# Default location for the packaged component database.
DEFAULT_DATABASE_PATH = Path(__file__).resolve().parent.parent / "database" / "database.h5"
R_J_PER_MOL_K = 8.314462618  # J·mol⁻¹·K⁻¹
REF_T_K = 298.15
REF_P_BAR = 1.0


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

    def ideal_gas_enthalpy(self, temperature: float) -> float:
        """Ideal-gas enthalpy relative to ``REF_T_K`` (J/mol)."""

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

        T = float(temperature)
        T0 = REF_T_K

        term_a = A * (T - T0)
        term_b = 0.5 * B * (T**2 - T0**2)
        term_c = (1.0 / 3.0) * C * (T**3 - T0**3)
        term_d = -D * (1.0 / T - 1.0 / T0)
        term_e = 0.25 * E * (T**4 - T0**4)

        delta_h_over_r = term_a + term_b + term_c + term_d + term_e
        return delta_h_over_r * R_J_PER_MOL_K

    def ideal_gas_entropy(self, temperature: float, pressure: float = REF_P_BAR) -> float:
        """Ideal-gas entropy relative to ``REF_T_K`` and ``REF_P_BAR`` (J/mol/K)."""

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

        T = float(temperature)
        T0 = REF_T_K
        P = float(pressure)
        P0 = REF_P_BAR

        term_a = A * math.log(T / T0)
        term_b = B * (T - T0)
        term_c = 0.5 * C * (T**2 - T0**2)
        term_d = -0.5 * D * (1.0 / T**2 - 1.0 / T0**2)
        term_e = (1.0 / 3.0) * E * (T**3 - T0**3)

        delta_s_over_r = term_a + term_b + term_c + term_d + term_e - math.log(P / P0)
        return delta_s_over_r * R_J_PER_MOL_K

    def ideal_gas_internal_energy(self, temperature: float) -> float:
        """Ideal-gas internal energy relative to ``REF_T_K`` (J/mol)."""

        h = self.ideal_gas_enthalpy(temperature)
        return h - R_J_PER_MOL_K * float(temperature)


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

    def ideal_enthalpy(self, temperature: float) -> float:
        """Return mixture ideal-gas enthalpy (J/mol) relative to ``REF_T_K``."""

        h = 0.0
        for comp, frac in zip(self.components, self.mole_fractions):
            h += frac * comp.ideal_gas_enthalpy(temperature)
        return h

    def ideal_entropy(self, temperature: float, pressure: float) -> float:
        """Return mixture ideal-gas entropy (J/mol/K) relative to ``REF_T_K`` and ``REF_P_BAR``."""

        s = 0.0
        for comp, frac in zip(self.components, self.mole_fractions):
            s += frac * comp.ideal_gas_entropy(temperature, pressure)
        return s - R_J_PER_MOL_K * sum(
            frac * math.log(frac) for frac in self.mole_fractions if frac > 0
        )

    def ideal_internal_energy(self, temperature: float) -> float:
        """Return mixture ideal-gas internal energy (J/mol) relative to ``REF_T_K``."""

        u = 0.0
        for comp, frac in zip(self.components, self.mole_fractions):
            u += frac * comp.ideal_gas_internal_energy(temperature)
        return u


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
    committing to a particular thermodynamic model. With an EOS attached, a
    simple phi–phi flash is performed; otherwise, the feed is echoed back as a
    single vapor phase.
    """

    eos_model: Optional[EquationOfState] = None
    activity_model: Optional[ActivityModel] = None

    def _wilson_K(self, mixture: Mixture, temperature: float, pressure: float) -> np.ndarray:
        """Wilson correlation for initial K-values."""

        Ks = []
        for comp in mixture.components:
            try:
                Tc = float(comp.metadata["Tc[K]"])
                Pc = float(comp.metadata["Pc[bar]"])
                omega = float(comp.metadata.get("omega", 0.0))
            except KeyError as exc:
                raise ValueError(
                    f"Critical properties required for Wilson correlation are missing for {comp.name}."
                ) from exc

            ln_K = math.log(Pc / pressure) + 5.373 * (1 + omega) * (1 - Tc / temperature)
            Ks.append(math.exp(ln_K))
        return np.array(Ks, dtype=float)

    def _rachford_rice(self, z: np.ndarray, K: np.ndarray, tol: float = 1e-8) -> Optional[float]:
        """Solve Rachford–Rice for vapor fraction in [0, 1]."""

        def f(v):
            return np.sum(z * (K - 1) / (1 + v * (K - 1)))

        f0 = f(0.0)
        f1 = f(1.0)
        if f0 * f1 > 0:
            return None

        lower, upper = 0.0, 1.0
        for _ in range(100):
            mid = 0.5 * (lower + upper)
            val = f(mid)
            if abs(val) < tol:
                return mid
            if val * f0 > 0:
                lower = mid
                f0 = val
            else:
                upper = mid
        return mid

    def _mixture_properties(
        self,
        mixture: Mixture,
        composition: np.ndarray,
        temperature: float,
        pressure: float,
        phase_label: str,
    ) -> Dict[str, float]:
        """Return Z, Vm, rho, h, s, u, ln_phi for a phase."""

        if self.eos_model is None:
            raise ValueError("An EOS model is required for property calculations.")

        phi, Z, Vm, rho, HR, SR = self.eos_model.fugacity_coefficients(
            mixture=mixture,
            composition=composition,
            temperature=temperature,
            pressure=pressure,
            phase_root="vapor" if phase_label == "vapor" else "liquid",
        )

        phase_mix = Mixture(components=mixture.components, mole_fractions=list(composition))
        h_ideal = phase_mix.ideal_enthalpy(temperature)
        s_ideal = phase_mix.ideal_entropy(temperature, pressure)
        u_ideal = phase_mix.ideal_internal_energy(temperature)

        h = h_ideal + HR
        s = s_ideal + SR
        u = h - pressure * 1e5 * Vm  # P in bar -> Pa

        return {
            "Z": Z,
            "Vm[m3/mol]": Vm,
            "rho[kg/m3]": rho,
            "h[J/mol]": h,
            "s[J/mol/K]": s,
            "u[J/mol]": u,
            "phi": phi,
        }

    def _single_phase_result(
        self, mixture: Mixture, temperature: float, pressure: float, phase: str
    ) -> FlashResult:
        """Fallback single-phase result (ideal gas if no EOS)."""

        df = mixture.as_dataframe().copy()
        df["phase"] = phase
        if self.eos_model is not None:
            props = self._mixture_properties(
                mixture=mixture,
                composition=np.array(mixture.mole_fractions, dtype=float),
                temperature=temperature,
                pressure=pressure,
                phase_label=phase,
            )
            df["Z"] = props["Z"]
            df["Vm[m3/mol]"] = props["Vm[m3/mol]"]
            df["rho[kg/m3]"] = props["rho[kg/m3]"]
            df["h[J/mol]"] = props["h[J/mol]"]
            df["s[J/mol/K]"] = props["s[J/mol/K]"]
            df["u[J/mol]"] = props["u[J/mol]"]
        else:
            h = mixture.ideal_enthalpy(temperature)
            s = mixture.ideal_entropy(temperature, pressure)
            u = mixture.ideal_internal_energy(temperature)
            df["h[J/mol]"] = h
            df["s[J/mol/K]"] = s
            df["u[J/mol]"] = u

        # expose compositions explicitly
        if phase.lower().startswith("vapor"):
            df["y"] = df["z"]
            df["x"] = df["z"]
        else:
            df["x"] = df["z"]
            df["y"] = df["z"]

        model_details: Dict[str, str] = {}
        if self.eos_model:
            model_details["eos"] = self.eos_model.name
        if self.activity_model:
            model_details["activity_model"] = self.activity_model.name

        return FlashResult(
            temperature=temperature,
            pressure=pressure,
            phases={phase: df},
            model_details=model_details,
        )

    def flash(
        self, mixture: Mixture, temperature: float, pressure: float
    ) -> FlashResult:
        if self.eos_model is None:
            return self._single_phase_result(mixture, temperature, pressure, "vapor")

        z = np.array(mixture.mole_fractions, dtype=float)
        if np.allclose(z, 0):
            return self._single_phase_result(mixture, temperature, pressure, "vapor")

        K = self._wilson_K(mixture, temperature, pressure)

        # Detect single-phase conditions
        if np.all(K <= 1):
            return self._single_phase_result(mixture, temperature, pressure, "liquid")
        if np.all(K >= 1):
            return self._single_phase_result(mixture, temperature, pressure, "vapor")

        vapor_fraction = self._rachford_rice(z, K)
        if vapor_fraction is None:
            # No split possible; revert to vapor phase
            return self._single_phase_result(mixture, temperature, pressure, "vapor")

        for _ in range(50):
            denom = 1 + vapor_fraction * (K - 1)
            x = z / denom
            x /= x.sum()
            y = K * x
            y /= y.sum()

            phi_v, Z_v, Vm_v, rho_v, HR_v, SR_v = self.eos_model.fugacity_coefficients(
                mixture=mixture,
                composition=y,
                temperature=temperature,
                pressure=pressure,
                phase_root="vapor",
            )
            phi_l, Z_l, Vm_l, rho_l, HR_l, SR_l = self.eos_model.fugacity_coefficients(
                mixture=mixture,
                composition=x,
                temperature=temperature,
                pressure=pressure,
                phase_root="liquid",
            )

            K_new = phi_l / phi_v
            if np.max(np.abs(K_new - K)) < 1e-6:
                K = K_new
                break
            K = K_new
            rr = self._rachford_rice(z, K)
            if rr is None:
                break
            vapor_fraction = rr

        # Build phase dataframes
        liquid_df = mixture.as_dataframe().copy()
        liquid_df["phase"] = "liquid"
        vapor_df = mixture.as_dataframe().copy()
        vapor_df["phase"] = "vapor"

        liquid_df["x"] = x
        liquid_df["y"] = y
        vapor_df["y"] = y
        vapor_df["x"] = x
        liquid_df["Z"] = Z_l
        vapor_df["Z"] = Z_v
        liquid_df["Vm[m3/mol]"] = Vm_l
        vapor_df["Vm[m3/mol]"] = Vm_v
        liquid_df["rho[kg/m3]"] = rho_l
        vapor_df["rho[kg/m3]"] = rho_v

        liquid_mix = Mixture(components=mixture.components, mole_fractions=list(x))
        vapor_mix = Mixture(components=mixture.components, mole_fractions=list(y))

        h_l_ideal = liquid_mix.ideal_enthalpy(temperature)
        s_l_ideal = liquid_mix.ideal_entropy(temperature, pressure)
        u_l_ideal = liquid_mix.ideal_internal_energy(temperature)

        h_v_ideal = vapor_mix.ideal_enthalpy(temperature)
        s_v_ideal = vapor_mix.ideal_entropy(temperature, pressure)
        u_v_ideal = vapor_mix.ideal_internal_energy(temperature)

        h_l = h_l_ideal + HR_l
        h_v = h_v_ideal + HR_v
        s_l = s_l_ideal + SR_l
        s_v = s_v_ideal + SR_v
        u_l = h_l - pressure * 1e5 * Vm_l
        u_v = h_v - pressure * 1e5 * Vm_v

        liquid_df["h[J/mol]"] = h_l
        vapor_df["h[J/mol]"] = h_v
        liquid_df["s[J/mol/K]"] = s_l
        vapor_df["s[J/mol/K]"] = s_v
        liquid_df["u[J/mol]"] = u_l
        vapor_df["u[J/mol]"] = u_v

        model_details: Dict[str, str] = {"vapor_fraction": f"{vapor_fraction:.5f}"}
        if self.eos_model:
            model_details["eos"] = self.eos_model.name
        if self.activity_model:
            model_details["activity_model"] = self.activity_model.name

        return FlashResult(
            temperature=temperature,
            pressure=pressure,
            phases={"liquid": liquid_df, "vapor": vapor_df},
            model_details=model_details,
        )

    def pure_PT_envelope(
        self, component: Component, temperatures: Iterable[float]
    ) -> pd.DataFrame:
        """Compute a saturation P–T curve for a pure component using phi–phi."""

        if self.eos_model is None:
            raise ValueError("EOS model required for P–T envelope calculation.")

        rows = []
        for T in temperatures:
            # Antoine as initial guess; fall back to simple estimate
            try:
                P = component.saturation_pressure(T, units="bar")
            except Exception:
                P = float(component.metadata.get("Pc[bar]", 1.0)) * 0.5

            mix = Mixture(components=[component], mole_fractions=[1.0])
            for _ in range(20):
                phi_v, Z_v, Vm_v, rho_v, HR_v, SR_v = self.eos_model.fugacity_coefficients(
                    mixture=mix, composition=[1.0], temperature=T, pressure=P, phase_root="vapor"
                )
                phi_l, Z_l, Vm_l, rho_l, HR_l, SR_l = self.eos_model.fugacity_coefficients(
                    mixture=mix, composition=[1.0], temperature=T, pressure=P, phase_root="liquid"
                )
                ratio = phi_l[0] / phi_v[0]
                if abs(ratio - 1.0) < 1e-6:
                    break
                P *= ratio

            rows.append(
                {
                    "T[K]": T,
                    "P[bar]": P,
                    "Z_v": Z_v,
                    "Z_l": Z_l,
                    "Vm_v[m3/mol]": Vm_v,
                    "Vm_l[m3/mol]": Vm_l,
                    "rho_v[kg/m3]": rho_v,
                    "rho_l[kg/m3]": rho_l,
                }
            )

        return pd.DataFrame(rows)

    def binary_Tx_diagram(
        self, components: List[Component], pressure: float, x1_grid: Iterable[float]
    ) -> pd.DataFrame:
        """Generate bubble/dew temperatures for a binary at fixed pressure."""

        if len(components) != 2:
            raise ValueError("binary_Tx_diagram expects exactly two components.")
        if self.eos_model is None:
            raise ValueError("EOS model required for T–x diagram calculation.")

        rows = []
        compA, compB = components
        mix = Mixture(components=components, mole_fractions=[0.5, 0.5])

        def bubble_temp(x):
            Tc_vals = [float(compA.metadata["Tc[K]"]), float(compB.metadata["Tc[K]"])]
            T_low = 0.5 * min(Tc_vals)
            T_high = 0.99 * max(Tc_vals)
            for _ in range(60):
                T_mid = 0.5 * (T_low + T_high)
                K = self._wilson_K(mix, T_mid, pressure)
                for _ in range(20):
                    y = K * x
                    y /= y.sum()
                    phi_v, _, _, _, _, _ = self.eos_model.fugacity_coefficients(
                        mixture=mix, composition=y, temperature=T_mid, pressure=pressure, phase_root="vapor"
                    )
                    phi_l, _, _, _, _, _ = self.eos_model.fugacity_coefficients(
                        mixture=mix, composition=x, temperature=T_mid, pressure=pressure, phase_root="liquid"
                    )
                    K_new = phi_l / phi_v
                    if np.max(np.abs(K_new - K)) < 1e-6:
                        K = K_new
                        break
                    K = K_new
                fval = np.sum(x * K) - 1.0
                if abs(fval) < 1e-6:
                    return T_mid, K
                if fval > 0:
                    T_low = T_mid
                else:
                    T_high = T_mid
            return T_mid, K

        def dew_temp(y):
            Tc_vals = [float(compA.metadata["Tc[K]"]), float(compB.metadata["Tc[K]"])]
            T_low = 0.5 * min(Tc_vals)
            T_high = 0.99 * max(Tc_vals)
            for _ in range(60):
                T_mid = 0.5 * (T_low + T_high)
                K = self._wilson_K(mix, T_mid, pressure)
                for _ in range(20):
                    x = y / K
                    x /= x.sum()
                    phi_v, _, _, _, _, _ = self.eos_model.fugacity_coefficients(
                        mixture=mix, composition=y, temperature=T_mid, pressure=pressure, phase_root="vapor"
                    )
                    phi_l, _, _, _, _, _ = self.eos_model.fugacity_coefficients(
                        mixture=mix, composition=x, temperature=T_mid, pressure=pressure, phase_root="liquid"
                    )
                    K_new = phi_l / phi_v
                    if np.max(np.abs(K_new - K)) < 1e-6:
                        K = K_new
                        break
                    K = K_new
                fval = np.sum(y / K) - 1.0
                if abs(fval) < 1e-6:
                    return T_mid, K
                if fval > 0:
                    T_high = T_mid
                else:
                    T_low = T_mid
            return T_mid, K

        for x1 in x1_grid:
            x = np.array([x1, 1 - x1], dtype=float)
            Tb, Kb = bubble_temp(x)
            yb = (Kb * x) / np.sum(Kb * x)

            y = np.array([x1, 1 - x1], dtype=float)
            Td, Kd = dew_temp(y)
            xd = (y / Kd) / np.sum(y / Kd)

            rows.append(
                {
                    "x1": x1,
                    "Tbubble[K]": Tb,
                    "y1_bubble": yb[0],
                    "Tdew[K]": Td,
                    "x1_dew": xd[0],
                }
            )

        return pd.DataFrame(rows)


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
