"""Abstractions for equilibrium calculations.

This module defines high-level classes that model mixtures, components,
thermodynamic models, and the system-level helper that orchestrates flash
calculations. The goal is to provide clear interfaces that later
implementations can build upon.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd


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
    def from_database(cls, name: str, database_path: Path) -> "Component":
        """Load component metadata from the HDF database if available.

        The current implementation is intentionally lightweight and forgiving:
        missing data results in an empty metadata dictionary rather than an
        exception. Downstream code can make stricter assertions once the schema
        stabilizes.
        """

        store_path = database_path.expanduser().resolve()
        if not store_path.exists():
            return cls(name=name)

        metadata: Dict[str, float] = {}
        with pd.HDFStore(store_path) as store:
            if "all" in store:
                table = store["all"]
                match = table[table["Name"].str.match(name, case=False)]
                if not match.empty:
                    # For now keep the entire row except for string identifiers.
                    row = match.iloc[0]
                    metadata = {
                        key: float(value)
                        for key, value in row.items()
                        if isinstance(value, (int, float))
                    }
        return cls(name=name, metadata=metadata)


@dataclass
class Mixture:
    """A collection of components with associated mole fractions."""

    components: List[Component]
    mole_fractions: List[float]

    @classmethod
    def from_names(
        cls, names: Iterable[str], mole_fractions: Iterable[float], database_path: Path
    ) -> "Mixture":
        """Build a mixture by looking up each component in the database."""

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
