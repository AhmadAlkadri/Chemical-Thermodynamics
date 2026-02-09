"""Chemical engineering thermodynamics utilities."""

from __future__ import annotations

from .core.component import Component
from .core.composition import Composition
from .core.mixture import Mixture

from .data import list_component_names
from .exceptions import CompositionError, ConvergenceError, InputRangeError, ModelError, PropertyNotFoundError, ThermoError as ChemThermoError
from .eos import get_eos, list_eos
from .flash import FlashResult, FlashSettings, PhaseResult, flash_tp
from .models import ActivityModel, EquationOfState, NRTL, PengRobinsonEOS
from .parameters import ActivityParameters, NRTLParameters, PCSAFTParameterError
from .units import PA_PER_BAR, STANDARD_P_PA
from .validation import validate_fractions, validate_pressure, validate_temperature

__version__ = "0.0.0"

__all__ = [
    "ChemThermoError",
    "Component",
    "Mixture",

    "InputRangeError",
    "PropertyNotFoundError",
    "flash_tp",
    "get_eos",
    "list_eos",
    "ActivityModel",
    "EquationOfState",
    "NRTL",
    "PengRobinsonEOS",
    "ActivityParameters",
    "NRTLParameters",
    "CompositionError",
    "PA_PER_BAR",
    "STANDARD_P_PA",
    "validate_fractions",
    "Composition",
    "PCSAFTParameterError",
    "ConvergenceError",
    "ModelError",
    "FlashResult",
    "FlashSettings",
    "PhaseResult",
    "validate_pressure",
    "validate_temperature",
    "list_component_names",
    "cite",
]


def cite(component_name: str, property_name: str) -> str:
    """Return the citation for a specific component property.
    
    Args:
        component_name: The name of the component (e.g., 'Methane').
        property_name: The property to cite (e.g., 'Tc', 'MW', 'antoine').
        
    Returns:
        The citation text from the bibliography.
    """
    try:
        comp = Component.from_database(component_name)
        return comp.get_citation(property_name)
    except Exception as exc:
        return f"Could not retrieve citation: {exc}"
