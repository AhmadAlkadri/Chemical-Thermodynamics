"""Chemical engineering thermodynamics utilities."""

from __future__ import annotations

from .core.component import Component
from .core.mixture import Mixture

from .data import list_component_names
from .exceptions import PropertyNotFoundError, ThermoError as ChemThermoError
from .flash import flash_tp
from .models import NRTL, PengRobinsonEOS

__version__ = "0.0.0"

__all__ = [
    "ChemThermoError",
    "Component",
    "Mixture",

    "PropertyNotFoundError",
    "flash_tp",
    "NRTL",
    "PengRobinsonEOS",
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
