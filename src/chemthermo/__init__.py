"""Chemical engineering thermodynamics package (work in progress)."""

from .core import Component, Composition, CompositionBasis, Mixture
from .exceptions import (
    CompositionError,
    ConvergenceError,
    InputRangeError,
    ModelError,
    PropertyNotFoundError,
    ThermoError,
)
from .flash import FlashResult, FlashSettings, PhaseResult, flash_tp
from .models import ActivityModel, EquationOfState, PengRobinsonEOS
from .units import (
    BAR_PER_PA,
    PA_PER_BAR,
    PA_PER_KPA,
    PA_PER_MPA,
    R_J_PER_MOL_K,
    STANDARD_P_PA,
    STANDARD_T_K,
)
from .validation import (
    COMPOSITION_SUM_TOL,
    validate_fractions,
    validate_pressure,
    validate_temperature,
)

__all__ = [
    "__version__",
    "BAR_PER_PA",
    "COMPOSITION_SUM_TOL",
    "Component",
    "CompositionError",
    "Composition",
    "CompositionBasis",
    "ConvergenceError",
    "EquationOfState",
    "FlashResult",
    "FlashSettings",
    "PhaseResult",
    "InputRangeError",
    "Mixture",
    "ModelError",
    "PA_PER_BAR",
    "PA_PER_KPA",
    "PA_PER_MPA",
    "PropertyNotFoundError",
    "PengRobinsonEOS",
    "R_J_PER_MOL_K",
    "STANDARD_P_PA",
    "STANDARD_T_K",
    "ThermoError",
    "ActivityModel",
    "flash_tp",
    "validate_fractions",
    "validate_pressure",
    "validate_temperature",
]

__version__ = "0.0.0"
