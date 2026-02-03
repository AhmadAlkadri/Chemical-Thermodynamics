"""Chemical engineering thermodynamics package (work in progress)."""

from .core import Component, Composition, CompositionBasis, Mixture
from .eos import PCSAFTEOS, EOSProtocol, get_eos, list_eos, register_eos
from .exceptions import (
    CompositionError,
    ConvergenceError,
    InputRangeError,
    ModelError,
    PropertyNotFoundError,
    ThermoError,
)
from .flash import FlashResult, FlashSettings, PhaseResult, flash_tp
from .models import NRTL, ActivityModel, EquationOfState, PengRobinsonEOS
from .parameters import (
    ActivityParameters,
    NRTLParameters,
    PCSAFTParameterError,
    PCSAFTParameterRegistry,
    get_pcsaft_parameters,
)
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
    "EOSProtocol",
    "FlashResult",
    "FlashSettings",
    "PhaseResult",
    "InputRangeError",
    "Mixture",
    "ModelError",
    "NRTL",
    "NRTLParameters",
    "ActivityParameters",
    "PA_PER_BAR",
    "PA_PER_KPA",
    "PA_PER_MPA",
    "PCSAFTEOS",
    "PCSAFTParameterError",
    "PCSAFTParameterRegistry",
    "PropertyNotFoundError",
    "PengRobinsonEOS",
    "R_J_PER_MOL_K",
    "STANDARD_P_PA",
    "STANDARD_T_K",
    "ThermoError",
    "ActivityModel",
    "flash_tp",
    "get_eos",
    "get_pcsaft_parameters",
    "list_eos",
    "register_eos",
    "validate_fractions",
    "validate_pressure",
    "validate_temperature",
]

__version__ = "0.0.0"
