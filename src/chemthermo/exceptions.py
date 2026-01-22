"""Custom exception types for chemthermo."""


class ThermoError(Exception):
    """Base class for chemthermo errors."""


class InputRangeError(ThermoError):
    """Raised when inputs are outside physical or supported ranges."""


class CompositionError(ThermoError):
    """Raised when a composition is invalid or inconsistent."""


class PropertyNotFoundError(ThermoError):
    """Raised when required properties are missing from data sources."""


class ModelError(ThermoError):
    """Raised when a thermodynamic model fails or is misused."""


class ConvergenceError(ThermoError):
    """Raised when a numerical procedure fails to converge."""
