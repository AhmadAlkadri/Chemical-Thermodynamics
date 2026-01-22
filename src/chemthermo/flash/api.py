"""User-facing flash calculation entry points."""

from __future__ import annotations

from ..core import Mixture
from ..exceptions import ModelError
from ..models import EquationOfState
from ..validation import validate_pressure, validate_temperature
from .results import FlashResult
from .settings import FlashSettings


def flash_tp(
    mixture: Mixture,
    *,
    temperature_K: float,
    pressure_Pa: float,
    eos: EquationOfState | None,
    settings: FlashSettings | None = None,
) -> FlashResult:
    """Perform a TP flash calculation.

    This entry point validates inputs and enforces the error contract. The
    numerical solver implementation will be added in a later milestone.
    """

    validate_temperature(temperature_K)
    validate_pressure(pressure_Pa)

    if eos is None:
        raise ModelError("An equation-of-state model is required for flash_tp.")

    _ = mixture
    _ = settings or FlashSettings()

    raise ModelError("flash_tp solver is not implemented yet.")
