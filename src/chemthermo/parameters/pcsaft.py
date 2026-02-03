"""PC-SAFT parameter registry stub."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Mapping, Sequence

from ..data import normalize_name
from ..exceptions import ModelError


class PCSAFTParameterError(ModelError):
    """Raised when PC-SAFT parameters are missing or invalid."""


@dataclass(frozen=True)
class PCSAFTParameterRegistry:
    """Registry for PC-SAFT parameters.

    This stub stores per-component parameter dictionaries keyed by normalized
    component name. The open-source build ships without PC-SAFT parameters.
    """

    parameters: Mapping[str, Mapping[str, float]] = field(default_factory=dict)

    def for_components(self, components: Sequence[str]) -> Mapping[str, Mapping[str, float]]:
        if not components:
            raise PCSAFTParameterError("PC-SAFT requires at least one component.")
        if not self.parameters:
            raise PCSAFTParameterError(
                "PC-SAFT parameters are not available in the open-source build."
            )

        canonical = [normalize_name(name) for name in components]
        missing = [name for name in canonical if name not in self.parameters]
        if missing:
            missing_list = ", ".join(missing)
            raise PCSAFTParameterError(
                f"Missing PC-SAFT parameters for: {missing_list}."
            )

        return {name: self.parameters[name] for name in canonical}


@lru_cache(maxsize=1)
def get_pcsaft_parameters() -> PCSAFTParameterRegistry:
    """Return the default PC-SAFT parameter registry."""

    return PCSAFTParameterRegistry()
