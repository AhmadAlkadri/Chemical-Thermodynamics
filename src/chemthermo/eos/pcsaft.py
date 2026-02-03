"""PC-SAFT EOS scaffolding (placeholder implementation)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ..exceptions import ModelError
from ..parameters import PCSAFTParameterRegistry, get_pcsaft_parameters
from .api import EOSProtocol
from .registry import register_eos


@dataclass(frozen=True)
class PCSAFTEOS(EOSProtocol):
    """Placeholder PC-SAFT EOS interface.

    The open-source build ships without PC-SAFT parameters or physics.
    """

    components: tuple[str, ...]
    parameters: PCSAFTParameterRegistry | None = None
    name: str = "PC-SAFT"

    def __post_init__(self) -> None:
        if not self.components:
            raise ModelError("PC-SAFT requires at least one component.")
        object.__setattr__(self, "components", tuple(self.components))

    def num_components(self) -> int:
        return len(self.components)

    def residual_helmholtz(
        self,
        *,
        temperature_K: float,
        volume_m3: float,
        composition: Sequence[float],
    ) -> float:
        registry = self.parameters or get_pcsaft_parameters()
        registry.for_components(self.components)
        raise NotImplementedError(
            "PC-SAFT residual Helmholtz energy is not implemented in the open-source build."
        )


register_eos("pcsaft", PCSAFTEOS)
