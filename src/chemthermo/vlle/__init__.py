"""VLLE plugin boundary for chemthermo (DEPRECATED).

This subpackage predates the in-tree multiphase equilibrium engine
(ADR-0011): ``flash_tp(..., flash_mode="modified-raoult")`` with
``FlashSettings(max_phases=3)`` (the default) now *discovers* vapor-liquid-
liquid equilibrium directly - no external engine is required. Keeping this
boundary as the documented way to reach "VLLE support" is therefore
misleading, so it is deprecated (ADR-0013): importing it emits a
``DeprecationWarning`` naming the in-tree replacement. The names below remain
importable for one deprecation cycle; removal needs its own ADR once no
internal user remains.
"""

import warnings

from .api import VLLEEngine
from .errors import VLLEError, VLLEPluginError, VLLEPluginNotInstalledError
from .loader import get_vlle_engine
from .types import VLLEInputs, VLLEPhase, VLLEResult

warnings.warn(
    "chemthermo.vlle is deprecated: three-phase (vapor-liquid-liquid) "
    "equilibrium is now discovered automatically in-tree by "
    "flash_tp(..., flash_mode='modified-raoult') with "
    "FlashSettings(max_phases=3) (the default). See ADR-0011 and ADR-0013. "
    "This subpackage remains importable for one deprecation cycle.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "VLLEEngine",
    "VLLEError",
    "VLLEInputs",
    "VLLEPhase",
    "VLLEResult",
    "VLLEPluginError",
    "VLLEPluginNotInstalledError",
    "get_vlle_engine",
]
