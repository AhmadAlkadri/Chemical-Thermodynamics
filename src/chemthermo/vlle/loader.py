"""Plugin loader for the VLLE engine (DEPRECATED, see ADR-0013)."""

from __future__ import annotations

import warnings
from importlib import import_module
from typing import Any

from .api import VLLEEngine
from .errors import VLLEPluginError, VLLEPluginNotInstalledError


def get_vlle_engine() -> VLLEEngine:
    """Return the VLLE engine from the optional chemthermo_vlle plugin.

    .. deprecated::
        This plugin boundary is deprecated (ADR-0013). Three-phase
        (vapor-liquid-liquid) equilibrium is now discovered automatically
        in-tree by ``flash_tp(..., flash_mode="modified-raoult")`` with
        ``FlashSettings(max_phases=3)`` (the default); see ADR-0011. Calling
        this function emits a ``DeprecationWarning`` and, absent an installed
        ``chemthermo_vlle`` plugin, still raises
        :class:`VLLEPluginNotInstalledError`.

    The plugin must expose a callable ``get_engine()`` that returns an object
    implementing the ``VLLEEngine`` protocol.
    """

    warnings.warn(
        "chemthermo.vlle.get_vlle_engine() is deprecated: three-phase "
        "(vapor-liquid-liquid) equilibrium is now discovered automatically "
        "in-tree by flash_tp(..., flash_mode='modified-raoult') with "
        "FlashSettings(max_phases=3) (the default). See ADR-0011 and "
        "ADR-0013.",
        DeprecationWarning,
        stacklevel=2,
    )

    try:
        module = import_module("chemthermo_vlle")
    except ModuleNotFoundError as exc:
        raise VLLEPluginNotInstalledError(
            "Install chemthermo_vlle to enable VLLE support."
        ) from exc

    get_engine = getattr(module, "get_engine", None)
    if not callable(get_engine):
        raise VLLEPluginError("chemthermo_vlle must provide a callable get_engine() factory.")

    engine: Any = get_engine()
    return engine
