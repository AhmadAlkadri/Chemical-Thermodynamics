"""Plugin loader for the VLLE engine."""

from __future__ import annotations

from importlib import import_module
from typing import Any

from .api import VLLEEngine
from .errors import VLLEPluginError, VLLEPluginNotInstalledError


def get_vlle_engine() -> VLLEEngine:
    """Return the VLLE engine from the optional chemthermo_vlle plugin.

    The plugin must expose a callable ``get_engine()`` that returns an object
    implementing the ``VLLEEngine`` protocol.
    """

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
