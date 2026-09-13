"""The `chemthermo.vlle` plugin boundary is deprecated (ADR-0013).

`flash_tp(..., flash_mode="modified-raoult")` with `FlashSettings(max_phases=3)`
(the default) now discovers vapor-liquid-liquid equilibrium in-tree (ADR-0011),
so this out-of-tree plugin boundary is no longer the way to reach "VLLE
support" and is deprecated rather than removed. Importing `chemthermo.vlle`,
or calling `get_vlle_engine()`, emits a `DeprecationWarning` naming the
in-tree replacement; the exported names and the loader's
`VLLEPluginNotInstalledError` behavior are otherwise unchanged for one
deprecation cycle.
"""

from __future__ import annotations

import importlib

import pytest

import chemthermo.vlle as vlle
from chemthermo.vlle import VLLEInputs, VLLEPhase, VLLEResult
from chemthermo.vlle.errors import VLLEPluginNotInstalledError
from chemthermo.vlle.loader import get_vlle_engine


def test_importing_vlle_warns_deprecation() -> None:
    """The module-level warning fires on import.

    The module was already imported once at collection time (the module
    imports above), so it is re-imported here via ``importlib.reload`` to
    re-run its module body - including the ``warnings.warn`` call - inside
    this test's ``pytest.warns`` capture.
    """
    with pytest.warns(DeprecationWarning, match="chemthermo.vlle is deprecated"):
        importlib.reload(vlle)


def test_vlle_plugin_missing_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    def _raise(name: str) -> None:  # type: ignore[return-value]
        raise ModuleNotFoundError(name)

    monkeypatch.setattr("chemthermo.vlle.loader.import_module", _raise)

    with pytest.warns(DeprecationWarning, match="get_vlle_engine.*is deprecated"):
        with pytest.raises(
            VLLEPluginNotInstalledError,
            match="Install chemthermo_vlle to enable VLLE support.",
        ):
            get_vlle_engine()


def test_vlle_types_construct() -> None:
    inputs = VLLEInputs(temperature_K=300.0, pressure_Pa=101325.0, z=[0.4, 0.6])
    phases = [
        VLLEPhase(name="liquid1", fraction=0.6, composition=[0.45, 0.55]),
        VLLEPhase(name="liquid2", fraction=0.4, composition=[0.3, 0.7]),
    ]
    result = VLLEResult(phases=phases, diagnostics={"phase_count": 2})

    assert inputs.z[0] == pytest.approx(0.4)
    assert result.phases[0].name == "liquid1"
