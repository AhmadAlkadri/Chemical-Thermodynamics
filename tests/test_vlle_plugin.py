from __future__ import annotations

import pytest

from chemthermo.vlle import VLLEInputs, VLLEPhase, VLLEResult
from chemthermo.vlle.errors import VLLEPluginNotInstalledError
from chemthermo.vlle.loader import get_vlle_engine


def test_vlle_plugin_missing_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    def _raise(name: str) -> None:  # type: ignore[return-value]
        raise ModuleNotFoundError(name)

    monkeypatch.setattr("chemthermo.vlle.loader.import_module", _raise)

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
