import pytest

import chemthermo as ct


class DummyEOS(ct.EquationOfState):
    def fugacity_coefficients(
        self,
        *,
        mixture: ct.Mixture,
        temperature_K: float,
        pressure_Pa: float,
        composition: list[float],
        phase: str,
    ) -> list[float]:
        return [1.0 for _ in composition]


def test_flash_tp_requires_eos() -> None:
    mix = ct.Mixture.from_database(["Methane"], [1.0])
    with pytest.raises(ct.ModelError):
        ct.flash_tp(mix, temperature_K=300.0, pressure_Pa=101325.0, eos=None)


def test_flash_tp_validates_inputs() -> None:
    mix = ct.Mixture.from_database(["Methane"], [1.0])
    eos = DummyEOS()
    with pytest.raises(ct.InputRangeError):
        ct.flash_tp(mix, temperature_K=0.0, pressure_Pa=101325.0, eos=eos)
    with pytest.raises(ct.InputRangeError):
        ct.flash_tp(mix, temperature_K=300.0, pressure_Pa=0.0, eos=eos)


def test_flash_tp_unimplemented_model_error() -> None:
    mix = ct.Mixture.from_database(["Methane"], [1.0])
    eos = DummyEOS()
    result = ct.flash_tp(mix, temperature_K=300.0, pressure_Pa=101325.0, eos=eos)
    assert result.phase_names() == ["vapor"]


def test_flash_tp_vlle_mode_points_to_modified_raoult(recwarn: pytest.WarningsRecorder) -> None:
    """`flash_mode="vlle"` is refused; the message names the in-tree replacement.

    Three-phase equilibrium is discovered automatically by `"modified-raoult"`
    with `FlashSettings(max_phases=3)` (ADR-0011); the out-of-tree
    `chemthermo_vlle` plugin is no longer the way to get it (ADR-0013), and
    the `chemthermo.vlle` package itself was removed in 0.4.0 (ADR-0037).
    The call emits no `DeprecationWarning`.
    """
    mix = ct.Mixture.from_database(["Methane"], [1.0])
    with pytest.raises(
        ct.ModelError,
        match="flash_mode='vlle' is not a supported mode",
    ) as excinfo:
        ct.flash_tp(mix, temperature_K=300.0, pressure_Pa=101325.0, flash_mode="vlle")

    message = str(excinfo.value)
    assert "modified-raoult" in message
    assert "FlashSettings(max_phases=3)" in message
    assert "ADR-0011" in message
    assert "chemthermo.vlle plugin package was removed in 0.4.0" in message
    assert "ADR-0037" in message
    assert "chemthermo_vlle plugin" not in message
    assert "Install chemthermo_vlle" not in message
    assert not [w for w in recwarn.list if issubclass(w.category, DeprecationWarning)]


def test_the_removed_vlle_package_is_gone() -> None:
    """ADR-0037: `chemthermo.vlle` (deprecated by ADR-0013) was removed in 0.4.0."""
    import importlib

    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("chemthermo.vlle")
