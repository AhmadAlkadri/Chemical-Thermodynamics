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
