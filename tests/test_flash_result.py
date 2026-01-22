import pytest

import chemthermo as ct


def test_flash_result_basic_fields() -> None:
    composition = ct.Composition([1.0])
    phase = ct.PhaseResult(name="vapor", composition=composition)
    result = ct.FlashResult(
        temperature_K=300.0,
        pressure_Pa=101325.0,
        phases={"vapor": phase},
        vapor_fraction=1.0,
    )
    assert result.phase_names() == ["vapor"]


def test_flash_result_requires_phase_key_match() -> None:
    composition = ct.Composition([1.0])
    phase = ct.PhaseResult(name="liquid", composition=composition)
    with pytest.raises(ValueError):
        ct.FlashResult(temperature_K=300.0, pressure_Pa=101325.0, phases={"vapor": phase})


def test_flash_result_rejects_out_of_range_vapor_fraction() -> None:
    composition = ct.Composition([1.0])
    phase = ct.PhaseResult(name="vapor", composition=composition)
    with pytest.raises(ct.InputRangeError):
        ct.FlashResult(
            temperature_K=300.0,
            pressure_Pa=101325.0,
            phases={"vapor": phase},
            vapor_fraction=1.5,
        )
