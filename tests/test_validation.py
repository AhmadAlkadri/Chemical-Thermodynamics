import pytest

import chemthermo as ct


def test_validate_temperature_rejects_nonpositive() -> None:
    with pytest.raises(ct.InputRangeError):
        ct.validate_temperature(0.0)
    with pytest.raises(ct.InputRangeError):
        ct.validate_temperature(-10.0)


def test_validate_pressure_rejects_nonpositive() -> None:
    with pytest.raises(ct.InputRangeError):
        ct.validate_pressure(0.0)
    with pytest.raises(ct.InputRangeError):
        ct.validate_pressure(-2.0)


def test_validate_fractions_requires_sum_to_one_by_default() -> None:
    with pytest.raises(ct.CompositionError):
        ct.validate_fractions([0.25, 0.25])


def test_validate_fractions_normalize_explicit() -> None:
    values = ct.validate_fractions([1.0, 1.0, 2.0], normalize=True)
    assert values == pytest.approx([0.25, 0.25, 0.5])


def test_validate_fractions_rejects_negative() -> None:
    with pytest.raises(ct.CompositionError):
        ct.validate_fractions([0.8, -0.2], normalize=False)


def test_validate_fractions_empty_raises() -> None:
    with pytest.raises(ct.CompositionError):
        ct.validate_fractions([])


def test_validate_fractions_sum_zero_raises() -> None:
    with pytest.raises(ct.CompositionError):
        ct.validate_fractions([0.0, 0.0], normalize=True)


def test_units_pressure_constants() -> None:
    assert ct.PA_PER_BAR == 100000.0
    assert ct.STANDARD_P_PA == 101325.0
    assert ct.PA_PER_BAR != ct.STANDARD_P_PA
