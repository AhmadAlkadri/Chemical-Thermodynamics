import pytest

import chemthermo as ct


def test_composition_requires_sum_to_one_by_default() -> None:
    with pytest.raises(ct.CompositionError):
        ct.Composition([0.2, 0.2])


def test_composition_normalize_explicit() -> None:
    comp = ct.Composition([1.0, 3.0], normalize=True)
    assert comp.fractions == pytest.approx((0.25, 0.75))


def test_composition_invalid_basis_raises() -> None:
    with pytest.raises(ct.CompositionError):
        ct.Composition([1.0], basis="volume")  # type: ignore[arg-type]
