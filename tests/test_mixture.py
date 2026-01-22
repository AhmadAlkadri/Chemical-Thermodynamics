import pytest

import chemthermo as ct


def test_mixture_requires_matching_lengths() -> None:
    composition = ct.Composition([1.0])
    with pytest.raises(ct.CompositionError):
        ct.Mixture(components=(), composition=composition)


def test_mixture_from_database_normalize_explicit() -> None:
    mixture = ct.Mixture.from_database(["Methane", "Ethane"], [1.0, 3.0], normalize=True)
    assert mixture.fractions == pytest.approx((0.25, 0.75))
    assert mixture.basis == "mole"


def test_mixture_from_database_requires_sum_to_one_by_default() -> None:
    with pytest.raises(ct.CompositionError):
        ct.Mixture.from_database(["Methane", "Ethane"], [0.2, 0.2])
