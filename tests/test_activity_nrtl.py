import pytest

import chemthermo as ct


def test_nrtl_gamma_unity_with_zero_tau() -> None:
    components = tuple(ct.Component.from_database(name) for name in ("Methane", "Ethane"))
    mixture = ct.Mixture(
        components=components,
        composition=ct.Composition(fractions=(0.4, 0.6), basis="mole", normalize=False),
    )

    params = ct.NRTLParameters.from_pairs([("Methane", "Ethane", 0.0, 0.0, 0.3, 0.3)])
    model = ct.NRTL(parameters=params)

    gamma = model.activity_coefficients(
        mixture=mixture, temperature_K=300.0, composition=mixture.fractions
    )

    assert gamma == pytest.approx((1.0, 1.0))


def test_nrtl_symmetric_pair_yields_equal_gamma() -> None:
    components = tuple(ct.Component.from_database(name) for name in ("Methane", "Ethane"))
    mixture = ct.Mixture(
        components=components,
        composition=ct.Composition(fractions=(0.5, 0.5), basis="mole", normalize=False),
    )

    params = ct.NRTLParameters.from_pairs([("Methane", "Ethane", 0.2, 0.2, 0.3, 0.3)])
    model = ct.NRTL(parameters=params)

    gamma = model.activity_coefficients(
        mixture=mixture, temperature_K=300.0, composition=mixture.fractions
    )

    assert gamma[0] == pytest.approx(gamma[1])
    assert gamma[0] > 0.0


def test_nrtl_parameter_loader_matches_database_pair() -> None:
    params = ct.NRTLParameters.load()
    tau, alpha = params.for_components(["Methane", "Ethane"])

    assert tau[0, 1] == pytest.approx(0.2)
    assert tau[1, 0] == pytest.approx(0.1)
    assert alpha[0, 1] == pytest.approx(0.3)
    assert alpha[1, 0] == pytest.approx(0.3)
