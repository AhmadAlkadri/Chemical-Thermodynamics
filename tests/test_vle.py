import pytest

import chemthermo as ct
from chemthermo.vle import bubble_pressure, bubble_temperature, dew_pressure, dew_temperature


@pytest.fixture
def binary_mixture():
    names = ("Propane", "n-Butane")
    return tuple(ct.Component.from_database(name) for name in names)


@pytest.fixture
def eos():
    return ct.PengRobinsonEOS()


def test_bubble_temperature_sanity(binary_mixture, eos):
    # Pure Propane (x1=1) at 1 atm
    # Propane Tb ~ 231 K
    mix = ct.Mixture(components=binary_mixture, composition=ct.Composition(fractions=(1.0, 0.0)))
    T_bub, y = bubble_temperature(mix, 101325.0, eos)

    assert 225.0 < T_bub < 235.0
    assert y[0] > 0.99  # Vapor should be pure propane


def test_dew_temperature_sanity(binary_mixture, eos):
    # Pure n-Butane (x2=1) at 1 atm
    # n-Butane Tb ~ 272 K
    mix = ct.Mixture(components=binary_mixture, composition=ct.Composition(fractions=(0.0, 1.0)))
    T_dew, x = dew_temperature(mix, 101325.0, eos)

    assert 268.0 < T_dew < 278.0
    assert x[1] > 0.99  # Liquid should be pure n-Butane


def test_bubble_pressure_sanity(binary_mixture, eos):
    # Equimolar at 300 K
    mix = ct.Mixture(components=binary_mixture, composition=ct.Composition(fractions=(0.5, 0.5)))
    P_bub, y = bubble_pressure(mix, 300.0, eos)

    # 300K is above Tb for both. So P > 1 atm.
    assert P_bub > 101325.0
    # Propane is more volatile, so vapor should be richer in propane (y[0] > 0.5)
    assert y[0] > 0.5


def test_dew_pressure_sanity(binary_mixture, eos):
    # Equimolar at 300 K
    mix = ct.Mixture(components=binary_mixture, composition=ct.Composition(fractions=(0.5, 0.5)))
    P_dew, x = dew_pressure(mix, 300.0, eos)

    assert P_dew > 101325.0
    assert P_dew < 100e5  # Reasonable upper bound

    # Check monotonicity with bubble pressure
    # For same T, Bubble P > Dew P usually?
    # Bubble P is pressure required to keep all liquid.
    # Dew P is pressure where first drop liquid forms.
    # At fixed T, compression: Vapor -> Dew P (Liquid forms) -> Bubble P (All liquid).
    # So Bubble P should be HIGHER than Dew P (for system of same z at same T). NO.
    # Bubble P (Start of boiling, from Liquid) vs Dew P (Start of condensing, from Vapor).
    # If we have z=0.5.
    # Bubble P is the P at the liquid boundary.
    # Dew P is the P at the vapor boundary.
    # At fixed T, liquid is stable at high P. Vapor at low P.
    # So transition Liquid -> Vapor happens as we lower P.
    # First bubble (Bubble Point) at P_bub.
    # Last liquid (Dew Point) at P_dew.
    # So P_bub > P_dew.

    # Wait, my manual validation showed:
    # x1=0.50: P_bub=3.48 bar
    # x1=0.50: P_dew=2.21 bar
    # So P_bub > P_dew. Correct.

    P_bub, _ = bubble_pressure(mix, 300.0, eos)
    assert P_bub > P_dew
