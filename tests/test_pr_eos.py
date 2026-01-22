import math

import pytest

import chemthermo as ct


def test_pr_eos_pure_methane_phi_regression() -> None:
    methane = ct.Component.from_database("Methane")
    mixture = ct.Mixture.from_components([methane], [1.0])
    eos = ct.PengRobinsonEOS()

    phi = eos.fugacity_coefficients(
        mixture=mixture,
        temperature_K=300.0,
        pressure_Pa=101325.0,
        composition=[1.0],
        phase="vapor",
    )

    assert phi == pytest.approx([0.9977890849], rel=1e-10)


def test_pr_eos_fugacity_coefficients_shape() -> None:
    mixture = ct.Mixture.from_database(["Methane", "Ethane"], [0.25, 0.75], normalize=True)
    eos = ct.PengRobinsonEOS()

    phi = eos.fugacity_coefficients(
        mixture=mixture,
        temperature_K=280.0,
        pressure_Pa=5.0e6,
        composition=[0.25, 0.75],
        phase="vapor",
    )

    assert len(phi) == 2
    assert all(math.isfinite(value) and value > 0.0 for value in phi)
