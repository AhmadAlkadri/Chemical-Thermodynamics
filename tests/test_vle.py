import numpy as np
import pytest

import chemthermo as ct
from chemthermo.vle import bubble_pressure, bubble_temperature, dew_pressure, dew_temperature


@pytest.fixture
def binary_components() -> tuple[ct.Component, ct.Component]:
    return (
        ct.Component.from_database("Methane"),
        ct.Component.from_database("Ethane"),
    )


@pytest.fixture
def eos() -> ct.PengRobinsonEOS:
    return ct.PengRobinsonEOS()


def _k_values(
    mixture: ct.Mixture,
    eos: ct.PengRobinsonEOS,
    temperature_K: float,
    pressure_Pa: float,
    liquid_composition: tuple[float, ...],
    vapor_composition: tuple[float, ...],
) -> np.ndarray:
    phi_v = np.array(
        eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=vapor_composition,
            phase="vapor",
        ),
        dtype=float,
    )
    phi_l = np.array(
        eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=liquid_composition,
            phase="liquid",
        ),
        dtype=float,
    )
    return phi_l / phi_v


def test_bubble_temperature_residual_closure_pr_binary(binary_components, eos):
    tol = 1e-7
    x = (0.4, 0.6)
    mixture = ct.Mixture(components=binary_components, composition=ct.Composition(fractions=x))

    temperature_K, y = bubble_temperature(
        mixture,
        pressure_Pa=2.0e6,
        eos=eos,
        settings={"tol": tol, "max_iter": 120},
    )

    K = _k_values(
        mixture,
        eos,
        temperature_K,
        2.0e6,
        liquid_composition=x,
        vapor_composition=tuple(y),
    )
    residual = float(np.dot(K, np.array(x, dtype=float)) - 1.0)
    assert abs(residual) < tol


def test_dew_temperature_residual_closure_pr_binary(binary_components, eos):
    tol = 1e-7
    y = (0.7, 0.3)
    mixture = ct.Mixture(components=binary_components, composition=ct.Composition(fractions=y))

    temperature_K, x = dew_temperature(
        mixture,
        pressure_Pa=2.0e6,
        eos=eos,
        settings={"tol": tol, "max_iter": 120},
    )

    K = _k_values(
        mixture,
        eos,
        temperature_K,
        2.0e6,
        liquid_composition=tuple(x),
        vapor_composition=y,
    )
    residual = float(np.sum(np.array(y, dtype=float) / K) - 1.0)
    assert abs(residual) < tol


def test_bubble_pressure_residual_closure_pr_binary(binary_components, eos):
    tol = 1e-7
    x = (0.5, 0.5)
    mixture = ct.Mixture(components=binary_components, composition=ct.Composition(fractions=x))

    pressure_Pa, y = bubble_pressure(
        mixture,
        temperature_K=220.0,
        eos=eos,
        settings={"tol": tol, "max_iter": 120},
    )

    K = _k_values(
        mixture,
        eos,
        220.0,
        pressure_Pa,
        liquid_composition=x,
        vapor_composition=tuple(y),
    )
    residual = float(np.dot(K, np.array(x, dtype=float)) - 1.0)
    assert abs(residual) < tol


def test_dew_pressure_residual_closure_pr_binary(binary_components, eos):
    tol = 1e-7
    y = (0.6, 0.4)
    mixture = ct.Mixture(components=binary_components, composition=ct.Composition(fractions=y))

    pressure_Pa, x = dew_pressure(
        mixture,
        temperature_K=220.0,
        eos=eos,
        settings={"tol": tol, "max_iter": 120},
    )

    K = _k_values(
        mixture,
        eos,
        220.0,
        pressure_Pa,
        liquid_composition=tuple(x),
        vapor_composition=y,
    )
    residual = float(np.sum(np.array(y, dtype=float) / K) - 1.0)
    assert abs(residual) < tol


def test_vle_raises_convergence_error_when_tolerance_not_met(binary_components, eos):
    mixture = ct.Mixture(
        components=binary_components,
        composition=ct.Composition(fractions=(0.5, 0.5)),
    )

    with pytest.raises(ct.ConvergenceError):
        bubble_pressure(
            mixture,
            temperature_K=220.0,
            eos=eos,
            settings={
                "tol": 1e-12,
                "max_iter": 1,
                "inner_max_iter": 2,
                "inner_tol": 1e-14,
            },
        )
