"""PC-SAFT density roots and the ``(T, P, x)`` fugacity interface (ADR-0015).

Validation Case P-3. Nothing here needs an optional dependency; the numbers
that came from ``teqp`` are written out as literals with their provenance and
are re-derived from ``teqp`` in ``tests/validation/test_pcsaft_flash_vs_teqp.py``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.eos._pcsaft_density import (
    _ETA_GRID,
    ETA_MAX,
    ETA_MIN,
    ETA_UNIFORM_STEP,
    solve_density_roots,
)
from chemthermo.eos.pcsaft import PCSAFTEOS, R_J_PER_MOL_K

#: n-hexane saturation at 300 K, from ``teqp.pure_VLE_T`` with the packaged
#: Gross & Sadowski (2001) parameters. Recomputed against teqp in
#: ``tests/validation/test_pcsaft_flash_vs_teqp.py`` (Case P-2/P-3); repeated
#: here as literals so this file needs no optional dependency.
HEXANE_300K_PSAT_PA = 21858.084278856164
HEXANE_300K_RHO_LIQUID = 7518.498733715526
HEXANE_300K_RHO_VAPOR = 8.868596301925571


def _hexane() -> PCSAFTEOS:
    return PCSAFTEOS(components=("n-Hexane",))


def _isotherm(eos: PCSAFTEOS, names: tuple[str, ...], temperature: float, x: list[float]):
    return eos._isotherm(names=names, temperature_K=temperature, composition=x)


# ---------------------------------------------------------------------------
# The one-variable rewrite is the same model
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("components", "x", "temperature"),
    [
        (("n-Hexane",), [1.0], 300.0),
        (("n-Hexane",), [1.0], 450.0),
        (("Methane", "n-Hexane"), [0.3, 0.7], 300.0),
        (("Methane", "n-Hexane"), [0.9, 0.1], 250.0),
        (("Methane", "n-Decane"), [0.4, 0.6], 350.0),
    ],
)
def test_the_eta_isotherm_reproduces_the_rho_model(
    components: tuple[str, ...], x: list[float], temperature: float
) -> None:
    """``_pcsaft_density`` re-derives the model in ``eta``; it must not drift.

    The complexity receipt in that module's docstring names this test as what
    keeps the two derivations honest.
    """
    eos = PCSAFTEOS(components=components)
    isotherm = _isotherm(eos, components, temperature, x)

    etas = np.array([1e-6, 1e-4, 1e-2, 0.05, 0.15, 0.3, 0.45, 0.6])
    pressures = isotherm.pressure(etas)
    for eta, pressure in zip(etas, pressures):
        density = isotherm.density_per_eta * float(eta)
        reference = eos.pressure_Pa(
            temperature_K=temperature, density_mol_m3=density, composition=x
        )
        assert float(pressure) == pytest.approx(reference, rel=1e-12, abs=1e-6)

        derivatives = isotherm.derivatives(np.array([float(eta)]), second=True)
        assert float(derivatives.a[0]) == pytest.approx(
            eos.residual_helmholtz(
                temperature_K=temperature, volume_m3=1.0 / density, composition=x
            ),
            rel=1e-12,
            abs=1e-12,
        )
        assert 1.0 + float(eta) * float(derivatives.a1[0]) == pytest.approx(
            eos.compressibility_factor(
                temperature_K=temperature, density_mol_m3=density, composition=x
            ),
            rel=1e-12,
            abs=1e-12,
        )


@pytest.mark.parametrize(
    ("components", "x", "temperature", "density"),
    [
        (("n-Hexane",), [1.0], 300.0, 10.0),
        (("n-Hexane",), [1.0], 300.0, 7518.5),
        (("n-Hexane",), [1.0], 300.0, 3000.0),
        (("Methane", "n-Hexane"), [0.3, 0.7], 300.0, 500.0),
        (("Methane", "n-Hexane"), [0.3, 0.7], 300.0, 9000.0),
        (("Methane", "n-Decane"), [0.4, 0.6], 350.0, 5000.0),
    ],
)
def test_dp_drho_matches_a_central_difference_of_the_other_module(
    components: tuple[str, ...], x: list[float], temperature: float, density: float
) -> None:
    """The analytic ``dP/drho`` against a difference of ``PCSAFTEOS.pressure_Pa``.

    The finite difference is taken on the *other* implementation, so this tests
    the new second derivative rather than differencing it against itself.
    """
    eos = PCSAFTEOS(components=components)
    isotherm = _isotherm(eos, components, temperature, x)
    eta = density / isotherm.density_per_eta
    _pressure, analytic = isotherm.pressure_and_slope(eta)

    step = density * 1e-6
    forward = eos.pressure_Pa(
        temperature_K=temperature, density_mol_m3=density + step, composition=x
    )
    backward = eos.pressure_Pa(
        temperature_K=temperature, density_mol_m3=density - step, composition=x
    )
    numeric = (forward - backward) / (2.0 * step)
    assert analytic == pytest.approx(numeric, rel=1e-6)


def test_the_scan_grid_is_deterministic_and_spans_the_documented_range() -> None:
    assert _ETA_GRID[0] == pytest.approx(ETA_MIN, rel=1e-15)
    assert _ETA_GRID[-1] == pytest.approx(ETA_MAX, rel=1e-15)
    assert np.all(np.diff(_ETA_GRID) > 0.0)
    uniform = _ETA_GRID[_ETA_GRID >= 1e-3]
    assert np.max(np.diff(uniform)) <= ETA_UNIFORM_STEP * (1.0 + 1e-12)


# ---------------------------------------------------------------------------
# Case P-3: the roots themselves
# ---------------------------------------------------------------------------


def test_pure_hexane_at_saturation_has_exactly_two_admissible_roots() -> None:
    """Both saturation densities, and the spinodal root discarded, not returned."""
    eos = _hexane()
    isotherm = _isotherm(eos, ("n-Hexane",), 300.0, [1.0])
    result = solve_density_roots(isotherm, HEXANE_300K_PSAT_PA)

    assert len(result.densities) == 2
    vapor, liquid = result.densities
    assert vapor == pytest.approx(HEXANE_300K_RHO_VAPOR, rel=1e-8)
    assert liquid == pytest.approx(HEXANE_300K_RHO_LIQUID, rel=1e-8)

    # Three sign changes were bracketed; the middle one is the spinodal branch
    # and is not among the two returned.
    assert result.bracket_count == 3
    for density in result.densities:
        eta = density / isotherm.density_per_eta
        _pressure, slope = isotherm.pressure_and_slope(eta)
        assert slope > 0.0


def test_every_returned_root_is_mechanically_stable_and_solves_the_pressure() -> None:
    eos = _hexane()
    isotherm = _isotherm(eos, ("n-Hexane",), 300.0, [1.0])
    for pressure in (0.5 * HEXANE_300K_PSAT_PA, HEXANE_300K_PSAT_PA, 6.0e5, 1.0e7):
        result = solve_density_roots(isotherm, pressure)
        assert result.densities
        assert result.max_relative_residual < 1e-9
        for density in result.densities:
            eta = density / isotherm.density_per_eta
            value, slope = isotherm.pressure_and_slope(eta)
            assert slope > 0.0
            assert abs(value - pressure) / pressure < 1e-9
            assert value == pytest.approx(
                eos.pressure_Pa(temperature_K=300.0, density_mol_m3=density, composition=[1.0]),
                rel=1e-10,
                abs=1e-6,
            )


def test_below_saturation_the_metastable_liquid_root_still_exists() -> None:
    """A root set is mechanical, not thermodynamic.

    At half the saturation pressure n-hexane at 300 K still has *two*
    mechanically stable roots: the vapour, and a stretched (metastable) liquid,
    because the liquid spinodal of this model sits near -40 MPa, far below any
    positive pressure. A cubic behaves identically. Which of the two is the
    real phase is decided by Gibbs energy, not by the root finder - see
    :func:`test_min_gibbs_selection_picks_the_vapor_below_saturation`.
    """
    eos = _hexane()
    isotherm = _isotherm(eos, ("n-Hexane",), 300.0, [1.0])
    result = solve_density_roots(isotherm, 0.5 * HEXANE_300K_PSAT_PA)
    assert len(result.densities) == 2
    assert result.densities[0] < HEXANE_300K_RHO_VAPOR
    assert result.densities[1] < HEXANE_300K_RHO_LIQUID

    spinodal_pressure = float(np.min(isotherm.pressure(_ETA_GRID)))
    assert spinodal_pressure < -1e7


@pytest.mark.parametrize(
    ("components", "x", "temperature", "pressure", "label"),
    [
        (("n-Hexane",), [1.0], 300.0, 1.0e7, "compressed liquid"),
        (("Methane",), [1.0], 300.0, 5.0e6, "supercritical methane"),
        (("Methane", "n-Hexane"), [0.3, 0.7], 300.0, 8.0e6, "compressed mixture liquid"),
    ],
)
def test_single_root_states_return_one_density(
    components: tuple[str, ...],
    x: list[float],
    temperature: float,
    pressure: float,
    label: str,
) -> None:
    eos = PCSAFTEOS(components=components)
    roots = eos.density_roots(temperature_K=temperature, pressure_Pa=pressure, composition=x)
    assert len(roots) == 1, f"{label}: expected one root, got {roots}"


def test_no_root_raises_a_model_error_naming_the_state() -> None:
    """Above the pressure the close-packing limit can reach there is no root."""
    eos = _hexane()
    with pytest.raises(ct.ModelError, match="no density root"):
        eos.density_roots(temperature_K=300.0, pressure_Pa=1.0e30, composition=[1.0])


def test_density_roots_are_deterministic() -> None:
    eos = PCSAFTEOS(components=("Methane", "n-Hexane"))
    first = eos.density_roots(temperature_K=300.0, pressure_Pa=3.0e6, composition=[0.3, 0.7])
    for _ in range(3):
        assert (
            eos.density_roots(temperature_K=300.0, pressure_Pa=3.0e6, composition=[0.3, 0.7])
            == first
        )


# ---------------------------------------------------------------------------
# Case P-3: the fugacity interface on top of the roots
# ---------------------------------------------------------------------------


def test_fugacity_coefficients_are_the_root_values_exponentiated() -> None:
    mixture = ct.Mixture.from_database(["Methane", "n-Hexane"], [0.3, 0.7])
    eos = PCSAFTEOS()
    temperature, pressure = 300.0, 9.0e5
    x = [0.3, 0.7]
    roots = eos.density_roots(
        temperature_K=temperature, pressure_Pa=pressure, composition=x, mixture=mixture
    )
    assert len(roots) == 2

    for phase, density in (("vapor", roots[0]), ("liquid", roots[-1])):
        expected = np.exp(
            eos._ln_fugacity_coefficients(
                names=("Methane", "n-Hexane"),
                temperature_K=temperature,
                density_mol_m3=density,
                composition=x,
            )
        )
        actual = eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=temperature,
            pressure_Pa=pressure,
            composition=x,
            phase=phase,
        )
        np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0.0)
        assert eos.molar_volume(
            temperature_K=temperature,
            pressure_Pa=pressure,
            composition=x,
            phase=phase,
            mixture=mixture,
        ) == pytest.approx(1.0 / density, rel=1e-15)


def test_a_single_root_state_returns_the_same_values_for_both_labels() -> None:
    mixture = ct.Mixture.from_database(["Methane", "n-Hexane"], [0.3, 0.7])
    eos = PCSAFTEOS()
    assert (
        len(
            eos.density_roots(
                temperature_K=300.0, pressure_Pa=8.0e6, composition=[0.3, 0.7], mixture=mixture
            )
        )
        == 1
    )
    vapor = eos.fugacity_coefficients(
        mixture=mixture,
        temperature_K=300.0,
        pressure_Pa=8.0e6,
        composition=[0.3, 0.7],
        phase="vapor",
    )
    liquid = eos.fugacity_coefficients(
        mixture=mixture,
        temperature_K=300.0,
        pressure_Pa=8.0e6,
        composition=[0.3, 0.7],
        phase="liquid",
    )
    assert vapor == liquid


def test_min_gibbs_selection_picks_the_vapor_below_saturation() -> None:
    """``_ln_phi_min_gibbs`` is the rule that turns two roots into one phase.

    Below the saturation pressure the vapour root has the lower
    ``sum_i x_i ln phi_i``; above it the liquid does; at saturation the two are
    equal, which for a pure fluid *is* the saturation condition.
    """
    from chemthermo.stability.tp import _ln_phi_min_gibbs

    eos = _hexane()
    mixture = ct.Mixture.from_database(["n-Hexane"], [1.0])
    composition = np.array([1.0])

    def ln_phi(pressure: float, phase: str) -> float:
        values = eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=300.0,
            pressure_Pa=pressure,
            composition=[1.0],
            phase=phase,
        )
        return math.log(values[0])

    below = 0.5 * HEXANE_300K_PSAT_PA
    assert ln_phi(below, "vapor") < ln_phi(below, "liquid")
    _terms, branch = _ln_phi_min_gibbs(
        eos, mixture=mixture, temperature=300.0, pressure=below, composition=composition
    )
    assert branch == "vapor"

    above = 1.5 * HEXANE_300K_PSAT_PA
    assert ln_phi(above, "liquid") < ln_phi(above, "vapor")
    _terms, branch = _ln_phi_min_gibbs(
        eos, mixture=mixture, temperature=300.0, pressure=above, composition=composition
    )
    assert branch == "liquid"

    tie = abs(ln_phi(HEXANE_300K_PSAT_PA, "liquid") - ln_phi(HEXANE_300K_PSAT_PA, "vapor"))
    assert tie < 1e-9


def test_saturation_pressure_is_where_the_two_roots_tie() -> None:
    """Bisect on the fugacity difference and land on the teqp saturation state."""
    eos = _hexane()
    mixture = ct.Mixture.from_database(["n-Hexane"], [1.0])

    def gap(pressure: float) -> float:
        liquid = eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=300.0,
            pressure_Pa=pressure,
            composition=[1.0],
            phase="liquid",
        )[0]
        vapor = eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=300.0,
            pressure_Pa=pressure,
            composition=[1.0],
            phase="vapor",
        )[0]
        return math.log(liquid) - math.log(vapor)

    low, high = 1.0e3, 1.0e5
    for _ in range(200):
        mid = 0.5 * (low + high)
        if gap(low) * gap(mid) <= 0.0:
            high = mid
        else:
            low = mid
    assert 0.5 * (low + high) == pytest.approx(HEXANE_300K_PSAT_PA, rel=1e-10)


# ---------------------------------------------------------------------------
# Component resolution (ADR-0015)
# ---------------------------------------------------------------------------


def test_components_may_come_from_the_mixture_or_be_pinned() -> None:
    mixture = ct.Mixture.from_database(["Methane", "n-Hexane"], [0.3, 0.7])

    def phi(eos: PCSAFTEOS) -> list[float]:
        return list(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=300.0,
                pressure_Pa=3.0e6,
                composition=[0.3, 0.7],
                phase="liquid",
            )
        )

    pinned = phi(PCSAFTEOS(components=("Methane", "n-Hexane")))
    assert phi(PCSAFTEOS()) == pinned
    assert phi(PCSAFTEOS(components=("METHANE", "n-hexane"))) == pinned


def test_a_component_order_that_disagrees_with_the_mixture_is_refused() -> None:
    mixture = ct.Mixture.from_database(["Methane", "n-Hexane"], [0.3, 0.7])
    swapped = PCSAFTEOS(components=("n-Hexane", "Methane"))
    with pytest.raises(ct.ModelError, match="do not match the mixture"):
        swapped.fugacity_coefficients(
            mixture=mixture,
            temperature_K=300.0,
            pressure_Pa=3.0e6,
            composition=[0.3, 0.7],
            phase="liquid",
        )


def test_the_phase_label_is_validated() -> None:
    mixture = ct.Mixture.from_database(["Methane", "n-Hexane"], [0.3, 0.7])
    with pytest.raises(ValueError, match="'vapor' or 'liquid'"):
        PCSAFTEOS().fugacity_coefficients(
            mixture=mixture,
            temperature_K=300.0,
            pressure_Pa=3.0e6,
            composition=[0.3, 0.7],
            phase="gas",
        )


def test_pcsaft_is_an_equation_of_state_and_is_exported() -> None:
    assert isinstance(PCSAFTEOS(), ct.EquationOfState)
    assert ct.PCSAFTEOS is PCSAFTEOS
    assert "PCSAFTEOS" in ct.__all__


def test_the_pressure_of_a_root_uses_the_exact_si_gas_constant() -> None:
    eos = _hexane()
    isotherm = _isotherm(eos, ("n-Hexane",), 300.0, [1.0])
    result = solve_density_roots(isotherm, HEXANE_300K_PSAT_PA)
    density = result.densities[0]
    z_factor = eos.compressibility_factor(
        temperature_K=300.0, density_mol_m3=density, composition=[1.0]
    )
    assert z_factor * density * R_J_PER_MOL_K * 300.0 == pytest.approx(
        HEXANE_300K_PSAT_PA, rel=1e-9
    )
