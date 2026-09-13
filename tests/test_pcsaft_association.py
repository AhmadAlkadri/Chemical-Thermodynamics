"""PC-SAFT association: parameters, site fractions, derivatives, API (ADR-0018).

Everything here runs without an optional dependency. The external cross-check
against FeOs - the one that decides whether the *equations* are right rather
than only self-consistent - lives in
``tests/validation/test_pcsaft_association_vs_feos.py`` (validation Cases P-6
and P-7).

What this file is responsible for:

* the packaged Gross & Sadowski (2002) records and the rules that guard them;
* the site-fraction solve, against the closed form that exists for a pure 2B
  fluid and against the mass-action equations themselves;
* **derivative discipline** - every analytic derivative of the association
  term against a finite difference of the quantity it differentiates, the
  Euler identity, and the Michelsen-Hendriks stationarity that the whole
  derivative route rests on;
* the ADR-0014 guard: with no associating component, nothing changes.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.eos import _pcsaft_association as assoc
from chemthermo.eos.pcsaft import _evaluate
from chemthermo.parameters.pcsaft import (
    PCSAFTAssociationRecord,
    PCSAFTParameterError,
    PCSAFTParameters,
)

# --------------------------------------------------------------------------
# States. Every one has Z > 0 (they are density roots of a real pressure), so
# ln phi exists and the finite-difference checks have something to compare.
# --------------------------------------------------------------------------

#: ``(components, x, T/K, rho/(mol/m^3))``.
STATES: list[tuple[tuple[str, ...], list[float], float, float]] = [
    (("Water",), [1.0], 300.0, 51180.849),
    (("Water",), [1.0], 373.15, 33.2726),
    (("Water",), [1.0], 373.15, 48755.5163),
    (("Ethanol",), [1.0], 450.0, 691.1054),
    (("Ethanol",), [1.0], 350.0, 15928.6188),
    (("Water", "Ethanol"), [0.5, 0.5], 320.0, 25401.3411),
    (("Water", "Ethanol"), [0.2, 0.8], 351.0, 18574.9573),
    (("Water", "n-Hexane"), [0.3, 0.7], 298.15, 10229.0998),
    (("Water", "n-Hexane"), [0.9, 0.1], 298.15, 33472.7843),
    (("Methanol", "Water", "n-Hexane"), [0.3, 0.4, 0.3], 320.0, 16684.835),
]
STATE_IDS = [f"{'/'.join(n)} {t:g}K {r:g}" for n, _, t, r in STATES]


def _state_id(index: int) -> str:
    return STATE_IDS[index]


# --------------------------------------------------------------------------
# Packaged parameters
# --------------------------------------------------------------------------

#: Gross & Sadowski (2002) Table 1, as it must appear in the packaged JSON:
#: ``(m, sigma/A, eps/k in K, kappa^AB, eps^AB/k in K)``. Written out here so
#: the assertion compares the package against this file, not against itself.
#: Verified against two independent secondary sources that both cite
#: DOI 10.1021/ie010954d - FeOs ``parameters/pcsaft/gross2002.json`` and
#: Clapeyron.jl ``PCSAFT_like.csv`` / ``PCSAFT_assoc.csv``. The paper itself
#: is paywalled and was not read; see validation Case P-6.
PUBLISHED_2002: dict[str, tuple[float, float, float, float, float]] = {
    "Water": (1.0656, 3.0007, 366.51, 0.034868, 2500.7),
    "Methanol": (1.5255, 3.2300, 188.90, 0.035176, 2899.5),
    "Ethanol": (2.3827, 3.1771, 198.24, 0.032384, 2653.4),
    "1-Propanol": (2.9997, 3.2522, 233.40, 0.015268, 2276.8),
    "n-Butanol": (2.7515, 3.6139, 259.59, 0.006692, 2544.6),
}


@pytest.mark.parametrize(("name", "expected"), sorted(PUBLISHED_2002.items()))
def test_packaged_association_parameter_values(
    name: str, expected: tuple[float, float, float, float, float]
) -> None:
    parameters = PCSAFTParameters.load()
    m, sigma_A, epsilon_k_K = parameters.for_components([name])
    (record,) = parameters.association_for_components([name])
    assert record is not None
    actual = (
        float(m[0]),
        float(sigma_A[0]),
        float(epsilon_k_K[0]),
        record.kappa_ab,
        record.epsilon_ab_k_K,
    )
    assert actual == expected
    assert (record.na, record.nb) == (1.0, 1.0)
    assert record.scheme == "2B"


def test_every_associating_record_names_a_databank_component() -> None:
    for name in PUBLISHED_2002:
        assert ct.Component.from_database(name).name


def test_the_non_associating_records_stay_non_associating() -> None:
    """The ADR-0014 set must not have grown an association block by accident."""
    parameters = PCSAFTParameters.load()
    non_associating = [
        name for name in parameters.names() if name not in {n.casefold() for n in PUBLISHED_2002}
    ]
    records = parameters.association_for_components(non_associating)
    assert all(record is None for record in records)
    assert len(non_associating) == 11


def test_user_supplied_association_parameters_work() -> None:
    """A component the package does not ship can be given sites by the caller."""
    parameters = PCSAFTParameters.from_records(
        [
            {
                "name": "Acetic acid",
                "m": 1.3403,
                "sigma_A": 3.8582,
                "epsilon_k_K": 211.59,
                "source": "made up for this test - not a published parameter set",
                "association": {
                    "scheme": "2B",
                    "kappa_ab": 0.075550,
                    "epsilon_ab_k_K": 3044.4,
                },
            }
        ]
    )
    eos = PCSAFTEOS(components=("Acetic acid",), parameters=parameters)
    assert eos.associates()
    value = eos.residual_helmholtz(temperature_K=350.0, volume_m3=1.0 / 15000.0, composition=[1.0])
    assert math.isfinite(value) and value < 0.0
    # Dropping the sites must change the answer, or the block was ignored.
    bare = PCSAFTParameters.from_records(
        [{"name": "Acetic acid", "m": 1.3403, "sigma_A": 3.8582, "epsilon_k_K": 211.59}]
    )
    without = PCSAFTEOS(components=("Acetic acid",), parameters=bare)
    assert not without.associates()
    assert (
        abs(
            value
            - without.residual_helmholtz(
                temperature_K=350.0, volume_m3=1.0 / 15000.0, composition=[1.0]
            )
        )
        > 1.0
    )


def test_an_association_record_object_round_trips_through_from_records() -> None:
    record = PCSAFTAssociationRecord(kappa_ab=0.01, epsilon_ab_k_K=2000.0, source="test")
    parameters = PCSAFTParameters.from_records(
        [
            {
                "name": "Water",
                "m": 1.0656,
                "sigma_A": 3.0007,
                "epsilon_k_K": 366.51,
                "association": record,
            }
        ]
    )
    (stored,) = parameters.association_for_components(["Water"])
    assert stored == record


@pytest.mark.parametrize(
    "payload",
    [
        {"kappa_ab": 0.0, "epsilon_ab_k_K": 2500.0},
        {"kappa_ab": -1e-3, "epsilon_ab_k_K": 2500.0},
        {"kappa_ab": 0.03, "epsilon_ab_k_K": 0.0},
        {"kappa_ab": 0.03, "epsilon_ab_k_K": 2500.0, "na": 0, "nb": 1},
        {"kappa_ab": 0.03, "epsilon_ab_k_K": 2500.0, "na": 1, "nb": 0},
        {"kappa_ab": 0.03, "epsilon_ab_k_K": 2500.0, "na": 2, "nb": 1, "scheme": "2B"},
        {"epsilon_ab_k_K": 2500.0},
        {"kappa_ab": "wet", "epsilon_ab_k_K": 2500.0},
    ],
)
def test_invalid_association_blocks_are_rejected(payload: dict[str, object]) -> None:
    with pytest.raises(PCSAFTParameterError):
        PCSAFTParameters.from_records(
            [
                {
                    "name": "Water",
                    "m": 1.0656,
                    "sigma_A": 3.0007,
                    "epsilon_k_K": 366.51,
                    "association": payload,
                }
            ]
        )


def test_a_scheme_label_that_agrees_with_the_counts_is_accepted() -> None:
    record = PCSAFTAssociationRecord(
        kappa_ab=0.03, epsilon_ab_k_K=2500.0, na=2.0, nb=2.0, scheme="4C"
    )
    assert (record.na, record.nb) == (2.0, 2.0)


# --------------------------------------------------------------------------
# Bit-identity of the non-associating path (the ADR-0014 guard)
# --------------------------------------------------------------------------

#: Values pinned in validation Cases P-1 and P-3 for the non-associating
#: model. They must be reproduced **exactly**, not to a tolerance: the
#: association term is supposed not to run at all here.
PINNED_NON_ASSOCIATING = {
    "hexane_300_7700": (-5.783742760059239, 0.661534529144653, -5.709015132378622),
    "hexane_300_100": (-0.13257462366272255, 0.8693321635961244, -0.12321246990538429),
}
PINNED_HEXANE_ROOTS = (8.868596301913758, 7518.498733715524)


@pytest.mark.parametrize(
    ("key", "density"), [("hexane_300_7700", 7700.0), ("hexane_300_100", 100.0)]
)
def test_non_associating_values_are_bit_identical(key: str, density: float) -> None:
    eos = PCSAFTEOS(components=("n-Hexane",))
    a_res, z_factor, ln_phi = PINNED_NON_ASSOCIATING[key]
    assert (
        eos.residual_helmholtz(temperature_K=300.0, volume_m3=1.0 / density, composition=[1.0])
        == a_res
    )
    assert (
        eos.compressibility_factor(temperature_K=300.0, density_mol_m3=density, composition=[1.0])
        == z_factor
    )
    assert (
        eos.ln_fugacity_coefficients(
            temperature_K=300.0, density_mol_m3=density, composition=[1.0]
        )[0]
        == ln_phi
    )


def test_non_associating_density_roots_are_bit_identical() -> None:
    eos = PCSAFTEOS(components=("n-Hexane",))
    roots = eos.density_roots(
        temperature_K=300.0, pressure_Pa=21858.084278856164, composition=[1.0]
    )
    assert roots == PINNED_HEXANE_ROOTS


def test_an_all_none_association_sequence_is_the_same_object_as_no_association() -> None:
    """The guard is on *sites present*, not on the argument being ``None``."""
    kwargs = dict(
        temperature_K=300.0,
        density_mol_m3=7700.0,
        x=np.array([1.0]),
        m=np.array([3.0576]),
        sigma_A=np.array([3.7983]),
        epsilon_k_K=np.array([236.77]),
        kij=np.zeros((1, 1)),
    )
    without = _evaluate(**kwargs)  # type: ignore[arg-type]
    with_empty = _evaluate(association=(None,), **kwargs)  # type: ignore[arg-type]
    assert without.a_res == with_empty.a_res
    assert without.z_minus_one == with_empty.z_minus_one
    assert with_empty.a_assoc == 0.0
    assert with_empty.site_fractions is None


def test_a_non_associating_instance_reports_so() -> None:
    eos = PCSAFTEOS(components=("Methane", "n-Hexane"))
    assert not eos.associates()
    assert eos.association_parameters() == (None, None)


def test_associates_needs_a_component_order() -> None:
    with pytest.raises(ct.ModelError, match="at least one component"):
        PCSAFTEOS().associates()


# --------------------------------------------------------------------------
# Site fractions
# --------------------------------------------------------------------------


def _setup_for(names: tuple[str, ...], x: list[float], temperature: float):
    """Return ``(setup, weights, delta, rho_a3, zeta_2, eta)`` at one state.

    Rebuilt from the public parameter accessors rather than reaching into the
    model, so the test's arithmetic is its own.
    """
    eos = PCSAFTEOS(components=names)
    m, sigma_A, epsilon_k_K = eos.component_parameters()
    d = sigma_A * (1.0 - 0.12 * np.exp(-3.0 * epsilon_k_K / temperature))
    setup = assoc.build_setup(
        temperature_K=temperature,
        sigma_A=sigma_A,
        d=d,
        association=eos.association_parameters(),
    )
    assert setup is not None
    return setup, np.asarray(x, dtype=float), m, d


def _delta_at(setup, x, m, d, rho_a3: float) -> tuple[np.ndarray, float, float]:
    moments = (math.pi / 6.0) * ((x * m)[None, :] * d[None, :] ** np.arange(4)[:, None]).sum(axis=1)
    zeta = moments * rho_a3
    zeta_2, eta = float(zeta[2]), float(zeta[3])
    g, _, _ = assoc.contact_values(pair_c=setup.pair_c, zeta_2=zeta_2, eta=eta)
    return setup.pair_constant * g, zeta_2, eta


def test_pure_2b_site_fractions_match_the_closed_form() -> None:
    """A pure 2B fluid has one equation with a quadratic solution."""
    temperature, density = 300.0, 40000.0
    setup, x, m, d = _setup_for(("Water",), [1.0], temperature)
    rho_a3 = density * 6.02214076e23 * 1e-30
    delta, _, _ = _delta_at(setup, x, m, d, rho_a3)
    weights = setup.weights(x)
    solved = assoc.solve_site_fractions(rho=rho_a3, weights=weights, delta=delta)

    strength = rho_a3 * float(delta[0, 1])
    closed_form = (-1.0 + math.sqrt(1.0 + 4.0 * strength)) / (2.0 * strength)
    assert solved.shape == (2,)
    assert solved[0] == pytest.approx(solved[1], abs=0.0, rel=1e-15)
    assert float(solved[0]) == pytest.approx(closed_form, rel=1e-14)
    assert 0.0 < float(solved[0]) < 1.0


@pytest.mark.parametrize("index", range(len(STATES)), ids=_state_id)
def test_site_fractions_satisfy_the_mass_action_equations(index: int) -> None:
    names, x_list, temperature, density = STATES[index]
    setup, x, m, d = _setup_for(names, x_list, temperature)
    rho_a3 = density * 6.02214076e23 * 1e-30
    delta, _, _ = _delta_at(setup, x, m, d, rho_a3)
    weights = setup.weights(x)
    x_sites = assoc.solve_site_fractions(rho=rho_a3, weights=weights, delta=delta)
    residual = assoc.mass_action_residual(rho=rho_a3, weights=weights, delta=delta, x_sites=x_sites)
    assert float(np.max(np.abs(residual))) < 1e-13
    assert np.all(x_sites > 0.0) and np.all(x_sites <= 1.0)


def test_site_fractions_go_to_one_in_the_ideal_gas_limit() -> None:
    setup, x, m, d = _setup_for(("Water",), [1.0], 500.0)
    rho_a3 = 1e-6 * 6.02214076e23 * 1e-30
    delta, _, _ = _delta_at(setup, x, m, d, rho_a3)
    x_sites = assoc.solve_site_fractions(rho=rho_a3, weights=setup.weights(x), delta=delta)
    assert float(np.min(x_sites)) > 1.0 - 1e-10


def test_the_model_of_the_eos_reproduces_those_site_fractions() -> None:
    """The site fractions the EOS carries are the ones this file solves for."""
    names, x_list, temperature, density = STATES[5]
    eos = PCSAFTEOS(components=names)
    carried = eos._state(temperature, density, x_list).site_fractions
    setup, x, m, d = _setup_for(names, x_list, temperature)
    rho_a3 = density * 6.02214076e23 * 1e-30
    delta, _, _ = _delta_at(setup, x, m, d, rho_a3)
    expected = assoc.solve_site_fractions(rho=rho_a3, weights=setup.weights(x), delta=delta)
    assert carried is not None
    np.testing.assert_array_equal(carried, expected)


def test_cross_association_moves_the_site_fractions() -> None:
    """Water in ethanol is not water on its own: the cross term has to act."""
    names, x_list, temperature, density = ("Water", "Ethanol"), [0.5, 0.5], 320.0, 25401.3411
    mixture_sites = PCSAFTEOS(components=names)._state(temperature, density, x_list).site_fractions
    pure_sites = PCSAFTEOS(components=("Water",))._state(temperature, density, [1.0]).site_fractions
    assert mixture_sites is not None and pure_sites is not None
    assert abs(float(mixture_sites[0]) - float(pure_sites[0])) > 1e-3


def test_the_wolbach_sandler_cross_rules_are_the_ones_implemented() -> None:
    """Cross ``kappa`` and ``epsilon`` against the rules written out here.

    Clapeyron.jl stores the water/ethanol cross pair explicitly as
    ``epsilon = 2577.05`` K and ``bondvol = 0.03356196748232913`` with the
    source DOI 10.1021/ie010954d, which is what these rules give - an
    independent confirmation of the rules, not only of the numbers.
    """
    temperature = 300.0
    setup, _, _, _ = _setup_for(("Water", "Ethanol"), [0.5, 0.5], temperature)
    kappa_w, kappa_e = 0.034868, 0.032384
    sigma_w, sigma_e = 3.0007, 3.1771
    eps_w, eps_e = 2500.7, 2653.4
    sigma_ij = 0.5 * (sigma_w + sigma_e)
    kappa_ij = math.sqrt(kappa_w * kappa_e) * (math.sqrt(sigma_w * sigma_e) / sigma_ij) ** 3
    eps_ij = 0.5 * (eps_w + eps_e)
    expected = sigma_ij**3 * kappa_ij * math.expm1(eps_ij / temperature)

    assert kappa_ij == pytest.approx(0.03356196748232913, rel=1e-15)
    assert eps_ij == pytest.approx(2577.05, rel=1e-15)
    # Site order is A(water), B(water), A(ethanol), B(ethanol): the A(water) /
    # B(ethanol) entry is the cross pair, and the like-type pairs are zero.
    assert float(setup.pair_constant[0, 3]) == pytest.approx(expected, rel=1e-14)
    assert float(setup.pair_constant[0, 2]) == 0.0
    assert float(setup.pair_constant[1, 3]) == 0.0


# --------------------------------------------------------------------------
# Derivative discipline
# --------------------------------------------------------------------------


@pytest.mark.parametrize("index", range(len(STATES)), ids=_state_id)
def test_compressibility_matches_a_finite_difference_in_density(index: int) -> None:
    names, x, temperature, density = STATES[index]
    eos = PCSAFTEOS(components=names)

    def a_res(rho: float) -> float:
        return eos.residual_helmholtz(temperature_K=temperature, volume_m3=1.0 / rho, composition=x)

    step = 1e-6 * density
    expected = 1.0 + density * (a_res(density + step) - a_res(density - step)) / (2.0 * step)
    actual = eos.compressibility_factor(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    assert actual == pytest.approx(expected, rel=1e-8, abs=1e-8)


@pytest.mark.parametrize("index", range(len(STATES)), ids=_state_id)
def test_ln_phi_matches_a_finite_difference_in_mole_numbers(index: int) -> None:
    """``ln phi_i = d(n a_res)/dn_i|_{T,V} - ln Z``, differenced directly."""
    names, x, temperature, density = STATES[index]
    eos = PCSAFTEOS(components=names)
    moles = np.asarray(x, dtype=float)
    volume = float(moles.sum()) / density

    def n_a_res(counts: np.ndarray) -> float:
        total = float(counts.sum())
        return total * eos.residual_helmholtz(
            temperature_K=temperature,
            volume_m3=volume / total,
            composition=(counts / total).tolist(),
        )

    z_factor = eos.compressibility_factor(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    expected = []
    for component in range(moles.size):
        step = 1e-6 * moles[component]
        up = moles.copy()
        up[component] += step
        down = moles.copy()
        down[component] -= step
        expected.append((n_a_res(up) - n_a_res(down)) / (2.0 * step) - math.log(z_factor))

    actual = eos.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-7, atol=1e-7)


@pytest.mark.parametrize("index", range(len(STATES)), ids=_state_id)
def test_the_euler_identity_holds(index: int) -> None:
    """``sum_i x_i ln phi_i = A^res/RT + Z - 1 - ln Z``, to round-off."""
    names, x, temperature, density = STATES[index]
    eos = PCSAFTEOS(components=names)
    ln_phi = np.asarray(
        eos.ln_fugacity_coefficients(
            temperature_K=temperature, density_mol_m3=density, composition=x
        )
    )
    a_res = eos.residual_helmholtz(
        temperature_K=temperature, volume_m3=1.0 / density, composition=x
    )
    z_factor = eos.compressibility_factor(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    left = float(np.asarray(x) @ ln_phi)
    right = a_res + z_factor - 1.0 - math.log(z_factor)
    assert left == pytest.approx(right, abs=1e-13, rel=0.0)


@pytest.mark.parametrize("index", range(len(STATES)), ids=_state_id)
def test_the_eta_path_agrees_with_the_density_path(index: int) -> None:
    names, x, temperature, density = STATES[index]
    eos = PCSAFTEOS(components=names)
    isotherm = eos._isotherm(names=names, temperature_K=temperature, composition=x)
    eta = density / isotherm.density_per_eta
    derivatives = isotherm.derivatives(np.array([eta]), second=True)
    a_res = eos.residual_helmholtz(
        temperature_K=temperature, volume_m3=1.0 / density, composition=x
    )
    assert float(derivatives.a[0]) == pytest.approx(a_res, rel=1e-12, abs=1e-12)
    z_factor = eos.compressibility_factor(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    assert 1.0 + eta * float(derivatives.a1[0]) == pytest.approx(z_factor, rel=1e-12, abs=1e-12)


@pytest.mark.parametrize("index", range(len(STATES)), ids=_state_id)
def test_the_eta_derivatives_match_finite_differences(index: int) -> None:
    """``a'`` and ``a''`` of the packing-fraction path, including association.

    ``a''`` is the one derivative Michelsen-Hendriks stationarity does *not*
    supply for free (it needs the site-fraction sensitivity, Eq. 7 of
    ``_pcsaft_association``), and it is what decides ``dP/drho`` and therefore
    which density roots survive the mechanical-stability filter.
    """
    names, x, temperature, density = STATES[index]
    eos = PCSAFTEOS(components=names)
    isotherm = eos._isotherm(names=names, temperature_K=temperature, composition=x)
    eta = density / isotherm.density_per_eta
    step = 1e-6 * eta

    def derivative(point: float):
        return isotherm.derivatives(np.array([point]), second=True)

    here = derivative(eta)
    a1_fd = (float(derivative(eta + step).a[0]) - float(derivative(eta - step).a[0])) / (2.0 * step)
    a2_fd = (float(derivative(eta + step).a1[0]) - float(derivative(eta - step).a1[0])) / (
        2.0 * step
    )
    assert float(here.a1[0]) == pytest.approx(a1_fd, rel=1e-7, abs=1e-7)
    assert float(here.a2[0]) == pytest.approx(a2_fd, rel=1e-7, abs=1e-7)


@pytest.mark.parametrize("index", range(len(STATES)), ids=_state_id)
def test_the_pressure_slope_matches_a_finite_difference(index: int) -> None:
    names, x, temperature, density = STATES[index]
    eos = PCSAFTEOS(components=names)
    isotherm = eos._isotherm(names=names, temperature_K=temperature, composition=x)
    eta = density / isotherm.density_per_eta
    _, slope = isotherm.pressure_and_slope(eta)
    step = 1e-6 * density

    def pressure(rho: float) -> float:
        return float(isotherm.pressure(np.array([rho / isotherm.density_per_eta]))[0])

    expected = (pressure(density + step) - pressure(density - step)) / (2.0 * step)
    assert slope == pytest.approx(expected, rel=1e-7)
    assert slope > 0.0


@pytest.mark.parametrize("index", range(len(STATES)), ids=_state_id)
def test_michelsen_hendriks_stationarity(index: int) -> None:
    """``dQ/dX = 0`` at the solution - the fact the whole derivative route uses.

    ``Q`` is written out here from its definition (Eq. 5 of
    ``_pcsaft_association``) rather than called, so this checks the module's
    claim and not its arithmetic. Two things are asserted: the gradient
    vanishes, and ``Q`` at the solution equals the ``a_assoc`` the model
    reports - the substitution that makes the stationary form usable.
    """
    names, x_list, temperature, density = STATES[index]
    setup, x, m, d = _setup_for(names, x_list, temperature)
    rho_a3 = density * 6.02214076e23 * 1e-30
    delta, _, _ = _delta_at(setup, x, m, d, rho_a3)
    weights = setup.weights(x)
    x_sites = assoc.solve_site_fractions(rho=rho_a3, weights=weights, delta=delta)

    def q_value(sites: np.ndarray) -> float:
        explicit = float((weights * (np.log(sites) - sites + 1.0)).sum())
        bonded = float((weights * sites) @ delta @ (weights * sites))
        return explicit - 0.5 * rho_a3 * bonded

    gradient = weights * (1.0 / x_sites - 1.0) - rho_a3 * weights * (delta @ (weights * x_sites))
    assert float(np.max(np.abs(gradient))) < 1e-12

    # Finite-difference gradient too: the analytic one above could be wrong in
    # exactly the way the module is.
    for site in range(x_sites.size):
        step = 1e-7 * float(x_sites[site])
        up = x_sites.copy()
        up[site] += step
        down = x_sites.copy()
        down[site] -= step
        assert abs((q_value(up) - q_value(down)) / (2.0 * step)) < 1e-6

    carried = PCSAFTEOS(components=names)._state(temperature, density, x_list)
    assert q_value(x_sites) == pytest.approx(carried.a_assoc, rel=1e-12, abs=1e-14)
    assert float(assoc.helmholtz(weights, x_sites)) == pytest.approx(
        carried.a_assoc, rel=1e-14, abs=1e-15
    )
    assert float(assoc.q_site_term(weights, x_sites)) != pytest.approx(carried.a_assoc)


# --------------------------------------------------------------------------
# API behaviour
# --------------------------------------------------------------------------


def test_results_are_deterministic() -> None:
    names, x, temperature, density = STATES[5]
    eos = PCSAFTEOS(components=names)
    first = eos.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    second = eos.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    assert first == second


def test_results_are_permutation_invariant() -> None:
    temperature, density = 320.0, 16684.835
    forward = PCSAFTEOS(components=("Methanol", "Water", "n-Hexane"))
    reversed_order = PCSAFTEOS(components=("n-Hexane", "Water", "Methanol"))
    x = [0.3, 0.4, 0.3]
    a_forward = forward.residual_helmholtz(
        temperature_K=temperature, volume_m3=1.0 / density, composition=x
    )
    a_reversed = reversed_order.residual_helmholtz(
        temperature_K=temperature, volume_m3=1.0 / density, composition=list(reversed(x))
    )
    assert a_forward == pytest.approx(a_reversed, rel=1e-14, abs=1e-14)

    ln_phi_forward = forward.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    ln_phi_reversed = reversed_order.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=density, composition=list(reversed(x))
    )
    # Measured worst relative deviation with association: 3.8e-13 (the site
    # sums are accumulated in component order, so reordering changes the
    # summation order). Without association it is at round-off.
    np.testing.assert_allclose(ln_phi_forward, list(reversed(ln_phi_reversed)), rtol=1e-11)


def test_a_non_associating_component_may_share_a_mixture_with_an_associating_one() -> None:
    eos = PCSAFTEOS(components=("Water", "n-Hexane"))
    assert eos.associates()
    assert [record is None for record in eos.association_parameters()] == [False, True]
    ln_phi = eos.ln_fugacity_coefficients(
        temperature_K=298.15, density_mol_m3=10229.0998, composition=[0.3, 0.7]
    )
    assert all(math.isfinite(value) for value in ln_phi)


def test_pure_hexane_in_a_water_mixture_is_unchanged_at_zero_water() -> None:
    """A vanishing associating component must not perturb the rest.

    At ``x_water = 0`` the association term is identically zero, so the answer
    has to be the pure n-hexane one - including the site-fraction solve, which
    is where a division by a zero weight would show up.
    """
    mixed = PCSAFTEOS(components=("Water", "n-Hexane"))
    pure = PCSAFTEOS(components=("n-Hexane",))
    a_mixed = mixed.residual_helmholtz(
        temperature_K=300.0, volume_m3=1.0 / 7700.0, composition=[0.0, 1.0]
    )
    a_pure = pure.residual_helmholtz(temperature_K=300.0, volume_m3=1.0 / 7700.0, composition=[1.0])
    assert a_mixed == pytest.approx(a_pure, rel=0.0, abs=1e-15)


def test_density_roots_and_phase_identity_work_for_an_associating_fluid() -> None:
    eos = PCSAFTEOS(components=("Water",))
    mixture = ct.Mixture.from_database(["Water"], [1.0])
    roots = eos.density_roots(temperature_K=373.15, pressure_Pa=101325.0, composition=[1.0])
    assert len(roots) == 2
    assert roots[0] < 100.0 < roots[-1]
    assert (
        eos.phase_identity(
            mixture=mixture,
            temperature_K=373.15,
            pressure_Pa=101325.0,
            composition=[1.0],
            phase="vapor",
        )
        == "vapor"
    )
    assert (
        eos.phase_identity(
            mixture=mixture,
            temperature_K=373.15,
            pressure_Pa=101325.0,
            composition=[1.0],
            phase="liquid",
        )
        == "liquid"
    )


def test_pressure_on_a_density_root_reproduces_the_target() -> None:
    eos = PCSAFTEOS(components=("Water", "Ethanol"))
    target = 101325.0
    for root in eos.density_roots(temperature_K=351.0, pressure_Pa=target, composition=[0.2, 0.8]):
        actual = eos.pressure_Pa(temperature_K=351.0, density_mol_m3=root, composition=[0.2, 0.8])
        assert actual == pytest.approx(target, rel=1e-9)


def test_no_reachable_state_produces_a_nan() -> None:
    """A coarse sweep: finite answers, or the pre-existing ``eta`` refusal.

    1000 states over five associating systems. The only acceptable failure is
    the ADR-0014 packing-fraction guard, which fires **before** the association
    term is ever built.
    """
    systems: list[tuple[tuple[str, ...], list[float]]] = [
        (("Water",), [1.0]),
        (("Water", "Ethanol"), [0.5, 0.5]),
        (("Water", "n-Hexane"), [0.5, 0.5]),
        (("Methanol", "Water", "n-Hexane"), [0.3, 0.4, 0.3]),
        (("1-Propanol", "n-Butanol", "Water"), [0.3, 0.3, 0.4]),
    ]
    checked = 0
    for names, x in systems:
        eos = PCSAFTEOS(components=names)
        for temperature in (250.0, 298.15, 350.0, 500.0, 800.0):
            for density in np.geomspace(1e-6, 60000.0, 40):
                try:
                    a_res = eos.residual_helmholtz(
                        temperature_K=temperature, volume_m3=1.0 / float(density), composition=x
                    )
                    z_factor = eos.compressibility_factor(
                        temperature_K=temperature, density_mol_m3=float(density), composition=x
                    )
                    sites = eos.site_fractions(
                        temperature_K=temperature, density_mol_m3=float(density), composition=x
                    )
                except ct.ModelError as error:
                    assert "outside (0, 1)" in str(error), str(error)
                    continue
                checked += 1
                assert math.isfinite(a_res) and math.isfinite(z_factor)
                assert all(0.0 < value <= 1.0 and math.isfinite(value) for value in sites)
    assert checked > 800, checked


def test_stability_tp_accepts_an_associating_mixture() -> None:
    mixture = ct.Mixture.from_database(["Water", "n-Hexane"], [0.5, 0.5])
    result = ct.stability_tp(mixture, temperature_K=298.15, pressure_Pa=1.0e6, eos=PCSAFTEOS())
    assert result.status == "unstable"
