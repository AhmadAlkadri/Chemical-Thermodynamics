"""Modified-Raoult `flash_tp`: low-pressure VLE and LLE from one tangent plane.

Covers the mode contract, agreement with the existing `gamma-gamma` path where
the vapor candidate cannot win, the vapor-liquid invariants (modified Raoult's
law, mass balance, `delta_g_split_rt < 0`, the lever rule), bubble and dew
temperatures found by bisecting `flash_tp` *verdicts* and checked against an
independent scalar equation, the azeotrope of 1-propanol / water, and the
three-phase neighbourhood of water / 1-butanol, where a state the solver cannot
answer must raise rather than return.

Every reference number in this module is produced either by an independent
solve written inside this file (no `chemthermo.flash` code) or by a scalar
identity of the model. The NRTL parameters are the binary sub-systems of
Tessier, Brennecke and Stadtherr, Chem. Eng. Sci. 55 (2000) 1785, Table 1;
they are LLE-fitted and temperature independent, which is stated wherever a
comparison with literature is made.

See validation Cases R-1, R-2 and R-3.
"""

from __future__ import annotations

from typing import Any, Callable, Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.models._antoine import antoine_saturation_pressures

PRESSURE_PA = 101325.0

# Tessier (2000) Table 1 component order: 1 n-propanol, 2 n-butanol, 3 water.
_PROPANOL, _BUTANOL, _WATER = 0, 1, 2

#: A stability tolerance far below the default, used wherever a *boundary*
#: (bubble or dew temperature) is located by bisecting flash verdicts. The
#: default `tpd_tol = 1e-8` puts the verdict boundary at
#: `sum_i x_i gamma_i Psat_i / P = exp(tpd_tol)`, which is a 1e-8 relative
#: offset from the true bubble point - see
#: `test_the_verdict_boundary_sits_exactly_at_the_stability_tolerance`.
SHARP = ct.FlashSettings(stability_settings=ct.StabilitySettings(tpd_tol=1e-14))

LnGamma = Callable[[np.ndarray], np.ndarray]


# --------------------------------------------------------------------------
# Fixtures built from the cited Tessier Problem 1 parameters
# --------------------------------------------------------------------------


def _binary(payload: dict[str, Any], first: int, second: int) -> tuple[tuple[str, str], ct.NRTL]:
    """An (names, model) pair for one binary sub-system of Problem 1."""
    names = [str(entry["chemthermo_name"]) for entry in payload["components"]]
    tau = payload["tau"]
    alpha = payload["alpha"]
    pair = (names[first], names[second])
    model = ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [
                (
                    pair[0],
                    pair[1],
                    float(tau[first][second]),
                    float(tau[second][first]),
                    float(alpha[first][second]),
                    float(alpha[second][first]),
                )
            ]
        )
    )
    return pair, model


@pytest.fixture(scope="module")
def propanol_water(tessier2000_payload: dict[str, Any]) -> tuple[tuple[str, str], ct.NRTL]:
    """1-Propanol(1) / Water(2): fully miscible, minimum-boiling azeotrope."""
    return _binary(tessier2000_payload, _PROPANOL, _WATER)


@pytest.fixture(scope="module")
def butanol_water(tessier2000_payload: dict[str, Any]) -> tuple[tuple[str, str], ct.NRTL]:
    """n-Butanol(1) / Water(2): partially miscible, heteroazeotrope."""
    return _binary(tessier2000_payload, _BUTANOL, _WATER)


def _mixture(names: Sequence[str], z: Sequence[float]) -> ct.Mixture:
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


def _flash(
    names: Sequence[str],
    model: ct.NRTL,
    z: Sequence[float],
    temperature_K: float,
    settings: ct.FlashSettings | None = None,
) -> ct.FlashResult:
    return ct.flash_tp(
        _mixture(names, z),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        flash_mode="modified-raoult",
        settings=settings,
    )


def _ln_gamma(names: Sequence[str], model: ct.NRTL) -> LnGamma:
    """`ln gamma(x)` as a plain callable, written here and not taken from flash."""
    mixture = _mixture(names, (0.5, 0.5))

    def ln_gamma(x: np.ndarray) -> np.ndarray:
        values = np.asarray(x, dtype=float)
        values = values / float(np.sum(values))
        gamma = model.activity_coefficients(
            mixture=mixture,
            temperature_K=298.15,
            composition=[float(value) for value in values],
        )
        return np.log(np.asarray(gamma, dtype=float))

    return ln_gamma


def _psat(names: Sequence[str], temperature_K: float) -> np.ndarray:
    return antoine_saturation_pressures(_mixture(names, (0.5, 0.5)), temperature_K)


def _phase_set(result: ct.FlashResult) -> list[tuple[float, ...]]:
    return sorted(tuple(result.phases[name].composition.fractions) for name in result.phase_names())


def _bisect_verdict(
    names: Sequence[str],
    model: ct.NRTL,
    z: Sequence[float],
    low_K: float,
    high_K: float,
    *,
    iterations: int = 60,
) -> float:
    """Temperature where the `flash_tp` phase *count* changes between the bounds."""
    low_is_two = len(_flash(names, model, z, low_K, SHARP).phase_names()) == 2
    assert low_is_two != (len(_flash(names, model, z, high_K, SHARP).phase_names()) == 2), (
        "bracket does not straddle a verdict change"
    )
    for _ in range(iterations):
        middle = 0.5 * (low_K + high_K)
        if (len(_flash(names, model, z, middle, SHARP).phase_names()) == 2) == low_is_two:
            low_K = middle
        else:
            high_K = middle
    return 0.5 * (low_K + high_K)


def _bubble_temperature(ln_gamma: LnGamma, names: Sequence[str], x: np.ndarray) -> float:
    """Independent scalar solve of `sum_i x_i gamma_i Psat_i(T) = P`."""

    def residual(temperature: float) -> float:
        return float(np.sum(x * np.exp(ln_gamma(x)) * _psat(names, temperature))) - PRESSURE_PA

    low, high = 300.0, 399.0
    for _ in range(200):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if residual(middle) < 0.0 else (low, middle)
    return 0.5 * (low + high)


def _bubble_vapor(ln_gamma: LnGamma, names: Sequence[str], x: np.ndarray) -> tuple[float, float]:
    """Bubble temperature and the incipient vapor's first mole fraction."""
    temperature = _bubble_temperature(ln_gamma, names, x)
    y = x * np.exp(ln_gamma(x)) * _psat(names, temperature) / PRESSURE_PA
    return temperature, float(y[0])


def _dew_liquid(ln_gamma: LnGamma, names: Sequence[str], y: np.ndarray, temperature: float):
    """Liquid in equilibrium with vapor `y` at `temperature`, by fixed point."""
    x = y.copy()
    for _ in range(20000):
        updated = y * PRESSURE_PA / (np.exp(ln_gamma(x)) * _psat(names, temperature))
        updated = updated / float(np.sum(updated))
        if float(np.max(np.abs(updated - x))) < 1e-16:
            return updated
        x = updated
    raise AssertionError("dew-point fixed point did not converge")


# --------------------------------------------------------------------------
# Mode contract
# --------------------------------------------------------------------------


def test_modified_raoult_requires_an_activity_model(propanol_water) -> None:
    names, _model = propanol_water
    with pytest.raises(ct.ModelError, match="activity model is required"):
        ct.flash_tp(
            _mixture(names, (0.5, 0.5)),
            temperature_K=361.0,
            pressure_Pa=PRESSURE_PA,
            eos=ct.PengRobinsonEOS(),
            flash_mode="modified-raoult",
        )


def test_modified_raoult_forbids_an_eos(propanol_water) -> None:
    names, model = propanol_water
    with pytest.raises(ct.ModelError, match="must not be given"):
        ct.flash_tp(
            _mixture(names, (0.5, 0.5)),
            temperature_K=361.0,
            pressure_Pa=PRESSURE_PA,
            eos=ct.PengRobinsonEOS(),
            activity_model=model,
            flash_mode="modified-raoult",
        )


def test_the_mode_is_never_inferred(propanol_water) -> None:
    """An activity model on its own still means gamma-gamma, not modified Raoult.

    Which vapor model applies at a given pressure is the caller's physical
    judgement. Inferring it would silently change every existing
    liquid-liquid call.
    """
    names, model = propanol_water
    result = ct.flash_tp(
        _mixture(names, (0.5, 0.5)),
        temperature_K=361.0,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
    )
    assert result.diagnostics["flash_mode"] == "gamma-gamma"


def test_antoine_range_is_enforced_not_extrapolated(propanol_water) -> None:
    """1-Propanol's Antoine fit stops at 400 K; 420 K must raise, not guess."""
    names, model = propanol_water
    with pytest.raises(ct.InputRangeError, match=r"\[285.00, 400.00\] K"):
        _flash(names, model, (0.5, 0.5), 420.0)


def test_every_result_records_the_antoine_validity_window(propanol_water) -> None:
    names, model = propanol_water
    diagnostics = _flash(names, model, (0.5, 0.5), 330.0).diagnostics
    # Intersection of 1-Propanol [285, 400] and Water [284, 441].
    assert diagnostics["antoine_valid_Tmin_K"] == 285.0
    assert diagnostics["antoine_valid_Tmax_K"] == 400.0


# --------------------------------------------------------------------------
# Case R-1: consistency with the gamma-gamma path
# --------------------------------------------------------------------------


@pytest.mark.parametrize("z1", [0.05, 0.10, 0.20, 0.30])
def test_liquid_liquid_matches_the_gamma_gamma_path_exactly(butanol_water, z1: float) -> None:
    """At 330 K the vapor candidate is never the lower-Gibbs one, so the two
    paths must return the *same* tie-line.

    The liquid candidate's reference offset `ln(Psat_i / P)` cancels between two
    liquid phases, so `K_i = gamma_i^I / gamma_i^II` is literally the
    gamma-gamma update; anything but agreement at round-off would mean the
    reference had leaked into the split.
    """
    names, model = butanol_water
    mixture = _mixture(names, (z1, 1.0 - z1))
    liquid_liquid = ct.flash_tp(
        mixture,
        temperature_K=330.0,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        flash_mode="gamma-gamma",
    )
    modified_raoult = _flash(names, model, (z1, 1.0 - z1), 330.0)

    assert modified_raoult.diagnostics["phase_regime"] == "LLE"
    assert modified_raoult.vapor_fraction is None
    assert sorted(modified_raoult.phase_names()) == ["liquid1", "liquid2"]

    first, second = _phase_set(liquid_liquid), _phase_set(modified_raoult)
    worst = max(abs(a - b) for pa, pb in zip(first, second) for a, b in zip(pa, pb))
    assert worst < 1e-10, worst

    fractions_ll = sorted(liquid_liquid.phase_fractions.values())
    fractions_mr = sorted(modified_raoult.phase_fractions.values())
    assert max(abs(a - b) for a, b in zip(fractions_ll, fractions_mr)) < 1e-10


def test_a_feed_outside_the_gap_is_one_liquid(butanol_water) -> None:
    names, model = butanol_water
    result = _flash(names, model, (0.45, 0.55), 330.0)
    assert result.phase_names() == ["liquid"]
    assert result.diagnostics["phase_state"] == "liquid"
    assert result.diagnostics["feed_branch"] == "liquid"
    assert float(result.diagnostics["tpd_min"]) > 0.0
    assert result.vapor_fraction == 0.0


def test_both_liquid_phases_are_stable_against_the_vapor_candidate(butanol_water) -> None:
    """Case R-1's post-split half: at 330 K neither liquid boils."""
    names, model = butanol_water
    diagnostics = _flash(names, model, (0.20, 0.80), 330.0).diagnostics
    assert diagnostics["post_split_checked"] is True
    assert diagnostics["post_split_stable"] is True
    assert diagnostics["phase_stability_liquid1"] == "stable"
    assert diagnostics["phase_stability_liquid2"] == "stable"


# --------------------------------------------------------------------------
# Case R-2: vapor-liquid correctness
# --------------------------------------------------------------------------


@pytest.mark.parametrize("x1", [0.3, 0.5])
def test_bubble_temperature_from_flash_verdicts_satisfies_the_scalar_equation(
    propanol_water, x1: float
) -> None:
    """The single-liquid / two-phase boundary is the bubble point.

    `flash_tp` is treated as a black box returning a phase *count*; the
    temperature where that count changes is bisected, and the claim under test
    is the scalar identity `sum_i x_i gamma_i(x) Psat_i(T) = P`, evaluated by
    code in this file.
    """
    names, model = propanol_water
    ln_gamma = _ln_gamma(names, model)
    x = np.array([x1, 1.0 - x1])

    temperature = _bisect_verdict(names, model, (x1, 1.0 - x1), 350.0, 361.0)
    ratio = float(np.sum(x * np.exp(ln_gamma(x)) * _psat(names, temperature))) / PRESSURE_PA
    assert abs(ratio - 1.0) < 1e-8, (temperature, ratio)

    independent = _bubble_temperature(ln_gamma, names, x)
    assert abs(temperature - independent) < 1e-6


@pytest.mark.parametrize("y1, low_K, high_K", [(0.3, 364.0, 365.0), (0.5, 361.5, 362.0)])
def test_dew_temperature_from_flash_verdicts_satisfies_the_scalar_equation(
    propanol_water, y1: float, low_K: float, high_K: float
) -> None:
    """The two-phase / single-vapor boundary is the dew point.

    The check is `sum_i y_i P / (gamma_i(x) Psat_i(T)) = 1` with the implicit
    `x` solved here by this file's own fixed point.
    """
    names, model = propanol_water
    ln_gamma = _ln_gamma(names, model)
    y = np.array([y1, 1.0 - y1])

    temperature = _bisect_verdict(names, model, (y1, 1.0 - y1), low_K, high_K)
    x = _dew_liquid(ln_gamma, names, y, temperature)
    residual = float(np.sum(y * PRESSURE_PA / (np.exp(ln_gamma(x)) * _psat(names, temperature))))
    assert abs(residual - 1.0) < 1e-8, (temperature, residual)


def test_the_verdict_boundary_sits_exactly_at_the_stability_tolerance(propanol_water) -> None:
    """Why `SHARP` exists, stated as an identity rather than a fudge.

    At a stationary point `tpd = -ln(sum_i W_i)`, and the vapor trial's first
    substitution from a liquid feed gives `W_i = x_i gamma_i Psat_i / P`. So a
    feed is declared unstable exactly when
    `sum_i x_i gamma_i Psat_i / P > exp(tpd_tol)`: the verdict boundary is the
    bubble point displaced by the tolerance, and nothing else.
    """
    names, model = propanol_water
    ln_gamma = _ln_gamma(names, model)
    x = np.array([0.5, 0.5])

    settings = ct.FlashSettings()
    low, high = 350.0, 361.0
    for _ in range(60):
        middle = 0.5 * (low + high)
        if len(_flash(names, model, (0.5, 0.5), middle, settings).phase_names()) == 1:
            low = middle
        else:
            high = middle
    boundary = 0.5 * (low + high)
    ratio = float(np.sum(x * np.exp(ln_gamma(x)) * _psat(names, boundary))) / PRESSURE_PA
    assert abs(ratio - np.exp(1e-8)) < 1e-11, ratio


def test_vapor_liquid_split_satisfies_modified_raoults_law(propanol_water) -> None:
    """`y_i P = x_i gamma_i Psat_i` to 1e-10, plus mass balance and `dG < 0`."""
    names, model = propanol_water
    ln_gamma = _ln_gamma(names, model)
    z = (0.5, 0.5)
    result = _flash(names, model, z, 361.0)

    assert result.phase_names() == ["liquid", "vapor"]
    assert result.diagnostics["phase_regime"] == "VLE"
    assert result.diagnostics["incipient_phase"] == "vapor"
    assert result.diagnostics["feed_branch"] == "liquid"

    x = np.array(result.phases["liquid"].composition.fractions)
    y = np.array(result.phases["vapor"].composition.fractions)
    beta = result.vapor_fraction
    assert beta is not None

    psat = _psat(names, 361.0)
    residual = float(np.max(np.abs(y * PRESSURE_PA - x * np.exp(ln_gamma(x)) * psat)))
    assert residual / PRESSURE_PA < 1e-10, residual

    assert float(np.max(np.abs(np.array(z) - (beta * y + (1.0 - beta) * x)))) < 1e-12
    assert float(result.diagnostics["mass_balance_residual"]) < 1e-12
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert result.diagnostics["post_split_stable"] is True

    # Lever rule on component 1, independent of the reported phase fractions.
    lever = (z[0] - x[0]) / (y[0] - x[0])
    assert abs(lever - beta) < 1e-12
    assert abs(result.phase_fractions["vapor"] - beta) < 1e-15
    assert abs(result.phase_fractions["liquid"] - (1.0 - beta)) < 1e-15


def test_a_superheated_feed_is_a_single_vapor(propanol_water) -> None:
    names, model = propanol_water
    result = _flash(names, model, (0.5, 0.5), 380.0)
    assert result.phase_names() == ["vapor"]
    assert result.vapor_fraction == 1.0
    assert result.diagnostics["feed_branch"] == "vapor"


def test_the_azeotrope_matches_an_independent_solve_and_is_a_vanishing_band(
    propanol_water,
) -> None:
    """Locate `y = x` on the bubble curve here, then check `flash_tp` agrees.

    The independent route is a 1-D bisection on `y_1(x_1) - x_1` along the
    bubble curve solved in this file. The package-side signature of an
    azeotrope is that the two-phase band closes: at the azeotropic composition
    `flash_tp` goes straight from one liquid to one vapor.

    Literature (commonly tabulated, **not** read from a primary source in this
    work): the 1-propanol / water minimum-boiling azeotrope at 1 atm is quoted
    near 87.7 C and x(1-propanol) ~ 0.43. The NRTL parameters used here were
    fitted to *liquid-liquid* data for the ternary and are temperature
    independent, so any agreement is partly fortuitous; the deviation is
    reported, not tuned away.
    """
    names, model = propanol_water
    ln_gamma = _ln_gamma(names, model)

    def excess(x1: float) -> float:
        return _bubble_vapor(ln_gamma, names, np.array([x1, 1.0 - x1]))[1] - x1

    low, high = 0.2, 0.8
    assert excess(low) > 0.0 and excess(high) < 0.0
    for _ in range(100):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if excess(middle) > 0.0 else (low, middle)
    x_azeotrope = 0.5 * (low + high)
    t_azeotrope = _bubble_temperature(ln_gamma, names, np.array([x_azeotrope, 1.0 - x_azeotrope]))

    assert abs(x_azeotrope - 0.419874) < 1e-5
    assert abs(t_azeotrope - 360.917975) < 1e-4
    # Deviation from the commonly tabulated values, recorded not asserted away.
    assert abs((t_azeotrope - 273.15) - 87.7) < 0.2
    assert abs(x_azeotrope - 0.43) < 0.02

    # Package side: the incipient vapor from `stability_tp` straddles it.
    below = _incipient_vapor(names, model, ln_gamma, x_azeotrope - 0.05)
    above = _incipient_vapor(names, model, ln_gamma, x_azeotrope + 0.05)
    assert below > x_azeotrope - 0.05
    assert above < x_azeotrope + 0.05

    # ... and the two-phase band has closed at the azeotrope itself.
    z = (x_azeotrope, 1.0 - x_azeotrope)
    assert _flash(names, model, z, t_azeotrope - 0.01, SHARP).phase_names() == ["liquid"]
    assert _flash(names, model, z, t_azeotrope + 0.01, SHARP).phase_names() == ["vapor"]


def _incipient_vapor(names: Sequence[str], model: ct.NRTL, ln_gamma: LnGamma, x1: float) -> float:
    """The package's incipient vapor at the bubble point of `x1`.

    At the bubble temperature the tangent-plane minimizer *is* the equilibrium
    vapor, so this reads `y(x)` straight off `stability_tp`.
    """
    x = np.array([x1, 1.0 - x1])
    temperature = _bubble_temperature(ln_gamma, names, x)
    result = ct.stability_tp(
        _mixture(names, (x1, 1.0 - x1)),
        temperature_K=temperature + 1e-8,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        vapor="ideal",
        settings=ct.StabilitySettings(tpd_tol=1e-14),
    )
    assert result.status == "unstable"
    assert result.phase_branch == "vapor"
    assert result.trial_composition is not None
    expected = _bubble_vapor(ln_gamma, names, x)[1]
    assert abs(result.trial_composition[0] - expected) < 1e-8
    return float(result.trial_composition[0])


# --------------------------------------------------------------------------
# Case R-3: the three-phase neighbourhood of water / 1-butanol
# --------------------------------------------------------------------------

#: Three-phase temperature at 101325 Pa, from the independent solve in
#: `test_three_phase_temperature_is_where_both_liquids_boil`.
T3_K = 366.2137741


def test_three_phase_temperature_is_where_both_liquids_boil(butanol_water) -> None:
    """T3 and the implied vapor, computed here from Antoine + NRTL only.

    At the three-phase point both conjugate liquids are at their bubble point
    simultaneously: `sum_i x_i gamma_i Psat_i(T3) = P` for *both*. Solving that
    on the water-rich branch and checking it on the butanol-rich branch is a
    non-trivial identity, because the binodal came from equal activities and
    knows nothing about Antoine.
    """
    names, model = butanol_water
    ln_gamma = _ln_gamma(names, model)
    first, second = _binodal(ln_gamma)
    assert abs(first - 0.019998419467) < 1e-10
    assert abs(second - 0.359999661508) < 1e-10

    x_i = np.array([first, 1.0 - first])
    x_ii = np.array([second, 1.0 - second])

    def residual(temperature: float) -> float:
        return float(np.sum(x_i * np.exp(ln_gamma(x_i)) * _psat(names, temperature))) - PRESSURE_PA

    low, high = 340.0, 399.0
    for _ in range(200):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if residual(middle) < 0.0 else (low, middle)
    t3 = 0.5 * (low + high)
    assert abs(t3 - T3_K) < 1e-4

    for x in (x_i, x_ii):
        ratio = float(np.sum(x * np.exp(ln_gamma(x)) * _psat(names, t3))) / PRESSURE_PA
        assert abs(ratio - 1.0) < 1e-10, ratio

    y = x_i * np.exp(ln_gamma(x_i)) * _psat(names, t3) / PRESSURE_PA
    assert abs(y[0] - 0.23406) < 1e-5
    assert abs(y[1] - 0.76594) < 1e-5

    gamma_i = np.exp(ln_gamma(x_i))
    gamma_ii = np.exp(ln_gamma(x_ii))
    assert abs(gamma_i[0] - 30.558) < 1e-3 and abs(gamma_i[1] - 1.00658) < 1e-5
    assert abs(gamma_ii[0] - 1.69755) < 1e-5 and abs(gamma_ii[1] - 1.54133) < 1e-5


def _binodal(ln_gamma: LnGamma) -> tuple[float, float]:
    """Conjugate liquid compositions from equal activities (Newton, FD Jacobian)."""

    def activities(x1: float) -> np.ndarray:
        x = np.array([x1, 1.0 - x1], dtype=float)
        return x * np.exp(ln_gamma(x))

    def mismatch(u: np.ndarray) -> np.ndarray:
        return activities(float(u[0])) - activities(float(u[1]))

    u = np.array([0.02, 0.50], dtype=float)
    for _ in range(200):
        residual = mismatch(u)
        if float(np.max(np.abs(residual))) < 1e-15:
            break
        jacobian = np.zeros((2, 2), dtype=float)
        for column in range(2):
            plus, minus = u.copy(), u.copy()
            plus[column] += 1e-7
            minus[column] -= 1e-7
            jacobian[:, column] = (mismatch(plus) - mismatch(minus)) / 2e-7
        u = u + np.linalg.solve(jacobian, -residual)
    assert float(np.max(np.abs(mismatch(u)))) < 1e-14
    return float(u[0]), float(u[1])


def test_two_kelvin_below_t3_is_a_stable_liquid_liquid_split(butanol_water) -> None:
    """Control (i): both liquids present and neither boils."""
    names, model = butanol_water
    result = _flash(names, model, (0.20, 0.80), T3_K - 2.0)
    assert sorted(result.phase_names()) == ["liquid1", "liquid2"]
    assert result.vapor_fraction is None
    diagnostics = result.diagnostics
    assert diagnostics["post_split_stable"] is True
    assert diagnostics["post_split_status"] == "stable"
    assert float(diagnostics["post_split_tpd_min"]) > -1e-8
    compositions = _phase_set(result)
    assert abs(compositions[0][0] - 0.019998419467) < 1e-9
    assert abs(compositions[1][0] - 0.359999661508) < 1e-9


def test_two_kelvin_above_t3_the_feed_has_evaporated(butanol_water) -> None:
    """Control (ii): a single vapor, verified against the dew-point equation.

    The brief allowed either a vapor-liquid answer or a fully evaporated feed;
    the solver returns the latter, so it is checked as such: `z` lies above its
    dew point, i.e. `sum_i z_i P / (gamma_i(x) Psat_i) < 1` for the liquid `x`
    that would be in equilibrium with it.
    """
    names, model = butanol_water
    ln_gamma = _ln_gamma(names, model)
    temperature = T3_K + 2.0
    result = _flash(names, model, (0.20, 0.80), temperature)
    assert result.phase_names() == ["vapor"]
    assert result.vapor_fraction == 1.0
    assert result.diagnostics["feed_branch"] == "vapor"
    assert float(result.diagnostics["tpd_min"]) > 0.0

    z = np.array([0.20, 0.80])
    x = _dew_liquid(ln_gamma, names, z, temperature)
    dew_residual = float(
        np.sum(z * PRESSURE_PA / (np.exp(ln_gamma(x)) * _psat(names, temperature)))
    )
    assert dew_residual < 1.0, dew_residual


def test_at_t3_the_two_phase_answer_is_refused(butanol_water) -> None:
    """Control (iii): the knife edge at T3, and the refusal window with `max_phases=2`.

    At T3 the solver converges a vapor-liquid pair whose liquid sits *on* the
    binodal. Which side of the knife-edge the post-split test lands on is set
    by how exactly T3 is known, so both outcomes are accepted here and the
    numbers are recorded: at T3 to 1e-7 K the pair is marginal
    (`post_split_tpd_min` at round-off), and 0.01 K below it the same pair is
    reported unstable.

    Below T3 the answer is now *resolved* by phase addition and removal
    (ADR-0011, Case V-3); this test pins the pre-ADR-0011 behavior, which
    `FlashSettings(max_phases=2)` still reproduces exactly.
    """
    names, model = butanol_water
    two_phase = ct.FlashSettings(max_phases=2)

    try:
        result = _flash(names, model, (0.20, 0.80), T3_K, two_phase)
    except ct.ConvergenceError as error:
        assert "third phase is required" in str(error)
    else:
        assert sorted(result.phase_names()) == ["liquid", "vapor"]
        y = np.array(result.phases["vapor"].composition.fractions)
        x = np.array(result.phases["liquid"].composition.fractions)
        assert abs(y[0] - 0.23406) < 1e-5
        assert abs(x[0] - 0.019998419467) < 1e-6
        assert abs(float(result.diagnostics["post_split_tpd_min"])) < 1e-10

    # Just below T3 the same vapor-liquid pair is genuinely not the answer: a
    # third stationary point lies below its tangent plane, and both phases see
    # it (two coexisting phases share one plane, so the tpd is the same from
    # either). This is the documented negative control.
    with pytest.raises(ct.ConvergenceError, match="third phase is required"):
        _flash(names, model, (0.20, 0.80), T3_K - 0.01, two_phase)

    unchecked = _flash(
        names,
        model,
        (0.20, 0.80),
        T3_K - 0.01,
        ct.FlashSettings(post_split_stability=False),
    )
    diagnostics = unchecked.diagnostics
    assert diagnostics["post_split_status"] == "unstable"
    assert diagnostics["phase_stability_vapor"] == "unstable"
    assert diagnostics["phase_stability_liquid"] == "unstable"
    assert float(diagnostics["post_split_tpd_min"]) < -1e-4
    assert float(diagnostics["phase_stability_tpd_min_vapor"]) == pytest.approx(
        float(diagnostics["phase_stability_tpd_min_liquid"]), abs=1e-12
    )
