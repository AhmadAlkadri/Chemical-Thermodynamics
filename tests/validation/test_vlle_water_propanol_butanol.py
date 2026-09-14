"""Validation Cases V-1 and V-2: the ternary vapor-liquid-liquid tie-triangle.

System
------
1-Propanol(1) / n-Butanol(2) / Water(3) at P = 101325 Pa, modified Raoult:
an NRTL liquid with Antoine pure-liquid reference fugacities against an ideal
vapor (ADR-0010). NRTL parameters are Table 1 of S. R. Tessier, J. F. Brennecke
and M. A. Stadtherr, Chem. Eng. Sci. 55 (2000) 1785-1796 (attributed there to
McDonald and Floudas, AIChE J. 41 (1995) 1798), loaded from
`tests/fixtures/nrtl/tessier2000_problem1.json`; Antoine coefficients from the
packaged databank (Koretsky 2012).

**These parameters were fitted to liquid-liquid data and are temperature
independent, and no experimental ternary VLLE data is used here.** Nothing in
this module is a comparison against measurement. What is validated is that
`flash_tp` finds the three-phase state that *this model* has, and finds it to
the accuracy of an independent solve.

The independent route
---------------------
At fixed ``(T, P)`` a three-phase state of a ternary is determined by six
unknowns - the two liquid compositions - and six equations, written and solved
in this module by a damped Newton iteration with a finite-difference Jacobian
that shares no code with `chemthermo.flash`:

    ln x_i^I + ln gamma_i(x^I) - ln x_i^II - ln gamma_i(x^II) = 0     (3 eq)
    sum_i x_i^I - 1 = 0,   sum_i x_i^II - 1 = 0                       (2 eq)
    sum_i y_i - 1 = 0,     y_i = x_i^I gamma_i(x^I) Psat_i(T) / P     (1 eq)

The last equation is the bubble condition on liquid I; that liquid II is *also*
at its bubble point is then a consequence, not an input, and is checked. The
three vertices ``(x^I, x^II, y)`` are the tie-triangle; a feed is inside it
exactly when the weights solving ``z = beta_I x^I + beta_II x^II + beta_V y``
are all positive, which is how the feeds below are chosen.

Gibbs ordering
--------------
For each three-phase answer this module also computes, from its own equations,

    G/RT = sum_j beta_j sum_i x_i^j [ ln x_i^j + t_i^j ]

with ``t_i = ln gamma_i + ln(Psat_i/P)`` for a liquid and ``0`` for the vapor,
and checks ``G(three phases) < G(the two-phase candidate) < G(one phase)``. A
lower Gibbs energy is what makes the three-phase answer the answer; equal
fugacities alone would be satisfied by the two-phase state as well.

No external reference
---------------------
`thermo` 0.6.0 cannot be made to split two liquids over a single excess-Gibbs
model (validation Case L-2), so there is no external three-phase reference for
this system and none is claimed. The independent Newton solve above and the
Gibbs ordering are the whole of the evidence.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pytest

import chemthermo as ct

PRESSURE_PA = 101325.0
FIXTURE_PATH = (
    Path(__file__).resolve().parents[1] / "fixtures" / "nrtl" / "tessier2000_problem1.json"
)

#: Newton tolerance of the independent tie-triangle solve (achieved <= 8.9e-16).
TRIANGLE_TOL = 1e-14

LnGamma = Callable[[np.ndarray], np.ndarray]


# --------------------------------------------------------------------------
# The model, written here from the fixture, not taken from chemthermo.models
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def system() -> tuple[list[str], ct.NRTL, LnGamma, Callable[[float], np.ndarray]]:
    with FIXTURE_PATH.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    names = [str(entry["chemthermo_name"]) for entry in payload["components"]]
    tau = np.array(payload["tau"], dtype=float)
    alpha = np.array(payload["alpha"], dtype=float)
    capital_g = np.exp(-alpha * tau)

    model = ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [
                (
                    names[i],
                    names[j],
                    float(tau[i][j]),
                    float(tau[j][i]),
                    float(alpha[i][j]),
                    float(alpha[j][i]),
                )
                for i in range(len(names))
                for j in range(i + 1, len(names))
            ]
        )
    )

    def ln_gamma(x: np.ndarray) -> np.ndarray:
        """Renon-Prausnitz NRTL with column sums, written out here."""
        values = np.asarray(x, dtype=float)
        values = values / float(np.sum(values))
        s = capital_g.T @ values
        c = (tau * capital_g).T @ values
        return np.array(
            [
                c[i] / s[i] + float(np.sum(values * capital_g[i, :] / s * (tau[i, :] - c / s)))
                for i in range(len(names))
            ]
        )

    antoine = []
    for name in names:
        record = ct.Component.from_database(name).antoine
        assert record is not None, name
        antoine.append(record)

    def psat(temperature_K: float) -> np.ndarray:
        """``ln(P^sat / bar) = A - B / (T + C)`` from the packaged databank."""
        return np.array(
            [np.exp(a.A - a.B / (temperature_K + a.C)) * 1.0e5 for a in antoine],
            dtype=float,
        )

    return names, model, ln_gamma, psat


# --------------------------------------------------------------------------
# The independent tie-triangle solve
# --------------------------------------------------------------------------


def _triangle_residual(
    u: np.ndarray, temperature_K: float, ln_gamma: LnGamma, psat: Callable[[float], np.ndarray]
) -> np.ndarray:
    x_i, x_ii = u[:3], u[3:]
    activity_i = np.log(np.abs(x_i)) + ln_gamma(x_i)
    activity_ii = np.log(np.abs(x_ii)) + ln_gamma(x_ii)
    y = np.exp(activity_i) * psat(temperature_K) / PRESSURE_PA
    return np.concatenate(
        [
            activity_i - activity_ii,
            [float(np.sum(x_i)) - 1.0, float(np.sum(x_ii)) - 1.0, float(np.sum(y)) - 1.0],
        ]
    )


def _tie_triangle(
    temperature_K: float,
    start: Sequence[float],
    ln_gamma: LnGamma,
    psat: Callable[[float], np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Damped Newton on the six equations of the module docstring."""
    u = np.array(start, dtype=float)
    residual = float("inf")
    for _ in range(200):
        f = _triangle_residual(u, temperature_K, ln_gamma, psat)
        residual = float(np.max(np.abs(f)))
        if residual < TRIANGLE_TOL:
            break
        jacobian = np.zeros((6, 6), dtype=float)
        step = 1e-8
        for column in range(6):
            forward = u.copy()
            forward[column] += step
            jacobian[:, column] = (
                _triangle_residual(forward, temperature_K, ln_gamma, psat) - f
            ) / step
        direction = np.linalg.solve(jacobian, -f)
        scale = 1.0
        while np.any(u + scale * direction <= 0.0):
            scale *= 0.5
        u = u + scale * direction

    x_i, x_ii = u[:3], u[3:]
    y = np.exp(np.log(x_i) + ln_gamma(x_i)) * psat(temperature_K) / PRESSURE_PA
    return x_i, x_ii, y, residual


#: Newton starting points, one per temperature (any point near the triangle works).
TRIANGLE_SEEDS = {
    365.0: (0.023616, 0.024407, 0.951977, 0.098737, 0.201955, 0.699308),
    364.0: (0.052144, 0.029172, 0.918684, 0.139639, 0.123478, 0.736883),
    363.0: (0.102828, 0.035390, 0.861782, 0.156303, 0.064227, 0.779470),
}

#: Barycentric weights of the feeds tested inside each triangle. All positive,
#: so every feed is inside by construction; the test asserts that too.
INSIDE_WEIGHTS = ((1 / 3, 1 / 3, 1 / 3), (0.2, 0.3, 0.5))


def _reduced_g(x: np.ndarray, is_vapor: bool, temperature_K: float, ln_gamma, psat) -> float:
    """``sum_i x_i (ln x_i + t_i)``, the composition-dependent part of ``G/RT``."""
    values = np.asarray(x, dtype=float)
    terms = (
        np.zeros_like(values)
        if is_vapor
        else ln_gamma(values) + np.log(psat(temperature_K) / PRESSURE_PA)
    )
    mask = values > 0.0
    return float(np.sum(values[mask] * (np.log(values[mask]) + terms[mask])))


def _flash(names, model, z, temperature_K, settings=None) -> ct.FlashResult:
    return ct.flash_tp(
        ct.Mixture.from_database(list(names), [float(value) for value in z], normalize=True),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        flash_mode="modified-raoult",
        settings=settings,
    )


def _phase_set(result: ct.FlashResult) -> list[tuple[float, ...]]:
    return sorted(tuple(result.phases[name].composition.fractions) for name in result.phase_names())


# --------------------------------------------------------------------------
# Case V-1: the tie-triangle
# --------------------------------------------------------------------------


@pytest.mark.parametrize("temperature_K", sorted(TRIANGLE_SEEDS))
def test_the_independent_tie_triangle_is_a_genuine_three_phase_state(
    system, temperature_K: float
) -> None:
    """Both liquids boil at once and share one vapor - the input to Case V-1.

    The Newton solve imposes the bubble condition on liquid I only. That liquid
    II is simultaneously at its bubble point, and that the vapor computed from
    either liquid is the same vapor, are consequences of the six equations and
    are checked here before any of these numbers is used as a reference.
    """
    _names, _model, ln_gamma, psat = system
    x_i, x_ii, y, residual = _tie_triangle(
        temperature_K, TRIANGLE_SEEDS[temperature_K], ln_gamma, psat
    )

    assert residual < TRIANGLE_TOL, residual
    # The equation that was *not* solved for: liquid II is also at its bubble
    # point. Achieved <= 2e-15 at all three temperatures.
    bubble_ii = float(np.sum(x_ii * np.exp(ln_gamma(x_ii)) * psat(temperature_K) / PRESSURE_PA))
    assert abs(bubble_ii - 1.0) < 1e-10, (temperature_K, bubble_ii)
    # And the vapor computed from liquid II is the vapor computed from liquid I.
    y_from_ii = x_ii * np.exp(ln_gamma(x_ii)) * psat(temperature_K) / PRESSURE_PA
    assert np.max(np.abs(y_from_ii - y)) < 1e-10
    # Three distinct phases, each a composition.
    for phase in (x_i, x_ii, y):
        assert abs(float(np.sum(phase)) - 1.0) < 1e-12
        assert np.all(phase > 0.0)
    assert np.max(np.abs(x_i - x_ii)) > 1e-3


@pytest.mark.parametrize("temperature_K", sorted(TRIANGLE_SEEDS))
@pytest.mark.parametrize(
    "weights",
    # ADR-0028 runtime trim: a second feed inside the same tie triangle says
    # the same thing about the same state - the compositions are the triangle's
    # and only the amounts move, which is what
    # `test_every_three_phase_answer_is_the_same_tie_triangle` asserts.
    [INSIDE_WEIGHTS[0], *(pytest.param(w, marks=pytest.mark.slow) for w in INSIDE_WEIGHTS[1:])],
)
def test_a_feed_inside_the_tie_triangle_returns_three_verified_phases(
    system, temperature_K: float, weights
) -> None:
    """Case V-1. Compositions and phase fractions against the independent solve.

    Achieved over the six (temperature, feed) combinations: worst composition
    deviation 2.2e-14, worst phase-fraction deviation 1.6e-13, worst
    equilibrium residual 1.8e-15, worst mass-balance residual 6.9e-18, worst
    post-split ``tpd_min`` -4.1e-16.
    """
    names, model, ln_gamma, psat = system
    x_i, x_ii, y, _ = _tie_triangle(temperature_K, TRIANGLE_SEEDS[temperature_K], ln_gamma, psat)
    vertices = np.column_stack([x_i, x_ii, y])

    z = vertices @ np.array(weights, dtype=float)
    z = z / float(np.sum(z))
    exact = np.linalg.solve(vertices, z)
    assert np.all(exact > 0.0), (temperature_K, exact)

    result = _flash(names, model, z, temperature_K)
    diagnostics = result.diagnostics

    assert sorted(result.phase_names()) == ["liquid1", "liquid2", "vapor"]
    assert diagnostics["phase_count"] == 3
    assert diagnostics["phase_regime"] == "VLLE"
    assert diagnostics["phase_state"] == "three_phase"
    assert str(diagnostics["phase_set_history"]).endswith("LLV")
    assert diagnostics["phases_added"] == 1
    assert diagnostics["phases_removed"] == 0

    # Compositions: the phase set, matched to the reference by nearest vertex.
    computed = {
        name: np.array(result.phases[name].composition.fractions) for name in result.phase_names()
    }
    for reference, label in ((x_i, "I"), (x_ii, "II"), (y, "V")):
        nearest = min(computed, key=lambda name: float(np.max(np.abs(computed[name] - reference))))
        deviation = float(np.max(np.abs(computed[nearest] - reference)))
        assert deviation < 1e-6, (temperature_K, label, deviation)
        assert (nearest == "vapor") == (label == "V")

    # Phase fractions, matched the same way.
    for index, reference in enumerate((x_i, x_ii, y)):
        nearest = min(computed, key=lambda name: float(np.max(np.abs(computed[name] - reference))))
        assert result.phase_fractions[nearest] == pytest.approx(exact[index], abs=1e-6)

    assert float(diagnostics["equilibrium_residual"]) < 1e-10
    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert diagnostics["post_split_stable"] is True
    assert diagnostics["post_split_status"] == "stable"
    for name in result.phase_names():
        assert diagnostics[f"phase_stability_{name}"] == "stable"
    assert float(diagnostics["post_split_tpd_min"]) > -1e-10

    assert result.vapor_fraction == pytest.approx(result.phase_fractions["vapor"], abs=0.0)


@pytest.mark.parametrize("temperature_K", sorted(TRIANGLE_SEEDS))
@pytest.mark.parametrize(
    "weights",
    # ADR-0028 runtime trim, as just above.
    [INSIDE_WEIGHTS[0], *(pytest.param(w, marks=pytest.mark.slow) for w in INSIDE_WEIGHTS[1:])],
)
def test_the_three_phase_state_has_the_lowest_gibbs_energy(
    system, temperature_K: float, weights
) -> None:
    """``G(3 phases) < G(2-phase candidate) < G(feed)``, all computed here.

    The two-phase candidate is the split `flash_tp` converged *before* the
    third phase was added, recovered by re-running with
    ``FlashSettings(post_split_stability=False)``. It satisfies equal
    fugacities, which is why the ordering rather than the residual is what
    identifies the answer.
    """
    names, model, ln_gamma, psat = system
    x_i, x_ii, y, _ = _tie_triangle(temperature_K, TRIANGLE_SEEDS[temperature_K], ln_gamma, psat)
    vertices = np.column_stack([x_i, x_ii, y])
    z = vertices @ np.array(weights, dtype=float)
    z = z / float(np.sum(z))

    three = _flash(names, model, z, temperature_K)
    two = _flash(names, model, z, temperature_K, ct.FlashSettings(post_split_stability=False))
    assert len(two.phase_names()) == 2

    def energy(result: ct.FlashResult) -> float:
        return float(
            sum(
                result.phase_fractions[name]
                * _reduced_g(
                    np.array(result.phases[name].composition.fractions),
                    name == "vapor",
                    temperature_K,
                    ln_gamma,
                    psat,
                )
                for name in result.phase_names()
            )
        )

    g_three = energy(three)
    g_two = energy(two)
    g_feed = _reduced_g(z, False, temperature_K, ln_gamma, psat)

    assert g_three < g_two < g_feed, (temperature_K, g_three, g_two, g_feed)
    # The package's own numbers agree with this independent computation.
    assert float(three.diagnostics["delta_g_split_rt"]) == pytest.approx(
        g_three - g_feed, abs=1e-12
    )
    assert float(three.diagnostics["delta_g_vs_two_phase_rt"]) == pytest.approx(
        g_three - g_two, abs=1e-9
    )
    assert float(three.diagnostics["delta_g_vs_two_phase_rt"]) < 0.0


def test_three_phase_results_are_deterministic_and_permutation_invariant(system) -> None:
    """Case V-1 invariants: the phase *set* survives reordering the components.

    The `liquid1` / `liquid2` labels are roles and may swap, exactly as for a
    two-phase liquid-liquid result, so the comparison is on the sorted set of
    compositions. `vapor` is the one name with a model-level meaning and must
    stay on the vapor. Achieved over all six orderings: worst composition
    difference after undoing the permutation 4.3e-15, worst phase-fraction
    difference 1.5e-14.
    """
    names, model, ln_gamma, psat = system
    temperature_K = 364.0
    x_i, x_ii, y, _ = _tie_triangle(temperature_K, TRIANGLE_SEEDS[temperature_K], ln_gamma, psat)
    z = np.column_stack([x_i, x_ii, y]) @ np.array([1 / 3, 1 / 3, 1 / 3])
    z = z / float(np.sum(z))

    base = _flash(names, model, z, temperature_K)
    repeat = _flash(names, model, z, temperature_K)
    assert dict(repeat.diagnostics) == dict(base.diagnostics)
    assert _phase_set(repeat) == _phase_set(base)

    reference = np.array(sorted(tuple(value) for value in _phase_set(base)))
    for order in itertools.permutations(range(len(names))):
        permuted_names = [names[i] for i in order]
        permuted_model = ct.NRTL(parameters=model.parameters)
        permuted = _flash(permuted_names, permuted_model, [z[i] for i in order], temperature_K)

        assert len(permuted.phase_names()) == 3
        assert sorted(permuted.phase_names()) == ["liquid1", "liquid2", "vapor"]

        undo = np.argsort(order)
        restored = np.array(
            sorted(
                tuple(np.array(permuted.phases[name].composition.fractions)[undo].tolist())
                for name in permuted.phase_names()
            )
        )
        assert np.allclose(restored, reference, rtol=0.0, atol=1e-9)

        # The vapor is the vapor whatever order the components came in.
        vapor = np.array(permuted.phases["vapor"].composition.fractions)[undo]
        assert np.max(np.abs(vapor - np.array(base.phases["vapor"].composition.fractions))) < 1e-9
        assert permuted.vapor_fraction == pytest.approx(base.vapor_fraction, abs=1e-9)

        for name in permuted.phase_names():
            composition = np.array(permuted.phases[name].composition.fractions)[undo]
            partner = min(
                base.phase_names(),
                key=lambda other: float(
                    np.max(np.abs(np.array(base.phases[other].composition.fractions) - composition))
                ),
            )
            assert permuted.phase_fractions[name] == pytest.approx(
                base.phase_fractions[partner], abs=1e-9
            )


# --------------------------------------------------------------------------
# Case V-2: negative controls
# --------------------------------------------------------------------------


def test_a_feed_in_the_vapor_liquid_region_returns_two_phases(system) -> None:
    """Outside the triangle past the ``x^II - y`` edge: vapor plus one liquid.

    Verified independently of the phase count: the liquid returned is *stable*
    against a second liquid (`stability_tp` with `vapor="none"` on it), so no
    third phase exists at this feed, and the barycentric weight of ``x^I`` in
    the tie-triangle is negative, which is what "outside" means.
    """
    names, model, ln_gamma, psat = system
    temperature_K = 364.0
    x_i, x_ii, y, _ = _tie_triangle(temperature_K, TRIANGLE_SEEDS[temperature_K], ln_gamma, psat)
    vertices = np.column_stack([x_i, x_ii, y])

    z = vertices @ np.array([-0.25, 0.45, 0.80])
    z = z / float(np.sum(z))
    weights = np.linalg.solve(vertices, z)
    assert weights[0] < 0.0, weights

    result = _flash(names, model, z, temperature_K)
    assert sorted(result.phase_names()) == ["liquid", "vapor"]
    assert result.diagnostics["phase_regime"] == "VLE"
    assert "phase_set_history" not in result.diagnostics
    assert result.diagnostics["post_split_stable"] is True

    liquid = ct.Mixture.from_database(
        list(names), list(result.phases["liquid"].composition.fractions), normalize=True
    )
    liquid_liquid = ct.stability_tp(
        liquid, temperature_K=temperature_K, pressure_Pa=PRESSURE_PA, activity_model=model
    )
    assert liquid_liquid.status == "stable", liquid_liquid.tpd_min
    assert float(liquid_liquid.tpd_min) >= -1e-8


def test_a_feed_in_the_liquid_liquid_region_returns_two_liquids(system) -> None:
    """Outside the triangle past the ``x^I - x^II`` edge: two liquids, no vapor.

    Verified independently: both converged liquids are below their bubble
    point, ``sum_i x_i gamma_i Psat_i / P < 1``, computed here from the model
    equations, so neither can boil and the vapor really is absent.
    """
    names, model, ln_gamma, psat = system
    temperature_K = 364.0
    x_i, x_ii, y, _ = _tie_triangle(temperature_K, TRIANGLE_SEEDS[temperature_K], ln_gamma, psat)
    vertices = np.column_stack([x_i, x_ii, y])

    z = vertices @ np.array([0.55, 0.60, -0.15])
    z = z / float(np.sum(z))
    weights = np.linalg.solve(vertices, z)
    assert weights[2] < 0.0, weights

    result = _flash(names, model, z, temperature_K)
    assert sorted(result.phase_names()) == ["liquid1", "liquid2"]
    assert result.diagnostics["phase_regime"] == "LLE"
    assert result.vapor_fraction is None
    assert result.diagnostics["post_split_stable"] is True

    for name in result.phase_names():
        x = np.array(result.phases[name].composition.fractions)
        bubble = float(np.sum(x * np.exp(ln_gamma(x)) * psat(temperature_K) / PRESSURE_PA))
        assert bubble < 1.0, (name, bubble)


def test_the_water_rich_corner_is_a_single_liquid(system) -> None:
    """Far from both the miscibility gap and the bubble surface."""
    names, model, ln_gamma, psat = system
    temperature_K = 364.0
    z = np.array([0.01, 0.005, 0.985])

    result = _flash(names, model, z, temperature_K)
    assert result.phase_names() == ["liquid"]
    assert result.diagnostics["phase_regime"] == "single-phase"
    assert result.diagnostics["feed_branch"] == "liquid"
    bubble = float(np.sum(z * np.exp(ln_gamma(z)) * psat(temperature_K) / PRESSURE_PA))
    assert bubble < 1.0, bubble


def test_a_superheated_feed_is_a_single_vapor(system) -> None:
    """At 380 K the feed is above its dew point, checked here independently."""
    names, model, ln_gamma, psat = system
    temperature_K = 380.0
    z = np.array([0.20, 0.15, 0.65])

    result = _flash(names, model, z, temperature_K)
    assert result.phase_names() == ["vapor"]
    assert result.vapor_fraction == 1.0
    assert result.diagnostics["feed_branch"] == "vapor"

    # Dew point: sum_i z_i P / (gamma_i(x) Psat_i) = 1 with x the incipient
    # liquid. Below 1 means the feed is above its dew point.
    x = z.copy()
    for _ in range(500):
        x = z * PRESSURE_PA / (np.exp(ln_gamma(x)) * psat(temperature_K))
        x = x / float(np.sum(x))
    ratio = float(np.sum(z * PRESSURE_PA / (np.exp(ln_gamma(x)) * psat(temperature_K))))
    assert ratio < 1.0, ratio


def test_the_thin_tie_triangle_at_363_k_is_found(system) -> None:
    """Case V-2 amended: the 363 K near-plait miss is fixed by ADR-0012.

    At 363 K the tie-triangle's two liquid vertices are close together (the
    system is near its plait point there). For a feed weighted towards those
    two vertices, *every* trial of the deterministic stability set used to
    collapse onto the trivial solution and `flash_tp` returned a single liquid,
    which is wrong for this model. The failure was in the stability test, not
    in the phase-addition search: the search was never entered.

    The cause was candidate switching *inside* a trial. From the Raoult-vapor
    start the liquid candidate has the lower Gibbs energy at the intermediate
    compositions, so the successive-substitution update used the liquid terms
    and the iterate was dragged onto the liquid surface and onto the trivial
    solution - even though the tangent-plane distance at the equilibrium vapor
    is -9.92e-03. ADR-0012 pins each trial to one candidate surface, and the
    vapor-surface trial then reaches its stationary point in a single
    substitution (the ideal-gas term is zero, so ``ln W_i = d_i``).

    The whole verdict map around this feed is validation Case V-5,
    `tests/validation/test_vlle_verdict_map.py`.
    """
    names, model, ln_gamma, psat = system
    temperature_K = 363.0
    x_i, x_ii, y, _ = _tie_triangle(temperature_K, TRIANGLE_SEEDS[temperature_K], ln_gamma, psat)
    vertices = np.column_stack([x_i, x_ii, y])
    # The two liquid vertices differ by only 0.053 in x_1 at this temperature.
    assert float(np.max(np.abs(x_i - x_ii))) < 0.09

    weights = np.array([0.5, 0.3, 0.2])
    z = vertices @ weights
    z = z / float(np.sum(z))
    assert np.all(np.linalg.solve(vertices, z) > 0.0)

    stability = ct.stability_tp(
        ct.Mixture.from_database(list(names), [float(value) for value in z], normalize=True),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        vapor="ideal",
    )
    assert stability.status == "unstable"
    assert stability.feed_branch == "liquid"
    assert stability.phase_branch == "vapor"
    assert stability.diagnostics["minimizing_trial"] == "raoult-vapor"
    assert stability.diagnostics["minimizing_trial_surface"] == "vapor"
    assert abs(stability.tpd_min - (-0.011680)) < 1e-5, stability.tpd_min

    result = _flash(names, model, z, temperature_K)
    assert sorted(result.phase_names()) == ["liquid1", "liquid2", "vapor"]
    assert result.diagnostics["phase_regime"] == "VLLE"
    for name, reference, weight in (
        ("liquid1", x_i, weights[0]),
        ("liquid2", x_ii, weights[1]),
        ("vapor", y, weights[2]),
    ):
        composition = np.array(result.phases[name].composition.fractions, dtype=float)
        assert float(np.max(np.abs(composition - reference))) < 1e-8, name
        assert abs(result.phase_fractions[name] - float(weight)) < 1e-8, name


# --------------------------------------------------------------------------
# The external cross-check that is not available, checked rather than assumed
# --------------------------------------------------------------------------


def test_thermo_cannot_hold_two_distinct_excess_gibbs_liquids(system) -> None:
    """`thermo` 0.6.0 gives no three-phase reference here, and this records why.

    Validation Case L-2 found that `thermo`'s `FlashVLN` collapses two
    `GibbsExcessLiquid` phases built on one excess-Gibbs model into one
    (`unique_liquid_count == 1`). That was measured for a *binary* liquid-liquid
    split; before claiming "no external reference exists" for a three-phase
    state it is checked again here, directly, with two distinct
    `GibbsExcessLiquid` objects.

    Measured with `thermo` 0.6.0, identical NRTL parameters and chemthermo's own
    Antoine records pushed into `thermo` (so `Psat` is shared):

    - `liquids=[liquid]`: `unique_liquid_count == 1`, and the flash returns a
      **two-phase** vapor-liquid answer at a feed that is inside the
      tie-triangle - vapor (0.21027032, 0.09978016, 0.68994952), liquid
      (0.09928882, 0.07716435, 0.82354683), betas (0.31446, 0.68554). That
      liquid is neither conjugate liquid; it is the single-liquid answer.
    - `liquids=[liquid, liquid]`: still `unique_liquid_count == 1`, and the
      flash raises `TypeError: 'NoneType' object is not subscriptable`.

    So there is no external three-phase reference for this system, and none is
    claimed anywhere in this module. The test also *adjudicates*: chemthermo's
    three-phase answer has the lower Gibbs energy of the two, computed here from
    this module's own equations.

    The assertions are deliberately loose - they pin the *failure*, so that a
    future `thermo` which fixes it is noticed rather than silently ignored.
    """
    pytest.importorskip("thermo")
    from thermo import NRTL as ThermoNRTL  # noqa: PLC0415
    from thermo import (  # noqa: PLC0415 - optional dependency
        ChemicalConstantsPackage,
        FlashVLN,
        GibbsExcessLiquid,
        IdealGas,
    )

    names, model, ln_gamma, psat = system
    temperature_K = 364.0
    x_i, x_ii, y, _ = _tie_triangle(temperature_K, TRIANGLE_SEEDS[temperature_K], ln_gamma, psat)
    z = np.column_stack([x_i, x_ii, y]) @ np.array([1 / 3, 1 / 3, 1 / 3])
    z = z / float(np.sum(z))

    constants, correlations = ChemicalConstantsPackage.from_IDs(["1-propanol", "butanol", "water"])
    # Share Psat: chemthermo stores ln(P/bar), thermo's Antoine is in Pa.
    for index, name in enumerate(names):
        record = ct.Component.from_database(name).antoine
        assert record is not None
        correlations.VaporPressures[index].add_correlation(
            name="chemthermo",
            model="Antoine",
            Tmin=record.Tmin_K,
            Tmax=record.Tmax_K,
            A=record.A + float(np.log(1e5)),
            B=record.B,
            C=record.C,
            base=float(np.e),
        )

    tau = [[0.0, -0.61259, -0.07149], [0.7164, 0.0, 0.90047], [2.7425, 3.51307, 0.0]]
    alpha = [[0.0, 0.3, 0.3], [0.3, 0.0, 0.48], [0.3, 0.48, 0.0]]
    excess = ThermoNRTL(T=temperature_K, xs=[1 / 3] * 3, tau_as=tau, alpha_cs=alpha)

    def liquid() -> object:
        return GibbsExcessLiquid(
            VaporPressures=correlations.VaporPressures,
            GibbsExcessModel=excess,
            equilibrium_basis=None,
            eos_pure_instances=None,
            use_Poynting=False,
            use_phis_sat=False,
            HeatCapacityGases=correlations.HeatCapacityGases,
            T=temperature_K,
            P=PRESSURE_PA,
            zs=[1 / 3] * 3,
        )

    gas = IdealGas(
        HeatCapacityGases=correlations.HeatCapacityGases,
        T=temperature_K,
        P=PRESSURE_PA,
        zs=[1 / 3] * 3,
    )

    single = FlashVLN(constants, correlations, liquids=[liquid()], gas=gas)
    assert single.unique_liquid_count == 1
    external = single.flash(T=temperature_K, P=PRESSURE_PA, zs=[float(v) for v in z])
    assert external.phase_count == 2, (
        "thermo now returns more than two phases here; Case V-1's 'no external "
        "reference' note needs revisiting"
    )

    doubled = FlashVLN(constants, correlations, liquids=[liquid(), liquid()], gas=gas)
    assert doubled.unique_liquid_count == 1
    with pytest.raises((TypeError, ValueError, AttributeError)):
        doubled.flash(T=temperature_K, P=PRESSURE_PA, zs=[float(v) for v in z])

    # Adjudication: chemthermo's three-phase answer is the lower-Gibbs one.
    ours = _flash(names, model, z, temperature_K)
    assert len(ours.phase_names()) == 3

    def energy_of(betas, compositions, vapors) -> float:
        return float(
            sum(
                float(beta)
                * _reduced_g(np.array(composition), is_vapor, temperature_K, ln_gamma, psat)
                for beta, composition, is_vapor in zip(betas, compositions, vapors)
            )
        )

    theirs_g = energy_of(
        external.betas,
        [phase.zs for phase in external.phases],
        [phase.__class__.__name__ == "IdealGas" for phase in external.phases],
    )
    ours_g = energy_of(
        [ours.phase_fractions[name] for name in ours.phase_names()],
        [ours.phases[name].composition.fractions for name in ours.phase_names()],
        [name == "vapor" for name in ours.phase_names()],
    )
    # Achieved with thermo 0.6.0: ours -0.693912758, theirs -0.693543321.
    assert ours_g < theirs_g, (ours_g, theirs_g)

    # And the intermediate *is* externally confirmed: thermo's two-phase answer
    # has the same Gibbs energy as the two-phase candidate chemthermo converged
    # before adding the third phase, to 2.4e-09. So the disagreement is about
    # the phase count, not about the two-phase thermodynamics.
    two = _flash(names, model, z, temperature_K, ct.FlashSettings(post_split_stability=False))
    ours_two_g = energy_of(
        [two.phase_fractions[name] for name in two.phase_names()],
        [two.phases[name].composition.fractions for name in two.phase_names()],
        [name == "vapor" for name in two.phase_names()],
    )
    assert ours_two_g == pytest.approx(theirs_g, abs=1e-7)
