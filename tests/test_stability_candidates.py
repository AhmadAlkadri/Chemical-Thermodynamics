"""Phase candidates in the tangent-plane evaluator, and `stability_tp(vapor=...)`.

Three things are under test here (ADR-0010):

1. the Antoine pure-liquid reference fugacity itself - its equation form, its
   units, and the refusal to extrapolate outside the fit's stated range;
2. the generalized evaluator - one minimum-Gibbs selection rule serving cubic
   roots, a single activity liquid, and the modified-Raoult pair, with the two
   pre-existing families behaving exactly as before;
3. the public `stability_tp(..., activity_model=..., vapor="ideal")` surface:
   the verdict and the *candidate label* of the feed and of the incipient
   phase, permutation invariance, determinism, and the error contract.

See validation Cases R-2 and R-3.
"""

from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.models._antoine import antoine_saturation_pressures, antoine_temperature_range
from chemthermo.stability._evaluator import (
    _ActivityTangentPlane,
    _EOSTangentPlane,
    _ModifiedRaoultTangentPlane,
    _select_surface,
    modified_raoult_candidates,
)

PRESSURE_PA = 101325.0
_PROPANOL, _BUTANOL, _WATER = 0, 1, 2


def _binary(payload: dict[str, Any], first: int, second: int) -> tuple[tuple[str, str], ct.NRTL]:
    names = [str(entry["chemthermo_name"]) for entry in payload["components"]]
    tau = payload["tau"]
    alpha = payload["alpha"]
    pair = (names[first], names[second])
    return pair, ct.NRTL(
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


@pytest.fixture(scope="module")
def propanol_water(tessier2000_payload: dict[str, Any]) -> tuple[tuple[str, str], ct.NRTL]:
    return _binary(tessier2000_payload, _PROPANOL, _WATER)


@pytest.fixture(scope="module")
def butanol_water(tessier2000_payload: dict[str, Any]) -> tuple[tuple[str, str], ct.NRTL]:
    return _binary(tessier2000_payload, _BUTANOL, _WATER)


def _mixture(names: Sequence[str], z: Sequence[float]) -> ct.Mixture:
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


def _stability(
    names: Sequence[str],
    model: ct.NRTL,
    z: Sequence[float],
    temperature_K: float,
    **kwargs: Any,
) -> ct.StabilityResult:
    return ct.stability_tp(
        _mixture(names, z),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        vapor="ideal",
        **kwargs,
    )


# --------------------------------------------------------------------------
# The Antoine reference fugacity
# --------------------------------------------------------------------------


def test_antoine_is_the_base_e_bar_form_and_boils_water_at_one_atmosphere() -> None:
    """`ln(P/bar) = A - B/(T + C)` - base e, not base 10, and T in K.

    The packaged water record (A = 11.6834, B = 3816.44, C = -46.13) is a
    Koretsky (2012) Appendix A.1 entry. Evaluated at the normal boiling point
    it must return one atmosphere; a base-10 reading, a Celsius reading, or a
    sign error on C would all fail this by orders of magnitude.
    """
    mixture = _mixture(("Water",), (1.0,))
    psat = antoine_saturation_pressures(mixture, 373.15)
    assert psat.shape == (1,)
    assert abs(float(psat[0]) - 101325.0) / 101325.0 < 2e-4

    antoine = mixture.components[0].antoine
    assert antoine is not None
    assert antoine.units == "bar"
    expected = math.exp(antoine.A - antoine.B / (373.15 + antoine.C)) * 1.0e5
    assert float(psat[0]) == expected


def test_antoine_ranges_are_the_intersection_of_the_component_windows() -> None:
    # 1-Propanol [285, 400] K and Water [284, 441] K -> [285, 400] K.
    mixture = _mixture(("1-Propanol", "Water"), (0.5, 0.5))
    assert antoine_temperature_range(mixture) == (285.0, 400.0)


@pytest.mark.parametrize("temperature_K", [280.0, 420.0])
def test_out_of_range_temperatures_raise_rather_than_extrapolate(
    propanol_water, temperature_K: float
) -> None:
    names, model = propanol_water
    with pytest.raises(ct.InputRangeError, match="Antoine correlation"):
        _stability(names, model, (0.5, 0.5), temperature_K)


def test_the_activity_only_path_has_no_antoine_requirement(butanol_water) -> None:
    """`vapor="none"` never touches Antoine, so its temperature range is free."""
    names, model = butanol_water
    result = ct.stability_tp(
        _mixture(names, (0.20, 0.80)),
        temperature_K=250.0,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
    )
    assert result.status == "unstable"


# --------------------------------------------------------------------------
# The `vapor` keyword contract
# --------------------------------------------------------------------------


def test_vapor_ideal_is_rejected_with_an_eos() -> None:
    with pytest.raises(ct.ModelError, match="only valid with 'activity_model'"):
        ct.stability_tp(
            _mixture(("Methane", "Ethane"), (0.5, 0.5)),
            temperature_K=240.0,
            pressure_Pa=3.0e6,
            eos=ct.PengRobinsonEOS(),
            vapor="ideal",
        )


def test_an_unknown_vapor_value_is_rejected(propanol_water) -> None:
    names, model = propanol_water
    with pytest.raises(ct.ModelError, match="vapor must be one of"):
        ct.stability_tp(
            _mixture(names, (0.5, 0.5)),
            temperature_K=361.0,
            pressure_Pa=PRESSURE_PA,
            activity_model=model,
            vapor="peng-robinson",  # type: ignore[arg-type]
        )


def test_the_default_leaves_the_activity_only_path_untouched(butanol_water) -> None:
    names, model = butanol_water
    explicit = ct.stability_tp(
        _mixture(names, (0.20, 0.80)),
        temperature_K=298.15,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        vapor="none",
    )
    implicit = ct.stability_tp(
        _mixture(names, (0.20, 0.80)),
        temperature_K=298.15,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
    )
    assert explicit.tpd_min == implicit.tpd_min
    assert explicit.feed_branch is None and implicit.feed_branch is None
    assert explicit.diagnostics["model_family"] == "activity"


# --------------------------------------------------------------------------
# The evaluator generalization
# --------------------------------------------------------------------------


def test_one_candidate_reports_no_label_and_two_report_the_winner(
    butanol_water, propanol_water
) -> None:
    """The single-candidate case has nothing to select, so its label is None."""
    ll_names, ll_model = butanol_water
    activity = _ActivityTangentPlane(
        ll_model, mixture=_mixture(ll_names, (0.5, 0.5)), temperature=298.15
    )
    terms, label = activity.ln_fugacity_terms(np.array([0.5, 0.5]))
    assert label is None
    assert terms.shape == (2,)

    names, model = propanol_water
    raoult = _ModifiedRaoultTangentPlane(
        model,
        mixture=_mixture(names, (0.5, 0.5)),
        temperature=361.0,
        pressure=PRESSURE_PA,
    )
    _terms, raoult_label = raoult.ln_fugacity_terms(np.array([0.5, 0.5]))
    assert raoult_label in ("liquid", "vapor")

    eos = _EOSTangentPlane(
        ct.PengRobinsonEOS(),
        mixture=_mixture(("Methane", "Ethane"), (0.5, 0.5)),
        temperature=240.0,
        pressure=3.0e6,
    )
    _eos_terms, eos_label = eos.ln_fugacity_terms(np.array([0.5, 0.5]))
    assert eos_label in ("vapor", "liquid")


def test_the_selected_candidate_is_the_one_with_the_lower_gibbs_energy(propanol_water) -> None:
    """`sum_i w_i term_i` is the only candidate-dependent part of `G/RT`.

    Checked against the *full* reduced Gibbs energy
    `sum_i w_i (ln w_i + term_i)`, so the test does not assume the
    ideal-mixing part cancels - it verifies that it does.
    """
    names, model = propanol_water
    mixture = _mixture(names, (0.5, 0.5))
    for temperature in (330.0, 355.0, 361.0, 365.0, 380.0):
        evaluator = _ModifiedRaoultTangentPlane(
            model, mixture=mixture, temperature=temperature, pressure=PRESSURE_PA
        )
        liquid, vapor = modified_raoult_candidates(
            model, mixture=mixture, temperature=temperature, pressure=PRESSURE_PA
        )
        for w1 in (0.05, 0.3, 0.5, 0.8, 0.95):
            w = np.array([w1, 1.0 - w1])
            terms, label = evaluator.ln_fugacity_terms(w)
            energies = {
                candidate.label: float(np.sum(w * (np.log(w) + candidate.ln_fugacity_terms(w))))
                for candidate in (liquid, vapor)
            }
            assert label == min(energies, key=lambda key: energies[key])
            chosen = liquid if label == "liquid" else vapor
            assert np.array_equal(terms, chosen.ln_fugacity_terms(w))


def test_the_ideal_vapor_candidate_contributes_nothing(propanol_water) -> None:
    names, model = propanol_water
    _liquid, vapor = modified_raoult_candidates(
        model,
        mixture=_mixture(names, (0.5, 0.5)),
        temperature=361.0,
        pressure=PRESSURE_PA,
    )
    w = np.array([0.4, 0.6])
    assert np.array_equal(vapor.ln_fugacity_terms(w), np.zeros(2))


def test_the_liquid_candidate_is_gamma_times_the_raoult_k_value(propanol_water) -> None:
    names, model = propanol_water
    mixture = _mixture(names, (0.5, 0.5))
    liquid, _vapor = modified_raoult_candidates(
        model, mixture=mixture, temperature=361.0, pressure=PRESSURE_PA
    )
    w = np.array([0.4, 0.6])
    gamma = np.asarray(
        model.activity_coefficients(mixture=mixture, temperature_K=361.0, composition=w.tolist()),
        dtype=float,
    )
    psat = antoine_saturation_pressures(mixture, 361.0)
    expected = np.log(gamma) + np.log(psat / PRESSURE_PA)
    assert float(np.max(np.abs(liquid.ln_fugacity_terms(w) - expected))) < 1e-15


def test_a_liquid_liquid_tangent_plane_is_unchanged_by_the_reference_offset(
    butanol_water,
) -> None:
    """The offset `ln(Psat_i / P)` cancels when feed and trial are both liquids.

    This is why the modified-Raoult path reproduces `gamma-gamma` tie-lines
    exactly (Case R-1): at a temperature where the vapor candidate never wins,
    the two evaluators define the *same* tangent-plane distance.
    """
    names, model = butanol_water
    with_vapor = _stability(names, model, (0.20, 0.80), 330.0)
    without_vapor = ct.stability_tp(
        _mixture(names, (0.20, 0.80)),
        temperature_K=330.0,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
    )
    assert with_vapor.status == "unstable"
    assert with_vapor.feed_branch == "liquid"
    assert with_vapor.phase_branch == "liquid"
    assert abs(with_vapor.tpd_min - without_vapor.tpd_min) < 1e-14
    assert with_vapor.trial_composition is not None
    assert without_vapor.trial_composition is not None
    # The two evaluators reach the same stationary point from different trial
    # sets (the modified-Raoult one also runs the two Raoult estimates), so the
    # converged compositions agree to the solver's own tolerance, not to
    # round-off. Measured: 5.1e-12 against `StabilitySettings.tol = 1e-10`.
    assert (
        max(
            abs(a - b)
            for a, b in zip(with_vapor.trial_composition, without_vapor.trial_composition)
        )
        < 1e-9
    )


def test_the_stationary_point_identity_holds_for_the_candidate_pair(propanol_water) -> None:
    """`tpd = -ln(sum_i W_i)` and `tm* = 1 - sum_i W_i` still hold (equations 6-7)."""
    names, model = propanol_water
    result = _stability(names, model, (0.5, 0.5), 361.0)
    assert result.status == "unstable"
    sum_w = float(result.diagnostics["sum_W"])
    assert abs(-math.log(sum_w) - result.tpd_min) < 1e-12
    assert abs(float(result.diagnostics["tm_at_stationary_point"]) - (1.0 - sum_w)) < 1e-15


# --------------------------------------------------------------------------
# The public `stability_tp(..., vapor="ideal")` surface
# --------------------------------------------------------------------------


def test_a_superheated_feed_is_stable_and_labelled_a_vapor(propanol_water) -> None:
    names, model = propanol_water
    result = _stability(names, model, (0.5, 0.5), 380.0)
    assert result.status == "stable"
    assert result.stable is True
    assert result.feed_branch == "vapor"
    assert result.diagnostics["model_family"] == "modified-raoult"
    assert result.diagnostics["pressure_dependent"] is True


def test_a_subcooled_feed_is_stable_and_labelled_a_liquid(propanol_water) -> None:
    names, model = propanol_water
    result = _stability(names, model, (0.5, 0.5), 330.0)
    assert result.status == "stable"
    assert result.feed_branch == "liquid"


def test_a_feed_inside_the_band_is_unstable_with_a_vapor_incipient_phase(
    propanol_water,
) -> None:
    names, model = propanol_water
    result = _stability(names, model, (0.5, 0.5), 361.0)
    assert result.status == "unstable"
    assert result.feed_branch == "liquid"
    assert result.phase_branch == "vapor"
    assert result.trial_composition is not None
    # The incipient vapor is richer in the more volatile 1-propanol than the
    # equilibrium liquid, and K = w / z is reported feed -> incipient.
    assert result.k_values is not None
    assert result.k_values[0] < 1.0 < result.k_values[1]


def test_the_antoine_window_is_reported_in_diagnostics(propanol_water) -> None:
    names, model = propanol_water
    diagnostics = _stability(names, model, (0.5, 0.5), 361.0).diagnostics
    assert diagnostics["antoine_valid_Tmin_K"] == 285.0
    assert diagnostics["antoine_valid_Tmax_K"] == 400.0


def test_results_are_deterministic(propanol_water) -> None:
    names, model = propanol_water
    first = _stability(names, model, (0.4, 0.6), 361.0)
    second = _stability(names, model, (0.4, 0.6), 361.0)
    assert first.tpd_min == second.tpd_min
    assert first.trial_composition == second.trial_composition
    assert first.feed_branch == second.feed_branch


@pytest.mark.parametrize(
    "z, temperature_K", [((0.4, 0.6), 361.0), ((0.20, 0.80), 330.0), ((0.5, 0.5), 380.0)]
)
def test_results_are_invariant_under_component_reordering(
    tessier2000_payload: dict[str, Any], z: tuple[float, float], temperature_K: float
) -> None:
    """Reordering the components must permute the answer, not change it.

    Compared trial by trial: each start is named by its label in both orders,
    and every one must agree to 1e-12 (measured <= 7.1e-15). The *reported*
    minimizer is compared to 1e-12 too when both orders report the same trial.
    At ``(0.4, 0.6)`` three trials reach one stationary point with ``tpd``
    equal to 1e-17, and which of them is reported is decided by the last bit
    of the summation order - it differs between machines (ADR-0032). Then the
    tie itself is asserted and the two reported compositions only have to agree
    to the stationarity tolerance, because one of the tied trials stopped at
    residual 3.3e-11 against ``StabilitySettings.tol = 1e-10``.
    """
    forward_names, forward_model = _binary(tessier2000_payload, _PROPANOL, _WATER)
    reverse_names, reverse_model = _binary(tessier2000_payload, _WATER, _PROPANOL)

    forward = _stability(forward_names, forward_model, z, temperature_K)
    reverse = _stability(reverse_names, reverse_model, tuple(reversed(z)), temperature_K)

    assert forward.status == reverse.status
    assert forward.feed_branch == reverse.feed_branch
    assert abs(forward.tpd_min - reverse.tpd_min) < 1e-12

    reverse_trials = {trial.label: trial for trial in reverse.trials}
    # Pure-component starts follow component order, so compare label sets.
    assert sorted(trial.label for trial in forward.trials) == sorted(reverse_trials)
    for trial in forward.trials:
        mirrored = reverse_trials[trial.label]
        assert (trial.converged, trial.trivial) == (mirrored.converged, mirrored.trivial)
        assert abs(trial.tpd - mirrored.tpd) < 1e-12, trial.label
        assert (trial.composition is None) == (mirrored.composition is None)
        if trial.composition is not None and mirrored.composition is not None:
            assert (
                max(abs(a - b) for a, b in zip(trial.composition, reversed(mirrored.composition)))
                < 1e-12
            ), trial.label

    if forward.trial_composition is not None:
        assert reverse.trial_composition is not None
        reported_move = max(
            abs(a - b)
            for a, b in zip(forward.trial_composition, reversed(reverse.trial_composition))
        )
        forward_label = forward.diagnostics["minimizing_trial"]
        reverse_label = reverse.diagnostics["minimizing_trial"]
        if forward_label == reverse_label:
            assert reported_move < 1e-12
        else:
            tied = reverse_trials[str(forward_label)]
            assert abs(tied.tpd - reverse.tpd_min) < 1e-12
            assert reported_move < 1e-9


def test_a_tpd_tie_is_reported_from_the_better_converged_trial(
    tessier2000_payload: dict[str, Any],
) -> None:
    """ADR-0035: at (0.4, 0.6), 361 K, three trials tie on one stationary point.

    One of them stopped at residual 3.3e-11 (inside ``tol = 1e-10``), the
    others at ~1e-16. Before ADR-0035 whichever had the lowest ``tpd`` in the
    last bit was reported, so the reported composition carried 8.4e-11 of
    iteration error in one component order and not the other. Now a tied
    trial converged 1000x better is reported instead, and both orders give
    the same composition to 1e-12 (measured 4.4e-16).
    """
    forward_names, forward_model = _binary(tessier2000_payload, _PROPANOL, _WATER)
    reverse_names, reverse_model = _binary(tessier2000_payload, _WATER, _PROPANOL)
    forward = _stability(forward_names, forward_model, (0.4, 0.6), 361.0)
    reverse = _stability(reverse_names, reverse_model, (0.6, 0.4), 361.0)

    for result in (forward, reverse):
        assert result.status == "unstable"
        assert float(result.diagnostics["minimizing_trial_residual"]) < 1e-13
    assert forward.trial_composition is not None and reverse.trial_composition is not None
    move = max(
        abs(a - b) for a, b in zip(forward.trial_composition, reversed(reverse.trial_composition))
    )
    assert move < 1e-12
    # The key appears only where the rule acted, and names the rule.
    for result in (forward, reverse):
        assert result.diagnostics.get("minimizing_trial_tie_break", "residual") == "residual"


def test_the_water_butanol_feed_sees_a_liquid_incipient_phase_from_a_vapor_feed(
    butanol_water,
) -> None:
    """Just at the heteroazeotrope the feed is a vapor whose incipient phase is a liquid."""
    names, model = butanol_water
    result = _stability(names, model, (0.20, 0.80), 366.2137741)
    assert result.status == "unstable"
    assert result.feed_branch == "vapor"
    assert result.phase_branch == "liquid"


# --------------------------------------------------------------------------
# Trial candidate surfaces (ADR-0012)
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def near_plait_ternary(
    tessier2000_payload: dict[str, Any],
) -> tuple[list[str], ct.NRTL, np.ndarray]:
    """The 363 K feed of validation Case V-2, and its model.

    `z` is the 0.5 / 0.3 / 0.2 barycentric mix of the tie-triangle vertices of
    validation Case V-1, pinned here rather than solved (it is solved in
    `tests/validation/test_vlle_water_propanol_butanol.py`).
    """
    names = [str(entry["chemthermo_name"]) for entry in tessier2000_payload["components"]]
    tau = tessier2000_payload["tau"]
    alpha = tessier2000_payload["alpha"]
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
                for i in range(3)
                for j in range(i + 1, 3)
            ]
        )
    )
    vertices = np.column_stack(
        [
            np.array([0.10282779, 0.03539032, 0.86178189]),
            np.array([0.15630287, 0.06422717, 0.77946996]),
            np.array([0.28312929, 0.06046985, 0.65640086]),
        ]
    )
    z = vertices @ np.array([0.5, 0.3, 0.2])
    return names, model, z / float(np.sum(z))


def test_the_modified_raoult_trial_set_names_one_surface_per_trial(propanol_water) -> None:
    """One vapor-surface trial, one liquid-surface trial per liquid estimate."""
    names, model = propanol_water
    evaluator = _ModifiedRaoultTangentPlane(
        model, mixture=_mixture(names, (0.5, 0.5)), temperature=361.0, pressure=PRESSURE_PA
    )
    z = np.array([0.5, 0.5])
    estimates = evaluator.initial_estimates(z, z > 0.0)
    assert [(estimate.label, estimate.surface) for estimate in estimates] == [
        ("raoult-vapor", "vapor"),
        ("raoult-liquid", "liquid"),
        ("pure-1-Propanol", "liquid"),
        ("pure-Water", "liquid"),
    ]


def test_the_activity_only_family_names_no_surface(butanol_water) -> None:
    """The bit-identity guard at its source: no surface means the old iteration.

    An estimate that carries no surface makes the solver call
    `ln_fugacity_terms` exactly as it did before ADR-0012, so the activity-only
    path is untouched. The whole-result version of this statement is
    `tests/test_flash_refactor_bit_identity.py`.

    The equation-of-state family *did* take surfaces, in ADR-0021; its
    counterpart to this test is
    `tests/test_stability_eos_surfaces.py::test_every_multicomponent_eos_trial_names_the_root_its_start_estimates`.
    """
    ll_names, ll_model = butanol_water
    activity = _ActivityTangentPlane(
        ll_model, mixture=_mixture(ll_names, (0.2, 0.8)), temperature=298.15
    )
    z = np.array([0.2, 0.8])
    assert all(estimate.surface is None for estimate in activity.initial_estimates(z, z > 0.0))

    result = ct.stability_tp(
        _mixture(ll_names, (0.2, 0.8)),
        temperature_K=298.15,
        pressure_Pa=PRESSURE_PA,
        activity_model=ll_model,
    )
    assert all(trial.surface is None for trial in result.trials)
    assert all(trial.surface_fallback is False for trial in result.trials)
    assert "trial_surfaces" not in result.diagnostics


def test_the_vapor_surface_is_reached_in_one_substitution(near_plait_ternary) -> None:
    """`ln W_i = d_i` on the ideal-gas surface: a constant map, so one step.

    This is why the vapor surface needs exactly one trial and why its starting
    point does not matter (see `_ModifiedRaoultTangentPlane.initial_estimates`).
    """
    names, model, z = near_plait_ternary
    result = ct.stability_tp(
        _mixture(names, z),
        temperature_K=363.0,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        vapor="ideal",
    )
    vapor_trials = [trial for trial in result.trials if trial.surface == "vapor"]
    assert len(vapor_trials) == 1
    trial = vapor_trials[0]
    assert trial.label == "raoult-vapor"
    assert trial.converged is True
    assert trial.ssi_iterations == 2
    assert trial.second_order_iterations == 0
    assert trial.residual == 0.0


def test_the_near_plait_feed_is_unstable_on_the_vapor_surface(near_plait_ternary) -> None:
    """Validation Case V-2, repaired: the numbers of the fixed trial."""
    names, model, z = near_plait_ternary
    result = ct.stability_tp(
        _mixture(names, z),
        temperature_K=363.0,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        vapor="ideal",
    )
    assert result.status == "unstable"
    assert result.feed_branch == "liquid"
    assert result.phase_branch == "vapor"
    assert result.tpd_min == pytest.approx(-0.011680295426, abs=1e-6)
    assert result.diagnostics["minimizing_trial"] == "raoult-vapor"
    assert result.diagnostics["minimizing_trial_surface"] == "vapor"
    assert result.diagnostics["trial_surfaces"] == "vapor:1,liquid:4"
    assert result.diagnostics["surface_fallback_trial_count"] == 0
    assert result.trial_composition is not None
    assert np.allclose(
        np.array(result.trial_composition),
        [0.303242, 0.050098, 0.646659],
        rtol=0.0,
        atol=1e-6,
    )


def test_without_a_fixed_surface_the_same_feed_is_reported_stable(
    near_plait_ternary, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The defect ADR-0012 fixes, reproduced by putting the old rule back.

    Stripping the surface labels from the trial set restores the pre-ADR-0012
    iteration - the lowest-Gibbs candidate re-selected at every iterate - and
    nothing else. Every trial then collapses onto the trivial solution and the
    feed is reported *stable*, although the tangent-plane distance at the
    equilibrium vapor is -9.92e-03. That is the Case V-2 miss, and it is caused
    by the candidate switching, not by the initial estimates: the estimates
    here are exactly the ones the repaired trial set uses.
    """
    names, model, z = near_plait_ternary
    original = _ModifiedRaoultTangentPlane.initial_estimates

    def without_surfaces(self, z_local: np.ndarray, active: np.ndarray):
        return [estimate._replace(surface=None) for estimate in original(self, z_local, active)]

    monkeypatch.setattr(_ModifiedRaoultTangentPlane, "initial_estimates", without_surfaces)
    result = ct.stability_tp(
        _mixture(names, z),
        temperature_K=363.0,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        vapor="ideal",
    )
    assert result.status == "stable"
    assert result.tpd_min == 0.0
    assert all(trial.trivial for trial in result.trials if trial.converged)
    assert all(trial.surface is None for trial in result.trials)


def test_the_reported_distance_is_the_lowest_gibbs_one_at_the_converged_point(
    near_plait_ternary,
) -> None:
    """`tpd` uses the min-Gibbs candidate; `surface` records what was iterated on.

    The distance to the tangent plane is the distance from the *lower envelope*
    of the candidates, so it is never taken from the pinned surface alone.
    """
    names, model, z = near_plait_ternary
    mixture = _mixture(names, z)
    evaluator = _ModifiedRaoultTangentPlane(
        model, mixture=mixture, temperature=363.0, pressure=PRESSURE_PA
    )
    result = ct.stability_tp(
        mixture,
        temperature_K=363.0,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        vapor="ideal",
    )
    feed_terms, _feed_label = evaluator.ln_fugacity_terms(np.asarray(z, dtype=float))
    d = np.log(np.asarray(z, dtype=float)) + feed_terms

    for trial in result.trials:
        if not trial.converged or trial.composition is None:
            continue
        w = np.array(trial.composition, dtype=float)
        terms, label = evaluator.ln_fugacity_terms(w)
        assert trial.phase_branch == label
        assert trial.tpd == pytest.approx(float(np.sum(w * (np.log(w) + terms - d))), abs=1e-14)


def test_an_unknown_surface_is_a_model_error(propanol_water) -> None:
    names, model = propanol_water
    evaluator = _ModifiedRaoultTangentPlane(
        model, mixture=_mixture(names, (0.5, 0.5)), temperature=361.0, pressure=PRESSURE_PA
    )
    with pytest.raises(ct.ModelError, match="Unknown phase-candidate surface"):
        evaluator.ln_terms_on_surface(np.array([0.5, 0.5]), "solid")


def test_an_unavailable_optional_surface_falls_back_and_says_so() -> None:
    """Point 3 of ADR-0012: no surface to walk means the lowest-Gibbs candidate.

    Not reachable through the modified-Raoult pair, whose two candidates are
    both evaluable everywhere, so it is exercised here on two stub candidates.
    """

    class _Missing:
        label = "vapor"
        optional = True

        def ln_fugacity_terms(self, composition: np.ndarray) -> np.ndarray:
            raise ct.ModelError("no root here")

    class _Present:
        label = "liquid"
        optional = True

        def ln_fugacity_terms(self, composition: np.ndarray) -> np.ndarray:
            return np.full(composition.shape, -0.25)

    w = np.array([0.4, 0.6])
    evaluated = _select_surface(
        (_Missing(), _Present()), w, "vapor", failure_message="no candidate"
    )
    assert evaluated.fell_back is True
    assert evaluated.label == "liquid"
    assert np.array_equal(evaluated.terms, np.full(2, -0.25))

    evaluated = _select_surface(
        (_Missing(), _Present()), w, "liquid", failure_message="no candidate"
    )
    assert evaluated.fell_back is False
    assert evaluated.label == "liquid"


def test_the_pure_water_trial_stalls_near_the_plait_point(near_plait_ternary) -> None:
    """A recorded limit, not a numerical failure: a near-singular Jacobian.

    At the 363 K Case V-2 feed the `pure-Water` liquid-surface trial does not
    converge. It is not an overflow and nothing in the iteration is non-finite:
    the iterate walks into the near-plait region, where the stationarity
    Jacobian ``dg/d(ln W)`` has an eigenvalue of about 3.3e-09 (condition number
    ~8.2e+08 at the stalling point w = (0.13762, 0.04163, 0.82075)), so the
    Newton step is dominated by the near-null direction and the line search
    cannot reduce the residual below about 5.4e-04. Successive substitution on
    the same surface *does* converge, to the trivial solution, but only after
    ~1.2e+04 iterations - far beyond `StabilitySettings.max_iter = 300`. The
    `nan` tangent-plane distance is the documented sentinel for a trial that did
    not converge, not a number that blew up.

    The verdict does not depend on this trial: four others converge, and the
    instability is found on the vapor surface.
    """
    names, model, z = near_plait_ternary
    result = ct.stability_tp(
        _mixture(names, z),
        temperature_K=363.0,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        vapor="ideal",
    )
    stalled = [trial for trial in result.trials if trial.label == "pure-Water"]
    assert len(stalled) == 1
    trial = stalled[0]
    assert trial.converged is False
    assert trial.surface == "liquid"
    assert trial.termination_reason == "second_order_no_progress"
    assert math.isnan(trial.tpd)
    assert math.isfinite(trial.residual) and trial.residual < 1e-3
    assert trial.composition is not None
    assert all(math.isfinite(value) for value in trial.composition)
    assert result.status == "unstable"
    assert int(result.diagnostics["converged_trial_count"]) == 4
