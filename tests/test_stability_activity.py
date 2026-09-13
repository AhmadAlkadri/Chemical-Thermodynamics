"""Tangent-plane phase stability driven by an activity-coefficient model.

Equation numbers refer to the module docstring of
``src/chemthermo/stability/tp.py``. For an activity model the only change is
``ln gamma_i`` in place of ``ln phi_i``; every identity below is therefore the
same identity the Peng-Robinson tests assert, re-derived here from
``ln gamma`` so that a shared-machinery regression cannot hide behind the EOS
path.

Parameters come from the cited fixture
``tests/fixtures/nrtl/tessier2000_problem1.json`` (Tessier, Brennecke &
Stadtherr, Chem. Eng. Sci. 55 (2000) 1785, Table 1), never from the packaged
synthetic defaults.
"""

from __future__ import annotations

import math
from typing import Any, Callable, Sequence

import numpy as np
import pytest

import chemthermo as ct

LnGamma = Callable[[np.ndarray], np.ndarray]

TEMPERATURE_K = 298.15  # Immaterial: the fixture gives dimensionless tau.
PRESSURE_PA = 101325.0  # Validated, but an activity model does not use it.

# Near-plait-point feed of Table 2: successive substitution alone cannot solve
# it, which is what makes it the regression guard for the second-order stage.
PLAIT_FEED = (0.148, 0.052, 0.80)

# n-butanol / water pair of the same fixture (pair 2-3 of Table 1), used as a
# self-contained partially miscible binary.
BUTANOL_WATER_TAU_12 = 0.90047
BUTANOL_WATER_TAU_21 = 3.51307
BUTANOL_WATER_ALPHA = 0.48
BUTANOL_WATER_NAMES = ("n-Butanol", "Water")


def _mixture(names: Sequence[str], z: Sequence[float]) -> ct.Mixture:
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


def _butanol_water_model() -> ct.NRTL:
    return ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [
                (
                    BUTANOL_WATER_NAMES[0],
                    BUTANOL_WATER_NAMES[1],
                    BUTANOL_WATER_TAU_12,
                    BUTANOL_WATER_TAU_21,
                    BUTANOL_WATER_ALPHA,
                    BUTANOL_WATER_ALPHA,
                )
            ]
        )
    )


def _butanol_water_ln_gamma(model: ct.NRTL) -> LnGamma:
    mixture = _mixture(BUTANOL_WATER_NAMES, (0.5, 0.5))

    def ln_gamma(x: np.ndarray) -> np.ndarray:
        values = np.asarray(x, dtype=float)
        values = values / float(np.sum(values))
        gamma = model.activity_coefficients(
            mixture=mixture,
            temperature_K=TEMPERATURE_K,
            composition=[float(value) for value in values],
        )
        return np.log(np.asarray(gamma, dtype=float))

    return ln_gamma


def _d_vector(ln_gamma: LnGamma, z: np.ndarray) -> np.ndarray:
    """Feed tangent-plane intercept d_i = ln z_i + ln gamma_i(z), equation (1)."""
    return np.log(z) + ln_gamma(z)


def _tpd(ln_gamma: LnGamma, w: np.ndarray, d: np.ndarray) -> float:
    """Reduced tangent-plane distance, equation (2), from first principles."""
    w = np.asarray(w, dtype=float)
    w = w / float(np.sum(w))
    return float(np.sum(w * (np.log(w) + ln_gamma(w) - d)))


def _tm(ln_gamma: LnGamma, w_capital: np.ndarray, d: np.ndarray) -> float:
    """Modified (unnormalized) tangent-plane function, equation (3)."""
    w_capital = np.asarray(w_capital, dtype=float)
    w = w_capital / float(np.sum(w_capital))
    return 1.0 + float(np.sum(w_capital * (np.log(w_capital) + ln_gamma(w) - d - 1.0)))


# --------------------------------------------------------------------------
# 1) Model selection and the pressure argument
# --------------------------------------------------------------------------


def test_exactly_one_model_is_required(tessier2000_model: ct.NRTL, tessier2000_names: list[str]):
    mixture = _mixture(tessier2000_names, PLAIT_FEED)

    with pytest.raises(ct.ModelError):
        ct.stability_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=PRESSURE_PA)

    with pytest.raises(ct.ModelError) as excinfo:
        ct.stability_tp(
            mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            eos=ct.PengRobinsonEOS(),
            activity_model=tessier2000_model,
        )
    assert "gamma-phi" in str(excinfo.value)


def test_pressure_is_validated_but_does_not_change_the_activity_result(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> None:
    """`pressure_Pa` stays required for API uniformity; it is inert here."""
    mixture = _mixture(tessier2000_names, PLAIT_FEED)

    low = ct.stability_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=1.0e4,
        activity_model=tessier2000_model,
    )
    high = ct.stability_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=1.0e7,
        activity_model=tessier2000_model,
    )

    assert low.tpd_min == high.tpd_min
    assert low.trial_composition == high.trial_composition
    assert low.diagnostics["pressure_dependent"] is False
    assert low.diagnostics["model_family"] == "activity"

    with pytest.raises(ct.InputRangeError):
        ct.stability_tp(
            mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=0.0,
            activity_model=tessier2000_model,
        )
    with pytest.raises(ct.InputRangeError):
        ct.stability_tp(
            mixture,
            temperature_K=-1.0,
            pressure_Pa=PRESSURE_PA,
            activity_model=tessier2000_model,
        )


def test_mass_basis_is_rejected(tessier2000_model: ct.NRTL, tessier2000_names: list[str]) -> None:
    mass_basis = ct.Mixture.from_database(tessier2000_names, list(PLAIT_FEED), basis="mass")
    with pytest.raises(ct.ModelError):
        ct.stability_tp(
            mass_basis,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=tessier2000_model,
        )


# --------------------------------------------------------------------------
# 2) Identities and invariants
# --------------------------------------------------------------------------


def test_tpd_and_tm_vanish_at_the_feed_composition(tessier2000_ln_gamma: LnGamma) -> None:
    """tpd(z) = 0 and tm(W = z) = 0: the feed touches its own tangent plane."""
    z = np.asarray(PLAIT_FEED, dtype=float)
    d = _d_vector(tessier2000_ln_gamma, z)

    assert _tpd(tessier2000_ln_gamma, z, d) == pytest.approx(0.0, abs=1e-14)
    assert _tm(tessier2000_ln_gamma, z, d) == pytest.approx(0.0, abs=1e-14)


def test_converged_trials_are_stationary_points_of_tm(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str], tessier2000_ln_gamma: LnGamma
) -> None:
    """Equation (5) directly, plus a central finite-difference gradient of tm.

    Differentiating tm numerically is an independent check that the analytic
    gradient (4) used to derive the successive-substitution map and the Newton
    stage is the right one for an activity model too.
    """
    settings = ct.StabilitySettings()
    mixture = _mixture(tessier2000_names, PLAIT_FEED)
    result = ct.stability_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
        settings=settings,
    )
    d = _d_vector(tessier2000_ln_gamma, np.asarray(PLAIT_FEED, dtype=float))

    non_trivial = [t for t in result.trials if t.converged and not t.trivial]
    assert non_trivial, "expected at least one non-trivial stationary point"

    for trial in non_trivial:
        assert trial.composition is not None
        assert trial.residual < settings.tol

        w = np.asarray(trial.composition, dtype=float)
        w_capital = w * trial.sum_W

        residual = float(np.max(np.abs(np.log(w_capital) + tessier2000_ln_gamma(w) - d)))
        assert residual < 1e-9

        step = 1e-7
        for k in range(w_capital.size):
            plus = w_capital.copy()
            minus = w_capital.copy()
            plus[k] += step
            minus[k] -= step
            derivative = (
                _tm(tessier2000_ln_gamma, plus, d) - _tm(tessier2000_ln_gamma, minus, d)
            ) / (2.0 * step)
            assert abs(derivative) < 1e-7

        # Equations (6) and (7).
        tm_star = _tm(tessier2000_ln_gamma, w_capital, d)
        assert tm_star == pytest.approx(1.0 - trial.sum_W, abs=1e-10)
        assert trial.tpd == pytest.approx(-math.log(trial.sum_W), abs=1e-10)
        assert _tpd(tessier2000_ln_gamma, w, d) == pytest.approx(trial.tpd, abs=1e-12)
        assert tm_star == pytest.approx(1.0 - math.exp(-trial.tpd), abs=1e-10)


def test_component_reordering_permutes_the_result(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> None:
    """Relabelling components must not change the physics."""
    order = (2, 0, 1)
    base = ct.stability_tp(
        _mixture(tessier2000_names, PLAIT_FEED),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
    )
    other = ct.stability_tp(
        _mixture([tessier2000_names[i] for i in order], [PLAIT_FEED[i] for i in order]),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
    )

    assert base.status == other.status
    assert other.tpd_min == pytest.approx(base.tpd_min, rel=1e-9, abs=1e-14)
    assert base.trial_composition is not None
    assert other.trial_composition is not None
    assert base.k_values is not None
    assert other.k_values is not None

    expected_w = tuple(base.trial_composition[i] for i in order)
    expected_k = tuple(base.k_values[i] for i in order)
    assert np.allclose(other.trial_composition, expected_w, atol=1e-9)
    assert np.allclose(other.k_values, expected_k, atol=1e-9)


def test_results_are_deterministic(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> None:
    mixture = _mixture(tessier2000_names, PLAIT_FEED)
    kwargs: dict[str, Any] = {
        "temperature_K": TEMPERATURE_K,
        "pressure_Pa": PRESSURE_PA,
        "activity_model": tessier2000_model,
    }
    first = ct.stability_tp(mixture, **kwargs)
    second = ct.stability_tp(mixture, **kwargs)

    assert first.status == second.status
    assert first.tpd_min == second.tpd_min
    assert first.trial_composition == second.trial_composition
    assert first.k_values == second.k_values
    assert [t.iterations for t in first.trials] == [t.iterations for t in second.trials]


def test_every_trial_is_recorded_and_no_wilson_trial_is_used(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> None:
    """Wilson K-values are a vapor-liquid device; the activity path must not use them."""
    result = ct.stability_tp(
        _mixture(tessier2000_names, PLAIT_FEED),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
    )

    assert [trial.label for trial in result.trials] == [
        f"pure-{name}" for name in tessier2000_names
    ]
    assert result.diagnostics["trial_count"] == len(result.trials)
    assert all(trial.phase_branch is None for trial in result.trials)
    assert result.feed_branch is None
    assert result.phase_branch is None
    assert "feed_branch" not in result.diagnostics


# --------------------------------------------------------------------------
# 3) Negative controls
# --------------------------------------------------------------------------


def test_pure_component_feed_is_stable(tessier2000_model: ct.NRTL) -> None:
    """One component has no composition degree of freedom, so tpd == 0."""
    result = ct.stability_tp(
        _mixture(("Water",), (1.0,)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
    )

    assert result.status == "stable"
    assert result.tpd_min == pytest.approx(0.0, abs=1e-14)
    assert result.diagnostics["active_component_count"] == 1
    assert [trial.label for trial in result.trials] == ["pure-Water"]
    assert all(trial.trivial for trial in result.trials)


def test_component_at_zero_mole_fraction_is_treated_as_a_pure_feed(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> None:
    result = ct.stability_tp(
        _mixture(tessier2000_names, (0.0, 0.0, 1.0)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
    )
    assert result.diagnostics["active_component_count"] == 1
    assert result.status == "stable"


def test_ideal_solution_is_stable() -> None:
    """tau = 0 gives gamma = 1, so tpd(w) = sum_i w_i ln(w_i / z_i) >= 0."""
    model = ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs([("Benzene", "Water", 0.0, 0.0, 0.3, 0.3)])
    )
    for z in ((0.5, 0.5), (0.1, 0.9), (0.9, 0.1)):
        result = ct.stability_tp(
            _mixture(("Benzene", "Water"), z),
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=model,
        )
        assert result.status == "stable"
        assert result.tpd_min >= 0.0
        assert all(trial.converged for trial in result.trials)


def test_zero_fraction_component_stays_absent_from_the_trial_phase(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> None:
    result = ct.stability_tp(
        _mixture(tessier2000_names, (0.10, 0.0, 0.90)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
    )
    assert result.diagnostics["active_component_count"] == 2
    if result.trial_composition is not None:
        assert result.trial_composition[1] == 0.0


# --------------------------------------------------------------------------
# 4) A partially miscible binary (validation Case S-8)
# --------------------------------------------------------------------------


def test_butanol_water_binary_splits_and_both_phases_share_one_tangent_plane() -> None:
    """n-butanol / water at z1 = 0.10 is unstable; its LLE phases are marginal.

    The two liquid phases of an LLE split are the two points where one common
    hyperplane touches the Gibbs surface. So each of them, used as a feed, must
    show `tpd_min = 0` and must find *the other* phase as its stationary point.
    A strictly negative tpd at either would mean the split is not an
    equilibrium. This is the activity-model analogue of validation Case S-3.

    The binodal pair is computed here from the equal-activity condition
    ``x_i gamma_i`` equal in both phases, by a Newton solve that shares no code
    with `stability_tp`.
    """
    model = _butanol_water_model()
    ln_gamma = _butanol_water_ln_gamma(model)

    unstable = ct.stability_tp(
        _mixture(BUTANOL_WATER_NAMES, (0.10, 0.90)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
    )
    assert unstable.status == "unstable"
    assert unstable.tpd_min == pytest.approx(-0.029994488835, rel=1e-9)
    assert unstable.trial_composition is not None
    assert np.allclose(unstable.trial_composition, (0.419473294157, 0.580526705843), atol=1e-9)
    # tpd < 0 at the reported minimizer, recomputed from the definition.
    d_feed = _d_vector(ln_gamma, np.asarray((0.10, 0.90), dtype=float))
    assert _tpd(
        ln_gamma, np.asarray(unstable.trial_composition, dtype=float), d_feed
    ) == pytest.approx(unstable.tpd_min, abs=1e-12)

    # Independent binodal: activities equal in both phases.
    def activities(x1: float) -> np.ndarray:
        x = np.array([x1, 1.0 - x1], dtype=float)
        return x * np.exp(ln_gamma(x))

    def mismatch(u: np.ndarray) -> np.ndarray:
        return activities(float(u[0])) - activities(float(u[1]))

    u = np.array([0.02, 0.50], dtype=float)
    for _ in range(200):
        f = mismatch(u)
        if float(np.max(np.abs(f))) < 1e-15:
            break
        jacobian = np.zeros((2, 2), dtype=float)
        h = 1e-7
        for k in range(2):
            plus = u.copy()
            minus = u.copy()
            plus[k] += h
            minus[k] -= h
            jacobian[:, k] = (mismatch(plus) - mismatch(minus)) / (2.0 * h)
        u = u + np.linalg.solve(jacobian, -f)
    assert float(np.max(np.abs(mismatch(u)))) < 1e-14
    # The feed must lie strictly between the two conjugate phases.
    assert u[0] < 0.10 < u[1]

    for this_phase, other_phase in ((u[0], u[1]), (u[1], u[0])):
        result = ct.stability_tp(
            _mixture(BUTANOL_WATER_NAMES, (this_phase, 1.0 - this_phase)),
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=model,
        )
        assert result.status == "stable"
        assert result.tpd_min == pytest.approx(0.0, abs=1e-10)
        assert result.trial_composition is not None
        assert result.trial_composition[0] == pytest.approx(other_phase, abs=1e-9)


# --------------------------------------------------------------------------
# 5) The second-order stage
# --------------------------------------------------------------------------


def test_successive_substitution_alone_cannot_solve_the_near_plait_feed(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> None:
    """The reason the second-order stage exists, asserted rather than assumed.

    With `second_order=False` and a 1000-iteration budget, no trial at the
    near-plait-point feed of Table 2 reaches the stationarity tolerance; the
    fixed-point map (8) has a contraction ratio too close to one. With the
    default settings the same feed converges in a handful of Newton steps.
    """
    mixture = _mixture(tessier2000_names, PLAIT_FEED)

    ssi_only = ct.stability_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
        settings=ct.StabilitySettings(max_iter=1000, second_order=False),
    )
    assert ssi_only.status == "inconclusive"
    assert all(not trial.converged for trial in ssi_only.trials)
    assert all(trial.termination_reason == "max_iter" for trial in ssi_only.trials)

    default = ct.stability_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
    )
    assert default.status == "unstable"
    assert all(trial.converged for trial in default.trials)
    assert all(trial.converged_stage == "second-order" for trial in default.trials)
    assert all(0 < trial.second_order_iterations <= 10 for trial in default.trials)
    assert default.diagnostics["minimizing_trial_stage"] == "second-order"
    assert default.diagnostics["second_order_trial_count"] == len(default.trials)


def test_second_order_settings_are_validated() -> None:
    for kwargs in (
        {"ssi_iterations": 0},
        {"second_order_max_iter": 0},
        {"second_order_max_step": 0.0},
    ):
        with pytest.raises(ct.InputRangeError):
            ct.StabilitySettings(**kwargs)  # type: ignore[arg-type]


def test_iteration_counts_split_into_stages(
    tessier2000_model: ct.NRTL, tessier2000_names: list[str]
) -> None:
    result = ct.stability_tp(
        _mixture(tessier2000_names, PLAIT_FEED),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
    )
    for trial in result.trials:
        assert trial.iterations == trial.ssi_iterations + trial.second_order_iterations
        assert trial.ssi_iterations == 50  # the default ssi_iterations budget
