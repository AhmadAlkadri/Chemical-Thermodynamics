"""Tangent-plane phase stability tests (Michelsen).

Equation numbers refer to the module docstring of
``src/chemthermo/stability/tp.py``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.stability.tp import _ln_phi_min_gibbs

EOS = ct.PengRobinsonEOS()

TERNARY = ("Methane", "Ethane", "Propane")
TERNARY_Z = (0.50, 0.30, 0.20)
UNSTABLE_T_K = 240.0
UNSTABLE_P_PA = 3.0e6


def _mixture(names, z):
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


def _d_vector(mixture, temperature_K, pressure_Pa):
    """Feed tangent-plane intercept d_i = ln z_i + ln phi_i(z), equation (1)."""
    z = np.array(mixture.fractions, dtype=float)
    ln_phi, _ = _ln_phi_min_gibbs(
        EOS,
        mixture=mixture,
        temperature=temperature_K,
        pressure=pressure_Pa,
        composition=z,
    )
    return np.log(z) + ln_phi


def _tpd(mixture, temperature_K, pressure_Pa, w, d):
    """Reduced tangent-plane distance, equation (2), from first principles."""
    w = np.array(w, dtype=float)
    w = w / float(np.sum(w))
    ln_phi, _ = _ln_phi_min_gibbs(
        EOS,
        mixture=mixture,
        temperature=temperature_K,
        pressure=pressure_Pa,
        composition=w,
    )
    return float(np.sum(w * (np.log(w) + ln_phi - d)))


def _tm(mixture, temperature_K, pressure_Pa, w_capital, d):
    """Modified (unnormalized) tangent-plane function, equation (3)."""
    w_capital = np.array(w_capital, dtype=float)
    w = w_capital / float(np.sum(w_capital))
    ln_phi, _ = _ln_phi_min_gibbs(
        EOS,
        mixture=mixture,
        temperature=temperature_K,
        pressure=pressure_Pa,
        composition=w,
    )
    return 1.0 + float(np.sum(w_capital * (np.log(w_capital) + ln_phi - d - 1.0)))


# --------------------------------------------------------------------------
# 1) Identity at the feed
# --------------------------------------------------------------------------


def test_tpd_and_tm_vanish_at_the_feed_composition() -> None:
    """tpd(z) = 0 and tm(W = z) = 0 because the feed touches its own tangent plane."""
    mixture = _mixture(TERNARY, TERNARY_Z)
    d = _d_vector(mixture, UNSTABLE_T_K, UNSTABLE_P_PA)
    z = np.array(mixture.fractions, dtype=float)

    assert _tpd(mixture, UNSTABLE_T_K, UNSTABLE_P_PA, z, d) == pytest.approx(0.0, abs=1e-12)
    assert _tm(mixture, UNSTABLE_T_K, UNSTABLE_P_PA, z, d) == pytest.approx(0.0, abs=1e-12)


# --------------------------------------------------------------------------
# 2) Stationarity: residual, finite-difference gradient, and tm/tpd/sum_W
# --------------------------------------------------------------------------


def test_converged_non_trivial_trials_are_stationary_points_of_tm() -> None:
    """Check equation (5) directly and by central finite differences of tm.

    A converged trial must satisfy ln W_i + ln phi_i(w) - d_i = 0, which is the
    gradient of tm with respect to W_i (equation (4)). Differentiating tm
    numerically is an independent check that the analytic gradient used to
    derive the successive-substitution map is the right one.
    """
    mixture = _mixture(TERNARY, TERNARY_Z)
    settings = ct.StabilitySettings()
    result = ct.stability_tp(
        mixture,
        temperature_K=UNSTABLE_T_K,
        pressure_Pa=UNSTABLE_P_PA,
        eos=EOS,
        settings=settings,
    )
    d = _d_vector(mixture, UNSTABLE_T_K, UNSTABLE_P_PA)

    non_trivial = [t for t in result.trials if t.converged and not t.trivial]
    assert non_trivial, "expected at least one non-trivial stationary point"

    for trial in non_trivial:
        assert trial.composition is not None
        assert trial.residual < settings.tol

        w = np.array(trial.composition, dtype=float)
        w_capital = w * trial.sum_W

        # Analytic stationarity residual recomputed from first principles.
        ln_phi, _ = _ln_phi_min_gibbs(
            EOS,
            mixture=mixture,
            temperature=UNSTABLE_T_K,
            pressure=UNSTABLE_P_PA,
            composition=w,
        )
        residual = np.max(np.abs(np.log(w_capital) + ln_phi - d))
        assert residual < 1e-9

        # Central finite-difference gradient of tm at the stationary point.
        step = 1e-6
        for k in range(w_capital.size):
            plus = w_capital.copy()
            minus = w_capital.copy()
            plus[k] += step
            minus[k] -= step
            derivative = (
                _tm(mixture, UNSTABLE_T_K, UNSTABLE_P_PA, plus, d)
                - _tm(mixture, UNSTABLE_T_K, UNSTABLE_P_PA, minus, d)
            ) / (2.0 * step)
            assert abs(derivative) < 1e-7

        # Equations (6) and (7): tm* = 1 - sum_W and tpd = -ln(sum_W).
        tm_star = _tm(mixture, UNSTABLE_T_K, UNSTABLE_P_PA, w_capital, d)
        assert tm_star == pytest.approx(1.0 - trial.sum_W, abs=1e-10)
        assert trial.tpd == pytest.approx(-math.log(trial.sum_W), abs=1e-10)
        assert _tpd(mixture, UNSTABLE_T_K, UNSTABLE_P_PA, w, d) == pytest.approx(
            trial.tpd, abs=1e-12
        )
        # tm and tpd always agree in sign: tm* = 1 - exp(-tpd).
        assert tm_star == pytest.approx(1.0 - math.exp(-trial.tpd), abs=1e-10)


# --------------------------------------------------------------------------
# 3) Canonical Peng-Robinson states
# --------------------------------------------------------------------------


def test_canonical_ternary_is_unstable_and_matches_the_flash_split_direction() -> None:
    """240 K, 3 MPa methane/ethane/propane is a two-phase state, so tpd_min < 0."""
    mixture = _mixture(TERNARY, TERNARY_Z)
    result = ct.stability_tp(
        mixture, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS
    )

    assert result.status == "unstable"
    assert result.stable is False
    assert result.tpd_min < 0.0
    assert result.tpd_min == pytest.approx(-0.3492770207, abs=1e-6)
    assert result.trial_composition is not None
    assert result.k_values is not None
    assert result.phase_branch == "vapor"
    assert result.feed_branch == "liquid"
    assert sum(result.trial_composition) == pytest.approx(1.0, abs=1e-12)

    flash = ct.flash_tp(mixture, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS)
    x = np.array(flash.phases["liquid"].composition.fractions)
    y = np.array(flash.phases["vapor"].composition.fractions)
    flash_k = y / x

    stability_k = np.array(result.k_values)
    # Same split direction: the incipient phase is enriched in exactly the
    # components the equilibrium vapor is enriched in.
    assert np.all(np.sign(np.log(stability_k)) == np.sign(np.log(flash_k)))
    # Same volatility ordering.
    assert np.all(np.diff(stability_k) < 0.0)
    assert np.all(np.diff(flash_k) < 0.0)


def test_hot_dilute_binary_is_stable() -> None:
    """450 K, 1 bar methane/ethane is far above both critical temperatures."""
    mixture = _mixture(("Methane", "Ethane"), (0.5, 0.5))
    result = ct.stability_tp(mixture, temperature_K=450.0, pressure_Pa=1.0e5, eos=EOS)

    assert result.status == "stable"
    assert result.stable is True
    assert result.tpd_min >= 0.0
    assert result.feed_branch == "vapor"
    assert all(trial.converged for trial in result.trials)


@pytest.mark.parametrize(
    ("temperature_K", "pressure_Pa"),
    [(300.0, 5.0e7), (200.0, 1.0e8)],
)
def test_dense_high_pressure_states_are_stable(temperature_K: float, pressure_Pa: float) -> None:
    """Dense supercritical states of the canonical ternary show no negative tpd."""
    mixture = _mixture(TERNARY, TERNARY_Z)
    result = ct.stability_tp(mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=EOS)

    assert result.status == "stable"
    assert result.tpd_min >= 0.0


# --------------------------------------------------------------------------
# 4) Consistency with the converged flash
# --------------------------------------------------------------------------


def test_equilibrium_phases_are_marginally_stable() -> None:
    """Equilibrium phases share one tangent plane, so each is marginally stable.

    At a converged two-phase solution the equal-fugacity condition means the
    Gibbs surface is touched by the same hyperplane at x and at y. Running the
    stability test with x as the feed must therefore find its minimum tangent
    plane distance at zero, and the stationary composition it finds must be the
    equilibrium vapor y (the other point of tangency), and symmetrically with y
    as the feed. A strictly negative tpd here would mean the flash returned a
    non-equilibrium split.
    """
    mixture = _mixture(TERNARY, TERNARY_Z)
    flash = ct.flash_tp(mixture, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS)
    x = tuple(flash.phases["liquid"].composition.fractions)
    y = tuple(flash.phases["vapor"].composition.fractions)

    for feed, partner in ((x, y), (y, x)):
        feed_mixture = _mixture(TERNARY, feed)
        result = ct.stability_tp(
            feed_mixture, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS
        )

        assert result.tpd_min >= -1e-6
        assert result.status == "stable"

        if result.trial_composition is None:
            # Every trial collapsed onto the trivial solution; that is the other
            # admissible outcome of a marginally stable feed.
            continue
        assert np.allclose(np.array(result.trial_composition), np.array(partner), atol=1e-6)


# --------------------------------------------------------------------------
# 5) Permutation invariance
# --------------------------------------------------------------------------


def test_component_reordering_permutes_the_result() -> None:
    """Relabelling components must not change the physics."""
    order = (2, 0, 1)
    mixture = _mixture(TERNARY, TERNARY_Z)
    permuted = _mixture(tuple(TERNARY[i] for i in order), tuple(TERNARY_Z[i] for i in order))

    base = ct.stability_tp(mixture, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS)
    other = ct.stability_tp(
        permuted, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS
    )

    assert base.status == other.status
    assert other.tpd_min == pytest.approx(base.tpd_min, abs=1e-10)
    assert base.trial_composition is not None
    assert other.trial_composition is not None
    assert base.k_values is not None
    assert other.k_values is not None

    expected_w = tuple(base.trial_composition[i] for i in order)
    expected_k = tuple(base.k_values[i] for i in order)
    assert np.allclose(other.trial_composition, expected_w, atol=1e-10)
    assert np.allclose(other.k_values, expected_k, atol=1e-10)


def test_results_are_deterministic() -> None:
    mixture = _mixture(TERNARY, TERNARY_Z)
    first = ct.stability_tp(mixture, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS)
    second = ct.stability_tp(
        mixture, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS
    )
    assert first.status == second.status
    assert first.tpd_min == second.tpd_min
    assert first.trial_composition == second.trial_composition
    assert first.k_values == second.k_values


# --------------------------------------------------------------------------
# 6) Degenerate feeds
# --------------------------------------------------------------------------


def test_pure_component_feed_runs_and_reports_a_status() -> None:
    """A one-component feed has no composition degrees of freedom: tpd == 0."""
    mixture = _mixture(("Methane",), (1.0,))
    result = ct.stability_tp(mixture, temperature_K=300.0, pressure_Pa=1.0e5, eos=EOS)

    assert result.status in {"stable", "inconclusive"}
    assert result.tpd_min == pytest.approx(0.0, abs=1e-12)
    assert all(trial.trivial for trial in result.trials if trial.converged)
    assert result.diagnostics["active_component_count"] == 1


def test_near_pure_feed_runs_without_error() -> None:
    mixture = _mixture(("Methane", "Ethane"), (1.0 - 1e-6, 1e-6))
    result = ct.stability_tp(
        mixture, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS
    )

    assert result.status in {"stable", "unstable", "inconclusive"}
    assert math.isfinite(result.tpd_min)


def test_zero_fraction_component_stays_absent_from_the_trial_phase() -> None:
    mixture = _mixture(TERNARY, (0.5, 0.5, 0.0))
    result = ct.stability_tp(
        mixture, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS
    )

    assert result.diagnostics["active_component_count"] == 2
    if result.trial_composition is not None:
        assert result.trial_composition[2] == 0.0


# --------------------------------------------------------------------------
# Input validation and settings
# --------------------------------------------------------------------------


def test_invalid_inputs_are_rejected() -> None:
    mixture = _mixture(TERNARY, TERNARY_Z)

    with pytest.raises(ct.InputRangeError):
        ct.stability_tp(mixture, temperature_K=-1.0, pressure_Pa=UNSTABLE_P_PA, eos=EOS)
    with pytest.raises(ct.InputRangeError):
        ct.stability_tp(mixture, temperature_K=UNSTABLE_T_K, pressure_Pa=0.0, eos=EOS)
    with pytest.raises(ct.ModelError):
        ct.stability_tp(
            mixture,
            temperature_K=UNSTABLE_T_K,
            pressure_Pa=UNSTABLE_P_PA,
            eos=None,  # type: ignore[arg-type]
        )

    mass_basis = ct.Mixture.from_database(list(TERNARY), list(TERNARY_Z), basis="mass")
    with pytest.raises(ct.ModelError):
        ct.stability_tp(mass_basis, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS)


def test_stability_settings_validation() -> None:
    for kwargs in (
        {"max_iter": 0},
        {"tol": 0.0},
        {"trivial_tol": -1.0},
        {"tpd_tol": 0.0},
    ):
        with pytest.raises(ct.InputRangeError):
            ct.StabilitySettings(**kwargs)  # type: ignore[arg-type]


def test_every_trial_is_recorded() -> None:
    mixture = _mixture(TERNARY, TERNARY_Z)
    result = ct.stability_tp(
        mixture, temperature_K=UNSTABLE_T_K, pressure_Pa=UNSTABLE_P_PA, eos=EOS
    )
    labels = [trial.label for trial in result.trials]
    assert labels == [
        "wilson-vapor",
        "wilson-liquid",
        "pure-Methane",
        "pure-Ethane",
        "pure-Propane",
    ]
    assert result.diagnostics["trial_count"] == len(result.trials)
