"""Automatic phase detection in `flash_tp` (ADR-0008).

Covers the tangent-plane path's verdicts, its seeded split, the verification
residuals it reports, the legacy escape hatch, and the failure semantics. The
external cross-check against `thermo` lives in
`tests/validation/test_flash_phase_detection_vs_thermo.py`.
"""

from __future__ import annotations

import itertools
from typing import Iterator

import numpy as np
import pytest

import chemthermo as ct

EOS = ct.PengRobinsonEOS()
LEGACY = ct.FlashSettings(phase_detection="wilson-heuristic")

#: Binary and ternary hydrocarbon states used for the invariant sweep.
GRID_MIXTURES: tuple[tuple[tuple[str, ...], tuple[float, ...]], ...] = (
    (("Methane", "Ethane"), (0.5, 0.5)),
    (("Methane", "Propane"), (0.7, 0.3)),
    (("Ethane", "n-Heptane"), (0.7, 0.3)),
    (("Methane", "n-Pentane"), (0.6, 0.4)),
    (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2)),
    (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3)),
)
GRID_T_K = (170.0, 200.0, 240.0, 280.0, 320.0, 360.0)
GRID_P_PA = (2.0e5, 1.0e6, 3.0e6, 8.0e6)


def _mixture(names: tuple[str, ...], z: tuple[float, ...]) -> ct.Mixture:
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


def _grid_results(
    settings: ct.FlashSettings | None = None,
) -> Iterator[tuple[tuple[str, ...], tuple[float, ...], float, float, ct.FlashResult]]:
    """Yield every grid state the solver answers; skip the ones it refuses."""
    for names, z in GRID_MIXTURES:
        mixture = _mixture(names, z)
        for temperature_K in GRID_T_K:
            for pressure_Pa in GRID_P_PA:
                try:
                    result = ct.flash_tp(
                        mixture,
                        temperature_K=temperature_K,
                        pressure_Pa=pressure_Pa,
                        eos=EOS,
                        settings=settings,
                    )
                except ct.ConvergenceError:
                    continue
                yield names, z, temperature_K, pressure_Pa, result


# --------------------------------------------------------------------------
# Regression: the converged equilibrium must not move.
# --------------------------------------------------------------------------


def test_canonical_binary_split_is_preserved() -> None:
    """The pinned 240 K / 3 MPa binary vapor fraction survives the new seed."""
    mixture = _mixture(("Methane", "Ethane"), (0.5, 0.5))
    result = ct.flash_tp(mixture, temperature_K=240.0, pressure_Pa=3.0e6, eos=EOS)
    assert result.vapor_fraction == pytest.approx(0.67451818, rel=1e-6)


def test_legacy_path_reproduces_the_pre_slice_numbers_exactly() -> None:
    """`phase_detection="wilson-heuristic"` is the untouched old solver.

    The two values below were produced before this slice (they are the
    `result.vapor_fraction` entries of the CLI fixtures) and must still come
    back bit-identical through the legacy path.
    """
    binary = ct.flash_tp(
        _mixture(("Methane", "Ethane"), (0.5, 0.5)),
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=EOS,
        settings=LEGACY,
    )
    assert binary.vapor_fraction == 0.6745181801306899

    ternary = ct.flash_tp(
        _mixture(("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2)),
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=EOS,
        settings=LEGACY,
    )
    assert ternary.vapor_fraction == 0.46829043780053325


def test_both_paths_agree_on_every_state_where_both_find_two_phases() -> None:
    """Same equilibrium, different iteration path.

    The two seeds converge to the same fixed point, so the only difference is
    the width of the K-update tolerance (`tol = 1e-8` on `max_i |dK_i|`). Over
    the 47 states of this grid where both paths find two phases, the worst
    relative vapor-fraction difference is 8.60e-7.
    """
    compared = 0
    worst = 0.0
    for names, z, temperature_K, pressure_Pa, new in _grid_results():
        new_beta = new.vapor_fraction
        if new_beta is None or not (0.0 < new_beta < 1.0):
            continue
        try:
            old = ct.flash_tp(
                _mixture(names, z),
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                eos=EOS,
                settings=LEGACY,
            )
        except ct.ConvergenceError:
            continue
        old_beta = old.vapor_fraction
        if old_beta is None or not (0.0 < old_beta < 1.0):
            continue
        compared += 1
        worst = max(worst, abs(new_beta - old_beta) / old_beta)

    assert compared >= 20
    assert worst < 1e-6


def test_gamma_phi_stays_on_the_heuristic_path() -> None:
    """Gamma-phi numbers are untouched and say so in diagnostics (ADR-0008)."""
    components = tuple(ct.Component.from_database(name) for name in ("Methane", "Ethane"))
    mixture = ct.Mixture(
        components=components,
        composition=ct.Composition(fractions=(0.5, 0.5), basis="mole", normalize=False),
    )
    result = ct.flash_tp(
        mixture,
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=EOS,
        activity_model=ct.NRTL(),
        flash_mode="gamma-phi",
    )
    assert result.diagnostics["phase_detection"] == "wilson-heuristic"
    assert result.diagnostics["k_seed"] == "wilson"
    assert "stability_status" not in result.diagnostics
    assert result.vapor_fraction == 0.7648352545438684


# --------------------------------------------------------------------------
# The verdict is now a thermodynamic criterion.
# --------------------------------------------------------------------------


def test_single_phase_result_is_a_stability_verdict() -> None:
    mixture = _mixture(("Methane", "Ethane"), (0.5, 0.5))
    result = ct.flash_tp(mixture, temperature_K=450.0, pressure_Pa=1.0e5, eos=EOS)

    assert result.phase_names() == ["vapor"]
    assert result.vapor_fraction == pytest.approx(1.0)
    assert result.diagnostics["phase_detection"] == "tangent-plane"
    assert result.diagnostics["stability_status"] == "stable"
    assert result.diagnostics["termination_reason"] == "feed_stable_tangent_plane"
    assert result.diagnostics["phase_count"] == 1
    assert result.diagnostics["iterations"] == 0
    assert float(result.diagnostics["tpd_min"]) >= 0.0
    assert result.diagnostics["feed_branch"] == "vapor"
    assert int(result.diagnostics["stability_trials"]) >= 1


def test_two_phase_result_reports_the_seed_and_the_verification_residuals() -> None:
    mixture = _mixture(("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2))
    result = ct.flash_tp(mixture, temperature_K=240.0, pressure_Pa=3.0e6, eos=EOS)

    diagnostics = result.diagnostics
    assert diagnostics["phase_detection"] == "tangent-plane"
    assert diagnostics["stability_status"] == "unstable"
    assert diagnostics["k_seed"] == "stability"
    assert diagnostics["incipient_phase"] == "vapor"
    assert float(diagnostics["tpd_min"]) < 0.0
    assert float(diagnostics["mass_balance_residual"]) < 1e-10
    assert float(diagnostics["fugacity_residual"]) < 1e-6
    assert float(diagnostics["delta_g_split_rt"]) < 0.0


def test_heuristic_calls_a_two_phase_feed_single_phase_and_tangent_plane_does_not() -> None:
    """Validation Case F-2: the state where the two paths disagree.

    Methane(0.6) / n-Pentane(0.4) at 175 K and 1.778 MPa. The Wilson K-values
    straddle 1, so the K-bound test does not fire, but Rachford-Rice finds no
    root for them (`f(0) = -4.2e-2`, `f(1) = -1.7e4`, same sign) and the legacy
    path reports a single liquid. The feed is in fact unstable
    (`tpd_min = -4.266e-2`) and the split it seeds lowers the Gibbs energy.
    `thermo`'s `FlashVL` with the same constants gives VF = 0.0791276; see
    `tests/validation/test_flash_phase_detection_vs_thermo.py`.
    """
    mixture = _mixture(("Methane", "n-Pentane"), (0.6, 0.4))
    temperature_K, pressure_Pa = 175.0, 1.778e6

    legacy = ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=EOS,
        settings=LEGACY,
    )
    assert legacy.phase_names() == ["liquid"]
    assert legacy.diagnostics["termination_reason"] == "rr_no_root"
    assert float(legacy.diagnostics["k_min"]) < 1.0 < float(legacy.diagnostics["k_max"])

    result = ct.flash_tp(mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=EOS)
    assert set(result.phase_names()) == {"liquid", "vapor"}
    assert result.vapor_fraction == pytest.approx(0.0792360, abs=1e-6)
    assert float(result.diagnostics["tpd_min"]) == pytest.approx(-0.04265697, rel=1e-6)
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert result.diagnostics["k_seed"] == "stability"
    assert result.diagnostics["incipient_phase"] == "vapor"

    x = np.array(result.phases["liquid"].composition.fractions)
    y = np.array(result.phases["vapor"].composition.fractions)
    # The incipient phase is an almost pure methane vapor.
    assert y[0] > 0.999
    assert x[0] == pytest.approx(0.56558, abs=1e-4)


def test_legacy_nonconvergence_becomes_a_clean_single_phase_answer() -> None:
    """Validation Case F-1: a state the heuristic could not solve at all.

    At 300 K / 1e7 Pa the Wilson K-values straddle 1, so the legacy path starts
    a split that never converges. The feed is stable, so the tangent-plane path
    answers in one step.
    """
    mixture = _mixture(("Methane", "Ethane"), (0.5, 0.5))
    temperature_K, pressure_Pa = 300.0, 1.0e7

    with pytest.raises(ct.ConvergenceError):
        ct.flash_tp(
            mixture,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            eos=EOS,
            settings=LEGACY,
        )

    result = ct.flash_tp(mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=EOS)
    assert len(result.phase_names()) == 1
    assert result.diagnostics["stability_status"] == "stable"
    assert result.diagnostics["termination_reason"] == "feed_stable_tangent_plane"


# --------------------------------------------------------------------------
# Invariants over a grid.
# --------------------------------------------------------------------------


def test_grid_invariants() -> None:
    """Every two-phase result on the grid satisfies the verification identities.

    Achieved over the 47 two-phase states of this grid: mass balance
    <= 2.01e-13, fugacity residual <= 8.25e-9, largest (least negative)
    `delta_g_split_rt` = -1.458e-4. 97 further states are single-phase and every
    one of them reports `stability_status == "stable"`.
    """
    two_phase = 0
    single_phase = 0
    worst_mass = 0.0
    worst_fugacity = 0.0
    worst_delta_g = -np.inf

    for names, z, temperature_K, pressure_Pa, result in _grid_results():
        where = (names, temperature_K, pressure_Pa)
        beta = result.vapor_fraction
        assert beta is not None
        diagnostics = result.diagnostics
        assert diagnostics["phase_detection"] == "tangent-plane"

        if len(result.phase_names()) == 1:
            single_phase += 1
            assert diagnostics["stability_status"] == "stable", where
            assert beta in (0.0, 1.0), where
            continue

        two_phase += 1
        assert 0.0 < beta < 1.0, where
        assert diagnostics["stability_status"] == "unstable", where
        assert float(diagnostics["tpd_min"]) < 0.0, where

        mass = float(diagnostics["mass_balance_residual"])
        fugacity = float(diagnostics["fugacity_residual"])
        delta_g = float(diagnostics["delta_g_split_rt"])
        assert mass < 1e-10, where
        assert fugacity < 1e-6, where
        assert delta_g < 0.0, where
        worst_mass = max(worst_mass, mass)
        worst_fugacity = max(worst_fugacity, fugacity)
        worst_delta_g = max(worst_delta_g, delta_g)

        # Recompute the material balance from the returned phases rather than
        # trusting the residual the solver wrote into diagnostics.
        x = np.array(result.phases["liquid"].composition.fractions)
        y = np.array(result.phases["vapor"].composition.fractions)
        feed = np.array(_mixture(names, z).fractions)
        assert np.max(np.abs(feed - (beta * y + (1.0 - beta) * x))) < 1e-10, where

    assert two_phase >= 30, two_phase
    assert single_phase >= 1
    assert worst_mass < 1e-10
    assert worst_fugacity < 1e-6
    assert worst_delta_g < 0.0


def test_every_two_phase_grid_state_passes_the_post_split_stability_check() -> None:
    """Validation Case L-4: the post-split gate over the phi-phi grid (ADR-0009).

    Every converged phase of every two-phase state is fed back into
    `stability_tp`. Two coexisting phases share one tangent plane, so each
    stability test converges onto its *partner* and reports a tangent-plane
    distance of zero up to the split's own tolerance.

    Achieved over the 47 two-phase states (94 phases): every phase reports
    `"stable"`, none needed the `"marginal"` (converged-onto-the-partner)
    reclassification at the default `tol = 1e-8`, and the most negative
    post-split `tpd_min` seen anywhere on the grid is **-7.0055e-09**, at
    Ethane/n-Heptane (0.7, 0.3), 360 K, 1 MPa. That is inside the default
    `tpd_tol = 1e-8` by only a factor of 1.4, which is why the
    converged-onto-the-partner rule exists at all; `tests/test_flash_lle.py`
    exercises it directly by loosening `tol`.

    No state on this grid needs a third phase, so no state raises.
    """
    two_phase = 0
    worst = 0.0
    marginal = 0

    for names, _z, temperature_K, pressure_Pa, result in _grid_results():
        if len(result.phase_names()) != 2:
            assert "post_split_checked" not in result.diagnostics
            continue
        where = (names, temperature_K, pressure_Pa)
        two_phase += 1
        diagnostics = result.diagnostics
        assert diagnostics["post_split_checked"] is True, where
        assert diagnostics["post_split_stable"] is True, where
        assert diagnostics["post_split_status"] == "stable", where
        for phase in ("liquid", "vapor"):
            verdict = diagnostics[f"phase_stability_{phase}"]
            assert verdict in {"stable", "marginal"}, where
            marginal += verdict == "marginal"
        worst = min(worst, float(diagnostics["post_split_tpd_min"]))

    assert two_phase >= 30, two_phase
    assert worst > -1e-8
    assert worst == pytest.approx(-7.0055e-09, rel=1e-3)
    assert marginal == 0


def test_results_are_permutation_invariant_and_deterministic() -> None:
    cases = (
        (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2), 240.0, 3.0e6),
        (("Methane", "n-Pentane"), (0.6, 0.4), 175.0, 1.778e6),
        (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3), 300.0, 1.0e6),
        (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2), 450.0, 1.0e5),
    )
    for names, z, temperature_K, pressure_Pa in cases:
        base = ct.flash_tp(
            _mixture(names, z), temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=EOS
        )

        # Determinism: repeated calls are bit-identical, diagnostics included.
        repeat = ct.flash_tp(
            _mixture(names, z), temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=EOS
        )
        assert repeat.vapor_fraction == base.vapor_fraction
        assert dict(repeat.diagnostics) == dict(base.diagnostics)

        for order in itertools.permutations(range(len(names))):
            permuted = ct.flash_tp(
                _mixture(tuple(names[i] for i in order), tuple(z[i] for i in order)),
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                eos=EOS,
            )
            assert permuted.phase_names() == base.phase_names()
            assert permuted.vapor_fraction == pytest.approx(base.vapor_fraction, abs=1e-12)
            undo = np.argsort(order)
            for phase in base.phase_names():
                restored = np.array(permuted.phases[phase].composition.fractions)[undo]
                assert np.allclose(
                    restored,
                    np.array(base.phases[phase].composition.fractions),
                    rtol=0.0,
                    atol=1e-12,
                )


# --------------------------------------------------------------------------
# Failure semantics and settings.
# --------------------------------------------------------------------------


def test_inconclusive_stability_raises_and_the_legacy_path_still_answers() -> None:
    starved = ct.StabilitySettings(max_iter=1, second_order=False)
    mixture = _mixture(("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2))

    assert (
        ct.stability_tp(
            mixture, temperature_K=240.0, pressure_Pa=3.0e6, eos=EOS, settings=starved
        ).status
        == "inconclusive"
    )

    with pytest.raises(ct.ConvergenceError, match="inconclusive"):
        ct.flash_tp(
            mixture,
            temperature_K=240.0,
            pressure_Pa=3.0e6,
            eos=EOS,
            settings=ct.FlashSettings(stability_settings=starved),
        )

    fallback = ct.flash_tp(
        mixture,
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=EOS,
        settings=ct.FlashSettings(phase_detection="wilson-heuristic", stability_settings=starved),
    )
    assert fallback.vapor_fraction == 0.46829043780053325
    assert fallback.diagnostics["phase_detection"] == "wilson-heuristic"


def test_split_iteration_limit_still_raises() -> None:
    mixture = _mixture(("Methane", "Ethane"), (0.5, 0.5))
    with pytest.raises(ct.ConvergenceError, match="did not converge"):
        ct.flash_tp(
            mixture,
            temperature_K=240.0,
            pressure_Pa=3.0e6,
            eos=EOS,
            settings=ct.FlashSettings(max_iter=1, tol=1e-12),
        )


def test_phase_detection_setting_is_validated() -> None:
    with pytest.raises(ct.InputRangeError, match="phase_detection"):
        ct.FlashSettings(phase_detection="michelsen")
    assert ct.FlashSettings().phase_detection == "tangent-plane"
    assert ct.FlashSettings().stability_settings is None
    assert ct.FlashSettings(
        stability_settings=ct.StabilitySettings(tpd_tol=1e-6)
    ).stability_settings
