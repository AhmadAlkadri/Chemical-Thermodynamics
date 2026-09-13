"""Liquid-liquid `flash_tp` with an activity model only (ADR-0009).

Covers the mode contract (`eos=None` plus `activity_model=`), the one-versus-two
liquid verdict, the binary binodal and lever rule against an independently
written equal-activity solve, the `FlashResult` invariants of a gamma-gamma
answer, the post-split stability check on both paths, and the failure semantics.

The published multicomponent tie-lines live in
`tests/validation/test_flash_lle_tessier2000.py`.
"""

from __future__ import annotations

import itertools
from typing import Callable, Sequence

import numpy as np
import pytest

import chemthermo as ct

LnGamma = Callable[[np.ndarray], np.ndarray]

TEMPERATURE_K = 298.15
PRESSURE_PA = 101325.0

# n-butanol / water: the 2-3 pair of Tessier, Brennecke and Stadtherr,
# Chem. Eng. Sci. 55 (2000) 1785, Table 1 (alpha implied by the printed G).
# Same parameters and provenance as validation Case S-8.
NAMES = ("n-Butanol", "Water")
TAU_12 = 0.90047
TAU_21 = 3.51307
ALPHA = 0.48

EOS = ct.PengRobinsonEOS()


def _mixture(names: Sequence[str], z: Sequence[float]) -> ct.Mixture:
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


def _model() -> ct.NRTL:
    return ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [(NAMES[0], NAMES[1], TAU_12, TAU_21, ALPHA, ALPHA)]
        )
    )


def _ln_gamma(model: ct.NRTL) -> LnGamma:
    mixture = _mixture(NAMES, (0.5, 0.5))

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


def _binodal(ln_gamma: LnGamma) -> tuple[float, float]:
    """Conjugate compositions from equal activities, solved here from scratch.

    Unknowns are the two butanol mole fractions; equations are
    ``x_i gamma_i`` equal in both phases. Newton with a finite-difference
    Jacobian; shares no code with `chemthermo.flash`.
    """

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
        step = 1e-7
        for column in range(2):
            plus = u.copy()
            minus = u.copy()
            plus[column] += step
            minus[column] -= step
            jacobian[:, column] = (mismatch(plus) - mismatch(minus)) / (2.0 * step)
        u = u + np.linalg.solve(jacobian, -residual)
    assert float(np.max(np.abs(mismatch(u)))) < 1e-14
    return float(u[0]), float(u[1])


def _flash(z: Sequence[float], settings: ct.FlashSettings | None = None) -> ct.FlashResult:
    return ct.flash_tp(
        _mixture(NAMES, z),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=_model(),
        settings=settings,
    )


def _phase_set(result: ct.FlashResult) -> list[tuple[float, ...]]:
    """Sorted phase compositions: the label-independent content of a result."""
    return sorted(tuple(result.phases[name].composition.fractions) for name in result.phase_names())


# --------------------------------------------------------------------------
# Mode contract
# --------------------------------------------------------------------------


def test_activity_model_without_an_eos_infers_liquid_liquid_mode() -> None:
    result = _flash((0.10, 0.90))
    assert result.diagnostics["flash_mode"] == "gamma-gamma"
    assert result.diagnostics["phase_detection"] == "tangent-plane"
    assert set(result.phase_names()) == {"liquid1", "liquid2"}

    explicit = ct.flash_tp(
        _mixture(NAMES, (0.10, 0.90)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=_model(),
        flash_mode="gamma-gamma",
    )
    assert dict(explicit.diagnostics) == dict(result.diagnostics)
    assert _phase_set(explicit) == _phase_set(result)


def test_mode_and_model_combinations_are_validated() -> None:
    mixture = _mixture(NAMES, (0.10, 0.90))

    # phi-phi still requires an EOS, even when an activity model is present.
    with pytest.raises(ct.ModelError, match="equation-of-state"):
        ct.flash_tp(
            mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=_model(),
            flash_mode="phi-phi",
        )

    # gamma-gamma describes both phases with the activity model.
    with pytest.raises(ct.ModelError, match="must not be given"):
        ct.flash_tp(
            mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            eos=EOS,
            activity_model=_model(),
            flash_mode="gamma-gamma",
        )

    with pytest.raises(ct.ModelError, match="activity model is required"):
        ct.flash_tp(
            mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            flash_mode="gamma-gamma",
        )

    # No model at all is still an error, with or without the keyword.
    with pytest.raises(ct.ModelError):
        ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=PRESSURE_PA)
    with pytest.raises(ct.ModelError):
        ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=PRESSURE_PA, eos=None)

    with pytest.raises(ct.ModelError, match="Unsupported flash_mode"):
        ct.flash_tp(
            mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=_model(),
            flash_mode="gamma-gamma-gamma",
        )


# --------------------------------------------------------------------------
# The verdict, and the tie-line it produces
# --------------------------------------------------------------------------


def test_stable_feed_returns_one_liquid() -> None:
    result = _flash((0.45, 0.55))
    assert result.phase_names() == ["liquid"]
    assert result.vapor_fraction is None
    assert result.phase_fractions == {"liquid": 1.0}
    assert result.diagnostics["stability_status"] == "stable"
    assert result.diagnostics["termination_reason"] == "feed_stable_tangent_plane"
    assert result.diagnostics["phase_count"] == 1
    assert float(result.diagnostics["tpd_min"]) > 0.0
    assert result.phases["liquid"].composition.fractions == pytest.approx((0.45, 0.55))


def test_single_component_feed_is_one_liquid() -> None:
    result = ct.flash_tp(
        _mixture(("Water",), (1.0,)),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=_model(),
    )
    assert result.phase_names() == ["liquid"]
    assert result.vapor_fraction is None


def test_binary_binodal_is_independent_of_the_feed_and_obeys_the_lever_rule() -> None:
    """Validation Case L-3.

    Four feeds inside the miscibility gap must return the *same* pair of
    conjugate compositions - the tie-line is a property of (T, P, model), not of
    the feed - and the phase amounts must follow the lever rule exactly. The
    reference pair comes from an equal-activity solve written in this module.

    Achieved: binodal x(n-Butanol) = 0.0199984194669 and 0.3599996615084, which
    every feed reproduces to <= 1.7e-12; lever-rule error <= 2.8e-12.

    Note on the feed set: 0.30 is *inside* the gap (0.0199984 < 0.30 <
    0.3599997) and does split. Only a feed above the upper binodal branch, such
    as 0.45, is a single liquid.
    """
    model = _model()
    low, high = _binodal(_ln_gamma(model))
    assert low == pytest.approx(0.019998419467, abs=1e-9)
    assert high == pytest.approx(0.359999661508, abs=1e-9)

    worst_composition = 0.0
    worst_lever = 0.0
    for feed in (0.05, 0.10, 0.20, 0.30):
        assert low < feed < high
        result = _flash((feed, 1.0 - feed))
        assert set(result.phase_names()) == {"liquid1", "liquid2"}

        phases = _phase_set(result)
        # Sorted by the first component, so phases[0] is the water-rich phase.
        worst_composition = max(
            worst_composition,
            abs(phases[0][0] - low),
            abs(phases[1][0] - high),
            abs(phases[0][1] - (1.0 - low)),
            abs(phases[1][1] - (1.0 - high)),
        )

        rich_name = max(
            result.phase_names(), key=lambda name: result.phases[name].composition.fractions[0]
        )
        expected = (feed - low) / (high - low)
        worst_lever = max(worst_lever, abs(result.phase_fractions[rich_name] - expected))

        assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
        assert float(result.diagnostics["mass_balance_residual"]) < 1e-12
        assert float(result.diagnostics["equilibrium_residual"]) < 1e-10

    assert worst_composition < 1e-6
    assert worst_lever < 1e-6

    outside = _flash((0.45, 0.55))
    assert outside.phase_names() == ["liquid"]


def test_two_phase_result_reports_the_verification_and_post_split_diagnostics() -> None:
    result = _flash((0.10, 0.90))
    diagnostics = result.diagnostics

    assert diagnostics["stability_status"] == "unstable"
    assert float(diagnostics["tpd_min"]) < 0.0
    assert diagnostics["k_seed"] == "stability"
    assert diagnostics["phase_regime"] == "LLE"
    assert diagnostics["phase_count"] == 2
    assert diagnostics["converged"] is True
    assert diagnostics["converged_stage"] in {"successive-substitution", "second-order"}
    assert int(diagnostics["ssi_iterations"]) >= 1
    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["equilibrium_residual"]) < 1e-10
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert "fugacity_residual" not in diagnostics

    assert diagnostics["post_split_checked"] is True
    assert diagnostics["post_split_stable"] is True
    assert diagnostics["post_split_status"] == "stable"
    assert diagnostics["phase_stability_liquid1"] == "stable"
    assert diagnostics["phase_stability_liquid2"] == "stable"
    assert float(diagnostics["post_split_tpd_min"]) > -1e-8
    for name in ("liquid1", "liquid2"):
        assert abs(float(diagnostics[f"phase_stability_tpd_min_{name}"])) < 1e-8


def test_flash_result_invariants_of_a_liquid_liquid_answer() -> None:
    result = _flash((0.10, 0.90))
    assert result.vapor_fraction is None
    assert set(result.phase_fractions) == set(result.phases)
    assert sum(result.phase_fractions.values()) == pytest.approx(1.0, abs=1e-12)
    for name, phase in result.phases.items():
        assert phase.name == name
        assert sum(phase.composition.fractions) == pytest.approx(1.0, abs=1e-12)
        assert 0.0 < result.phase_fractions[name] < 1.0

    feed = np.array(_mixture(NAMES, (0.10, 0.90)).fractions)
    total = np.zeros_like(feed)
    for name, phase in result.phases.items():
        total += result.phase_fractions[name] * np.array(phase.composition.fractions)
    assert np.max(np.abs(total - feed)) < 1e-12


# --------------------------------------------------------------------------
# Invariants
# --------------------------------------------------------------------------


def test_results_are_deterministic_and_permutation_invariant() -> None:
    """The phase *set* is invariant; the labels are roles and may swap.

    Achieved over the two-phase cases below: worst composition difference after
    undoing the permutation 1.4e-13, worst phase-fraction difference 1.4e-13.
    """
    parameters = ct.NRTLParameters.from_pairs([(NAMES[0], NAMES[1], TAU_12, TAU_21, ALPHA, ALPHA)])
    cases = ((0.10, 0.90), (0.20, 0.80), (0.45, 0.55))
    for z in cases:
        base = _flash(z)
        repeat = _flash(z)
        assert dict(repeat.diagnostics) == dict(base.diagnostics)
        assert _phase_set(repeat) == _phase_set(base)

        for order in itertools.permutations(range(len(NAMES))):
            permuted_names = tuple(NAMES[i] for i in order)
            permuted = ct.flash_tp(
                _mixture(permuted_names, tuple(z[i] for i in order)),
                temperature_K=TEMPERATURE_K,
                pressure_Pa=PRESSURE_PA,
                activity_model=ct.NRTL(parameters=parameters),
            )
            assert len(permuted.phase_names()) == len(base.phase_names())

            undo = np.argsort(order)
            restored = sorted(
                tuple(np.array(permuted.phases[name].composition.fractions)[undo].tolist())
                for name in permuted.phase_names()
            )
            reference = _phase_set(base)
            assert np.allclose(np.array(restored), np.array(reference), rtol=0.0, atol=1e-9)

            # Phase fractions match once the phases are paired by composition.
            for name in permuted.phase_names():
                composition = np.array(permuted.phases[name].composition.fractions)[undo]
                partner = min(
                    base.phase_names(),
                    key=lambda other: float(
                        np.max(
                            np.abs(np.array(base.phases[other].composition.fractions) - composition)
                        )
                    ),
                )
                assert permuted.phase_fractions[name] == pytest.approx(
                    base.phase_fractions[partner], abs=1e-9
                )


# --------------------------------------------------------------------------
# The second-order stage
# --------------------------------------------------------------------------


def test_successive_substitution_alone_cannot_solve_the_near_plait_feed(
    tessier2000_names: list[str], tessier2000_model: ct.NRTL
) -> None:
    """The second-order stage is required, not decorative.

    The Tessier (2000) Problem 1 feed z = (0.148, 0.052, 0.80) sits close to a
    plait point. With `second_order=False` the split never converges inside the
    default 100-iteration budget; measured, successive substitution alone needs
    3922 iterations from the same seed. With the stage enabled the same feed
    converges in 50 substitutions plus 5 second-order steps to an equal-activity
    residual of 4.4e-16.
    """
    feed = (0.148, 0.052, 0.80)
    mixture = _mixture(tessier2000_names, feed)

    with pytest.raises(ct.ConvergenceError):
        ct.flash_tp(
            mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=tessier2000_model,
            settings=ct.FlashSettings(second_order=False),
        )

    converged = ct.flash_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
    )
    assert converged.diagnostics["converged_stage"] == "second-order"
    assert int(converged.diagnostics["second_order_iterations"]) > 0
    assert float(converged.diagnostics["equilibrium_residual"]) < 1e-12

    # Successive substitution does get there, but only with a budget two orders
    # of magnitude larger than the default.
    patient = ct.flash_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=tessier2000_model,
        settings=ct.FlashSettings(second_order=False, max_iter=6000, ssi_iterations=6000),
    )
    assert int(patient.diagnostics["ssi_iterations"]) > 1000
    patient_phases = np.array(_phase_set(patient))
    converged_phases = np.array(_phase_set(converged))
    assert np.max(np.abs(patient_phases - converged_phases)) < 1e-6


# --------------------------------------------------------------------------
# Failure semantics
# --------------------------------------------------------------------------


def test_inconclusive_stability_raises_in_gamma_gamma() -> None:
    starved = ct.StabilitySettings(max_iter=1, second_order=False)
    mixture = _mixture(NAMES, (0.10, 0.90))
    assert (
        ct.stability_tp(
            mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=_model(),
            settings=starved,
        ).status
        == "inconclusive"
    )
    with pytest.raises(ct.ConvergenceError, match="inconclusive"):
        _flash((0.10, 0.90), ct.FlashSettings(stability_settings=starved))


def test_post_split_failure_raises_and_post_split_stability_false_returns() -> None:
    """The post-split gate, demonstrated on a deliberately under-converged split.

    **This case is synthetic.** No feed in this repository's validated grids
    genuinely needs a third phase, so the failure is manufactured by loosening
    the split tolerance to `tol = 1e-3`: the two converged phases are then only
    accurate to ~1e-3, and the tangent-plane test run on one of them sees a
    negative `tpd` far larger than `tpd_tol` at a point that is not close enough
    to the partner phase to be classified `"marginal"`. That is the same code
    path a genuine third phase would take, and it is exercised here because it
    is the only way to reach it today. A real third-phase state would look
    identical from `flash_tp`'s side.
    """
    mixture = _mixture(("Ethane", "n-Heptane"), (0.7, 0.3))
    loose = ct.FlashSettings(tol=1e-3)

    with pytest.raises(ct.ConvergenceError, match="not a stable phase set"):
        ct.flash_tp(mixture, temperature_K=360.0, pressure_Pa=1.0e6, eos=EOS, settings=loose)

    downgraded = ct.flash_tp(
        mixture,
        temperature_K=360.0,
        pressure_Pa=1.0e6,
        eos=EOS,
        settings=ct.FlashSettings(tol=1e-3, post_split_stability=False),
    )
    assert set(downgraded.phase_names()) == {"liquid", "vapor"}
    assert downgraded.diagnostics["post_split_checked"] is True
    assert downgraded.diagnostics["post_split_stable"] is False
    assert downgraded.diagnostics["post_split_status"] == "unstable"
    assert float(downgraded.diagnostics["post_split_tpd_min"]) < -1e-8

    # At the default tolerance the same state passes.
    tight = ct.flash_tp(mixture, temperature_K=360.0, pressure_Pa=1.0e6, eos=EOS)
    assert tight.diagnostics["post_split_stable"] is True


def test_a_converged_phase_that_finds_its_partner_is_marginal_not_unstable() -> None:
    """Validation Case S-3 as a flash-level guard.

    Two coexisting phases share one tangent plane, so a stability test on either
    finds the other with `tpd = 0`. At `tol = 1e-4` that zero is numerically
    -2.0e-05 for the canonical ternary - well past `tpd_tol` - and the check
    must still pass, because the minimizer *is* the partner phase.
    """
    result = ct.flash_tp(
        _mixture(("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2)),
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=EOS,
        settings=ct.FlashSettings(tol=1e-4),
    )
    assert result.diagnostics["post_split_stable"] is True
    assert "marginal" in {
        result.diagnostics["phase_stability_liquid"],
        result.diagnostics["phase_stability_vapor"],
    }
    assert float(result.diagnostics["post_split_tpd_min"]) < -1e-8


def test_gamma_phi_and_legacy_paths_declare_that_they_were_not_checked() -> None:
    components = tuple(ct.Component.from_database(name) for name in ("Methane", "Ethane"))
    mixture = ct.Mixture(
        components=components,
        composition=ct.Composition(fractions=(0.5, 0.5), basis="mole", normalize=False),
    )
    gamma_phi = ct.flash_tp(
        mixture,
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=EOS,
        activity_model=ct.NRTL(),
        flash_mode="gamma-phi",
    )
    assert gamma_phi.diagnostics["post_split_checked"] is False
    assert gamma_phi.diagnostics["post_split_skipped_reason"] == "gamma_phi_stability_unsupported"

    legacy = ct.flash_tp(
        _mixture(("Methane", "Ethane"), (0.5, 0.5)),
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=EOS,
        settings=ct.FlashSettings(phase_detection="wilson-heuristic"),
    )
    assert legacy.diagnostics["post_split_checked"] is False
    assert legacy.diagnostics["post_split_skipped_reason"] == "legacy_wilson_heuristic_path"


def test_new_settings_are_validated() -> None:
    assert ct.FlashSettings().post_split_stability is True
    assert ct.FlashSettings().second_order is True
    assert ct.FlashSettings().ssi_iterations == 50
    assert ct.FlashSettings().second_order_tol == 1e-12
    with pytest.raises(ct.InputRangeError, match="ssi_iterations"):
        ct.FlashSettings(ssi_iterations=0)
    with pytest.raises(ct.InputRangeError, match="second_order_max_iter"):
        ct.FlashSettings(second_order_max_iter=0)
    with pytest.raises(ct.InputRangeError, match="second_order_tol"):
        ct.FlashSettings(second_order_tol=0.0)
