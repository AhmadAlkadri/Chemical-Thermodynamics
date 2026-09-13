"""The phi-phi negative-flash / second-order split of ADR-0016.

The defect this pins: PC-SAFT carbon dioxide / n-decane, ``z = (0.9, 0.1)``, at
240 K and 1.0 MPa. The tangent-plane test proves the feed unstable, the seeded
K-values bracket a Rachford-Rice root, and then successive substitution
oscillates: ``K_CO2`` crosses 1 repeatedly, the root of ``f`` leaves ``[0, 1]``
and at the ninth iterate there is no root at all. Before ADR-0016 the split
raised `ConvergenceError`.

**Nothing here is compared against a number copied from the implementation.**
The reference split is recomputed in this file by a damped Newton on the
equal-fugacity system in vapor mole numbers - a different formulation, with its
own Jacobian - and the density roots, the Gibbs-energy change and the
post-split stability verdict are recomputed from the public API.
"""

from __future__ import annotations

import numpy as np
import pytest

import chemthermo as ct

COMPONENTS = ("Carbon dioxide", "n-Decane")
FEED = (0.9, 0.1)

#: The three states the orchestrator's grid scan found failing, all on the
#: tangent-plane path and all on the legacy path.
FAILING_STATES = ((240.0, 1.0e6), (250.0, 1.0e6), (260.0, 1.5e6))

REFERENCE_STATE = (240.0, 1.0e6)


def _mixture(fractions: tuple[float, ...] = FEED) -> ct.Mixture:
    return ct.Mixture.from_database(list(COMPONENTS), list(fractions), normalize=True)


def _ln_phi(
    mixture: ct.Mixture,
    temperature_K: float,
    pressure_Pa: float,
    composition: np.ndarray,
    phase: str,
) -> np.ndarray:
    return np.log(
        np.array(
            ct.PCSAFTEOS().fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                composition=(composition / composition.sum()).tolist(),
                phase=phase,
            ),
            dtype=float,
        )
    )


def _independent_split(
    mixture: ct.Mixture,
    temperature_K: float,
    pressure_Pa: float,
    z: np.ndarray,
    vapor_seed: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Damped Newton on ``ln f_i^V - ln f_i^L = 0`` in vapor mole numbers.

    Unknowns are the vapor mole numbers ``v_i`` per mole of feed, with the
    liquid holding ``z_i - v_i``; the vapor fraction is ``sum_i v_i``. This is
    the equal-fugacity *system*, not the Gibbs minimization the solver under
    test performs, and the Jacobian here is a finite difference of this
    residual. Returns ``(x, y, beta, residual)``.
    """

    def residual(v: np.ndarray) -> np.ndarray:
        liquid = z - v
        x = liquid / liquid.sum()
        y = v / v.sum()
        return (
            np.log(y)
            + _ln_phi(mixture, temperature_K, pressure_Pa, y, "vapor")
            - np.log(x)
            - _ln_phi(mixture, temperature_K, pressure_Pa, x, "liquid")
        )

    v = np.clip(vapor_seed, 1e-14, z - 1e-14)
    size = z.size
    for _ in range(300):
        r = residual(v)
        if np.max(np.abs(r)) < 1e-13:
            break
        jacobian = np.zeros((size, size))
        for column in range(size):
            step = 1e-8 * max(abs(v[column]), 1e-8)
            plus, minus = v.copy(), v.copy()
            plus[column] += step
            minus[column] -= step
            jacobian[:, column] = (residual(plus) - residual(minus)) / (2.0 * step)
        direction = np.linalg.solve(jacobian, -r)
        scale = 1.0
        current = float(np.max(np.abs(r)))
        while scale > 1e-13:
            candidate = v + scale * direction
            if np.all(candidate > 0.0) and np.all(candidate < z):
                if float(np.max(np.abs(residual(candidate)))) < current:
                    break
            scale *= 0.5
        else:  # pragma: no cover - the seeds below all make progress
            break
        v = v + scale * direction

    liquid = z - v
    return (
        liquid / liquid.sum(),
        v / v.sum(),
        float(v.sum()),
        float(np.max(np.abs(residual(v)))),
    )


def _reduced_g(fractions: np.ndarray, ln_f: np.ndarray) -> float:
    mask = fractions > 0.0
    return float(np.sum(fractions[mask] * (np.log(fractions[mask]) + ln_f[mask])))


# --------------------------------------------------------------------------
# A. The three failing states converge, and to the right answer.
# --------------------------------------------------------------------------


@pytest.mark.parametrize(("temperature_K", "pressure_Pa"), FAILING_STATES)
def test_the_previously_failing_states_return_a_verified_two_phase_result(
    temperature_K: float, pressure_Pa: float
) -> None:
    mixture = _mixture()
    result = ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=ct.PCSAFTEOS(),
    )
    diagnostics = result.diagnostics

    assert sorted(result.phase_names()) == ["liquid", "vapor"]
    assert diagnostics["stability_status"] == "unstable"
    assert diagnostics["phase_count"] == 2
    beta = result.vapor_fraction
    assert beta is not None and 0.0 < beta < 1.0

    # The invariants every split must carry (brain.md: "a phase split must be
    # verified, not just converged").
    assert float(diagnostics["mass_balance_residual"]) < 1e-12
    assert float(diagnostics["fugacity_residual"]) < 1e-10
    assert float(diagnostics["delta_g_split_rt"]) < 0.0
    assert diagnostics["post_split_status"] == "stable"

    # The ADR-0016 stage keys are present exactly because the stage ran.
    assert diagnostics["converged_stage"] == "second-order"
    assert int(diagnostics["negative_flash_steps"]) > 0

    # Independently: mass balance and equal fugacities, recomputed here.
    z = np.array(mixture.composition.fractions, dtype=float)
    x = np.array(result.phases["liquid"].composition.fractions, dtype=float)
    y = np.array(result.phases["vapor"].composition.fractions, dtype=float)
    assert np.max(np.abs(z - (beta * y + (1.0 - beta) * x))) < 1e-12
    ln_phi_x = _ln_phi(mixture, temperature_K, pressure_Pa, x, "liquid")
    ln_phi_y = _ln_phi(mixture, temperature_K, pressure_Pa, y, "vapor")
    assert np.max(np.abs((np.log(x) + ln_phi_x) - (np.log(y) + ln_phi_y))) < 1e-10


def test_the_reference_state_matches_an_independent_newton_to_1e_8() -> None:
    """Validation Case F-4, the 240 K / 1.0 MPa reference split."""
    temperature_K, pressure_Pa = REFERENCE_STATE
    mixture = _mixture()
    z = np.array(mixture.composition.fractions, dtype=float)

    result = ct.flash_tp(
        mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=ct.PCSAFTEOS()
    )
    beta = result.vapor_fraction
    assert beta is not None
    x = np.array(result.phases["liquid"].composition.fractions, dtype=float)
    y = np.array(result.phases["vapor"].composition.fractions, dtype=float)

    # Two different seeds, so the reference is not a restatement of the
    # solver's own starting point.
    for seed_beta in (0.05, 0.30):
        reference = _independent_split(
            mixture,
            temperature_K,
            pressure_Pa,
            z,
            np.array([seed_beta * (1.0 - 1e-7), seed_beta * 1e-7]),
        )
        reference_x, reference_y, reference_beta, reference_residual = reference
        assert reference_residual < 1e-12
        assert 0.0 < reference_beta < 1.0
        assert beta == pytest.approx(reference_beta, abs=1e-8)
        assert x == pytest.approx(reference_x, abs=1e-8)
        assert y == pytest.approx(reference_y, abs=1e-8)

    # The Gibbs-energy reduction, recomputed from the public API.
    ln_phi_x = _ln_phi(mixture, temperature_K, pressure_Pa, x, "liquid")
    ln_phi_y = _ln_phi(mixture, temperature_K, pressure_Pa, y, "vapor")
    ln_phi_feed_liquid = _ln_phi(mixture, temperature_K, pressure_Pa, z, "liquid")
    ln_phi_feed_vapor = _ln_phi(mixture, temperature_K, pressure_Pa, z, "vapor")
    feed_g = min(_reduced_g(z, ln_phi_feed_liquid), _reduced_g(z, ln_phi_feed_vapor))
    delta_g = beta * _reduced_g(y, ln_phi_y) + (1.0 - beta) * _reduced_g(x, ln_phi_x) - feed_g
    assert delta_g < 0.0
    assert delta_g == pytest.approx(float(result.diagnostics["delta_g_split_rt"]), abs=1e-12)


def test_the_reference_state_phase_densities_are_the_pinned_roots() -> None:
    """``density_roots`` at the converged phases (Case F-4, ADR-0015).

    Both phases have two mechanically stable roots at this ``(T, P)``; the
    liquid phase is on the highest and the vapor phase on the lowest, which is
    the root each `phase=` label selects.
    """
    temperature_K, pressure_Pa = REFERENCE_STATE
    mixture = _mixture()
    eos = ct.PCSAFTEOS()
    result = ct.flash_tp(
        mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=ct.PCSAFTEOS()
    )

    liquid_roots = eos.density_roots(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        composition=list(result.phases["liquid"].composition.fractions),
        mixture=mixture,
    )
    vapor_roots = eos.density_roots(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        composition=list(result.phases["vapor"].composition.fractions),
        mixture=mixture,
    )
    assert len(liquid_roots) == 2 and len(vapor_roots) == 2
    assert liquid_roots[0] == pytest.approx(699.9787395, rel=1e-6)
    assert liquid_roots[1] == pytest.approx(17250.70876, rel=1e-6)
    assert vapor_roots[0] == pytest.approx(557.8191800, rel=1e-6)
    assert vapor_roots[1] == pytest.approx(24161.93786, rel=1e-6)

    # And the molar volumes the labels select are the reciprocals of the ends.
    for phase, roots in (("liquid", liquid_roots), ("vapor", vapor_roots)):
        volume = eos.molar_volume(
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            composition=list(result.phases[phase].composition.fractions),
            phase=phase,
            mixture=mixture,
        )
        expected = roots[-1] if phase == "liquid" else roots[0]
        assert volume == pytest.approx(1.0 / expected, rel=1e-12)


# --------------------------------------------------------------------------
# The two mechanisms, separately.
# --------------------------------------------------------------------------


def test_the_split_really_does_leave_the_unit_interval_on_the_way() -> None:
    """The negative flash is used, not merely available.

    Reproduced here by re-running the successive-substitution update from the
    stability seed: the Rachford-Rice root is negative at several iterates.
    """
    temperature_K, pressure_Pa = REFERENCE_STATE
    mixture = _mixture()
    z = np.array(mixture.composition.fractions, dtype=float)
    eos = ct.PCSAFTEOS()

    stability = ct.stability_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=ct.PCSAFTEOS(),
    )
    assert stability.status == "unstable"
    w = np.array(stability.trial_composition, dtype=float)
    K = (w * np.exp(-float(stability.tpd_min))) / z

    below_zero = 0
    none_at_all = 0
    beta = 0.3
    for _ in range(12):
        t = 1.0 + beta * (K - 1.0)
        x = z / t
        x = x / x.sum()
        y = K * x
        y = y / y.sum()
        K = np.array(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                composition=x.tolist(),
                phase="liquid",
            )
        ) / np.array(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                composition=y.tolist(),
                phase="vapor",
            )
        )
        if not (float(np.max(K)) > 1.0 > float(np.min(K))):
            none_at_all += 1
            break
        # The window's root, by bisection written here.
        lower, upper = 1.0 / (1.0 - float(np.max(K))), 1.0 / (1.0 - float(np.min(K)))
        low, high = lower + 1e-12 * (upper - lower), upper - 1e-12 * (upper - lower)
        for _ in range(300):
            beta = 0.5 * (low + high)
            if float(np.sum(z * (K - 1.0) / (1.0 + beta * (K - 1.0)))) > 0.0:
                low = beta
            else:
                high = beta
        if beta < 0.0:
            below_zero += 1

    assert below_zero >= 2, "the iterates must actually go negative"
    assert none_at_all == 1, "and then run out of a window entirely"


def test_disabling_the_second_order_stage_restores_the_old_failure() -> None:
    """``second_order=False`` is the pre-ADR-0016 phi-phi split."""
    mixture = _mixture()
    for temperature_K, pressure_Pa in FAILING_STATES:
        with pytest.raises(ct.ConvergenceError):
            ct.flash_tp(
                mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                eos=ct.PCSAFTEOS(),
                settings=ct.FlashSettings(second_order=False),
            )


def test_the_legacy_path_is_deliberately_left_failing() -> None:
    """ADR-0016 decision 8: `wilson-heuristic` reproduces the old behavior."""
    mixture = _mixture()
    for temperature_K, pressure_Pa in FAILING_STATES:
        with pytest.raises(ct.ConvergenceError, match="failed to bracket"):
            ct.flash_tp(
                mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                eos=ct.PCSAFTEOS(),
                settings=ct.FlashSettings(phase_detection="wilson-heuristic"),
            )


def test_a_converged_first_stage_carries_none_of_the_new_diagnostics_keys() -> None:
    """ADR-0016 decision 6: the keys are conditional, deliberately."""
    result = ct.flash_tp(
        ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True),
        temperature_K=240.0,
        pressure_Pa=3.0e6,
        eos=ct.PengRobinsonEOS(),
    )
    for key in (
        "ssi_iterations",
        "second_order_iterations",
        "converged_stage",
        "negative_flash_steps",
    ):
        assert key not in result.diagnostics


# --------------------------------------------------------------------------
# D. Error semantics: a converged beta outside (0, 1) is not an answer.
# --------------------------------------------------------------------------


def test_a_converged_vapor_fraction_outside_the_unit_interval_raises(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """SYNTHETIC state: no state in this repository reaches it naturally.

    A negative flash is a legitimate *iterate* and never a legitimate *answer*:
    it says the phase set collapsed to one phase, which contradicts the
    tangent-plane verdict that started the split. Every physical ingredient
    here is real - the feed is genuinely unstable under Peng-Robinson and the
    stability test runs untouched - and only the *converged split* is replaced,
    because a model that produces this outcome while still being unstable at
    the feed would itself be fiction. Marked synthetic for that reason.
    """
    from chemthermo.flash import _detect
    from chemthermo.flash._split import _SplitSolution

    mixture = ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True)
    z = np.array(mixture.composition.fractions, dtype=float)
    stability = ct.stability_tp(
        mixture, temperature_K=240.0, pressure_Pa=3.0e6, eos=ct.PengRobinsonEOS()
    )
    assert stability.status == "unstable", "the premise: the feed really is two-phase"

    def collapsed(*args: object, **kwargs: object) -> _SplitSolution:
        return _SplitSolution(
            x=z.copy(),
            y=z.copy(),
            vapor_fraction=-0.2,
            K=np.ones_like(z),
            ln_f_x=np.zeros_like(z),
            ln_f_y=np.zeros_like(z),
            iterations=3,
            max_delta=0.0,
            converged=True,
            negative_flash_steps=3,
        )

    monkeypatch.setattr(_detect, "_solve_k_loop", collapsed)

    with pytest.raises(ct.ConvergenceError, match="outside \\(0, 1\\)") as excinfo:
        ct.flash_tp(mixture, temperature_K=240.0, pressure_Pa=3.0e6, eos=ct.PengRobinsonEOS())
    message = str(excinfo.value)
    assert "beta=-2.000000e-01" in message
    assert "tpd_min" in message
    assert "not a single-phase state" in message
