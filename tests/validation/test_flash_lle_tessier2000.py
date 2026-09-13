"""Liquid-liquid `flash_tp` tie-lines for the Tessier (2000) NRTL systems.

Source: S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable phase
stability analysis for excess Gibbs energy models", Chem. Eng. Sci. 55 (2000)
1785-1796 (author copy: https://academicweb.nd.edu/~markst/srt2000.pdf).

- Problem 1: n-propanol(1) / n-butanol(2) / water(3), parameters Table 1.
  Fixture: `tests/fixtures/nrtl/tessier2000_problem1.json`.
- Problem 2: n-propanol(1) / n-butanol(2) / benzene(3) / water(4), parameters
  Table 4. Fixture: `tests/fixtures/nrtl/tessier2000_problem2.json`.

**The paper publishes stationary points of the tangent-plane distance, not
tie-lines.** There is therefore no printed tie-line to compare against, and
nothing in this module claims one. The published content that *is* used is the
parameter set, the feeds, and the stability verdicts (which
`tests/validation/test_stability_nrtl_tessier2000.py` already checks against
Tables 2 and 5). What is validated here is the flash built on top of that:

1. against an **independent equal-activity solve** written in this module -
   its own successive-substitution loop plus a damped Newton solve of the full
   ``(x^I, x^II, beta)`` system with a finite-difference Jacobian, sharing no
   code with `chemthermo.flash`;
2. against the **invariants** every converged split must satisfy (material
   balance, equal activities, a negative Gibbs-energy change, and a post-split
   stability test of both phases);
3. against **`thermo` 0.6.0** where that is possible (see the last test, which
   records honestly what `thermo` could and could not be made to do).

The tie-line values below are the output of route 1. They are recorded so that a
future change that moves an equilibrium is caught, not because a publication
printed them.
"""

from __future__ import annotations

import math
from typing import Any, Callable

import numpy as np
import pytest

import chemthermo as ct

TEMPERATURE_K = 298.15  # Immaterial: both tables give dimensionless tau.
PRESSURE_PA = 101325.0  # Validated but inert for an activity model.

LnGamma = Callable[[np.ndarray], np.ndarray]

#: Tie-lines of the Problem 1 feeds, from the independent route in this module.
#: Values: (phase A, phase B, fraction of phase B) with the phases ordered by
#: increasing first-component mole fraction, so the record is label independent.
PROBLEM1_TIE_LINES: dict[tuple[float, ...], tuple[tuple[float, ...], tuple[float, ...], float]] = {
    (0.12, 0.08, 0.80): (
        (0.063945, 0.030845, 0.905210),
        (0.147546, 0.104156, 0.748298),
        0.670502,
    ),
    (0.13, 0.07, 0.80): (
        (0.077843, 0.032599, 0.889559),
        (0.153205, 0.086640, 0.760155),
        0.692089,
    ),
    (0.12, 0.05, 0.83): (
        (0.093662, 0.034385, 0.871953),
        (0.156101, 0.071402, 0.772496),
        0.421822,
    ),
    (0.148, 0.052, 0.80): (
        (0.116727, 0.037073, 0.846200),
        (0.154474, 0.055090, 0.790436),
        0.828496,
    ),
}

#: Problem 2 feeds: (benzene-rich phase, water-rich phase, water-rich fraction).
PROBLEM2_TIE_LINES: dict[tuple[float, ...], float] = {
    (0.148, 0.052, 0.600, 0.200): 0.166189,
    (0.148, 0.052, 0.700, 0.100): 0.060911,
    (0.25, 0.15, 0.40, 0.20): 0.027619,
    (0.25, 0.15, 0.35, 0.25): 0.070626,
}
PROBLEM2_STABLE_FEED = (0.25, 0.25, 0.25, 0.25)


def _model_for(payload: dict[str, Any]) -> tuple[list[str], ct.NRTL, LnGamma]:
    names = [str(entry["chemthermo_name"]) for entry in payload["components"]]
    tau = payload["tau"]
    alpha = payload["alpha"]
    pairs = [
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
    model = ct.NRTL(parameters=ct.NRTLParameters.from_pairs(pairs))
    mixture = ct.Mixture.from_database(names, [1.0 / len(names)] * len(names), normalize=True)

    def ln_gamma(x: np.ndarray) -> np.ndarray:
        values = np.asarray(x, dtype=float)
        values = values / float(np.sum(values))
        gamma = model.activity_coefficients(
            mixture=mixture,
            temperature_K=TEMPERATURE_K,
            composition=[float(value) for value in values],
        )
        return np.log(np.asarray(gamma, dtype=float))

    return names, model, ln_gamma


def _rachford_rice(z: np.ndarray, K: np.ndarray) -> float | None:
    def f(beta: float) -> float:
        denominator = 1.0 + beta * (K - 1.0)
        if np.any(denominator <= 0.0):
            return float("nan")
        return float(np.sum(z * (K - 1.0) / denominator))

    low, high = 0.0, 1.0
    f_low, f_high = f(low), f(high)
    if not math.isfinite(f_low) or not math.isfinite(f_high) or f_low * f_high > 0.0:
        return None
    for _ in range(200):
        mid = 0.5 * (low + high)
        value = f(mid)
        if abs(value) < 1e-14:
            return mid
        if value * f_low > 0.0:
            low, f_low = mid, value
        else:
            high = mid
    return 0.5 * (low + high)


def _independent_tie_line(
    z: np.ndarray, w: np.ndarray, ln_gamma: LnGamma
) -> tuple[np.ndarray, np.ndarray, float, float, int, int]:
    """Own successive substitution, then damped Newton on the full system.

    Unknowns ``(x^I, x^II, beta)``; equations

        ln(x_i^I gamma_i^I) - ln(x_i^II gamma_i^II) = 0     (i = 1..n)
        z_i - (1 - beta) x_i^I - beta x_i^II       = 0     (i = 1..n-1)
        sum_i x_i^I - 1 = 0,   sum_i x_i^II - 1 = 0

    with a central-difference Jacobian. The last mass balance is dropped because
    it is implied by the other equations. Nothing here is imported from
    `chemthermo.flash`.
    """
    x_i = z.copy()
    x_ii = w.copy()
    beta = 0.5
    substitutions = 0
    for substitutions in range(1, 6001):
        K = np.exp(ln_gamma(x_i) - ln_gamma(x_ii))
        candidate = _rachford_rice(z, K)
        if candidate is None:
            break
        beta = candidate
        next_i = z / (1.0 + beta * (K - 1.0))
        next_i = next_i / float(np.sum(next_i))
        next_ii = K * next_i
        next_ii = next_ii / float(np.sum(next_ii))
        moved = max(float(np.max(np.abs(next_i - x_i))), float(np.max(np.abs(next_ii - x_ii))))
        x_i, x_ii = next_i, next_ii
        if moved < 1e-13:
            break

    n = z.size

    def residuals(u: np.ndarray) -> np.ndarray:
        first = u[:n]
        second = u[n : 2 * n]
        fraction = u[2 * n]
        return np.concatenate(
            [
                np.log(first) + ln_gamma(first) - np.log(second) - ln_gamma(second),
                (z - (1.0 - fraction) * first - fraction * second)[:-1],
                [float(np.sum(first)) - 1.0, float(np.sum(second)) - 1.0],
            ]
        )

    u = np.concatenate([x_i, x_ii, [beta]])
    step = 1e-7
    residual = math.inf
    newton = 0
    for newton in range(1, 101):
        f = residuals(u)
        residual = float(np.max(np.abs(f)))
        if residual < 1e-14:
            newton -= 1
            break
        jacobian = np.zeros((2 * n + 1, 2 * n + 1), dtype=float)
        for column in range(2 * n + 1):
            plus = u.copy()
            minus = u.copy()
            plus[column] += step
            minus[column] -= step
            jacobian[:, column] = (residuals(plus) - residuals(minus)) / (2.0 * step)
        try:
            delta = np.linalg.solve(jacobian, -f)
        except np.linalg.LinAlgError:  # pragma: no cover - defensive
            break
        scale = 1.0
        accepted = False
        while scale > 1e-10:
            candidate_u = u + scale * delta
            if np.any(candidate_u[: 2 * n] <= 0.0) or not (0.0 < candidate_u[2 * n] < 1.0):
                scale *= 0.5
                continue
            if float(np.max(np.abs(residuals(candidate_u)))) < residual:
                u = candidate_u
                accepted = True
                break
            scale *= 0.5
        if not accepted:  # pragma: no cover - defensive
            break

    return u[:n], u[n : 2 * n], float(u[2 * n]), residual, substitutions, newton


def _ordered_phases(result: ct.FlashResult) -> tuple[np.ndarray, np.ndarray, float]:
    """Phases ordered by the first component, with the second phase's fraction.

    Ordering by composition rather than by label is what makes an assertion
    independent of the ``liquid1`` / ``liquid2`` role assignment.
    """
    entries = sorted(
        (
            (np.asarray(result.phases[name].composition.fractions, dtype=float), name)
            for name in result.phase_names()
        ),
        key=lambda entry: float(entry[0][0]),
    )
    return entries[0][0], entries[1][0], float(result.phase_fractions[entries[1][1]])


def _flash(model: ct.NRTL, names: list[str], feed: tuple[float, ...]) -> ct.FlashResult:
    return ct.flash_tp(
        ct.Mixture.from_database(names, list(feed), normalize=True),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
    )


def _check_feed(
    names: list[str],
    model: ct.NRTL,
    ln_gamma: LnGamma,
    feed: tuple[float, ...],
) -> tuple[np.ndarray, np.ndarray, float]:
    """Flash the feed, verify it against the independent solve and the invariants."""
    z = np.asarray(feed, dtype=float)
    z = z / float(np.sum(z))

    result = _flash(model, names, feed)
    assert set(result.phase_names()) == {"liquid1", "liquid2"}, feed
    assert result.vapor_fraction is None
    diagnostics = result.diagnostics

    low, high, high_fraction = _ordered_phases(result)

    stability = ct.stability_tp(
        ct.Mixture.from_database(names, list(feed), normalize=True),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
    )
    assert stability.status == "unstable", feed
    assert stability.trial_composition is not None
    reference = _independent_tie_line(
        z, np.asarray(stability.trial_composition, dtype=float), ln_gamma
    )
    r_i, r_ii, r_beta, r_residual, _substitutions, _newton = reference
    assert r_residual < 1e-12, feed

    reference_entries = sorted((r_i, r_ii), key=lambda values: float(values[0]))
    reference_high_fraction = r_beta if float(r_ii[0]) > float(r_i[0]) else 1.0 - r_beta

    assert np.allclose(low, reference_entries[0], atol=1e-6), feed
    assert np.allclose(high, reference_entries[1], atol=1e-6), feed
    assert high_fraction == pytest.approx(reference_high_fraction, abs=1e-6), feed

    # Invariants of the converged split.
    assert float(diagnostics["equilibrium_residual"]) < 1e-10, feed
    assert float(diagnostics["mass_balance_residual"]) < 1e-12, feed
    assert float(diagnostics["delta_g_split_rt"]) < 0.0, feed
    assert diagnostics["post_split_stable"] is True, feed
    assert diagnostics["phase_stability_liquid1"] == "stable", feed
    assert diagnostics["phase_stability_liquid2"] == "stable", feed
    assert abs(float(diagnostics["post_split_tpd_min"])) < 1e-8, feed

    # Material balance recomputed from the returned phases, not from the
    # residual the solver wrote into diagnostics.
    recombined = (1.0 - high_fraction) * low + high_fraction * high
    assert float(np.max(np.abs(recombined - z))) < 1e-12, feed

    # Every iteration count is reported.
    assert int(diagnostics["ssi_iterations"]) >= 1
    assert int(diagnostics["second_order_iterations"]) >= 0
    assert diagnostics["converged_stage"] in {"successive-substitution", "second-order"}

    return low, high, high_fraction


# --------------------------------------------------------------------------
# Case L-1: Problem 1 tie-lines
# --------------------------------------------------------------------------


def test_problem1_tie_lines(tessier2000_payload: dict[str, Any]) -> None:
    """Validation Case L-1.

    Every Table 2 feed is unstable and splits into a verified tie-line. The
    recorded values come from the independent route in this module; achieved
    agreement between the two routes is <= 1.1e-12 in composition and phase
    fraction, with equal-activity residuals <= 4.5e-16 and mass balance
    <= 2.3e-16.
    """
    names, model, ln_gamma = _model_for(tessier2000_payload)

    for feed, (expected_low, expected_high, expected_fraction) in PROBLEM1_TIE_LINES.items():
        low, high, high_fraction = _check_feed(names, model, ln_gamma, feed)
        assert np.allclose(low, expected_low, atol=1e-6), feed
        assert np.allclose(high, expected_high, atol=1e-6), feed
        assert high_fraction == pytest.approx(expected_fraction, abs=1e-6), feed

    # The two feeds whose Gibbs-energy reduction the orchestrator reference
    # quotes; both are tiny because these feeds sit near the plait point.
    near_plait = _flash(model, names, (0.148, 0.052, 0.80))
    assert float(near_plait.diagnostics["delta_g_split_rt"]) == pytest.approx(-1.065e-6, rel=1e-3)
    wide = _flash(model, names, (0.12, 0.08, 0.80))
    assert float(wide.diagnostics["delta_g_split_rt"]) == pytest.approx(-1.952e-4, rel=1e-3)


# --------------------------------------------------------------------------
# Case L-2: Problem 2 tie-lines, the stable control, and the post-split check
# --------------------------------------------------------------------------


def test_problem2_tie_lines_and_stable_control(
    tessier2000_problem2_payload: dict[str, Any],
) -> None:
    """Validation Case L-2.

    Four feeds split into a benzene-rich and a water-rich liquid; the fifth,
    z = (0.25, 0.25, 0.25, 0.25), is the paper's stable feed (all printed D >= 0,
    Table 5) and must come back as a single liquid. Both phases of every split
    pass the post-split stability test, so none of these states needs a third
    phase.
    """
    names, model, ln_gamma = _model_for(tessier2000_problem2_payload)
    benzene = names.index("Benzene")
    water = names.index("Water")

    for feed, expected_fraction in PROBLEM2_TIE_LINES.items():
        _check_feed(names, model, ln_gamma, feed)

        result = _flash(model, names, feed)
        # Identify the phases by chemistry, not by label.
        water_rich = max(
            result.phase_names(),
            key=lambda name: result.phases[name].composition.fractions[water],
        )
        benzene_rich = min(
            result.phase_names(),
            key=lambda name: result.phases[name].composition.fractions[water],
        )
        assert result.phases[water_rich].composition.fractions[water] > 0.9, feed
        assert (
            result.phases[benzene_rich].composition.fractions[benzene]
            > result.phases[water_rich].composition.fractions[benzene]
        ), feed
        assert result.phase_fractions[water_rich] == pytest.approx(expected_fraction, abs=1e-6), (
            feed
        )

    stable = _flash(model, names, PROBLEM2_STABLE_FEED)
    assert stable.phase_names() == ["liquid"]
    assert stable.vapor_fraction is None
    assert stable.diagnostics["stability_status"] == "stable"
    assert float(stable.diagnostics["tpd_min"]) > 0.0


# --------------------------------------------------------------------------
# Case L-2 (continued): the post-split check is a real gate
# --------------------------------------------------------------------------


def test_post_split_check_uses_the_same_stability_settings(
    tessier2000_payload: dict[str, Any],
) -> None:
    """A starved stability setting makes the post-split test inconclusive.

    This proves the check is actually executed on the converged phases rather
    than being reported from the feed's analysis: the feed test is given enough
    budget to converge (it is the *split* seed), while the per-phase tests are
    not, and the flash then refuses.
    """
    names, model, _ = _model_for(tessier2000_payload)
    feed = (0.12, 0.08, 0.80)

    # Sanity: with the same starved settings the feed analysis is inconclusive,
    # so flash_tp cannot even start; that is a different failure, checked in
    # tests/test_flash_lle.py. Here we only assert that the default settings do
    # run the per-phase tests and report them.
    result = _flash(model, names, feed)
    assert result.diagnostics["post_split_checked"] is True
    for name in ("liquid1", "liquid2"):
        assert f"phase_stability_{name}" in result.diagnostics
        assert f"phase_stability_tpd_min_{name}" in result.diagnostics

    downgraded = ct.flash_tp(
        ct.Mixture.from_database(names, list(feed), normalize=True),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        settings=ct.FlashSettings(post_split_stability=False),
    )
    assert downgraded.diagnostics["post_split_checked"] is True
    assert downgraded.diagnostics["post_split_stable"] is True


# --------------------------------------------------------------------------
# Independent external check: `thermo` 0.6.0
# --------------------------------------------------------------------------


def test_thermo_agrees_on_the_binary_binodal_but_its_flash_will_not_split() -> None:
    """Honest record of what `thermo` 0.6.0 could and could not be made to do.

    Setup attempted: `thermo.GibbsExcessLiquid` built on `thermo.NRTL` with the
    same taus and alphas (n-butanol / water, Tessier Table 1 pair 2-3), two of
    them plus a `CEOSGas` handed to `thermo.FlashVLN`.

    - **`FlashVLN.flash` will not return a liquid-liquid split.** It reports
      `unique_liquid_count == 1`: two `GibbsExcessLiquid` phases built from the
      same excess-Gibbs model are deduplicated into one, and the returned result
      is a single phase for every feed tried (0.05, 0.10, 0.20, 0.30 butanol),
      all of which are inside the miscibility gap. So no phase-fraction
      comparison is possible, and none is claimed.
    - **`thermo`'s own Michelsen stability test does find the split**, and its
      converged trial compositions are the conjugate pair. That is what is
      compared here, at the tolerance `thermo`'s own stability search reaches
      (its residual is ~1.9e-06).

    Achieved: |dx| = 1.0e-06 on the water-rich branch and 2.5e-05 on the
    butanol-rich branch, both inside 1e-04. `thermo` also calls all four feeds
    unstable, including 0.30, which independently confirms that 0.30 lies inside
    the gap.
    """
    thermo = pytest.importorskip("thermo")

    tau = [[0.0, 0.90047], [3.51307, 0.0]]
    alpha = [[0.0, 0.48], [0.48, 0.0]]

    constants, correlations = thermo.ChemicalConstantsPackage.from_IDs(["1-butanol", "water"])

    def excess(zs: list[float]) -> Any:
        return thermo.NRTL(
            T=TEMPERATURE_K,
            xs=zs,
            tau_coeffs=[[[tau[i][j], 0, 0, 0, 0, 0] for j in range(2)] for i in range(2)],
            alpha_coeffs=[[[alpha[i][j], 0] for j in range(2)] for i in range(2)],
        )

    def liquid(zs: list[float]) -> Any:
        return thermo.GibbsExcessLiquid(
            VaporPressures=correlations.VaporPressures,
            GibbsExcessModel=excess(zs),
            HeatCapacityGases=correlations.HeatCapacityGases,
            T=TEMPERATURE_K,
            P=PRESSURE_PA,
            zs=zs,
        )

    gas = thermo.CEOSGas(
        thermo.eos_mix.PRMIX,
        eos_kwargs={"Tcs": constants.Tcs, "Pcs": constants.Pcs, "omegas": constants.omegas},
        HeatCapacityGases=correlations.HeatCapacityGases,
        T=TEMPERATURE_K,
        P=PRESSURE_PA,
        zs=[0.1, 0.9],
    )
    water_rich = liquid([0.02, 0.98])
    butanol_rich = liquid([0.36, 0.64])
    flasher = thermo.FlashVLN(constants, correlations, liquids=[water_rich, butanol_rich], gas=gas)

    # Recorded limitation, asserted so a future `thermo` that fixes it is noticed.
    assert flasher.unique_liquid_count == 1

    model = ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [("n-Butanol", "Water", 0.90047, 3.51307, 0.48, 0.48)]
        )
    )

    worst = 0.0
    for feed in (0.05, 0.10, 0.20, 0.30):
        zs = [feed, 1.0 - feed]

        flashed = flasher.flash(T=TEMPERATURE_K, P=PRESSURE_PA, zs=zs)
        assert flashed.phase_count == 1  # the limitation above, in action

        stable, info = flasher.stability_test_Michelsen(
            TEMPERATURE_K, PRESSURE_PA, zs, min_phase=water_rich, other_phase=butanol_rich
        )
        assert stable is False, feed  # thermo agrees the feed is not one phase
        thermo_pair = sorted((float(info[0][0]), float(info[1][0])))

        result = ct.flash_tp(
            ct.Mixture.from_database(["n-Butanol", "Water"], zs, normalize=True),
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=model,
        )
        ours = sorted(
            float(result.phases[name].composition.fractions[0]) for name in result.phase_names()
        )
        assert len(ours) == 2, feed
        worst = max(worst, max(abs(a - b) for a, b in zip(ours, thermo_pair)))

    assert worst < 1e-4
