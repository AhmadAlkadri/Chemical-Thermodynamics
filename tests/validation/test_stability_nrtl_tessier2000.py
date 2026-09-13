"""`stability_tp` against the published NRTL tangent-plane minima of Tessier (2000).

Source: S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable phase
stability analysis for excess Gibbs energy models", Chem. Eng. Sci. 55 (2000)
1785-1796 (author copy: https://academicweb.nd.edu/~markst/srt2000.pdf).

- Problem 1: n-propanol(1) / n-butanol(2) / water(3), parameters Table 1
  (attributed there to McDonald & Floudas, AIChE J. 41 (1995) 1798), stationary
  points Table 2. Fixture: `tests/fixtures/nrtl/tessier2000_problem1.json`.
- Problem 2: n-propanol(1) / n-butanol(2) / benzene(3) / water(4), parameters
  Table 4 (regressed from Gmehling et al., DECHEMA Chemistry Data Series
  1977-1990), stationary points Table 5. Fixture:
  `tests/fixtures/nrtl/tessier2000_problem2.json`.

Where `tests/validation/test_nrtl_tessier2000.py` checks the *model* (does our
`ln gamma` put the published stationary points where the paper says they are),
this module checks the *solver*: does `chemthermo.stability_tp`, given only the
feed and the model, find those points by itself and report the right verdict.

Two reference routes are used:

1. A from-scratch damped-Newton refinement of the printed stationary points,
   written here and sharing no code with `chemthermo.stability`. Its `D` values
   are the primary comparison target because the paper prints only five digits.
2. `thermo` 0.6.0 as an independent `ln gamma` implementation (skipped when the
   optional dependency is absent).
"""

from __future__ import annotations

from typing import Any, Callable

import numpy as np
import pytest

import chemthermo as ct

TEMPERATURE_K = 298.15  # Immaterial: both tables give dimensionless tau.
PRESSURE_PA = 101325.0  # Validated but inert for an activity model.

LnGamma = Callable[[np.ndarray], np.ndarray]

# Recomputed global minima of Problem 1, Table 2 (validation Case N-3). The
# paper's printed D is reproduced for two of these four and disputed for the
# other two; see the fixture's `dispute_note` fields.
PROBLEM1_GLOBAL_MINIMA: dict[tuple[float, ...], float] = {
    (0.148, 0.052, 0.80): -9.851037e-06,
    (0.12, 0.08, 0.80): -7.481797e-04,
    (0.13, 0.07, 0.80): -3.276225e-04,
    (0.12, 0.05, 0.83): -5.735988e-05,
}


def _ln_gamma_for(payload: dict[str, Any]) -> tuple[list[str], ct.NRTL, LnGamma]:
    """Return (names, model, ln_gamma) for a Tessier fixture payload."""
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


def _tangent_plane_distance(ln_gamma: LnGamma, z: np.ndarray, x: np.ndarray) -> float:
    d = np.log(z) + ln_gamma(z)
    return float(np.sum(x * (np.log(x) + ln_gamma(x) - d)))


def _refine(
    ln_gamma: LnGamma, z: np.ndarray, w0: np.ndarray, max_iter: int = 200
) -> tuple[np.ndarray, float, float]:
    """Damped Newton on the simplex-constrained stationarity system.

    Unknowns are (w, k); equations are ``ln w_i + ln gamma_i(w) - d_i - k = 0``
    plus ``sum_i w_i - 1 = 0``. The Jacobian is built by central differences, so
    the solve shares no derivative code with the model or with the package's
    stability solver. Returns (w, k, residual); at a stationary point D = k.
    """
    n = z.size
    d = np.log(z) + ln_gamma(z)

    def residuals(u: np.ndarray) -> np.ndarray:
        return np.concatenate(
            [np.log(u[:n]) + ln_gamma(u[:n]) - d - u[n], [float(np.sum(u[:n])) - 1.0]]
        )

    u = np.concatenate([np.asarray(w0, dtype=float) / float(np.sum(w0)), [0.0]])
    step = 1e-7
    for _ in range(max_iter):
        f = residuals(u)
        if float(np.max(np.abs(f))) < 1e-14:
            break
        jacobian = np.zeros((n + 1, n + 1), dtype=float)
        for column in range(n + 1):
            plus = u.copy()
            minus = u.copy()
            plus[column] += step
            minus[column] -= step
            jacobian[:, column] = (residuals(plus) - residuals(minus)) / (2.0 * step)
        delta = np.linalg.solve(jacobian, -f)
        scale = 1.0
        while np.any(u[:n] + scale * delta[:n] <= 0.0):
            scale *= 0.5
        u = u + scale * delta
    return u[:n], float(u[n]), float(np.max(np.abs(residuals(u))))


def _stability(model: ct.NRTL, names: list[str], z: tuple[float, ...]) -> ct.StabilityResult:
    mixture = ct.Mixture.from_database(names, list(z), normalize=True)
    return ct.stability_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
    )


# --------------------------------------------------------------------------
# Problem 2 fixture provenance
# --------------------------------------------------------------------------


def test_problem2_alpha_is_symmetric_and_reproduces_the_printed_G_matrix(
    tessier2000_problem2_payload: dict[str, Any],
) -> None:
    """Table 4 prints G_ij and tau_ij but not alpha_ij.

    alpha is implied by ``G_ij = exp(-alpha_ij tau_ij)``. The two directions of
    each pair must agree (alpha is physically symmetric) to the precision the
    five-decimal printed G supports, and re-exponentiating the rounded
    symmetric average must reproduce the printed G. Achieved: symmetry
    2.597e-05, G reconstruction 4.41e-06.
    """
    tau = np.asarray(tessier2000_problem2_payload["tau"], dtype=float)
    alpha = np.asarray(tessier2000_problem2_payload["alpha"], dtype=float)
    g_published = np.asarray(tessier2000_problem2_payload["G_published"], dtype=float)

    off_diagonal = ~np.eye(tau.shape[0], dtype=bool)
    recovered = np.zeros_like(tau)
    recovered[off_diagonal] = -np.log(g_published[off_diagonal]) / tau[off_diagonal]

    symmetry_error = float(np.max(np.abs(recovered - recovered.T)))
    assert symmetry_error < 1e-4
    assert tessier2000_problem2_payload["alpha_symmetry_max_abs_error"] == pytest.approx(
        symmetry_error, rel=1e-9
    )

    symmetric = 0.5 * (recovered + recovered.T)
    assert np.allclose(alpha[off_diagonal], np.round(symmetric, 3)[off_diagonal], atol=0.0)

    g_from_alpha = np.exp(-alpha * tau)
    np.fill_diagonal(g_from_alpha, 1.0)
    max_error = float(np.max(np.abs(g_from_alpha - g_published)))
    assert max_error < 1e-5
    assert tessier2000_problem2_payload["G_reconstruction_max_abs_error"] == pytest.approx(
        max_error, rel=1e-9
    )


def test_problem2_components_map_to_the_databank(
    tessier2000_problem2_names: list[str],
) -> None:
    assert tessier2000_problem2_names == ["1-Propanol", "n-Butanol", "Benzene", "Water"]
    for name in tessier2000_problem2_names:
        assert name in ct.list_component_names()


# --------------------------------------------------------------------------
# Acceptance A: Problem 1, Table 2
# --------------------------------------------------------------------------


def test_problem1_stability_tp_finds_the_table2_global_minima(
    tessier2000_payload: dict[str, Any],
) -> None:
    """Every Table 2 feed is unstable and `tpd_min` is that feed's lowest D.

    The near-plait-point feeds have |D| down to ~1e-05, which is three orders
    of magnitude above the default `tpd_tol` of 1e-08, so the verdict does not
    depend on the tolerance choice. The expected values are the recomputed
    minima of validation Case N-3, not the printed digit strings (two of which
    the ledger documents as typographical errors).
    """
    names, model, ln_gamma = _ln_gamma_for(tessier2000_payload)

    for feed, expected in PROBLEM1_GLOBAL_MINIMA.items():
        result = _stability(model, names, feed)

        assert result.status == "unstable", feed
        assert result.stable is False
        assert result.tpd_min == pytest.approx(expected, rel=1e-5), feed
        assert result.trial_composition is not None

        z = np.asarray(feed, dtype=float)
        w = np.asarray(result.trial_composition, dtype=float)

        # The reported point really is a stationary point of D, and D there is
        # the value reported.
        refined, k, residual = _refine(ln_gamma, z, w)
        assert residual < 1e-12, feed
        assert np.allclose(refined, w, atol=1e-6), feed
        assert _tangent_plane_distance(ln_gamma, z, w) == pytest.approx(result.tpd_min, abs=1e-12)
        assert k == pytest.approx(result.tpd_min, rel=1e-6, abs=1e-14)

        # Every trial is recorded and every one converged.
        assert all(trial.converged for trial in result.trials)
        assert len(result.trials) == len(names)


def test_problem1_reaches_every_printed_table2_point_except_the_plait_saddle(
    tessier2000_payload: dict[str, Any],
) -> None:
    """Honest record of which printed Table 2 points the trial set reaches.

    Seven of the eight non-trivial printed stationary points are reached by
    some trial. The one that is not is the near-plait-point *saddle* of the
    z1 = 0.148 feed (printed D = +4.5711e-08, recomputed +4.57105e-08): every
    trial that heads toward it is pulled into the trivial solution instead.
    That point has a positive D and therefore cannot change any verdict, but it
    is recorded here rather than quietly omitted.
    """
    names, model, ln_gamma = _ln_gamma_for(tessier2000_payload)

    reached: list[tuple[tuple[float, ...], float]] = []
    missed: list[tuple[tuple[float, ...], float]] = []

    for entry in tessier2000_payload["stationary_points"]:
        if entry["trivial"]:
            continue
        feed = tuple(float(value) for value in entry["feed"])
        z = np.asarray(feed, dtype=float)
        printed = np.asarray(entry["printed_composition"], dtype=float)
        printed = printed / float(np.sum(printed))
        target, _, residual = _refine(ln_gamma, z, printed)
        assert residual < 1e-12

        result = _stability(model, names, feed)
        found = any(
            trial.composition is not None
            and not trial.trivial
            and np.allclose(np.asarray(trial.composition, dtype=float), target, atol=1e-6)
            for trial in result.trials
        )
        (reached if found else missed).append((feed, float(entry["printed_D"])))

    assert len(reached) == 7
    assert missed == [((0.148, 0.052, 0.8), 4.5711e-08)]


# --------------------------------------------------------------------------
# Acceptance B: Problem 2, Table 5
# --------------------------------------------------------------------------


def test_problem2_stability_tp_reproduces_the_table5_verdicts_and_minima(
    tessier2000_problem2_payload: dict[str, Any],
) -> None:
    """Feed 2 is stable (all printed D >= 0); the other four are unstable.

    For each feed, `tpd_min` is compared against this module's own refinement
    of the corresponding printed point (relative 1e-6) and against the printed
    five-digit value. The printed comparison is asserted at rel 2.5e-04: five
    printed digits alone justify ~5e-05, and Table 4's alpha is only recoverable
    to the precision of the printed G matrix, which alone moves D by up to
    9.0e-05 relative. Achieved against the printed values:
    9.8e-06, 1.6e-05, 5.7e-05, 1.0e-04, 1.3e-04 relative. For the stable feed,
    `tpd_min` is the smallest *positive* stationary D, printed 0.03079.
    """
    names, model, ln_gamma = _ln_gamma_for(tessier2000_problem2_payload)

    by_feed: dict[tuple[float, ...], list[dict[str, Any]]] = {}
    for entry in tessier2000_problem2_payload["stationary_points"]:
        by_feed.setdefault(tuple(float(v) for v in entry["feed"]), []).append(entry)

    assert len(by_feed) == 5
    unstable_feeds = 0

    for feed, entries in by_feed.items():
        z = np.asarray(feed, dtype=float)
        target_entry = next(entry for entry in entries if entry["global_minimum"])
        printed_d = float(target_entry["printed_D"])
        printed_w = np.asarray(target_entry["printed_composition"], dtype=float)
        printed_w = printed_w / float(np.sum(printed_w))

        refined, k, residual = _refine(ln_gamma, z, printed_w)
        assert residual < 1e-12, feed
        refined_d = _tangent_plane_distance(ln_gamma, z, refined)
        assert refined_d == pytest.approx(k, rel=1e-9, abs=1e-15)
        # The refined point is still the printed one (three printed digits).
        assert np.allclose(refined, printed_w, atol=1e-3), feed
        # The refined D matches the printed D to the digits the paper prints.
        assert refined_d == pytest.approx(printed_d, rel=2.5e-4), feed

        result = _stability(model, names, feed)
        expected_status = "unstable" if printed_d < 0.0 else "stable"
        assert result.status == expected_status, feed
        unstable_feeds += expected_status == "unstable"

        assert result.tpd_min == pytest.approx(refined_d, rel=1e-6), feed
        assert result.trial_composition is not None
        assert np.allclose(np.asarray(result.trial_composition, dtype=float), refined, atol=1e-6), (
            feed
        )
        assert all(trial.converged for trial in result.trials)

    assert unstable_feeds == 4


def test_problem2_records_which_other_printed_points_are_reached(
    tessier2000_problem2_payload: dict[str, Any],
) -> None:
    """Honest record of the non-global printed Table 5 points.

    Reached: the benzene-rich points of the two z3-rich feeds (printed D
    -0.03365 and -3.1279e-03). Not reached: the three *positive* D points
    (+0.06532 of feed 2, +0.02268 of feed 4, +0.01066 of feed 5). Trials that
    head toward those are pulled into the trivial solution, exactly as for the
    Problem 1 saddle. All three have D > 0 and so cannot change a verdict.
    """
    names, model, ln_gamma = _ln_gamma_for(tessier2000_problem2_payload)

    reached: list[float] = []
    missed: list[float] = []

    for entry in tessier2000_problem2_payload["stationary_points"]:
        if entry["trivial"] or entry["global_minimum"]:
            continue
        feed = tuple(float(value) for value in entry["feed"])
        z = np.asarray(feed, dtype=float)
        printed = np.asarray(entry["printed_composition"], dtype=float)
        printed = printed / float(np.sum(printed))
        target, _, residual = _refine(ln_gamma, z, printed)
        assert residual < 1e-12

        result = _stability(model, names, feed)
        found = any(
            trial.composition is not None
            and not trial.trivial
            and np.allclose(np.asarray(trial.composition, dtype=float), target, atol=1e-6)
            for trial in result.trials
        )
        (reached if found else missed).append(float(entry["printed_D"]))

    assert sorted(reached) == [-0.03365, -0.0031279]
    assert sorted(missed) == [0.01066, 0.02268, 0.06532]


def test_problem2_disputed_printed_D_is_recorded_not_accommodated(
    tessier2000_problem2_payload: dict[str, Any],
) -> None:
    """One Table 5 D value does not reproduce; assert the mismatch, not the digits.

    Feed (0.25, 0.15, 0.40, 0.20), stationary point near
    (0.195, 0.0786, 0.114, 0.613): printed +0.02268, recomputed +0.0266799
    (17.6% apart) while the composition is reproduced to 3.4e-04 and every
    other point reproduces to 2.5e-04 relative or better. Treated as a
    typographical error; the assertion below fails if a future change ever
    "fixes" it, which would mean something else moved.
    """
    _, _, ln_gamma = _ln_gamma_for(tessier2000_problem2_payload)

    disputed = [
        entry
        for entry in tessier2000_problem2_payload["stationary_points"]
        if entry["printed_D_disputed"]
    ]
    assert len(disputed) == 1

    entry = disputed[0]
    z = np.asarray(entry["feed"], dtype=float)
    printed_w = np.asarray(entry["printed_composition"], dtype=float)
    printed_w = printed_w / float(np.sum(printed_w))
    refined, _, residual = _refine(ln_gamma, z, printed_w)
    assert residual < 1e-12
    assert np.allclose(refined, printed_w, atol=1e-3)

    refined_d = _tangent_plane_distance(ln_gamma, z, refined)
    assert refined_d == pytest.approx(0.0266799, rel=1e-5)
    assert abs(refined_d / float(entry["printed_D"]) - 1.0) > 1e-2


# --------------------------------------------------------------------------
# Acceptance E: independent cross-check against `thermo` 0.6.0
# --------------------------------------------------------------------------


@pytest.mark.parametrize("problem", ["problem1", "problem2"])
def test_stability_minima_agree_with_thermo_ln_gamma(
    problem: str,
    tessier2000_payload: dict[str, Any],
    tessier2000_problem2_payload: dict[str, Any],
) -> None:
    """Recompute D at chemthermo's minimizer using `thermo.NRTL` for ln gamma.

    Both sides evaluate the same closed-form Renon-Prausnitz equation, so the
    tangent-plane distance at the same composition must agree at round-off, not
    at a "physically reasonable" tolerance. Achieved: max |dD| = 4.7e-16 over
    the nine feeds.
    """
    thermo = pytest.importorskip("thermo")

    payload = tessier2000_payload if problem == "problem1" else tessier2000_problem2_payload
    names, model, _ = _ln_gamma_for(payload)
    tau = np.asarray(payload["tau"], dtype=float)
    alpha = np.asarray(payload["alpha"], dtype=float)
    size = tau.shape[0]

    def thermo_ln_gamma(x: np.ndarray) -> np.ndarray:
        values = np.asarray(x, dtype=float)
        values = values / float(np.sum(values))
        reference = thermo.NRTL(
            T=TEMPERATURE_K,
            xs=values.tolist(),
            tau_coeffs=[[[tau[i][j], 0, 0, 0, 0, 0] for j in range(size)] for i in range(size)],
            alpha_coeffs=[[[alpha[i][j], 0] for j in range(size)] for i in range(size)],
        )
        return np.asarray(reference.lngammas(), dtype=float)

    feeds = sorted({tuple(float(v) for v in e["feed"]) for e in payload["stationary_points"]})
    worst = 0.0
    for feed in feeds:
        result = _stability(model, names, feed)
        assert result.trial_composition is not None
        z = np.asarray(feed, dtype=float)
        w = np.asarray(result.trial_composition, dtype=float)
        reference_d = _tangent_plane_distance(thermo_ln_gamma, z, w)
        worst = max(worst, abs(reference_d - result.tpd_min))
        assert reference_d == pytest.approx(result.tpd_min, abs=1e-10), feed

    assert worst < 1e-12
