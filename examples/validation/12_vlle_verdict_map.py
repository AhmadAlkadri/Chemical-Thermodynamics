"""The verdict map of the ternary VLLE region: how many phases, feed by feed.

What is being validated
-----------------------
`flash_tp` is asked, one feed at a time, *how many phases* the model has at
1-propanol / n-butanol / water, 1 atm, 363 / 364 / 365 K - and every answer is
adjudicated against a state worked out independently in this script.

The independent route builds, for each feed, every state the model admits:

1. the **three-phase** state: the tie-triangle from a six-equation damped
   Newton solve (three equal activities, two normalizations, and the bubble
   condition on *one* liquid), admissible when the feed's barycentric weights
   in that triangle are all positive;
2. a **vapor-liquid** state from a four-equation Newton solve of
   ``z_i - (1 - b) x_i - b K_i(x) x_i = 0``, ``sum_i x_i = 1``, with
   ``K_i = gamma_i(x) Psat_i / P``;
3. a **liquid-liquid** state from a seven-equation Newton solve of the three
   equal activities, the two normalizations and two mass balances;
4. the **single-phase** state, the feed itself on whichever of the two
   candidates (NRTL liquid, ideal vapor) has the lower Gibbs energy.

Each admissible state is scored with the reduced molar Gibbs energy

    G/RT = sum_j beta_j sum_i x_i^j [ ln x_i^j + t_i^j ]

with ``t_i = ln gamma_i + ln(Psat_i/P)`` for a liquid and ``0`` for the vapor,
and the **lowest-Gibbs one is the expected answer**. Equal fugacities alone do
not pick a state out - more than one of these satisfies them - so the Gibbs
energy is what adjudicates.

The script then prints a confusion matrix of *expected* against *obtained*
phase counts over a deterministic grid of 75-76 feeds per temperature (21
strictly inside the triangle on a barycentric lattice, the rest on a 1/12
mole-fraction lattice outside it), lists every disagreement with its numbers,
and checks that

- every three-phase answer returns the *same* tie-triangle, whatever the feed,
  and phase fractions equal to the feed's barycentric weights;
- every two-phase answer has equal ``ln x_i + t_i`` across its phases and adds
  back to the feed;
- every single-phase answer is the feed itself, on the lower-Gibbs candidate.

Why this map exists
-------------------
Validation Case V-2 recorded a feed inside the 363 K tie-triangle that the
deterministic stability trial set missed, so `flash_tp` returned a single
liquid - wrong for this model. ADR-0012 fixed it by running each stability
trial on **one fixed phase-candidate surface** instead of re-selecting the
lowest-Gibbs candidate at every iteration. That changes which stationary points
are reachable, so the honest check is not the one repaired feed but the whole
map.

Parameters and their provenance
-------------------------------
NRTL Table 1 of S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, Chem. Eng.
Sci. 55 (2000) 1785-1796 (attributed there to McDonald and Floudas, AIChE J. 41
(1995) 1798); ``alpha`` implied by ``G_ij = exp(-alpha_ij tau_ij)``. Antoine
coefficients from the packaged databank (Koretsky 2012),
``ln(P^sat/bar) = A - B/(T + C)``.

**These parameters were fitted to liquid-liquid data and are temperature
independent, and no experimental ternary VLLE data is used anywhere in this
script.** Nothing printed is a comparison against measurement; `thermo` 0.6.0
cannot split two liquids over one excess-Gibbs model (validation Case L-2), so
there is no external three-phase reference for this system either.

Needs no optional dependency. All inputs are SI. See validation Case V-5.
"""

from __future__ import annotations

import argparse
import itertools
from typing import Callable, Sequence

import numpy as np

import chemthermo as ct

NAMES = ("1-Propanol", "n-Butanol", "Water")
PRESSURE_PA = 101325.0
TEMPERATURES = (363.0, 364.0, 365.0)

# Tessier et al. (2000) Table 1, Problem 1. tau[i][j] is tau_ij.
TAU = np.array(
    [
        [0.0, -0.61259, -0.07149],
        [0.7164, 0.0, 0.90047],
        [2.7425, 3.51307, 0.0],
    ]
)
ALPHA = np.array(
    [
        [0.0, 0.3, 0.3],
        [0.3, 0.0, 0.48],
        [0.3, 0.48, 0.0],
    ]
)
CAPITAL_G = np.exp(-ALPHA * TAU)

#: Newton starting points for the tie-triangle, one per temperature.
TRIANGLE_SEEDS = {
    365.0: (0.023616, 0.024407, 0.951977, 0.098737, 0.201955, 0.699308),
    364.0: (0.052144, 0.029172, 0.918684, 0.139639, 0.123478, 0.736883),
    363.0: (0.102828, 0.035390, 0.861782, 0.156303, 0.064227, 0.779470),
}

NEWTON_TOL = 1e-13
TRIANGLE_TOL = 1e-14
INSIDE_LATTICE = 8
OUTSIDE_LATTICE = 12
OUTSIDE_MARGIN = 1e-3

_FAILURES: list[str] = []


def record(label: str, passed: bool) -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {label}")
    if not passed:
        _FAILURES.append(label)


# --------------------------------------------------------------------------
# The model, written out here
# --------------------------------------------------------------------------


def ln_gamma(x: np.ndarray) -> np.ndarray:
    """Renon-Prausnitz NRTL with column sums."""
    values = np.asarray(x, dtype=float)
    values = values / float(np.sum(values))
    s = CAPITAL_G.T @ values
    c = (TAU * CAPITAL_G).T @ values
    return np.array(
        [
            c[i] / s[i] + float(np.sum(values * CAPITAL_G[i, :] / s * (TAU[i, :] - c / s)))
            for i in range(3)
        ]
    )


_ANTOINE = []
for _name in NAMES:
    _record = ct.Component.from_database(_name).antoine
    assert _record is not None, _name
    _ANTOINE.append(_record)


def psat(temperature_K: float) -> np.ndarray:
    """``ln(P^sat / bar) = A - B / (T + C)`` from the packaged databank."""
    return np.array(
        [np.exp(a.A - a.B / (temperature_K + a.C)) * 1.0e5 for a in _ANTOINE], dtype=float
    )


def model() -> ct.NRTL:
    return ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [
                (
                    NAMES[i],
                    NAMES[j],
                    float(TAU[i][j]),
                    float(TAU[j][i]),
                    float(ALPHA[i][j]),
                    float(ALPHA[j][i]),
                )
                for i in range(3)
                for j in range(i + 1, 3)
            ]
        )
    )


def liquid_terms(x: np.ndarray, temperature_K: float) -> np.ndarray:
    return ln_gamma(x) + np.log(psat(temperature_K) / PRESSURE_PA)


def reduced_g(x: np.ndarray, is_vapor: bool, temperature_K: float) -> float:
    values = np.asarray(x, dtype=float)
    terms = np.zeros_like(values) if is_vapor else liquid_terms(values, temperature_K)
    mask = values > 0.0
    return float(np.sum(values[mask] * (np.log(values[mask]) + terms[mask])))


# --------------------------------------------------------------------------
# Independent numerics
# --------------------------------------------------------------------------


def newton(
    residual: Callable[[np.ndarray], np.ndarray],
    start: np.ndarray,
    *,
    positive: Sequence[int],
    tol: float = NEWTON_TOL,
    max_iter: int = 200,
) -> tuple[np.ndarray, float]:
    """Damped Newton with a forward-difference Jacobian and a line search."""
    u = np.array(start, dtype=float)
    step = 1e-8
    for _ in range(max_iter):
        f = residual(u)
        norm = float(np.max(np.abs(f)))
        if not np.isfinite(norm):
            return u, float("inf")
        if norm < tol:
            return u, norm
        jacobian = np.zeros((u.size, u.size), dtype=float)
        for column in range(u.size):
            forward = u.copy()
            forward[column] += step
            jacobian[:, column] = (residual(forward) - f) / step
        try:
            direction = np.linalg.solve(jacobian, -f)
        except np.linalg.LinAlgError:
            return u, norm
        if not np.all(np.isfinite(direction)):
            return u, norm
        scale = 1.0
        for _ in range(60):
            candidate = u + scale * direction
            if np.any(candidate[list(positive)] <= 0.0):
                scale *= 0.5
                continue
            trial = float(np.max(np.abs(residual(candidate))))
            if np.isfinite(trial) and trial < norm:
                u = candidate
                break
            scale *= 0.5
        else:
            return u, norm
    return u, float(np.max(np.abs(residual(u))))


def tie_triangle(temperature_K: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Six unknowns, six equations; the vapor follows from liquid I."""

    def residual(u: np.ndarray) -> np.ndarray:
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

    u, norm = newton(
        residual,
        np.array(TRIANGLE_SEEDS[temperature_K], dtype=float),
        positive=range(6),
        tol=TRIANGLE_TOL,
    )
    x_i, x_ii = u[:3], u[3:]
    y = np.exp(np.log(x_i) + ln_gamma(x_i)) * psat(temperature_K) / PRESSURE_PA
    return x_i, x_ii, y, norm


def vapor_liquid(
    z: np.ndarray, temperature_K: float, seeds: Sequence[np.ndarray]
) -> tuple[np.ndarray, np.ndarray, float, float] | None:
    saturation = psat(temperature_K)

    def residual(u: np.ndarray) -> np.ndarray:
        x, beta = u[:3], u[3]
        k = np.exp(ln_gamma(x)) * saturation / PRESSURE_PA
        return np.concatenate([z - (1.0 - beta) * x - beta * k * x, [float(np.sum(x)) - 1.0]])

    for seed in seeds:
        u, norm = newton(residual, np.concatenate([seed, [0.5]]), positive=range(4))
        x, beta = u[:3], float(u[3])
        if norm >= 1e-12 or not (1e-8 < beta < 1.0 - 1e-8) or not np.all(x > 0.0):
            continue
        y = np.exp(ln_gamma(x)) * saturation / PRESSURE_PA * x
        y = y / float(np.sum(y))
        if float(np.max(np.abs(y - x))) > 1e-7:
            return x, y, beta, norm
    return None


def liquid_liquid(
    z: np.ndarray, temperature_K: float, seeds: Sequence[tuple[np.ndarray, np.ndarray]]
) -> tuple[np.ndarray, np.ndarray, float, float] | None:
    def residual(u: np.ndarray) -> np.ndarray:
        x_i, x_ii, beta = u[:3], u[3:6], u[6]
        activity_i = np.log(x_i) + liquid_terms(x_i, temperature_K)
        activity_ii = np.log(x_ii) + liquid_terms(x_ii, temperature_K)
        balance = z - (1.0 - beta) * x_i - beta * x_ii
        return np.concatenate(
            [
                activity_i - activity_ii,
                [float(np.sum(x_i)) - 1.0, float(np.sum(x_ii)) - 1.0],
                balance[:2],
            ]
        )

    for first, second in seeds:
        u, norm = newton(residual, np.concatenate([first, second, [0.5]]), positive=range(7))
        x_i, x_ii, beta = u[:3], u[3:6], float(u[6])
        if norm >= 1e-12 or not (1e-8 < beta < 1.0 - 1e-8):
            continue
        if not (np.all(x_i > 0.0) and np.all(x_ii > 0.0)):
            continue
        if float(np.max(np.abs(x_i - x_ii))) <= 1e-6:
            continue
        balance = float(np.max(np.abs(z - (1.0 - beta) * x_i - beta * x_ii)))
        if balance < 1e-10:
            return x_i, x_ii, beta, max(norm, balance)
    return None


# --------------------------------------------------------------------------
# The independent classifier and the grid
# --------------------------------------------------------------------------


def expected_state(
    z: np.ndarray,
    temperature_K: float,
    triangle: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[int, str, np.ndarray, list[tuple[float, int, str]]]:
    """The lowest-Gibbs admissible state of ``z``: ``(phases, kind, weights, all)``."""
    x_i, x_ii, y = triangle
    vertices = np.column_stack([x_i, x_ii, y])
    weights = np.linalg.solve(vertices, z)

    scored: list[tuple[float, int, str]] = []
    if np.all(weights > 0.0):
        scored.append(
            (
                float(
                    weights[0] * reduced_g(x_i, False, temperature_K)
                    + weights[1] * reduced_g(x_ii, False, temperature_K)
                    + weights[2] * reduced_g(y, True, temperature_K)
                ),
                3,
                "VLL",
            )
        )

    g_liquid = reduced_g(z, False, temperature_K)
    g_vapor = reduced_g(z, True, temperature_K)
    scored.append((min(g_liquid, g_vapor), 1, "V" if g_vapor < g_liquid else "L"))

    ideal_x = z * PRESSURE_PA / psat(temperature_K)
    ideal_x = ideal_x / float(np.sum(ideal_x))
    two_phase = vapor_liquid(z, temperature_K, (z, ideal_x, x_i, x_ii, y))
    if two_phase is not None:
        x, y_vl, beta, _ = two_phase
        scored.append(
            (
                (1.0 - beta) * reduced_g(x, False, temperature_K)
                + beta * reduced_g(y_vl, True, temperature_K),
                2,
                "VL",
            )
        )

    split = liquid_liquid(z, temperature_K, ((x_i, x_ii), (x_ii, x_i)))
    if split is not None:
        first, second, beta, _ = split
        scored.append(
            (
                (1.0 - beta) * reduced_g(first, False, temperature_K)
                + beta * reduced_g(second, False, temperature_K),
                2,
                "LL",
            )
        )

    scored.sort(key=lambda entry: entry[0])
    return scored[0][1], scored[0][2], weights, scored


def feed_grid(
    triangle: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> list[tuple[str, np.ndarray]]:
    vertices = np.column_stack(list(triangle))
    feeds: list[tuple[str, np.ndarray]] = []

    n = INSIDE_LATTICE
    for i in range(1, n):
        for j in range(1, n - i):
            k = n - i - j
            if k < 1:
                continue
            z = vertices @ (np.array([i, j, k], dtype=float) / n)
            feeds.append(("inside", z / float(np.sum(z))))

    m = OUTSIDE_LATTICE
    for i in range(1, m):
        for j in range(1, m - i):
            k = m - i - j
            if k < 1:
                continue
            z = np.array([i, j, k], dtype=float) / m
            if np.all(np.linalg.solve(vertices, z) > -OUTSIDE_MARGIN):
                continue
            feeds.append(("outside", z))
    return feeds


def flash(nrtl: ct.NRTL, z: np.ndarray, temperature_K: float) -> ct.FlashResult:
    return ct.flash_tp(
        ct.Mixture.from_database(list(NAMES), [float(value) for value in z], normalize=True),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=nrtl,
        flash_mode="modified-raoult",
    )


# --------------------------------------------------------------------------
# The map
# --------------------------------------------------------------------------


def section_map(nrtl: ct.NRTL, temperature_K: float) -> None:
    print(f"\n-- T = {temperature_K:.1f} K " + "-" * 52)
    x_i, x_ii, y, norm = tie_triangle(temperature_K)
    triangle = (x_i, x_ii, y)
    print(f"   independent tie-triangle, Newton residual = {norm:.2e}")
    for label, vertex in (("liquid I ", x_i), ("liquid II", x_ii), ("vapor    ", y)):
        print(f"     {label} = ({vertex[0]:.8f}, {vertex[1]:.8f}, {vertex[2]:.8f})")

    feeds = feed_grid(triangle)
    confusion: dict[tuple[int, int], int] = {}
    disagreements: list[str] = []
    worst_triangle = 0.0
    worst_fraction = 0.0
    worst_equilibrium = 0.0
    worst_balance = 0.0
    worst_single = 0.0
    counts = {1: 0, 2: 0, 3: 0}

    for _region, z in feeds:
        phases, _kind, weights, _scored = expected_state(z, temperature_K, triangle)
        result = flash(nrtl, z, temperature_K)
        obtained = len(result.phase_names())
        confusion[(phases, obtained)] = confusion.get((phases, obtained), 0) + 1
        counts[obtained] = counts.get(obtained, 0) + 1
        if phases != obtained:
            disagreements.append(
                f"z = ({z[0]:.6f}, {z[1]:.6f}, {z[2]:.6f}): expected {phases}, got {obtained}"
                f" ({result.phase_names()})"
            )
            continue

        names = result.phase_names()
        compositions = [
            np.array(result.phases[name].composition.fractions, dtype=float) for name in names
        ]
        fractions = [result.phase_fractions[name] for name in names]

        if obtained == 3:
            best_error = float("inf")
            best_order: tuple[int, ...] = (0, 1, 2)
            for order in itertools.permutations(range(3)):
                error = max(
                    float(np.max(np.abs(compositions[index] - triangle[order[index]])))
                    for index in range(3)
                )
                if error < best_error:
                    best_error, best_order = error, order
            worst_triangle = max(worst_triangle, best_error)
            worst_fraction = max(
                worst_fraction,
                max(abs(fractions[index] - weights[best_order[index]]) for index in range(3)),
            )
        elif obtained == 2:
            terms = [
                np.zeros(3)
                if name.startswith("vapor")
                else liquid_terms(composition, temperature_K)
                for name, composition in zip(names, compositions)
            ]
            worst_equilibrium = max(
                worst_equilibrium,
                float(
                    np.max(
                        np.abs(
                            np.log(compositions[0]) + terms[0] - np.log(compositions[1]) - terms[1]
                        )
                    )
                ),
            )
            worst_balance = max(
                worst_balance,
                float(
                    np.max(
                        np.abs(fractions[0] * compositions[0] + fractions[1] * compositions[1] - z)
                    )
                ),
            )
        else:
            worst_single = max(worst_single, float(np.max(np.abs(compositions[0] - z))))

    print(
        f"\n   {len(feeds)} feeds; obtained 1 / 2 / 3 phases: {counts[1]} / {counts[2]} / {counts[3]}"
    )
    print("   confusion matrix (rows = expected, columns = obtained)")
    print("            obtained 1   obtained 2   obtained 3")
    for expected in (1, 2, 3):
        cells = "".join(f"{confusion.get((expected, got), 0):>13d}" for got in (1, 2, 3))
        print(f"   exp {expected}{cells}")
    if disagreements:
        print("   disagreements:")
        for line in disagreements:
            print(f"     - {line}")

    record(f"{temperature_K:.0f} K: the independent tie-triangle converged", norm < TRIANGLE_TOL)
    record(f"{temperature_K:.0f} K: the grid has at least 60 feeds", len(feeds) >= 60)
    record(f"{temperature_K:.0f} K: no disagreement in the verdict map", not disagreements)
    print(
        f"     worst |dx| against the tie-triangle = {worst_triangle:.2e}; "
        f"worst |dbeta| against the barycentric weights = {worst_fraction:.2e}"
    )
    record(
        f"{temperature_K:.0f} K: every three-phase answer is the same tie-triangle",
        worst_triangle < 1e-8 and worst_fraction < 1e-8,
    )
    print(
        f"     worst two-phase equilibrium residual = {worst_equilibrium:.2e}; "
        f"worst mass-balance residual = {worst_balance:.2e}"
    )
    record(
        f"{temperature_K:.0f} K: every two-phase answer satisfies equal fugacity and mass balance",
        worst_equilibrium < 1e-8 and worst_balance < 1e-10,
    )
    record(
        f"{temperature_K:.0f} K: every single-phase answer is the feed itself",
        worst_single < 1e-12,
    )


def section_repaired_feed(nrtl: ct.NRTL) -> None:
    """The feed of validation Case V-2, which used to come back a single liquid."""
    print("\n-- The Case V-2 feed at 363 K " + "-" * 44)
    temperature_K = 363.0
    x_i, x_ii, y, _ = tie_triangle(temperature_K)
    z = np.column_stack([x_i, x_ii, y]) @ np.array([0.5, 0.3, 0.2])
    z = z / float(np.sum(z))
    print(f"   z = ({z[0]:.8f}, {z[1]:.8f}, {z[2]:.8f})")

    stability = ct.stability_tp(
        ct.Mixture.from_database(list(NAMES), [float(value) for value in z], normalize=True),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=nrtl,
        vapor="ideal",
    )
    print(f"   stability_tp: {stability.status}, tpd_min = {stability.tpd_min:.12e}")
    print(f"     feed candidate = {stability.feed_branch}, incipient = {stability.phase_branch}")
    print(f"     trial surfaces = {stability.diagnostics.get('trial_surfaces')}")
    for trial in stability.trials:
        print(
            f"       {trial.label:<18s} surface={str(trial.surface):<7s} "
            f"converged={str(trial.converged):<5s} trivial={str(trial.trivial):<5s} "
            f"tpd={trial.tpd:+.6e} branch={trial.phase_branch}"
        )
    record("Case V-2 feed: the tangent plane is now unstable", stability.status == "unstable")
    record(
        "Case V-2 feed: the minimizing trial ran on the vapor surface",
        stability.diagnostics.get("minimizing_trial_surface") == "vapor",
    )

    result = flash(nrtl, z, temperature_K)
    print(f"   flash_tp: {result.phase_names()}, regime = {result.diagnostics['phase_regime']}")
    for name in result.phase_names():
        composition = result.phases[name].composition.fractions
        print(
            f"     {name:<9s} = ({composition[0]:.8f}, {composition[1]:.8f}, "
            f"{composition[2]:.8f}), beta = {result.phase_fractions[name]:.8f}"
        )
    record(
        "Case V-2 feed: flash_tp returns three phases",
        sorted(result.phase_names()) == ["liquid1", "liquid2", "vapor"],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="map all three temperatures instead of the first one only",
    )
    # `parse_known_args`, not `parse_args`: `tests/test_examples.py` runs this
    # script via `runpy.run_path` with pytest's own `sys.argv` still in place.
    args, _unknown = parser.parse_known_args()

    print("=" * 78)
    print("Verdict map: 1-propanol / n-butanol / water at 101325 Pa")
    print("Modified Raoult (NRTL liquid + Antoine reference fugacity, ideal vapor)")
    print("Parameters: Tessier et al. (2000) Table 1; LLE-fitted, temperature independent.")
    print("No experimental data is used; nothing here is a comparison against measurement.")
    print("=" * 78)

    nrtl = model()
    temperatures = TEMPERATURES if args.full else TEMPERATURES[:1]
    for temperature_K in temperatures:
        section_map(nrtl, temperature_K)
    section_repaired_feed(nrtl)
    if not args.full:
        print(f"\n  (pass --full for all {len(TEMPERATURES)} temperatures)")

    print("\n" + "=" * 78)
    if _FAILURES:
        print(f"FAIL: {len(_FAILURES)} check(s) failed")
        for failure in _FAILURES:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
