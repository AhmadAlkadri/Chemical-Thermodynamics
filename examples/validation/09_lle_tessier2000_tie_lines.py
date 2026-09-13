"""Liquid-liquid tie-lines of the Tessier (2000) NRTL systems, verified two ways.

Source
------
S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable phase stability
analysis for excess Gibbs energy models", Chemical Engineering Science 55
(2000) 1785-1796.

* Problem 1 (section 4.1): n-propanol(1) / n-butanol(2) / water(3). Parameters
  Table 1 (attributed there to McDonald & Floudas, AIChE J. 41 (1995) 1798).
* Problem 2 (section 4.2): n-propanol(1) / n-butanol(2) / benzene(3) /
  water(4). Parameters Table 4 (regressed from Gmehling et al., DECHEMA
  Chemistry Data Series 1977-1990).

The paper publishes *stationary points of the tangent-plane distance*, not
tie-lines, so there is no printed tie-line to copy. What is validated here is
therefore the solver against an independent route, plus the invariants that any
converged liquid-liquid split must satisfy.

What this script does
---------------------
For each published feed:

1. runs ``chemthermo.flash_tp(mixture, ..., activity_model=NRTL(<fixture>))``
   with **no** equation of state - the mode is inferred as ``"gamma-gamma"`` -
   and prints the tie-line, the phase fractions and every verification residual
   the solver reports;
2. solves the same equilibrium again from scratch with code written here and
   shared with nothing in ``chemthermo.flash``: an own successive-substitution
   loop followed by a damped Newton solve of

       ln(x_i^I gamma_i^I) - ln(x_i^II gamma_i^II) = 0        (i = 1..n)
       z_i - (1 - beta) x_i^I - beta x_i^II       = 0        (i = 1..n-1)
       sum_i x_i^I - 1 = 0,   sum_i x_i^II - 1 = 0

   with a finite-difference Jacobian;
3. compares the two, **up to a swap of the two phase labels** (``liquid1`` and
   ``liquid2`` are roles assigned by the seed, not identities);
4. checks the invariants: material balance, equal activities,
   ``delta_g_split_rt < 0``, and that both converged phases pass the post-split
   stability test.

The stable control feed of Problem 2, z = (0.25, 0.25, 0.25, 0.25), must come
back as a single liquid.

Where ``examples/validation/08_stability_nrtl_tessier2000.py`` validates the
stability *test*, this script validates the *flash* built on it.

No optional dependency is required (in particular, not `thermo`).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

import chemthermo as ct

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "nrtl"

TEMPERATURE_K = 298.15  # Immaterial: both tables give dimensionless tau.
PRESSURE_PA = 101325.0  # Validated but inert for an activity model.

#: Feeds of Table 2 (Problem 1) and Table 5 (Problem 2), plus the stable control.
PROBLEM1_FEEDS = (
    (0.12, 0.08, 0.80),
    (0.13, 0.07, 0.80),
    (0.12, 0.05, 0.83),
    (0.148, 0.052, 0.80),
)
PROBLEM2_FEEDS = (
    (0.148, 0.052, 0.600, 0.200),
    (0.148, 0.052, 0.700, 0.100),
    (0.25, 0.15, 0.40, 0.20),
    (0.25, 0.15, 0.35, 0.25),
    (0.25, 0.25, 0.25, 0.25),  # stable control
)

#: max |x_solver - x_independent| and |beta_solver - beta_independent|.
TIE_LINE_TOL = 1e-6
#: Residual required of the solver's own equal-activity check.
EQUILIBRIUM_TOL = 1e-10
#: Material-balance residual required of a converged split.
MASS_BALANCE_TOL = 1e-12


def _load(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise SystemExit(f"Missing validation fixture {path}. Run from a repo checkout.")
    with path.open("r", encoding="utf-8") as handle:
        payload: dict[str, Any] = json.load(handle)
    return payload


def _build(payload: dict[str, Any]) -> tuple[list[str], ct.NRTL, Callable[[Any], Any]]:
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

    def ln_gamma(x: Any) -> Any:
        # NRTL ln gamma is homogeneous of degree zero, so rescaling is exact.
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
    """Bisection on the Rachford-Rice function over beta in (0, 1)."""

    def f(beta: float) -> float:
        denominator = 1.0 + beta * (K - 1.0)
        if np.any(denominator <= 0.0):
            return float("nan")
        return float(np.sum(z * (K - 1.0) / denominator))

    low, high = 0.0, 1.0
    f_low = f(low)
    if not math.isfinite(f_low) or not math.isfinite(f(high)) or f_low * f(high) > 0.0:
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
    z: np.ndarray, w: np.ndarray, ln_gamma: Callable[[Any], Any]
) -> tuple[np.ndarray, np.ndarray, float, float, int, int]:
    """Own successive substitution, then damped Newton on the full system."""
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
        except np.linalg.LinAlgError:
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
        if not accepted:
            break

    return u[:n], u[n : 2 * n], float(u[2 * n]), residual, substitutions, newton


def _vector(values: Any, width: int = 6) -> str:
    return "(" + ", ".join(f"{float(v):.{width}f}" for v in values) + ")"


def _compare(
    solver: tuple[np.ndarray, np.ndarray, float],
    reference: tuple[np.ndarray, np.ndarray, float],
) -> float:
    """Largest difference between two tie-lines, allowing a label swap."""
    x_i, x_ii, beta = solver
    r_i, r_ii, r_beta = reference
    direct = max(
        float(np.max(np.abs(x_i - r_i))),
        float(np.max(np.abs(x_ii - r_ii))),
        abs(beta - r_beta),
    )
    swapped = max(
        float(np.max(np.abs(x_i - r_ii))),
        float(np.max(np.abs(x_ii - r_i))),
        abs(beta - (1.0 - r_beta)),
    )
    return min(direct, swapped)


def _run_problem(title: str, path: Path, feeds: tuple[tuple[float, ...], ...]) -> int:
    payload = _load(path)
    names, model, ln_gamma = _build(payload)

    print("=" * 78)
    print(f"{title}: {payload['name']}")
    print(f"  Source: {payload['citation']['reference']}")
    print(f"  System: {' / '.join(names)}")
    print(
        f"  Tolerances: tie-line vs independent solve <= {TIE_LINE_TOL:g}, "
        f"equal-activity residual <= {EQUILIBRIUM_TOL:g}, "
        f"mass balance <= {MASS_BALANCE_TOL:g}"
    )

    failures = 0
    for feed in feeds:
        z = np.asarray(feed, dtype=float)
        z = z / float(np.sum(z))
        mixture = ct.Mixture.from_database(names, list(feed), normalize=True)
        result = ct.flash_tp(
            mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=model,
        )
        diagnostics = result.diagnostics

        print("-" * 78)
        print(f"feed z = {_vector(feed, 3)}")
        print(
            f"  stability: {diagnostics['stability_status']}, "
            f"tpd_min = {float(diagnostics['tpd_min']):+.6e}"
        )

        if len(result.phase_names()) == 1:
            ok = result.phase_names() == ["liquid"] and diagnostics["stability_status"] == "stable"
            print(f"  single liquid (feed found stable); vapor_fraction = {result.vapor_fraction}")
            print(f"  {'PASS' if ok else 'FAIL'}")
            failures += 0 if ok else 1
            continue

        x_i = np.asarray(result.phases["liquid1"].composition.fractions, dtype=float)
        x_ii = np.asarray(result.phases["liquid2"].composition.fractions, dtype=float)
        beta = float(result.phase_fractions["liquid2"])

        stability = ct.stability_tp(
            mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=model,
        )
        assert stability.trial_composition is not None
        reference = _independent_tie_line(
            z, np.asarray(stability.trial_composition, dtype=float), ln_gamma
        )
        r_i, r_ii, r_beta, r_residual, substitutions, newton = reference

        difference = _compare((x_i, x_ii, beta), (r_i, r_ii, r_beta))
        equilibrium = float(diagnostics["equilibrium_residual"])
        mass_balance = float(diagnostics["mass_balance_residual"])
        delta_g = float(diagnostics["delta_g_split_rt"])

        print(f"  flash_tp     liquid1 = {_vector(x_i)}  fraction = {1.0 - beta:.6f}")
        print(f"               liquid2 = {_vector(x_ii)}  fraction = {beta:.6f}")
        print(f"  independent  phase I = {_vector(r_i)}  phase II fraction = {r_beta:.6f}")
        print(f"               phase II= {_vector(r_ii)}  (residual {r_residual:.2e},")
        print(f"               {substitutions} substitutions + {newton} Newton steps)")
        print(
            f"  stages: {diagnostics['ssi_iterations']} ssi + "
            f"{diagnostics['second_order_iterations']} second-order "
            f"({diagnostics['converged_stage']})"
        )
        print(
            f"  |d tie-line| = {difference:.3e}; equal-activity residual = {equilibrium:.3e}; "
            f"mass balance = {mass_balance:.3e}; delta_g_split_rt = {delta_g:+.6e}"
        )
        print(
            f"  post-split: {diagnostics['post_split_status']} "
            f"(liquid1 {diagnostics['phase_stability_liquid1']}, "
            f"liquid2 {diagnostics['phase_stability_liquid2']}, "
            f"most negative tpd {float(diagnostics['post_split_tpd_min']):+.3e})"
        )

        ok = (
            r_residual < 1e-12
            and difference < TIE_LINE_TOL
            and equilibrium < EQUILIBRIUM_TOL
            and mass_balance < MASS_BALANCE_TOL
            and delta_g < 0.0
            and bool(diagnostics["post_split_stable"])
        )
        print(f"  {'PASS' if ok else 'FAIL'}")
        failures += 0 if ok else 1

    return failures


def main() -> int:
    failures = _run_problem(
        "Problem 1", FIXTURE_DIR / "tessier2000_problem1.json", PROBLEM1_FEEDS
    ) + _run_problem("Problem 2", FIXTURE_DIR / "tessier2000_problem2.json", PROBLEM2_FEEDS)

    print("=" * 78)
    if failures:
        print(f"FAIL: {failures} feed(s) outside tolerance.")
        return 1
    print(
        "PASS: every unstable feed splits into a verified liquid-liquid tie-line that an "
        "independent equal-activity solve reproduces, the stable control stays one phase, "
        "and both phases of every split pass the post-split stability test."
    )
    return 0


if __name__ == "__main__":
    exit_code = main()
    if exit_code != 0:
        raise SystemExit(exit_code)
