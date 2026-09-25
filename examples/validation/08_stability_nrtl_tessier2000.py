"""Reproduce the published NRTL tangent-plane global minima of Tessier (2000).

Source
------
S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable phase stability
analysis for excess Gibbs energy models", Chemical Engineering Science 55
(2000) 1785-1796.

* Problem 1 (section 4.1): n-propanol(1) / n-butanol(2) / water(3). Parameters
  Table 1 (attributed there to McDonald & Floudas, AIChE J. 41 (1995) 1798),
  stationary points and D values Table 2.
* Problem 2 (section 4.2): n-propanol(1) / n-butanol(2) / benzene(3) /
  water(4). Parameters Table 4 (regressed from Gmehling et al., DECHEMA
  Chemistry Data Series 1977-1990), stationary points and D values Table 5.

What this script does
---------------------
For each published feed:

1. runs `chemthermo.stability_tp(..., activity_model=NRTL(parameters=<fixture>))`
   -- i.e. the package's own Michelsen tangent-plane solver, given only the feed
   and the model;
2. independently refines every printed stationary composition of that feed to a
   true stationary point of

       d_i  = ln z_i + ln gamma_i(z)
       D(x) = sum_i x_i [ ln x_i + ln gamma_i(x) - d_i ]

   with a damped Newton solve of ``ln w_i + ln gamma_i(w) - d_i = k``,
   ``sum_i w_i = 1``, and takes the most negative (lowest) non-trivial D as the
   feed's reference global minimum;
3. compares the solver's `tpd_min`, its minimizing composition and its verdict
   against that reference and against the printed five-digit D.

Where `examples/validation/07_nrtl_tessier_stationary_points.py` validates the
NRTL *model* against Table 2, this script validates the *stability solver*.

No optional dependency is required (in particular, not `thermo`).

Printed D values that do not reproduce are reported as KNOWN-TYPO rather than
silently accommodated; see the fixtures' `dispute_note` fields and
`.agents/brain/validation-cases.md` Cases N-3, S-6 and S-7.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

import chemthermo as ct

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "nrtl"
PROBLEM_FIXTURES = (
    ("Problem 1", FIXTURE_DIR / "tessier2000_problem1.json"),
    ("Problem 2", FIXTURE_DIR / "tessier2000_problem2.json"),
)

TEMPERATURE_K = 298.15  # Immaterial: both tables give dimensionless tau.
PRESSURE_PA = 101325.0  # Validated but inert for an activity model.

# |tpd_min / D_refined - 1| for the solver against this script's own refinement.
SOLVER_REL_TOL = 1e-6
# max |w_solver - w_refined|.
COMPOSITION_TOL = 1e-6
# Against the printed five-digit D, as a pure relative tolerance. Five printed
# digits alone justify ~5e-05; the alpha implied by the printed G matrix moves D
# by up to 9e-05 relative for Problem 2, hence 2.5e-04.
PRINTED_REL_TOL = 2.5e-4
STATIONARITY_TOL = 1e-12


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


def _tangent_plane_distance(ln_gamma: Callable[[Any], Any], z: Any, x: Any) -> float:
    d = np.log(z) + ln_gamma(z)
    return float(np.sum(x * (np.log(x) + ln_gamma(x) - d)))


def _refine(ln_gamma: Callable[[Any], Any], z: Any, w0: Any) -> tuple[np.ndarray, float, float]:
    """Damped Newton on (w, k): ln w_i + ln gamma_i(w) - d_i - k = 0, sum w = 1."""
    n = z.size
    d = np.log(z) + ln_gamma(z)

    def residuals(u: np.ndarray) -> np.ndarray:
        return np.concatenate(
            [np.log(u[:n]) + ln_gamma(u[:n]) - d - u[n], [float(np.sum(u[:n])) - 1.0]]
        )

    u = np.concatenate([np.asarray(w0, dtype=float) / float(np.sum(w0)), [0.0]])
    step = 1e-7
    for _ in range(200):
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


def _feed_groups(payload: dict[str, Any]) -> list[tuple[tuple[float, ...], list[dict[str, Any]]]]:
    groups: dict[tuple[float, ...], list[dict[str, Any]]] = {}
    for entry in payload["stationary_points"]:
        groups.setdefault(tuple(float(v) for v in entry["feed"]), []).append(entry)
    return list(groups.items())


def _vector(values: Any, width: int = 8) -> str:
    return "(" + ", ".join(f"{float(v):.{width}f}" for v in values) + ")"


def _run_problem(title: str, path: Path) -> tuple[int, int]:
    payload = _load(path)
    names, model, ln_gamma = _build(payload)

    print("=" * 78)
    print(f"{title}: {payload['name']}")
    print(f"  Source: {payload['citation']['reference']}")
    print(f"  Location: {payload['citation']['location']}")
    print(f"  System: {' / '.join(names)}")
    print(
        f"  Tolerances: solver vs refined |dD/D| <= {SOLVER_REL_TOL:g}, "
        f"|dw| <= {COMPOSITION_TOL:g}; refined vs printed |dD/D| <= {PRINTED_REL_TOL:g}"
    )

    failures = 0
    disputed = 0

    for feed, entries in _feed_groups(payload):
        z = np.asarray(feed, dtype=float)
        print("-" * 78)
        print(f"feed z = {_vector(feed, 3)}")

        # Reference: refine every printed non-trivial point, keep the lowest D.
        best_entry: dict[str, Any] | None = None
        best_w: np.ndarray | None = None
        best_d = float("inf")
        worst_residual = 0.0
        for entry in entries:
            if entry["trivial"]:
                continue
            printed_w = np.asarray(entry["printed_composition"], dtype=float)
            printed_w = printed_w / float(np.sum(printed_w))
            w, _, residual = _refine(ln_gamma, z, printed_w)
            worst_residual = max(worst_residual, residual)
            value = _tangent_plane_distance(ln_gamma, z, w)
            if value < best_d:
                best_d, best_w, best_entry = value, w, entry

        if best_entry is None or best_w is None:
            print("  no non-trivial printed stationary point: skipped")
            continue

        printed_d = float(best_entry["printed_D"])
        expected_status = "unstable" if best_d < 0.0 else "stable"

        result = ct.stability_tp(
            ct.Mixture.from_database(names, list(feed), normalize=True),
            temperature_K=TEMPERATURE_K,
            pressure_Pa=PRESSURE_PA,
            activity_model=model,
        )

        print(
            f"  printed global minimum  D = {printed_d:+.6e}  at x = {_vector(best_entry['printed_composition'], 6)}"
        )
        print(f"  refined global minimum  D = {best_d:+.10e}  at w = {_vector(best_w, 8)}")
        print(
            f"  stability_tp     tpd_min = {result.tpd_min:+.10e}  status = {result.status.upper()}"
        )

        composition_error = float("nan")
        if result.trial_composition is not None:
            composition_error = float(
                np.max(np.abs(np.asarray(result.trial_composition, dtype=float) - best_w))
            )
            print(f"                         w = {_vector(result.trial_composition, 8)}")

        solver_rel = abs((result.tpd_min - best_d) / best_d) if best_d != 0.0 else float("nan")
        printed_rel = abs((best_d - printed_d) / printed_d) if printed_d != 0.0 else float("nan")
        stages = ", ".join(
            f"{trial.label}: {trial.ssi_iterations}+{trial.second_order_iterations}"
            f"{'' if trial.converged_stage is None else ' (' + trial.converged_stage + ')'}"
            for trial in result.trials
        )
        print(f"  refinement residual = {worst_residual:.2e}; trials (ssi+newton) = {stages}")
        print(
            f"  |dD/D| solver-vs-refined = {solver_rel:.3e}, "
            f"|dw| = {composition_error:.2e}, |dD/D| refined-vs-printed = {printed_rel:.3e}"
        )

        solver_ok = (
            worst_residual < STATIONARITY_TOL
            and result.status == expected_status
            and solver_rel < SOLVER_REL_TOL
            and composition_error < COMPOSITION_TOL
        )
        printed_ok = printed_rel < PRINTED_REL_TOL

        if best_entry["printed_D_disputed"]:
            disputed += 1
            status = "KNOWN-TYPO" if solver_ok and not printed_ok else "FAIL"
            note = best_entry.get("dispute_note", "printed D digits disputed")
            print(f"  {status}: {note}")
        else:
            status = "PASS" if solver_ok and printed_ok else "FAIL"
            print(f"  {status}")

        if status == "FAIL":
            failures += 1

    return failures, disputed


def main() -> int:
    failures = 0
    disputed = 0
    for title, path in PROBLEM_FIXTURES:
        problem_failures, problem_disputed = _run_problem(title, path)
        failures += problem_failures
        disputed += problem_disputed

    print("=" * 78)
    if failures:
        print(f"FAIL: {failures} feed(s) outside tolerance.")
        return 1
    print(
        f"PASS: every published global tangent-plane minimum reproduced by stability_tp; "
        f"{disputed} printed D value(s) flagged as typographical."
    )
    return 0


if __name__ == "__main__":
    exit_code = main()
    if exit_code != 0:
        raise SystemExit(exit_code)
