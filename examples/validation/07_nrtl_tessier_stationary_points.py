"""Reproduce the published NRTL tangent-plane stationary points of Tessier (2000).

Source
------
S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable phase stability
analysis for excess Gibbs energy models", Chemical Engineering Science 55
(2000) 1785-1796. Problem 1 is n-propanol(1) / n-butanol(2) / water(3) with the
NRTL parameters of Table 1 (attributed there to McDonald & Floudas, AIChE J. 41
(1995) 1798). Table 2 lists, for four feeds, every stationary point of the
tangent-plane distance D and the value of D at each.

What this script does
---------------------
For each feed z and each printed stationary composition x:

1. evaluates the activity-based tangent-plane distance

       d_i  = ln z_i + ln gamma_i(z)
       D(x) = sum_i x_i [ ln x_i + ln gamma_i(x) - d_i ]

   using `chemthermo.NRTL` for ln gamma;
2. refines the printed (3-digit) composition to a true stationary point of D on
   the composition simplex with a damped Newton solve of
   ``ln w_i + ln gamma_i(w) - d_i = k``, ``sum_i w_i = 1``;
3. compares the refined composition and D against the printed values.

No optional dependency is required (in particular, not `thermo`).

Two printed D values do not reproduce and are reported as KNOWN-TYPO rather
than silently accommodated; see the fixture's `dispute_note` fields and
`.agents/brain/validation-cases.md` Case N-3.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

import chemthermo as ct

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_PATH = REPO_ROOT / "tests" / "fixtures" / "nrtl" / "tessier2000_problem1.json"

# Tolerance on |D_recomputed / D_printed - 1| for the undisputed points.
RELATIVE_D_TOL = 2e-5
# Tolerance on the composition, which the paper prints to three digits.
COMPOSITION_TOL = 1e-3
STATIONARITY_TOL = 1e-10


def _load_fixture() -> dict[str, Any]:
    if not FIXTURE_PATH.is_file():
        raise SystemExit(
            f"Missing validation fixture {FIXTURE_PATH}. Run this script from a repo checkout."
        )
    with FIXTURE_PATH.open("r", encoding="utf-8") as handle:
        payload: dict[str, Any] = json.load(handle)
    return payload


def _ln_gamma(payload: dict[str, Any]) -> Callable[[np.ndarray], np.ndarray]:
    """Build ln gamma(x) for the Problem 1 system from the cited fixture."""
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
    mixture = ct.Mixture.from_database(names, [1 / 3] * len(names), normalize=True)

    def ln_gamma(x: np.ndarray) -> np.ndarray:
        # NRTL ln gamma is homogeneous of degree zero, so rescaling is exact.
        values = np.asarray(x, dtype=float)
        values = values / float(np.sum(values))
        gamma = model.activity_coefficients(
            mixture=mixture,
            temperature_K=298.15,
            composition=[float(value) for value in values],
        )
        return np.log(np.asarray(gamma, dtype=float))

    return ln_gamma


def _tangent_plane_distance(
    ln_gamma: Callable[[np.ndarray], np.ndarray], z: np.ndarray, x: np.ndarray
) -> float:
    d = np.log(z) + ln_gamma(z)
    return float(np.sum(x * (np.log(x) + ln_gamma(x) - d)))


def _stationary_point(
    ln_gamma: Callable[[np.ndarray], np.ndarray], z: np.ndarray, w0: np.ndarray
) -> tuple[np.ndarray, float, float]:
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
            up = u.copy()
            um = u.copy()
            up[column] += step
            um[column] -= step
            jacobian[:, column] = (residuals(up) - residuals(um)) / (2.0 * step)
        delta = np.linalg.solve(jacobian, -f)
        scale = 1.0
        while np.any(u[:n] + scale * delta[:n] <= 0.0):
            scale *= 0.5
        u = u + scale * delta
    return u[:n], float(u[n]), float(np.max(np.abs(residuals(u))))


def main() -> int:
    payload = _load_fixture()
    ln_gamma = _ln_gamma(payload)
    names = [str(entry["chemthermo_name"]) for entry in payload["components"]]

    print("NRTL tangent-plane stationary points, Tessier, Brennecke & Stadtherr (2000)")
    print(f"  Source: {payload['citation']['reference']}")
    print(f"  Location: {payload['citation']['location']}")
    print(f"  System: {' / '.join(names)}  (paper: n-propanol / n-butanol / water)")
    print(f"  Tolerances: |dD/D| <= {RELATIVE_D_TOL:g}, max |dx| <= {COMPOSITION_TOL:g}")
    print("  (*) printed compositions are renormalized to sum to one before use.")

    failures = 0
    disputed = 0
    current_feed: tuple[float, ...] | None = None

    for entry in payload["stationary_points"]:
        z = np.asarray(entry["feed"], dtype=float)
        feed_key = tuple(float(value) for value in entry["feed"])
        if feed_key != current_feed:
            current_feed = feed_key
            print("-" * 78)
            print(f"feed z = ({z[0]:.3f}, {z[1]:.3f}, {z[2]:.3f})")

        x_printed = np.asarray(entry["printed_composition"], dtype=float)
        x_printed = x_printed / float(np.sum(x_printed))
        printed_d = float(entry["printed_D"])

        if entry["trivial"]:
            d_value = _tangent_plane_distance(ln_gamma, z, z)
            status = "PASS" if abs(d_value) < 1e-12 else "FAIL"
            failures += status == "FAIL"
            print(
                f"  trivial point w = z                       "
                f"D printed = {printed_d:+.4e}  recomputed = {d_value:+.4e}   {status}"
            )
            continue

        w, k, residual = _stationary_point(ln_gamma, z, x_printed)
        d_value = _tangent_plane_distance(ln_gamma, z, w)
        comp_error = float(np.max(np.abs(w - x_printed)))
        relative = abs((d_value - printed_d) / printed_d)

        print(
            f"  printed* x = ({x_printed[0]:.6f}, {x_printed[1]:.6f}, {x_printed[2]:.6f})"
            f"  D printed    = {printed_d:+.6e}"
        )
        print(f"  refined  w = ({w[0]:.6f}, {w[1]:.6f}, {w[2]:.6f})  D recomputed = {d_value:+.6e}")
        print(
            f"           stationarity residual = {residual:.2e}, |dx|max = {comp_error:.2e}, "
            f"|dD/D| = {relative:.3e}, D - k = {d_value - k:+.2e}"
        )

        ok = residual < STATIONARITY_TOL and comp_error < COMPOSITION_TOL
        if entry["printed_D_disputed"]:
            disputed += 1
            status = "KNOWN-TYPO" if ok and relative > 1e-3 else "FAIL"
            print(f"           {status}: printed D digits disputed. {entry['dispute_note']}")
        else:
            status = "PASS" if ok and relative < RELATIVE_D_TOL else "FAIL"
            print(f"           {status}")
        if status == "FAIL":
            failures += 1

    print("=" * 78)
    if failures:
        print(f"FAIL: {failures} stationary point(s) outside tolerance.")
        return 1
    print(
        f"PASS: all stationary points reproduced; {disputed} printed D value(s) "
        "flagged as typographical (composition still reproduced)."
    )
    return 0


if __name__ == "__main__":
    exit_code = main()
    if exit_code != 0:
        raise SystemExit(exit_code)
