"""Published-anchor validation of NRTL against Tessier, Brennecke & Stadtherr (2000).

Two independent routes are exercised here:

1. `thermo.NRTL` (and `thermo.nrtl.NRTL_gammas`) as an external reference
   implementation of the same closed-form equation. Those tests skip when
   `thermo` is not installed.
2. A from-scratch tangent-plane distance and stationary-point solve written in
   this module, compared with the compositions and D values printed in the
   paper's Table 2. This route needs no optional dependency and always runs.

Source: S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable phase
stability analysis for excess Gibbs energy models", Chem. Eng. Sci. 55 (2000)
1785-1796 (author copy: https://academicweb.nd.edu/~markst/srt2000.pdf).
Parameters: Table 1, attributed to McDonald & Floudas, AIChE J. 41 (1995) 1798.
Stationary points and D: Table 2.

Definitions used below (paper section 2, equations 1-4):

    d_i  = ln z_i + ln gamma_i(z)
    D(x) = sum_i x_i [ ln x_i + ln gamma_i(x) - d_i ]

A stationary point of D on the composition simplex satisfies
``ln x_i + ln gamma_i(x) - d_i = k`` for all i with a single constant k, and
at such a point D = k. That identity is asserted directly.
"""

from __future__ import annotations

from typing import Any, Callable, Sequence

import numpy as np
import pytest

TEMPERATURE_K = 298.15  # Immaterial: Table 1 gives dimensionless tau directly.

LnGamma = Callable[[np.ndarray], np.ndarray]


def _tangent_plane_distance(
    ln_gamma: Callable[[np.ndarray], np.ndarray], z: np.ndarray, x: np.ndarray
) -> float:
    d = np.log(z) + ln_gamma(z)
    return float(np.sum(x * (np.log(x) + ln_gamma(x) - d)))


def _successive_substitution(
    ln_gamma: Callable[[np.ndarray], np.ndarray],
    z: np.ndarray,
    w0: np.ndarray,
    max_iter: int = 20000,
) -> tuple[np.ndarray, float]:
    """Iterate ln W_i = d_i - ln gamma_i(w), w = W / sum(W).

    Returns the converged composition and the stationarity residual
    max_i |ln W_i + ln gamma_i(w) - d_i|.
    """
    d = np.log(z) + ln_gamma(z)
    w = np.asarray(w0, dtype=float) / float(np.sum(w0))
    for _ in range(max_iter):
        big_w = np.exp(d - ln_gamma(w))
        w_next = big_w / float(np.sum(big_w))
        if float(np.max(np.abs(w_next - w))) < 1e-16:
            w = w_next
            break
        w = w_next
    big_w = np.exp(d - ln_gamma(w))
    residual = float(np.max(np.abs(np.log(big_w) + ln_gamma(w) - d)))
    return w, residual


def _stationary_point(
    ln_gamma: Callable[[np.ndarray], np.ndarray],
    z: np.ndarray,
    w0: np.ndarray,
    max_iter: int = 200,
) -> tuple[np.ndarray, float, float]:
    """Damped Newton solve of the simplex-constrained stationarity system.

    Unknowns are (w, k); equations are
    ``ln w_i + ln gamma_i(w) - d_i - k = 0`` for every i, plus
    ``sum_i w_i - 1 = 0``. The Jacobian is built by central differences, so
    the solve shares no analytic derivative code with the model.

    Returns (w, k, residual). Unlike successive substitution this converges to
    saddle points and maxima of D as well as to minima, which matters for the
    near-plait-point feed of Table 2.
    """
    n = z.size
    d = np.log(z) + ln_gamma(z)

    def residuals(u: np.ndarray) -> np.ndarray:
        w = u[:n]
        k = u[n]
        return np.concatenate([np.log(w) + ln_gamma(w) - d - k, [float(np.sum(w)) - 1.0]])

    u = np.concatenate([np.asarray(w0, dtype=float) / float(np.sum(w0)), [0.0]])
    step = 1e-7
    for _ in range(max_iter):
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


def test_published_G_matrix_is_reproduced_by_the_implied_alpha(
    tessier2000_payload: dict[str, Any],
) -> None:
    """Table 1 prints G_ij; alpha_ij is implied by G_ij = exp(-alpha_ij tau_ij).

    The quotients -ln(G_ij)/tau_ij round to 0.3 (pairs 1-2, 1-3) and 0.48
    (pair 2-3); re-exponentiating those rounded values must reproduce the
    printed G matrix, otherwise the fixture's alpha is not a faithful stand-in.
    """
    tau = np.asarray(tessier2000_payload["tau"], dtype=float)
    alpha = np.asarray(tessier2000_payload["alpha"], dtype=float)
    g_published = np.asarray(tessier2000_payload["G_published"], dtype=float)

    recovered = np.zeros_like(tau)
    off_diagonal = ~np.eye(tau.shape[0], dtype=bool)
    recovered[off_diagonal] = -np.log(g_published[off_diagonal]) / tau[off_diagonal]
    assert np.allclose(recovered[off_diagonal], alpha[off_diagonal], rtol=0.0, atol=1e-6)

    g_from_alpha = np.exp(-alpha * tau)
    max_error = float(np.max(np.abs(g_from_alpha - g_published)))
    # Achieved: 3.99e-08.
    assert max_error < 1e-7
    assert tessier2000_payload["G_reconstruction_max_abs_error"] == pytest.approx(
        max_error, rel=1e-9
    )


def test_nrtl_matches_thermo_for_asymmetric_parameters(
    tessier2000_payload: dict[str, Any], tessier2000_ln_gamma: LnGamma
) -> None:
    """External reference: `thermo` implements the same closed form.

    Both evaluate the standard Renon-Prausnitz equation, so agreement should
    be at round-off, not at a "physically reasonable" tolerance. Achieved max
    |d ln gamma| = 8.9e-16 over the compositions below.
    """
    thermo = pytest.importorskip("thermo")
    thermo_nrtl = thermo.NRTL
    nrtl_gammas = thermo.nrtl.NRTL_gammas

    tau = np.asarray(tessier2000_payload["tau"], dtype=float)
    alpha = np.asarray(tessier2000_payload["alpha"], dtype=float)
    size = tau.shape[0]
    ln_gamma = tessier2000_ln_gamma

    compositions: Sequence[Sequence[float]] = (
        (0.12, 0.08, 0.80),
        (0.50, 0.30, 0.20),
        (0.0597449, 0.0282358, 0.9120193),
        (0.20, 0.20, 0.60),
        (0.80, 0.10, 0.10),
        (0.05, 0.90, 0.05),
    )

    worst_class = 0.0
    worst_function = 0.0
    for composition in compositions:
        x = np.asarray(composition, dtype=float)
        x = x / float(np.sum(x))
        ours = ln_gamma(x)

        reference_class = thermo_nrtl(
            T=TEMPERATURE_K,
            xs=x.tolist(),
            tau_coeffs=[[[tau[i][j], 0, 0, 0, 0, 0] for j in range(size)] for i in range(size)],
            alpha_coeffs=[[[alpha[i][j], 0] for j in range(size)] for i in range(size)],
        )
        diff_class = float(np.max(np.abs(np.asarray(reference_class.lngammas()) - ours)))
        diff_function = float(
            np.max(
                np.abs(
                    np.log(np.asarray(nrtl_gammas(xs=x.tolist(), taus=tau, alphas=alpha))) - ours
                )
            )
        )
        worst_class = max(worst_class, diff_class)
        worst_function = max(worst_function, diff_function)

        assert diff_class < 1e-9
        assert diff_function < 1e-9

    assert worst_class < 1e-12
    assert worst_function < 1e-12


def test_tessier2000_table2_stationary_points_are_reproduced(
    tessier2000_payload: dict[str, Any], tessier2000_ln_gamma: LnGamma
) -> None:
    """The exam: reproduce every stationary point and D value in Table 2.

    Six of the seven non-trivial printed D values are reproduced to 1.2e-05
    relative or better. Two are not, and this test asserts against the
    recomputed values while documenting the mismatch rather than choosing
    whichever number passes -- see the fixture's `dispute_note` entries and the
    validation ledger (Case N-3).
    """
    ln_gamma = tessier2000_ln_gamma
    disputed_seen = 0

    for entry in tessier2000_payload["stationary_points"]:
        z = np.asarray(entry["feed"], dtype=float)
        x_printed = np.asarray(entry["printed_composition"], dtype=float)
        x_printed = x_printed / float(np.sum(x_printed))
        printed_d = float(entry["printed_D"])

        if entry["trivial"]:
            # The feed itself is always a stationary point with D = 0 exactly.
            assert _tangent_plane_distance(ln_gamma, z, z) == pytest.approx(0.0, abs=1e-15)
            assert printed_d == 0.0
            continue

        w, k, residual = _stationary_point(ln_gamma, z, x_printed)
        assert residual < 1e-10
        d_value = _tangent_plane_distance(ln_gamma, z, w)
        # At a stationary point the constant k equals D.
        assert d_value == pytest.approx(k, rel=1e-9, abs=1e-15)
        # The converged point must still be the printed one (3 printed digits).
        assert np.allclose(w, x_printed, rtol=0.0, atol=1e-3)

        relative = abs((d_value - printed_d) / printed_d)
        if entry["printed_D_disputed"]:
            disputed_seen += 1
            # Far outside the round-off of a 5-digit printed value, and the
            # composition is reproduced: the printed digits are the outlier.
            assert relative > 1e-3
        else:
            assert relative < 2e-5, f"feed {entry['feed']} point {entry['printed_composition']}"

    assert disputed_seen == 2


def test_successive_substitution_reaches_the_printed_minima(
    tessier2000_payload: dict[str, Any], tessier2000_ln_gamma: LnGamma
) -> None:
    """Successive substitution reproduces the D minima but not the saddle.

    Starting from each printed composition, ln W_i = d_i - ln gamma_i(w)
    converges to the printed root for every point except the near-plait-point
    root of the z1 = 0.148 feed, where it drifts to the trivial solution
    w = z. That is the initialization dependence Tessier et al. set out to
    eliminate, so it is asserted as observed behaviour rather than hidden.
    """
    ln_gamma = tessier2000_ln_gamma
    drifted_to_trivial = 0

    for entry in tessier2000_payload["stationary_points"]:
        if entry["trivial"]:
            continue
        z = np.asarray(entry["feed"], dtype=float)
        x_printed = np.asarray(entry["printed_composition"], dtype=float)
        x_printed = x_printed / float(np.sum(x_printed))

        w, residual = _successive_substitution(ln_gamma, z, x_printed)
        assert residual < 1e-10

        if np.allclose(w, z, rtol=0.0, atol=1e-6):
            drifted_to_trivial += 1
            continue
        assert np.allclose(w, x_printed, rtol=0.0, atol=1e-3)
        d_value = _tangent_plane_distance(ln_gamma, z, w)
        assert d_value < 0.0

    assert drifted_to_trivial == 1
