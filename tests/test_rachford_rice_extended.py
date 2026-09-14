"""The extended (negative-flash) Rachford-Rice solver of ADR-0016.

Three things are being pinned here, in increasing order of importance:

1. the Leibovici-Neoschil window is the right window - it is exactly the set of
   vapor fractions for which every phase mole fraction is non-negative;
2. a root outside ``[0, 1]`` is *found*, and it is the root - checked against a
   residual, against an independently derived closed form for a binary, and
   against a brute-force scan;
3. **bit-identity**: for every ``(z, K)`` on which the old in-window solver had
   an answer, the extended solver returns the same `float`, compared with
   ``==``. That is what makes ADR-0016 safe for
   ``tests/test_flash_refactor_bit_identity.py``, and it is asserted here on
   K-vectors taken from real Peng-Robinson flash iterations, not only on
   synthetic ones.
"""

from __future__ import annotations

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.flash._split import (
    _extended_rachford_rice,
    _rachford_rice,
    _rachford_rice_window,
)


def _f(z: np.ndarray, K: np.ndarray, beta: float) -> float:
    return float(np.sum(z * (K - 1.0) / (1.0 + beta * (K - 1.0))))


# --------------------------------------------------------------------------
# The window.
# --------------------------------------------------------------------------


def test_window_is_exactly_where_every_mole_fraction_is_non_negative() -> None:
    """``t_i > 0`` for all i, inside the window and nowhere outside it."""
    z = np.array([0.5, 0.3, 0.2])
    K = np.array([3.0, 0.9, 0.4])
    window = _rachford_rice_window(K)
    assert window is not None
    lower, upper = window
    assert lower == pytest.approx(1.0 / (1.0 - 3.0))
    assert upper == pytest.approx(1.0 / (1.0 - 0.4))
    # The window always contains the physical range.
    assert lower < 0.0 < 1.0 < upper

    for beta in np.linspace(lower + 1e-9, upper - 1e-9, 501):
        t = 1.0 + beta * (K - 1.0)
        assert np.all(t > 0.0)
        assert np.all(z / t >= 0.0)
        assert np.all(K * z / t >= 0.0)
    for beta in (lower - 1e-6, upper + 1e-6, lower - 1.0, upper + 1.0):
        assert np.any(1.0 + beta * (K - 1.0) <= 0.0)


def test_there_is_no_window_when_every_k_is_on_one_side_of_one() -> None:
    """All-vapor or all-liquid K: no root, and the callers' single-phase case."""
    z = np.array([0.5, 0.3, 0.2])
    for K in (np.array([1.5, 1.2, 1.1]), np.array([0.9, 0.8, 0.5])):
        assert _rachford_rice_window(K) is None
        beta, status = _extended_rachford_rice(z, K)
        assert beta is None
        assert status == "single-phase"
        # And indeed f has no zero anywhere the compositions are admissible:
        # the admissible set is a half line (only one of the two bounds of the
        # window exists), and f keeps one sign along the whole of it.
        bounds = 1.0 / (1.0 - K)
        if np.all(K > 1.0):
            grid = float(np.max(bounds)) + np.geomspace(1e-9, 1e9, 400)
        else:
            grid = float(np.min(bounds)) - np.geomspace(1e-9, 1e9, 400)
        signs = {np.sign(_f(z, K, float(b))) for b in grid}
        assert len(signs) == 1


def test_the_window_shrinks_to_the_physical_range_as_k_approaches_one() -> None:
    K = np.array([1.0 + 1e-12, 1.0 - 1e-12])
    window = _rachford_rice_window(K)
    assert window is not None
    lower, upper = window
    assert lower < -1e11 and upper > 1e11


# --------------------------------------------------------------------------
# The root outside [0, 1].
# --------------------------------------------------------------------------


def test_a_root_below_zero_is_found_and_is_the_root() -> None:
    """Constructed: f(0) < 0 and f(1) < 0, so the root is at negative beta."""
    z = np.array([0.5, 0.3, 0.2])
    K = np.array([1.2, 0.5, 0.3])
    assert _f(z, K, 0.0) < 0.0 and _f(z, K, 1.0) < 0.0
    assert _rachford_rice(z, K)[0] is None

    beta, status = _extended_rachford_rice(z, K)
    assert status == "negative-flash"
    assert beta is not None and beta < 0.0
    assert _f(z, K, beta) == pytest.approx(0.0, abs=1e-14)


def test_a_root_above_one_is_found_and_is_the_root() -> None:
    """Constructed: f(0) > 0 and f(1) > 0, so the root is above one."""
    z = np.array([0.5, 0.3, 0.2])
    K = np.array([2.0, 1.5, 0.99])
    assert _f(z, K, 0.0) > 0.0 and _f(z, K, 1.0) > 0.0
    assert _rachford_rice(z, K)[0] is None

    beta, status = _extended_rachford_rice(z, K)
    assert status == "negative-flash"
    assert beta is not None and beta > 1.0
    assert _f(z, K, beta) == pytest.approx(0.0, abs=1e-14)


def test_the_binary_root_matches_its_closed_form() -> None:
    """For two components the Rachford-Rice root is algebraic; use it.

    Write ``a_i = K_i - 1`` and ``t_i = 1 + beta a_i``. Clearing the two
    denominators of ``f(beta) = z1 a1 / t1 + z2 a2 / t2 = 0`` gives

        z1 a1 (1 + beta a2) + z2 a2 (1 + beta a1) = 0
        (z1 a1 + z2 a2) + beta a1 a2 (z1 + z2) = 0

    so, since ``z1 + z2 = 1``,

        beta = -(z1 a1 + z2 a2) / (a1 a2)

    which is a *closed form*, derived here, holding equally for a root inside
    ``[0, 1]`` and for one outside it.
    """
    for z1, k1, k2 in ((0.9, 3.0, 0.4), (0.2, 5.0, 0.95), (0.6, 1.02, 0.001)):
        z = np.array([z1, 1.0 - z1])
        K = np.array([k1, k2])
        a1, a2 = k1 - 1.0, k2 - 1.0
        numerator = z[0] * a1 + z[1] * a2
        denominator = -(z[0] * a1 * a2 + z[1] * a2 * a1)
        expected = numerator / denominator
        beta, _status = _extended_rachford_rice(z, K)
        assert beta is not None
        assert beta == pytest.approx(expected, rel=1e-12, abs=1e-12)


def test_the_root_agrees_with_a_brute_force_scan_of_the_window() -> None:
    """A dense sign-change scan of the window must bracket the same root."""
    z = np.array([0.4, 0.35, 0.25])
    K = np.array([4.0, 1.3, 0.97])
    window = _rachford_rice_window(K)
    assert window is not None
    lower, upper = window
    grid = np.linspace(lower + 1e-9, upper - 1e-9, 200_001)
    values = np.array([_f(z, K, float(b)) for b in grid])
    changes = np.flatnonzero(np.sign(values[:-1]) != np.sign(values[1:]))
    assert changes.size == 1, "f is monotone: exactly one sign change"
    low, high = float(grid[changes[0]]), float(grid[changes[0] + 1])

    beta, _status = _extended_rachford_rice(z, K)
    assert beta is not None
    assert low <= beta <= high


def test_f_is_strictly_decreasing_on_the_window() -> None:
    """The monotonicity that makes "at most one root" true."""
    z = np.array([0.5, 0.5])
    K = np.array([6.0, 0.2])
    window = _rachford_rice_window(K)
    assert window is not None
    lower, upper = window
    grid = np.linspace(lower + 1e-6, upper - 1e-6, 5_000)
    values = np.array([_f(z, K, float(b)) for b in grid])
    assert np.all(np.diff(values) < 0.0)


# --------------------------------------------------------------------------
# Bit-identity with the in-window solver.
# --------------------------------------------------------------------------


def test_an_in_window_root_is_the_old_solver_s_float_exactly() -> None:
    """Synthetic K-sets whose root is inside [0, 1]: same double, `==`."""
    rng = np.random.default_rng(20260913)
    compared = 0
    for _ in range(400):
        n = int(rng.integers(2, 5))
        z = rng.random(n) + 1e-3
        z = z / z.sum()
        K = np.exp(rng.normal(0.0, 1.5, size=n))
        old, _f0, _f1 = _rachford_rice(z, K)
        if old is None:
            continue
        compared += 1
        new, status = _extended_rachford_rice(z, K)
        assert status == "bracketed"
        assert new == old  # bit-for-bit, not approx
    assert compared >= 50, compared


def _iteration_k_vectors() -> list[tuple[np.ndarray, np.ndarray]]:
    """``(z, K)`` pairs taken from real Peng-Robinson split iterations.

    The pairs are collected by re-running the successive-substitution update
    here, in this test file, from the Wilson estimate - so they are the same
    kind of K-vector the loop actually forms, without reaching into the solver.
    """
    eos = ct.PengRobinsonEOS()
    states = (
        (("Methane", "Ethane"), (0.5, 0.5), 240.0, 3.0e6),
        (("Methane", "Propane"), (0.7, 0.3), 280.0, 8.0e6),
        (("Ethane", "n-Heptane"), (0.7, 0.3), 320.0, 1.0e6),
        (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2), 240.0, 3.0e6),
        (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3), 320.0, 1.0e6),
    )
    collected: list[tuple[np.ndarray, np.ndarray]] = []
    for names, fractions, temperature_K, pressure_Pa in states:
        mixture = ct.Mixture.from_database(list(names), list(fractions), normalize=True)
        z = np.array(mixture.composition.fractions, dtype=float)
        K = np.array(
            [
                (component.pc_pa / pressure_Pa)
                * np.exp(5.373 * (1.0 + component.omega) * (1.0 - component.tc_k / temperature_K))
                for component in mixture.components
            ]
        )
        beta = 0.5
        for _ in range(12):
            collected.append((z.copy(), K.copy()))
            t = 1.0 + beta * (K - 1.0)
            if np.any(t <= 0.0):
                break
            x = z / t
            x = x / x.sum()
            y = K * x
            y = y / y.sum()
            phi_l = np.array(
                eos.fugacity_coefficients(
                    mixture=mixture,
                    temperature_K=temperature_K,
                    pressure_Pa=pressure_Pa,
                    composition=x.tolist(),
                    phase="liquid",
                )
            )
            phi_v = np.array(
                eos.fugacity_coefficients(
                    mixture=mixture,
                    temperature_K=temperature_K,
                    pressure_Pa=pressure_Pa,
                    composition=y.tolist(),
                    phase="vapor",
                )
            )
            K = phi_l / phi_v
            updated, _f0, _f1 = _rachford_rice(z, K)
            if updated is None:
                break
            beta = updated
    return collected


def test_real_iteration_k_vectors_give_the_old_solver_s_float_exactly() -> None:
    pairs = _iteration_k_vectors()
    assert len(pairs) >= 40, len(pairs)
    compared = 0
    for z, K in pairs:
        old, _f0, _f1 = _rachford_rice(z, K)
        new, status = _extended_rachford_rice(z, K)
        if old is None:
            assert status in ("single-phase", "negative-flash", "no-root")
            continue
        compared += 1
        assert status == "bracketed"
        assert new == old
    assert compared >= 40, compared


def test_the_extended_solver_never_reports_bracketed_outside_the_unit_interval() -> None:
    """The status word and the value agree, whichever branch produced them."""
    rng = np.random.default_rng(7)
    for _ in range(300):
        n = int(rng.integers(2, 4))
        z = rng.random(n) + 1e-3
        z = z / z.sum()
        K = np.exp(rng.normal(0.0, 2.0, size=n))
        beta, status = _extended_rachford_rice(z, K)
        if beta is None:
            assert status in ("single-phase", "no-root")
            continue
        if status == "bracketed":
            assert 0.0 <= beta <= 1.0
        else:
            assert status == "negative-flash"
            assert not 0.0 < beta < 1.0
        window = _rachford_rice_window(K)
        assert window is not None
        assert window[0] < beta < window[1]
        assert np.all(1.0 + beta * (K - 1.0) > 0.0)


# --------------------------------------------------------------------------
# The underflowed-K cancellation (ADR-0026).
# --------------------------------------------------------------------------


#: The stability seed of polyethylene(53000) / n-pentane at 453 K and 3.0 MPa,
#: 5 wt% polymer: a trace polymer feed against an essentially pure solvent
#: incipient phase. Reproduced in
#: ``tests/test_pcsaft_polymer.py::test_the_underflowed_k_is_the_stability_seed_at_3_mpa``
#: from the model itself; kept here as plain numbers so this module stays a
#: test of the *equation* and needs no equation of state.
TRACE_Z = np.array([7.163935597e-05, 9.999283606e-01])
TRACE_K = np.array([1.11871119e-18, 1.00132622e00])


def test_the_naive_denominator_cancels_to_zero_at_beta_one() -> None:
    """The measured cause, stated as floating-point arithmetic.

    ``K - 1`` is an exact ``-1.0`` for any ``K`` below the spacing of doubles
    at one, so ``1 + 1 * (K - 1)`` is an exact ``0.0`` - and the positivity
    guard then reports "no admissible vapour fraction" for an equation that
    brackets a root perfectly well.
    """
    assert TRACE_K[0] - 1.0 == -1.0
    assert 1.0 + 1.0 * (TRACE_K[0] - 1.0) == 0.0
    # The convex combination is the same quantity and does not cancel.
    assert (1.0 - 1.0) + 1.0 * TRACE_K[0] == TRACE_K[0]

    # f(0) and f(1) straddle zero, so a root exists in [0, 1].
    f0 = float(np.sum(TRACE_Z * (TRACE_K - 1.0)))
    f1 = float(np.sum(TRACE_Z * (TRACE_K - 1.0) / TRACE_K))
    assert f0 > 0.0 > f1
    assert f1 < -1e13

    # And the default in-window solver - the one the phi-phi split consults
    # before it decides how to seed itself - reports none.
    assert _rachford_rice(TRACE_Z, TRACE_K)[0] is None


def test_the_convex_denominator_finds_that_root_and_it_is_the_root() -> None:
    """The repair is of the equation, not of the search.

    The wider (negative-flash) window happens to find a root here too - it
    never evaluates ``f`` at exactly ``beta = 1``, so the cancellation does not
    reach it - but that is luck, and the caller that refused this state
    (``chemthermo.flash._detect._flash_tp_tangent_plane``) asks the in-window
    solver. What ``convex_denominators`` does is give that solver back the
    ``f(1)`` it always had.
    """
    plain, plain_status = _extended_rachford_rice(TRACE_Z, TRACE_K)
    assert plain is not None and plain_status == "bracketed"

    beta, status = _extended_rachford_rice(TRACE_Z, TRACE_K, convex_denominators=True)
    assert beta is not None
    assert status == "bracketed"
    assert 0.0 < beta < 1.0

    # Independently written residual, in the cancellation-free form.
    def f(v: float) -> float:
        return float(np.sum(TRACE_Z * (TRACE_K - 1.0) / ((1.0 - v) + v * TRACE_K)))

    assert abs(f(beta)) < 1e-12
    # The closed form for a binary: clearing the two denominators of
    # ``z_1 a_1 / t_1 + z_2 a_2 / t_2 = 0`` with ``a_i = K_i - 1`` and
    # ``t_i = 1 + beta a_i`` leaves ``z_1 a_1 + z_2 a_2 + beta a_1 a_2 = 0``
    # for a normalized feed, which is linear in beta.
    z1, z2 = TRACE_Z
    a1, a2 = TRACE_K - 1.0
    closed = -(z1 * a1 + z2 * a2) / (a1 * a2)
    assert beta == pytest.approx(closed, rel=1e-9)

    # Both phase compositions are non-negative there.
    t = (1.0 - beta) + beta * TRACE_K
    assert np.all(t > 0.0)
    assert np.all(TRACE_Z / t >= 0.0)


def test_the_convex_form_is_off_by_default_and_moves_no_root_that_existed() -> None:
    """Bit-identity: the same double wherever the naive denominator worked.

    ``f(0)``'s denominator is an exact ``1.0`` either way and the bisection
    only ever evaluates ``f`` strictly inside ``(0, 1)``, where the naive form
    is admissible and is therefore what both spellings use. Asserted with
    ``==`` on synthetic K-sets and on K-vectors from real flash iterations.
    """
    rng = np.random.default_rng(20260914)
    compared = 0
    for _ in range(400):
        n = int(rng.integers(2, 5))
        z = rng.random(n) + 1e-3
        z = z / z.sum()
        K = np.exp(rng.normal(0.0, 1.5, size=n))
        plain, _f0, _f1 = _rachford_rice(z, K)
        convex, _g0, _g1 = _rachford_rice(z, K, convex_denominators=True)
        if plain is None:
            continue
        compared += 1
        assert convex == plain
        assert (
            _extended_rachford_rice(z, K, convex_denominators=True)[0]
            == (_extended_rachford_rice(z, K)[0])
        )
    assert compared >= 50, compared

    for z, K in _iteration_k_vectors():
        plain, _f0, _f1 = _rachford_rice(z, K)
        if plain is None:
            continue
        assert _rachford_rice(z, K, convex_denominators=True)[0] == plain


def test_the_convex_form_still_refuses_a_genuinely_rootless_k_set() -> None:
    """It repairs a rounding artefact, not the single-phase verdict itself."""
    z = np.array([0.4, 0.6])
    for K in (np.array([0.3, 0.5]), np.array([2.0, 4.0])):
        assert _rachford_rice(z, K, convex_denominators=True)[0] is None
        assert _extended_rachford_rice(z, K, convex_denominators=True) == (None, "single-phase")
