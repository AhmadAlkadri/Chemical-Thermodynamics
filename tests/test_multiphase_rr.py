"""Validation Case V-4: the multiphase Rachford-Rice solver on published data.

Source
------
R. Okuno, R. T. Johns and K. Sepehrnoori, "A new algorithm for Rachford-Rice
for multiphase compositional simulation", SPE Journal 15 (2010) 313-325,
Table 1 (four example overall compositions with constant K-values) and the
text, which prints the solution of Example 3 as
``(beta_1, beta_2) = (0.87, 2.2e-6)`` and gives Example 4 as the phase
compositions ``x_ij`` of a *negative* flash - an overall composition outside
the tie-triangle - from which the phase fractions follow uniquely from
``z = sum_j beta_j x^j``.

Why a constant-K test is the right unit test here
-------------------------------------------------
The Rachford-Rice stage of a multiphase flash takes ``z`` and ``K`` and returns
``beta``; it never calls a thermodynamic model. Fixing ``K`` therefore isolates
exactly this solver, and Table 1 supplies four ``(z, K)`` pairs chosen by the
authors to be hard: Example 2 is a case their reference root-finder cannot
solve at all, Example 3 sits next to a critical end point with a solution
2.2e-06 from a phase boundary, and Example 4 is a near-critical *negative*
flash.

What is asserted
----------------
- the Rachford-Rice equations ``f_j(beta) = sum_i z_i (K_i^j - 1) / t_i`` vanish;
- ``F(beta) = -sum_i z_i ln t_i`` is a minimum, checked both by the Hessian
  being positive definite and by direct comparison against feasible
  perturbations;
- the published solution of Example 3 is reproduced;
- Example 4's phase fractions equal the exact linear solve of
  ``z = sum_j beta_j x^j`` from the paper's printed compositions, including the
  negative one;
- a constructed negative flash returns the exact (negative) fractions, which is
  the signal `flash_tp` uses to *remove* a phase.
"""

from __future__ import annotations

import numpy as np
import pytest

from chemthermo.flash._multiphase_rr import (
    _feasible_region,
    _multiphase_rachford_rice,
    _NoMultiphaseSolution,
)

# --------------------------------------------------------------------------
# Okuno et al. (2010), Table 1. Reference phase is the last phase, so the two
# K columns are phases 1 and 2 against phase 3.
# --------------------------------------------------------------------------

EXAMPLE_1_Z = (
    0.204322076984,
    0.070970999150,
    0.267194323384,
    0.296291964579,
    0.067046080882,
    0.062489248292,
    0.031685306730,
)
EXAMPLE_1_K = (
    (1.23466988745, 1.52713341421),
    (0.89727701141, 0.02456487977),
    (2.29525708098, 1.46348240453),
    (1.58954899888, 1.16090546194),
    (0.23349348597, 0.24166289908),
    (0.02038108640, 0.14815282572),
    (1.40715641002, 14.3128010831),
)

EXAMPLE_2_Z = (
    0.132266176697,
    0.205357472415,
    0.170087543100,
    0.186151796211,
    0.111333894738,
    0.034955417168,
    0.159847699672,
)
EXAMPLE_2_K = (
    (26.3059904941, 66.7435876079),
    (1.91580344867, 1.26478653025),
    (1.42153325608, 0.94711004430),
    (3.21966622946, 3.94954222664),
    (0.22093634359, 0.35954341233),
    (0.01039336513, 0.09327536295),
    (19.4239894458, 12.0162990083),
)

EXAMPLE_3_Z = (
    0.896646630194,
    0.046757914522,
    0.000021572890,
    0.000026632729,
    0.016499094171,
    0.025646758089,
    0.014401397406,
)
EXAMPLE_3_K = (
    (1.64571122126, 1.61947897153),
    (1.91627717926, 2.65352105653),
    (0.71408616431, 0.68719907526),
    (0.28582415424, 0.18483049029),
    (0.04917567928, 0.01228448216),
    (0.00326226927, 0.00023212526),
    (0.00000570946, 0.00000003964),
)

EXAMPLE_4_Z = (0.08860, 0.81514, 0.09626)
EXAMPLE_4_K = (
    (0.112359551, 1.011235955),
    (13.72549020, 0.980392157),
    (3.389830508, 0.847457627),
)
#: Table 1's printed phase compositions for Example 4, columns = phases 1, 2, 3.
EXAMPLE_4_X = (
    (0.100, 0.900, 0.890),
    (0.700, 0.050, 0.051),
    (0.200, 0.050, 0.059),
)

EXAMPLES = {
    "example-1": (EXAMPLE_1_Z, EXAMPLE_1_K),
    "example-2": (EXAMPLE_2_Z, EXAMPLE_2_K),
    "example-3": (EXAMPLE_3_Z, EXAMPLE_3_K),
    "example-4": (EXAMPLE_4_Z, EXAMPLE_4_K),
}


def _rachford_rice_equations(z: np.ndarray, K: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """``f_j(beta) = sum_i z_i (K_i^j - 1) / t_i``, written here, not imported."""
    t = 1.0 + (K - 1.0) @ beta
    return np.array([float(np.sum(z * (K[:, j] - 1.0) / t)) for j in range(K.shape[1])])


def _objective(z: np.ndarray, K: np.ndarray, beta: np.ndarray) -> float:
    """``F(beta) = -sum_i z_i ln t_i``, written here, not imported."""
    t = 1.0 + (K - 1.0) @ beta
    if np.any(t <= 0.0):
        return float("inf")
    return float(-np.sum(z * np.log(t)))


def _hessian(z: np.ndarray, K: np.ndarray, beta: np.ndarray) -> np.ndarray:
    t = 1.0 + (K - 1.0) @ beta
    A = 1.0 - K
    return (A * (z / (t * t))[:, None]).T @ A


@pytest.mark.parametrize("label", sorted(EXAMPLES))
def test_okuno_table_1_rachford_rice_equations_vanish(label: str) -> None:
    """Every Table 1 example solves its Rachford-Rice equations to round-off."""
    z = np.array(EXAMPLES[label][0], dtype=float)
    K = np.array(EXAMPLES[label][1], dtype=float)

    solution = _multiphase_rachford_rice(z, K)

    equations = _rachford_rice_equations(z, K, solution.beta)
    assert np.max(np.abs(equations)) < 1e-11, (label, equations)
    # The phase fractions reproduce the feed exactly: x^r = z / t and
    # x^j = K^j z / t, so sum_j beta_j x^j must be z again.
    t = 1.0 + (K - 1.0) @ solution.beta
    phases = np.column_stack([z / t] + [K[:, j] * z / t for j in range(K.shape[1])])
    recombined = phases @ solution.phase_fractions
    assert np.max(np.abs(recombined - z)) < 1e-12, label
    # Every phase composition is a composition.
    # sum_i x_i^j = sum_i x_i^r is the Rachford-Rice equation itself, and
    # sum_i x_i^r = 1 then follows from sum_j beta_j = 1 and sum_i z_i = 1, so
    # the normalization is only as good as the residual (achieved <= 1.1e-12).
    for column in range(phases.shape[1]):
        assert np.all(phases[:, column] >= 0.0)
        assert abs(float(np.sum(phases[:, column])) - 1.0) < 1e-9


@pytest.mark.parametrize("label", sorted(EXAMPLES))
def test_okuno_table_1_solution_minimizes_the_convex_objective(label: str) -> None:
    """`F` is at a minimum: positive-definite Hessian and no feasible descent."""
    z = np.array(EXAMPLES[label][0], dtype=float)
    K = np.array(EXAMPLES[label][1], dtype=float)

    solution = _multiphase_rachford_rice(z, K)
    beta = solution.beta
    best = _objective(z, K, beta)

    eigenvalues = np.linalg.eigvalsh(_hessian(z, K, beta))
    assert float(np.min(eigenvalues)) > 0.0, (label, eigenvalues)

    # Deterministic feasible perturbations: the solution is interior to S, so a
    # small enough step in any direction stays inside it.
    A, b = _feasible_region(z, K)
    slack = float(np.min(b - A @ beta))
    assert slack > 0.0, label
    step = min(1e-4, 0.25 * slack / max(1.0, float(np.max(np.abs(A)))))
    directions = [
        np.array([1.0, 0.0]),
        np.array([0.0, 1.0]),
        np.array([1.0, 1.0]) / np.sqrt(2.0),
        np.array([1.0, -1.0]) / np.sqrt(2.0),
    ]
    for direction in directions:
        for sign in (1.0, -1.0):
            value = _objective(z, K, beta + sign * step * direction)
            assert value >= best - 1e-15, (label, direction, sign, value, best)


def test_okuno_example_3_reproduces_the_published_solution() -> None:
    """Example 3: the paper prints ``(beta_1, beta_2) = (0.87, 2.2e-6)``.

    The second phase is 2.2e-06 from vanishing, next to a critical end point.
    It is a real phase, which is why phase *removal* keys on a non-positive
    fraction rather than on a small one.
    """
    z = np.array(EXAMPLE_3_Z, dtype=float)
    K = np.array(EXAMPLE_3_K, dtype=float)

    solution = _multiphase_rachford_rice(z, K)

    assert solution.beta[0] == pytest.approx(0.87, abs=5e-3)
    assert solution.beta[1] == pytest.approx(2.2e-6, rel=5e-2)
    # Achieved values, pinned: (0.870163357, 2.18030300e-06), reference phase
    # 0.129834463, 4 Newton iterations, scaled residual 1.7e-15.
    assert solution.beta[0] == pytest.approx(0.8701633569, rel=1e-8)
    assert solution.beta[1] == pytest.approx(2.180303e-06, rel=1e-6)
    assert solution.reference_fraction == pytest.approx(0.1298344631, rel=1e-8)
    assert solution.iterations <= 8
    assert solution.residual < 1e-12
    assert np.all(solution.phase_fractions > 0.0)


def test_okuno_example_4_is_a_negative_flash_matching_the_printed_compositions() -> None:
    """Example 4: the overall composition is outside the tie-triangle.

    Table 1 prints the three phase compositions, so the phase fractions are the
    exact solution of the 3x3 linear system ``z = sum_j beta_j x^j``. One of
    them is negative - that is the paper's point, and it is the signal
    `flash_tp` uses to remove a phase.
    """
    z = np.array(EXAMPLE_4_Z, dtype=float)
    K = np.array(EXAMPLE_4_K, dtype=float)
    X = np.array(EXAMPLE_4_X, dtype=float)

    # The printed K-values are the printed compositions against phase 3.
    assert np.max(np.abs(X[:, :2] / X[:, 2:3] - K)) < 1e-8

    exact = np.linalg.solve(X, z)
    assert exact[2] < 0.0, exact

    solution = _multiphase_rachford_rice(z, K)
    fractions = solution.phase_fractions
    # The solver's own phase order is (reference, 1, 2) = (phase 3, 1, 2).
    assert fractions[1] == pytest.approx(exact[0], abs=1e-8)
    assert fractions[2] == pytest.approx(exact[1], abs=1e-6)
    assert fractions[0] == pytest.approx(exact[2], abs=1e-6)
    # Achieved: (1.2, 14.6599999, -14.859999903), against the exact
    # (1.2, 14.66, -14.86).
    assert fractions[1] == pytest.approx(1.2, abs=1e-8)
    assert fractions[2] == pytest.approx(14.66, abs=1e-6)
    assert fractions[0] == pytest.approx(-14.86, abs=1e-6)
    assert solution.residual < 1e-11


def test_a_constructed_negative_flash_returns_the_exact_negative_fraction() -> None:
    """A synthetic three-phase case whose feed lies outside the tie-triangle.

    Three phase compositions are chosen, ``K`` is built from them, and a feed is
    placed outside their triangle with known (one negative) weights. The solver
    must return exactly those weights: the Rachford-Rice system is linear in the
    phase amounts once ``K`` is fixed, so there is a unique answer to hit.
    """
    phases = np.array(
        [
            [0.70, 0.10, 0.20],
            [0.20, 0.70, 0.15],
            [0.10, 0.20, 0.65],
        ],
        dtype=float,
    )
    weights = np.array([0.8, 0.45, -0.25])
    z = phases @ weights
    assert abs(float(np.sum(z)) - 1.0) < 1e-12
    assert np.all(z > 0.0), z

    K = phases[:, :2] / phases[:, 2:3]
    solution = _multiphase_rachford_rice(z, K)

    assert solution.phase_fractions[1] == pytest.approx(weights[0], abs=1e-10)
    assert solution.phase_fractions[2] == pytest.approx(weights[1], abs=1e-10)
    assert solution.phase_fractions[0] == pytest.approx(weights[2], abs=1e-10)
    assert solution.phase_fractions[0] < 0.0


def test_the_feasible_region_is_the_okuno_set_and_contains_no_pole() -> None:
    """``S`` bounds ``t_i`` away from zero by construction (Okuno et al. eq. 10).

    The claim that makes the region usable is ``t_i >= max(z_i, max_j K_i^j z_i)``
    everywhere in ``S``, so ``F`` has no pole on it - not even on its boundary.
    """
    z = np.array(EXAMPLE_1_Z, dtype=float)
    K = np.array(EXAMPLE_1_K, dtype=float)
    A, b = _feasible_region(z, K)

    floor = np.maximum(z, np.max(K * z[:, None], axis=1))
    solution = _multiphase_rachford_rice(z, K)

    for beta in (solution.beta, np.zeros(2), solution.beta * 0.5):
        if np.any(A @ beta > b):
            continue
        t = 1.0 + (K - 1.0) @ beta
        assert np.all(t >= floor - 1e-12), (beta, t, floor)


def test_two_identical_phases_leave_the_objective_flat_rather_than_wrong() -> None:
    """``K^j = 1`` means phase ``j`` *is* the reference phase.

    ``F`` then does not depend on ``beta_j`` at all: that column of
    ``A = 1 - K`` is zero, the Hessian is singular, and every ``beta_j`` solves
    the system. The solver converges (there is nothing to solve) and the
    duplicate is detected upstream by
    :func:`chemthermo.flash._multiphase._collapsed_phase`, which uses the
    stability module's trivial-solution metric. This test pins that division of
    labour, so a future change that makes the solver guess here is noticed.
    """
    z = np.array([0.3, 0.4, 0.3], dtype=float)
    K = np.column_stack([np.array([3.0, 0.6, 0.4]), np.ones(3)])

    solution = _multiphase_rachford_rice(z, K)

    assert solution.residual < 1e-12
    eigenvalues = np.linalg.eigvalsh(_hessian(z, K, solution.beta))
    assert float(np.min(eigenvalues)) == pytest.approx(0.0, abs=1e-14)
    assert float(np.max(eigenvalues)) > 0.0


def test_a_phase_set_with_no_solution_is_reported_as_a_recession() -> None:
    """A binary with three phases away from its three-phase temperature.

    Gibbs' phase rule allows three phases in a binary at fixed pressure only at
    a single temperature, so at any other temperature the three-phase
    Rachford-Rice system has no finite solution: the feasible region recedes and
    ``F`` decreases for ever along it. These ``z`` and ``K`` are the n-Butanol /
    Water set the flash reaches at ``T3 - 0.05 K`` (Case V-3). The solver must
    report the recession, and name the phase whose amount runs to minus
    infinity, rather than raise a plain convergence error: that report is what
    becomes a phase *removal*.
    """
    z = np.array([0.20, 0.80], dtype=float)
    K = np.array([[0.086745, 1.56986996], [1.28152122, 0.82837416]], dtype=float)

    with pytest.raises(_NoMultiphaseSolution) as error:
        _multiphase_rachford_rice(z, K)

    rates = error.value.fraction_rates
    assert rates.shape == (3,)
    # The reference phase (the vapor, index 0) is the one that leaves.
    assert int(np.argmin(rates)) == 0
    assert float(rates[0]) < 0.0


def test_a_supplied_infeasible_starting_point_is_ignored() -> None:
    """``beta0`` is a hint, never a constraint: an infeasible one is dropped."""
    z = np.array(EXAMPLE_1_Z, dtype=float)
    K = np.array(EXAMPLE_1_K, dtype=float)

    reference = _multiphase_rachford_rice(z, K)
    hinted = _multiphase_rachford_rice(z, K, beta0=np.array([1e6, -1e6]))

    assert np.max(np.abs(hinted.beta - reference.beta)) < 1e-10


def test_components_absent_from_the_feed_do_not_affect_the_solution() -> None:
    """``z_i = 0`` contributes nothing to ``F`` and nothing to the constraints."""
    z = np.array(EXAMPLE_4_Z, dtype=float)
    K = np.array(EXAMPLE_4_K, dtype=float)

    padded_z = np.append(z, 0.0)
    padded_K = np.vstack([K, np.array([7.0, 0.05])])

    base = _multiphase_rachford_rice(z, K)
    padded = _multiphase_rachford_rice(padded_z, padded_K)

    assert np.max(np.abs(padded.beta - base.beta)) < 1e-12
    assert padded.t[-1] == 1.0
