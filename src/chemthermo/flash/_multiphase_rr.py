"""Multiphase Rachford-Rice as a constrained convex minimization.

References restated in our own notation; the implementation is written from
the mathematics below, not ported from other code.

- R. Okuno, R. T. Johns and K. Sepehrnoori, "A new algorithm for Rachford-Rice
  for multiphase compositional simulation", SPE Journal 15 (2010) 313-325.
- C. F. Leibovici and J. Neoschil, "A solution of Rachford-Rice equations for
  multiphase systems", Fluid Phase Equilibria 112 (1995) 217-221 (the
  ``t_i >= 0`` feasible region this module deliberately does *not* use).
- M. L. Michelsen, "The isothermal flash problem. Part II. Phase-split
  calculation", Fluid Phase Equilibria 9 (1982) 21-40, for the two-phase
  ancestor of the same objective.

The equations
-------------
Take one mole of feed ``z`` split into ``NP`` phases with mole fractions
``beta_j`` and compositions ``x^j``. Choose a **reference phase** ``r`` and
write the equilibrium ratios against it,

    K_i^j = x_i^j / x_i^r        (j != r)                                  (1)

Material balance ``z_i = sum_j beta_j x_i^j`` together with
``beta_r = 1 - sum_{j != r} beta_j`` gives

    z_i = x_i^r t_i,   t_i = 1 + sum_{j != r} beta_j (K_i^j - 1)           (2)

so ``x_i^r = z_i / t_i`` and ``x_i^j = K_i^j z_i / t_i``. Requiring each phase
to be normalized, ``sum_i (x_i^j - x_i^r) = 0``, gives the ``NP - 1``
Rachford-Rice equations

    f_j(beta) = sum_i z_i (K_i^j - 1) / t_i = 0        (j != r)            (3)

The Jacobian of (3), ``df_j / dbeta_k = -sum_i z_i (K_i^j - 1)(K_i^k - 1) /
t_i^2``, is **symmetric**, so (3) is the gradient of a scalar function. Its
antiderivative is (Okuno et al. eq. 9; Michelsen 1994)

    F(beta) = -sum_i z_i ln(t_i)                                           (4)

with ``dF / dbeta_j = -f_j(beta)`` and Hessian

    H_jk = sum_i z_i (1 - K_i^j)(1 - K_i^k) / t_i^2                        (5)

which is positive semidefinite (it is ``A^T diag(z / t^2) A`` with
``A_ij = 1 - K_i^j``), and positive definite whenever ``A`` has full column
rank. So **the multiphase Rachford-Rice equations are the stationarity
conditions of a convex function**, and solving them is a minimization rather
than a root find. That is what makes the solve unconditionally well posed:
there is one minimum or none, never a spurious root, and Newton's direction is
always a descent direction.

The feasible region
-------------------
``F`` has a pole wherever ``t_i = 0``. Leibovici and Neoschil (1995) take the
region ``t_i >= 0``, whose boundary *is* the set of poles: an iterate that
lands near it converges extremely slowly and the Jacobian becomes
catastrophically ill conditioned. Okuno et al. shrink the region instead, using
nothing more than the non-negativity of the phase mole fractions themselves.
From ``x_i^r = z_i / t_i >= 0`` and ``x_i^j = K_i^j z_i / t_i >= 0``, and
because ``sum_i x_i^j = 1`` bounds every mole fraction above by 1,

    t_i >= z_i        and      t_i >= K_i^j z_i   for every j != r          (6)

Writing ``a_i = (1 - K_i^j)_j`` so that ``t_i = 1 - a_i . beta``, (6) is the
polyhedron

    S = { beta : a_i . beta <= b_i },   b_i = min( 1 - z_i,
                                                  min_j (1 - K_i^j z_i) )   (7)

Every point of ``S`` has ``t_i >= max(z_i, max_j K_i^j z_i) > 0`` for an active
component, so ``S`` contains **no pole at all**, not even on its boundary, and
a full step to the boundary is safe.

Note what ``S`` does *not* constrain: the signs of the ``beta_j``. A converged
``beta_j <= 0`` is a meaningful answer - the "negative flash" - and it is
exactly the signal that phase ``j`` does not exist at this feed. Phase removal
in :mod:`chemthermo.flash._multiphase` reads that sign; it is not a failure
mode to be guarded against.

Initial estimate
----------------
Okuno et al. start from the equally weighted mean of the vertices of
``S`` intersected with ``P = { beta : beta_j >= 0, sum_j beta_j <= 1 }``, which
is both feasible and interior. The vertices are enumerated by solving every
``(NP - 1)``-subset of the constraint rows and keeping the feasible solutions -
cheap for the phase counts this package supports (15 subsets for a ternary
three-phase flash). If that intersection is empty - which happens exactly when
no all-positive split exists, i.e. for a negative flash - the enumeration falls
back to the vertices of ``S`` alone.

Solution
--------
Newton on (4) with the exact Hessian (5), the maximum feasible step size of
Okuno et al. step 4, and backtracking on ``F`` from there. Convergence is
measured on ``max_j |f_j|`` divided by the magnitude of the terms that sum
to it, which is the maximum Rachford-Rice residual and the maximum norm of
the gradient at the same time, made scale free.
"""

from __future__ import annotations

import itertools
import math

import numpy as np

from ..exceptions import ConvergenceError

#: Armijo constant for the backtracking line search on ``F``.
_ARMIJO_C = 1e-4
#: Fraction of the feasible step taken on the first trial of the line search.
#: A full step to the boundary of ``S`` is safe (no pole lies there), so this
#: is 1.0 and the search only shortens.
_MAX_FEASIBLE_FRACTION = 1.0
#: Smallest accepted backtracking scale.
_MIN_LINE_SEARCH_SCALE = 1e-16
#: Ridge added to the Hessian before the solve. The Hessian is only positive
#: *semi*definite when ``A`` loses column rank (two phases merging, a critical
#: end point), and this keeps the linear solve usable there.
_HESSIAN_RIDGE = 1e-14
#: Largest number of constraint subsets enumerated for the initial estimate.
_MAX_VERTEX_SUBSETS = 20000
#: Tolerance on the feasibility test used while enumerating vertices.
_VERTEX_FEASIBILITY_TOL = 1e-9
#: Largest phase-fraction magnitude the iteration may reach before the phase
#: set is declared to have no solution. A genuine negative flash stays small
#: (Okuno et al.'s Example 4 converges to beta = (1.2, 14.66, -14.86)); an
#: iterate running to 1e6 is travelling along a recession ray, on which the
#: residual tends to zero without a minimum ever being reached.
_RECESSION_MAGNITUDE = 1.0e6


class _NoMultiphaseSolution(Exception):
    """This phase set has no Rachford-Rice solution at this feed.

    Okuno et al. prove that when the feasible region is unbounded along a
    direction ``d`` with ``(K_i^j - 1) d_j >= 0`` for every component, ``F`` is
    non-increasing along ``d`` for ever: there is no minimum, so there is no
    ``NP``-phase split of this feed. That is not a numerical failure - it is
    the statement that one of the phases does not exist here, which is exactly
    what Gibbs' phase rule says for, say, a *binary* three-phase set at any
    temperature other than its single three-phase temperature.

    ``fraction_rates`` is the rate of change of every phase fraction along that
    recession direction, reference phase first. The phase whose rate is most
    negative is the one whose amount runs to minus infinity, and therefore the
    one to remove.
    """

    def __init__(self, message: str, *, fraction_rates: np.ndarray) -> None:
        super().__init__(message)
        self.fraction_rates = fraction_rates


class _RachfordRiceSolution:
    """Converged multiphase Rachford-Rice solution.

    Attributes:
        beta: Mole fractions of the ``NP - 1`` non-reference phases, in the
            column order of ``K``. May be negative (negative flash).
        reference_fraction: ``1 - sum_j beta_j``, the reference phase's mole
            fraction. May also be negative.
        t: The denominators ``t_i`` of equation (2), one per component.
        residual: Scaled ``max_j |f_j|`` of equation (3) at the solution.
        objective: ``F(beta)`` of equation (4).
        iterations: Newton iterations performed.
        converged: True when ``residual`` met the requested tolerance.
    """

    __slots__ = (
        "beta",
        "converged",
        "iterations",
        "objective",
        "reference_fraction",
        "residual",
        "t",
    )

    def __init__(
        self,
        *,
        beta: np.ndarray,
        t: np.ndarray,
        residual: float,
        objective: float,
        iterations: int,
        converged: bool,
    ) -> None:
        self.beta = beta
        self.reference_fraction = 1.0 - float(np.sum(beta))
        self.t = t
        self.residual = residual
        self.objective = objective
        self.iterations = iterations
        self.converged = converged

    @property
    def phase_fractions(self) -> np.ndarray:
        """All ``NP`` phase fractions, reference phase first."""
        return np.concatenate(([self.reference_fraction], self.beta))


def _feasible_region(z: np.ndarray, K: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(A, b)`` of the Okuno constraint set ``S`` (equation (7))."""
    A = 1.0 - K
    b = np.minimum(1.0 - z, np.min(1.0 - K * z[:, None], axis=1))
    return A, b


def _vertices(A: np.ndarray, b: np.ndarray, dimension: int) -> list[np.ndarray]:
    """Vertices of ``{beta : A beta <= b}`` by enumerating constraint subsets."""
    rows = A.shape[0]
    if rows < dimension:
        return []
    if math.comb(rows, dimension) > _MAX_VERTEX_SUBSETS:
        return []

    found: list[np.ndarray] = []
    for subset in itertools.combinations(range(rows), dimension):
        selection = list(subset)
        matrix = A[selection]
        try:
            if abs(float(np.linalg.det(matrix))) < 1e-14:
                continue
            vertex = np.linalg.solve(matrix, b[selection])
        except np.linalg.LinAlgError:  # pragma: no cover - guarded by the det test
            continue
        if not np.all(np.isfinite(vertex)):
            continue
        if np.all(A @ vertex <= b + _VERTEX_FEASIBILITY_TOL):
            found.append(vertex)
    return found


def _initial_beta(z: np.ndarray, K: np.ndarray) -> np.ndarray:
    """Equally weighted mean of the vertices of ``S`` intersected with ``P``.

    ``P = {beta_j >= 0, sum_j beta_j <= 1}`` is the physical region where every
    phase is present. Its intersection with ``S`` is empty exactly when no
    all-positive split exists, and the estimate then falls back to ``S`` alone
    so that a negative flash still starts from a feasible point.
    """
    A, b = _feasible_region(z, K)
    dimension = K.shape[1]
    positive = np.vstack([A, -np.eye(dimension), np.ones((1, dimension))])
    bounds = np.concatenate([b, np.zeros(dimension), [1.0]])

    found = _vertices(positive, bounds, dimension)
    if not found:
        found = _vertices(A, b, dimension)
    if not found:
        raise ConvergenceError(
            "The multiphase Rachford-Rice feasible region has no vertex, so this phase "
            "set admits no solution at this feed."
        )
    return np.mean(np.array(found), axis=0)


def _multiphase_rachford_rice(
    z: np.ndarray,
    K: np.ndarray,
    *,
    beta0: np.ndarray | None = None,
    tol: float = 1e-12,
    max_iter: int = 100,
) -> _RachfordRiceSolution:
    """Solve the ``NP``-phase Rachford-Rice equations (3) by minimizing (4).

    Args:
        z: Feed mole fractions, one per component. Components with ``z_i == 0``
            are dropped: they contribute nothing to ``F`` and their ``x_i^j``
            are zero in every phase.
        K: ``(NC, NP - 1)`` equilibrium ratios ``K_i^j = x_i^j / x_i^r`` of the
            non-reference phases against the reference phase ``r``.
        beta0: Optional starting point, normally the previous outer iterate.
            Used only when it is feasible; otherwise the vertex-mean estimate
            of :func:`_initial_beta` is used.
        tol: Target for the scaled Rachford-Rice residual (see
            ``scaled_residual`` in the body).
        max_iter: Maximum Newton iterations.

    Returns:
        The converged :class:`_RachfordRiceSolution`. ``t`` is returned at full
        component length, with ``t_i = 1`` for the dropped components so that
        ``x_i^j = K_i^j z_i / t_i = 0`` still holds.

    Raises:
        ConvergenceError: If the feasible region is empty or unbounded (no
            solution with this many phases exists at this feed), or if the line
            search can make no progress before ``tol`` is met.
    """
    z = np.asarray(z, dtype=float)
    K = np.asarray(K, dtype=float)
    if K.ndim != 2 or K.shape[0] != z.size:
        raise ValueError("K must have shape (n_components, n_phases - 1).")
    if np.any(~np.isfinite(K)) or np.any(K <= 0.0):
        raise ConvergenceError("Non-finite or non-positive K-values in the multiphase flash.")

    active = z > 0.0
    z_active = z[active]
    k_active = K[active]
    dimension = K.shape[1]

    A, b = _feasible_region(z_active, k_active)

    beta = None
    if beta0 is not None:
        candidate = np.asarray(beta0, dtype=float)
        if candidate.shape == (dimension,) and np.all(np.isfinite(candidate)):
            if np.all(A @ candidate <= b):
                beta = candidate
    if beta is None:
        beta = _initial_beta(z_active, k_active)

    def denominators(values: np.ndarray) -> np.ndarray:
        return 1.0 + (k_active - 1.0) @ values

    def objective(values: np.ndarray) -> float:
        t = denominators(values)
        if np.any(t <= 0.0):
            return math.inf
        return float(-np.sum(z_active * np.log(t)))

    def equations(t: np.ndarray) -> np.ndarray:
        return ((k_active - 1.0) * (z_active / t)[:, None]).sum(axis=0)

    def scaled_residual(t: np.ndarray, values: np.ndarray) -> float:
        """``max_j |f_j|`` divided by the magnitude of the terms summed in it.

        ``f_j`` is a sum of ``NC`` terms ``z_i (K_i^j - 1) / t_i`` that cancel
        at the solution, so its attainable accuracy is set by the size of those
        terms, not by 1. Dividing by ``sum_i z_i |K_i^j - 1| / t_i`` (floored at
        1 so the measure is never *loosened*) makes the stopping test scale
        free, which matters because a vapor/liquid ``K`` of 30 and a
        liquid/liquid ``K`` of 1.1 appear in the same phase set.
        """
        magnitude = (np.abs(k_active - 1.0) * (z_active / t)[:, None]).sum(axis=0)
        return float(np.max(np.abs(values) / np.maximum(magnitude, 1.0)))

    identity = np.eye(dimension)
    iterations = 0
    converged = False
    for iteration in range(1, max_iter + 1):
        t = denominators(beta)
        if np.any(t <= 0.0):  # pragma: no cover - S contains no pole
            raise ConvergenceError("Multiphase Rachford-Rice left its feasible region.")
        current_equations = equations(t)
        residual = scaled_residual(t, current_equations)
        if residual < tol:
            converged = True
            break
        iterations = iteration

        gradient = -current_equations
        weights = z_active / (t * t)
        hessian = (A * weights[:, None]).T @ A
        direction: np.ndarray
        try:
            direction = np.linalg.solve(hessian + _HESSIAN_RIDGE * identity, -gradient)
        except np.linalg.LinAlgError:  # pragma: no cover - the ridge keeps it solvable
            direction = -gradient
        if not np.all(np.isfinite(direction)) or float(gradient @ direction) >= 0.0:
            direction = -gradient

        # Okuno et al. step 4: the largest step that stays inside S.
        travel = A @ direction
        slack = b - A @ beta
        moving = travel > 0.0
        feasible_step = (
            float(np.min(slack[moving] / travel[moving])) if np.any(moving) else math.inf
        )
        if not math.isfinite(feasible_step):
            raise _NoMultiphaseSolution(
                "The multiphase Rachford-Rice feasible region is unbounded along a descent "
                "direction, so this phase set has no solution at this feed.",
                fraction_rates=np.concatenate(([-float(np.sum(direction))], direction)),
            )
        feasible_step = max(feasible_step, 0.0)

        slope = float(gradient @ direction)
        scale = min(_MAX_FEASIBLE_FRACTION, feasible_step)
        current = objective(beta)
        accepted = False
        while scale >= _MIN_LINE_SEARCH_SCALE:
            candidate = beta + scale * direction
            value = objective(candidate)
            if math.isfinite(value):
                candidate_t = denominators(candidate)
                candidate_residual = scaled_residual(candidate_t, equations(candidate_t))
                # Accept an Armijo decrease of F *or* a decrease of the residual.
                # Close to the minimum the change in F drops below the floating
                # point resolution of F itself and the Armijo test can no longer
                # be satisfied, while the residual still has digits left.
                if value <= current + _ARMIJO_C * scale * slope or candidate_residual < residual:
                    beta = candidate
                    accepted = True
                    break
            scale *= 0.5
        if not accepted:
            break

        if float(np.max(np.abs(beta))) > _RECESSION_MAGNITUDE:
            # F is still decreasing this far out, so the feasible region recedes
            # and there is no minimum: this phase set has no split of this feed.
            # beta ~ s d for large s, so beta itself is the recession direction.
            raise _NoMultiphaseSolution(
                "The multiphase Rachford-Rice objective still decreases at |beta| > "
                f"{_RECESSION_MAGNITUDE:.0e}, so the feasible region recedes and this "
                "phase set has no solution at this feed.",
                fraction_rates=np.concatenate(([-float(np.sum(beta))], beta)),
            )

    t = denominators(beta)
    final_equations = equations(t)
    residual = scaled_residual(t, final_equations)
    converged = converged or residual < tol
    if not converged:
        raise ConvergenceError(
            "Multiphase Rachford-Rice did not converge; max |f_j| = "
            f"{residual:.3e} after {iterations} Newton iterations."
        )

    full_t = np.ones_like(z)
    full_t[active] = t
    return _RachfordRiceSolution(
        beta=beta,
        t=full_t,
        residual=residual,
        objective=objective(beta),
        iterations=iterations,
        converged=True,
    )
