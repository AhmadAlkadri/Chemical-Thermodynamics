"""The liquid-liquid split's second-order stage: a damped Newton minimization.

Successive substitution on the equal-activity condition converges linearly
with a ratio close to one near a plait point, so a second-order stage is
required for the gamma-gamma (liquid-liquid) split (ADR-0009); see
:func:`_second_order_split` for the full derivation. This stage applies to
the liquid-liquid split only - the phi-phi and gamma-phi splits never call
it, which is what keeps every phi-phi and gamma-phi number unchanged.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np

from ..exceptions import ModelError
from .settings import FlashSettings

#: Central-difference step for the second-order stage's Hessian.
_HESSIAN_STEP = 1e-7
#: Smallest accepted backtracking scale in the second-order line search.
_MIN_LINE_SEARCH_SCALE = 1e-14
#: Armijo constant for the line search on the two-phase Gibbs energy.
_ARMIJO_C = 1e-4


class _SecondOrderSplit:
    """Outcome of the second-order liquid-liquid stage."""

    __slots__ = ("beta", "iterations", "residual", "x_i", "x_ii")

    def __init__(
        self,
        *,
        x_i: np.ndarray,
        x_ii: np.ndarray,
        beta: float,
        iterations: int,
        residual: float,
    ) -> None:
        self.x_i = x_i
        self.x_ii = x_ii
        self.beta = beta
        self.iterations = iterations
        self.residual = residual


def _second_order_split(
    *,
    z: np.ndarray,
    x_ii: np.ndarray,
    beta: float,
    ln_gamma: Callable[[np.ndarray], np.ndarray],
    settings: FlashSettings,
) -> _SecondOrderSplit:
    """Damped Newton *minimization* of the two-phase Gibbs energy.

    Derivation
    ----------
    Take one mole of feed and let ``n_i`` be the moles of component ``i`` in
    phase II, so phase I holds ``z_i - n_i``. Write ``L = sum_i (z_i - n_i)``,
    ``V = sum_i n_i`` (``V`` is ``beta``), ``x_i^I = (z_i - n_i) / L`` and
    ``x_i^II = n_i / V``. Both phases are liquids with the same pure-liquid
    reference, so the reference terms contribute the constant
    ``sum_i z_i mu_i^0`` and the *reduced* Gibbs energy that depends on the
    split is

        g(n) = sum_i (z_i - n_i) ln(x_i^I gamma_i^I)
             + sum_i n_i ln(x_i^II gamma_i^II)                            (1)

    Differentiating (1) with respect to ``n_k``, the terms in which the
    *logarithms* move cancel: for either phase, with mole numbers ``N_i`` and
    total ``N``,

        sum_i N_i d ln(x_i gamma_i) = sum_i N_i d ln x_i + sum_i N_i d ln gamma_i
                                    = (sum_i dN_i - dN) + 0 = 0           (2)

    the first bracket because ``d ln x_i = dN_i / N_i - dN / N``, the second by
    Gibbs-Duhem at fixed ``T, P``. What survives is

        dg / dn_k = ln(x_k^II gamma_k^II) - ln(x_k^I gamma_k^I)            (3)

    so **the gradient of the objective is exactly the equal-activity residual**
    that the result reports as ``equilibrium_residual``: a stationary point of
    (1) is equation (1) of the module docstring, and a *minimum* of (1) is the
    equilibrium rather than any other stationary point.

    The Hessian ``d^2 g / dn_j dn_k`` is built by central differences of (3),
    so no derivative of the activity model is needed and no analytic derivative
    code is shared with it. It is symmetrized, and when its smallest eigenvalue
    is not positive a multiple of the identity is added (a standard modified
    Newton step) so the step is a descent direction; if the solve still fails,
    the step falls back to steepest descent.

    Why minimize ``g`` instead of solving the equal-activity system directly
    -----------------------------------------------------------------------
    Newton on the residual system in ``(ln K, beta)`` was implemented and
    measured first, and it is **not** robust here: from the 50-iteration
    successive-substitution iterate of the Tessier et al. (2000) Problem 1 feed
    ``z = (0.12, 0.05, 0.83)`` it walks into the trivial branch
    (``beta -> -8``, the two phases merging) and stalls at a residual of
    1.4e-07; it only converges if given 200 or more substitutions first. A
    residual system cannot tell the equilibrium from the trivial solution -
    both are roots. Minimizing ``g`` can: the trivial solution
    (``n_i = beta z_i`` for any ``beta``) is a stationary *ridge* with
    ``g = g(feed)``, and an unstable feed has ``g < g(feed)`` at the true split,
    so a monotone descent from any iterate below the feed energy cannot reach
    it. With the same 50 substitutions the descent method converges on that feed
    in 7 iterations to a residual of 4.4e-16.

    Line search and box
    -------------------
    ``n`` is kept strictly inside ``0 < n_i < z_i`` (both phases present, no
    negative mole numbers) by backtracking; a step is accepted when it gives an
    Armijo decrease of ``g`` or when it decreases the residual. The stage stops
    at ``settings.second_order_tol``, at ``settings.second_order_max_iter``, or
    when no admissible step improves either measure.

    Args:
        z: Feed mole fractions (normalized).
        x_ii: Phase-II composition of the starting iterate.
        beta: Phase-II mole fraction of the starting iterate.
        ln_gamma: Callable returning ``ln gamma`` at a normalized composition.
        settings: Flash settings (``second_order_tol``,
            ``second_order_max_iter``).

    Returns:
        The refined split; the caller keeps it only if it improved on the
        starting iterate.
    """
    active = z > 0.0
    index = np.flatnonzero(active)

    def clamp(values: np.ndarray) -> np.ndarray:
        floor = 1e-300
        return np.minimum(np.maximum(values, floor), z - floor * np.ones_like(z))

    def energy_and_gradient(n: np.ndarray) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
        liquid_i = z - n
        total_i = float(np.sum(liquid_i))
        total_ii = float(np.sum(n))
        if not (total_i > 0.0 and total_ii > 0.0):
            raise ValueError("degenerate split")
        composition_i = liquid_i / total_i
        composition_ii = n / total_ii
        activity_i = np.zeros_like(z)
        activity_ii = np.zeros_like(z)
        activity_i[active] = np.log(composition_i[active]) + ln_gamma(composition_i)[active]
        activity_ii[active] = np.log(composition_ii[active]) + ln_gamma(composition_ii)[active]
        value = float(
            np.sum(liquid_i[active] * activity_i[active]) + np.sum(n[active] * activity_ii[active])
        )
        return value, (activity_ii - activity_i)[active], composition_i, composition_ii

    n = clamp(beta * x_ii)
    try:
        energy, gradient, composition_i, composition_ii = energy_and_gradient(n)
    except (ValueError, ModelError):
        return _SecondOrderSplit(
            x_i=np.array(z), x_ii=np.array(x_ii), beta=beta, iterations=0, residual=math.inf
        )
    residual = float(np.max(np.abs(gradient)))

    iterations = 0
    size = index.size
    identity = np.eye(size)
    for iteration in range(1, settings.second_order_max_iter + 1):
        if residual < settings.second_order_tol:
            break
        iterations = iteration

        hessian = np.zeros((size, size), dtype=float)
        try:
            for column, component in enumerate(index):
                step = min(_HESSIAN_STEP, 0.25 * n[component], 0.25 * (z[component] - n[component]))
                if step <= 0.0:
                    raise ValueError("degenerate finite-difference step")
                plus = n.copy()
                minus = n.copy()
                plus[component] += step
                minus[component] -= step
                _, gradient_plus, _, _ = energy_and_gradient(plus)
                _, gradient_minus, _, _ = energy_and_gradient(minus)
                hessian[:, column] = (gradient_plus - gradient_minus) / (2.0 * step)
        except (ValueError, ModelError):
            break

        hessian = 0.5 * (hessian + hessian.T)
        direction: np.ndarray
        try:
            smallest = float(np.min(np.linalg.eigvalsh(hessian)))
            shift = 0.0 if smallest > 1e-10 else (1e-10 - smallest)
            direction = np.linalg.solve(hessian + shift * identity, -gradient)
        except np.linalg.LinAlgError:
            direction = -gradient
        if not np.all(np.isfinite(direction)) or float(gradient @ direction) >= 0.0:
            direction = -gradient

        slope = float(gradient @ direction)
        scale = 1.0
        accepted = False
        while scale >= _MIN_LINE_SEARCH_SCALE:
            candidate = n.copy()
            candidate[index] = n[index] + scale * direction
            if np.any(candidate[index] <= 0.0) or np.any(candidate[index] >= z[index]):
                scale *= 0.5
                continue
            try:
                (
                    candidate_energy,
                    candidate_gradient,
                    candidate_i,
                    candidate_ii,
                ) = energy_and_gradient(candidate)
            except (ValueError, ModelError):
                scale *= 0.5
                continue
            candidate_residual = float(np.max(np.abs(candidate_gradient)))
            if (
                candidate_energy < energy + _ARMIJO_C * scale * slope
                or candidate_residual < residual
            ):
                n = candidate
                energy = candidate_energy
                gradient = candidate_gradient
                residual = candidate_residual
                composition_i, composition_ii = candidate_i, candidate_ii
                accepted = True
                break
            scale *= 0.5

        if not accepted:
            break

    return _SecondOrderSplit(
        x_i=composition_i,
        x_ii=composition_ii,
        beta=float(np.sum(n)),
        iterations=iterations,
        residual=residual,
    )
