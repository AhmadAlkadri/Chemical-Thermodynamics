"""The two-phase split written in **log mole numbers** (ADR-0024).

:mod:`chemthermo.flash._second_order` minimizes the two-phase Gibbs energy in
the *linear* mole numbers ``n`` of one phase, boxed by ``0 < n_i < z_i`` with a
floor of ``1e-300``. That is the right parametrization whenever the two phases
have compositions of comparable magnitude, and it is what every converging
state in this repository uses.

It is the wrong one for a polymer/solvent vapour-liquid split. The equilibrium
vapour above a polyethylene melt holds a polymer mole fraction of
``exp(-450)``: a mole number that a floor of ``1e-300`` cannot represent, that
a finite-difference step of ``1e-7`` cannot resolve, and whose Newton
correction spans hundreds of orders of magnitude. Successive substitution
cannot start either - the tangent-plane minimizer is an essentially pure
polymer melt, its K-values span ``1e+180``, and the Rachford-Rice window
``1/(1 - K_max) < beta < 1/(1 - K_min)`` collapses because every K falls on one
side of one (validation Case P-13 (vi), resolved as Case P-14).

What this module does instead is to carry the same minimization in
``u_i = ln n_i``:

- ``n_i = exp(u_i)`` is positive for every finite ``u``, so no floor and no box
  on the small side is needed; ``u ~ -410`` is an ordinary double.
- the phase-II mole fractions are formed as ``ln y_i = u_i - logsumexp(u)``, so
  a mole fraction below the exponential's range still has a finite, accurate
  logarithm. ``y_i`` itself may be an exact ``0.0``; see
  :func:`chemthermo.flash.tp.flash_tp` for what that means in a result.
- the objective, its gradient and the equal-fugacity residual are unchanged.
  Writing ``r_k`` for the residual of
  :func:`chemthermo.flash._second_order._second_order_split` equation (3),
  the gradient in the new variable is ``dg/du_k = n_k r_k`` by the chain rule,
  so a stationary point of one parametrization is a stationary point of the
  other and the two stages have the same fixed points.

The Newton step is taken on ``r(u) = 0`` rather than on ``dg/du = 0``. The two
systems differ by the diagonal factor ``diag(n)``: with ``J = dr/dn`` the true
Hessian is ``diag(n) J diag(n) + diag(n r)``, whose condition number is of
order ``n_max / n_min`` - ``1e+190`` for the polymer state, i.e. numerically
singular - while ``dr/du = J diag(n)`` has entries of order one there. Both
give the same Newton direction in exact arithmetic; only one of them survives
in doubles. Descent is still measured on ``g`` through
``(n r) . du``, and a direction that is not a descent direction is replaced by
``-r``, which is a successive substitution step in log space.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np

from ..exceptions import ModelError
from .settings import FlashSettings

#: Central-difference step in ``u`` for the log-space Hessian. It is a
#: *multiplicative* step on the mole numbers (``n -> n exp(+-h)``), which is
#: what makes one step size usable across mole numbers spanning ``1e+190``.
_HESSIAN_STEP = 1e-6
#: Smallest accepted backtracking scale.
_MIN_LINE_SEARCH_SCALE = 1e-14
#: Armijo constant on the two-phase Gibbs energy.
_ARMIJO_C = 1e-4
#: Phase-II mole fraction of the log-space seed (:func:`log_space_seed`). Any
#: value in ``(0, 1)`` gives an admissible seed - see that function - and the
#: stage is a descent method from wherever it starts; a half-and-half split is
#: the neutral choice between the two single-phase limits.
_SEED_PHASE_FRACTION = 0.5
#: Widest ``|ln K|`` the seed carries. A stationary point with a component at
#: an exact ``0.0`` gives ``ln K = -inf``; clamping keeps the seed finite
#: without changing any ``ln K`` that is representable at all.
_SEED_LN_K_CLAMP = 1.0e4

#: A mole fraction at or below this is "outside machine range" for the purposes
#: of choosing a parametrization (ADR-0024 decision 2). It is far below any
#: composition the linear stage has ever been asked to resolve - the smallest
#: mole fraction in a converged split anywhere in this repository is ``6.3e-10``
#: - and far above the smallest positive double, so a stationary point that
#: trips it is unambiguously the polymer/solvent geometry rather than a merely
#: dilute one.
TRACE_MOLE_FRACTION = 1e-30


class LogSpaceSplit:
    """Outcome of the log-space stage.

    Attributes:
        x_i: Phase-I composition (the reference phase, ``z - n``).
        x_ii: Phase-II composition (``n / sum n``). May contain exact zeros.
        ln_x_ii: ``ln x_ii``, finite for every component present in the feed
            even where ``x_ii`` itself underflowed to ``0.0``.
        beta: Phase-II mole fraction, ``sum_i n_i``.
        ln_f_i: Phase I's tangent-plane fugacity terms at ``x_i``.
        ln_f_ii: Phase II's terms at ``x_ii``.
        residual: ``max_i |ln(x_i^II f_i^II) - ln(x_i^I f_i^I)|``.
        iterations: Newton iterations actually taken.
        max_delta_k: Largest change in ``K = x^II / x^I`` over the final
            accepted step; ``inf`` when no step was taken.
    """

    __slots__ = (
        "beta",
        "iterations",
        "ln_f_i",
        "ln_f_ii",
        "ln_x_ii",
        "max_delta_k",
        "residual",
        "x_i",
        "x_ii",
    )

    def __init__(
        self,
        *,
        x_i: np.ndarray,
        x_ii: np.ndarray,
        ln_x_ii: np.ndarray,
        beta: float,
        ln_f_i: np.ndarray,
        ln_f_ii: np.ndarray,
        residual: float,
        iterations: int,
        max_delta_k: float,
    ) -> None:
        self.x_i = x_i
        self.x_ii = x_ii
        self.ln_x_ii = ln_x_ii
        self.beta = beta
        self.ln_f_i = ln_f_i
        self.ln_f_ii = ln_f_ii
        self.residual = residual
        self.iterations = iterations
        self.max_delta_k = max_delta_k


class _Iterate:
    """One admissible point of the log-space stage."""

    __slots__ = (
        "energy",
        "ln_f_i",
        "ln_f_ii",
        "ln_k",
        "ln_x_ii",
        "moles",
        "residual",
        "x_i",
        "x_ii",
    )

    def __init__(
        self,
        *,
        energy: float,
        residual: np.ndarray,
        x_i: np.ndarray,
        x_ii: np.ndarray,
        ln_x_ii: np.ndarray,
        ln_k: np.ndarray,
        moles: np.ndarray,
        ln_f_i: np.ndarray,
        ln_f_ii: np.ndarray,
    ) -> None:
        self.energy = energy
        self.residual = residual
        self.x_i = x_i
        self.x_ii = x_ii
        self.ln_x_ii = ln_x_ii
        self.ln_k = ln_k
        self.moles = moles
        self.ln_f_i = ln_f_i
        self.ln_f_ii = ln_f_ii


def _log_sum_exp(values: np.ndarray) -> float:
    """``ln sum_i exp(values_i)``, shifted by the maximum so it cannot overflow."""
    largest = float(np.max(values))
    if not math.isfinite(largest):
        raise ValueError("non-finite log mole numbers")
    with np.errstate(under="ignore"):
        total = float(np.sum(np.exp(values - largest)))
    if total <= 0.0:  # pragma: no cover - the shifted maximum contributes 1.0
        raise ValueError("empty phase")
    return largest + math.log(total)


def has_trace_component(composition: np.ndarray, active: np.ndarray) -> bool:
    """True when a component present in the feed is below :data:`TRACE_MOLE_FRACTION`."""
    values = np.asarray(composition, dtype=float)[active]
    return bool(values.size) and bool(np.min(values) <= TRACE_MOLE_FRACTION)


def log_space_seed(
    *,
    z: np.ndarray,
    w: np.ndarray,
    tpd_min: float,
    incipient_vapor: bool,
    ln_capital_w: np.ndarray | None = None,
) -> np.ndarray:
    """Phase-II log mole numbers from the tangent-plane stationary point.

    The stability test returns the normalized stationary composition ``w`` and
    the reduced tangent-plane distance there, from which Michelsen's
    unnormalized mole numbers are ``ln W_i = ln w_i - tpd``
    (:func:`chemthermo.flash._detect._stability_k_seed` uses the same identity
    in linear form). The K-values of that stationary point are ``K = W / z``
    when the incipient phase is the vapour-like one and ``K = z / W`` when it is
    the liquid-like one, and they are formed here **as logarithms**: for the
    polymer states this slice exists for they span ``e^{+-460}``, which is not a
    double.

    ``ln K`` alone is not a split. What turns it into one is the Rachford-Rice
    denominator ``t_i = (1 - beta) + beta K_i`` at a *fixed* ``beta``, rather
    than at the ``beta`` that solves Rachford-Rice:

        n_i = beta K_i z_i / t_i

    Two properties make this safe for any ``K > 0`` and any
    ``beta`` in ``(0, 1)``, and neither holds for the Rachford-Rice root:

    - ``t_i`` is a convex combination of ``1`` and ``K_i`` and so is strictly
      positive - there is no pole to step over and no window to bracket;
    - ``n_i / z_i = beta K_i / t_i`` lies strictly in ``(0, 1)``, so the seed
      satisfies the two-phase box componentwise by construction.

    Rachford-Rice itself has **no root at all** at such a stationary point: an
    essentially pure polymer melt has every ``K_i`` on one side of one, so the
    Leibovici-Neoschil window is empty and both the plain and the extended
    solver report "single phase" for a feed the same stability test has just
    proved unstable. That is why this seed exists rather than a better
    bracketing rule; see ADR-0024.

    Args:
        z: Feed mole fractions.
        w: Normalized tangent-plane stationary composition.
        tpd_min: Reduced tangent-plane distance at ``w``.
        incipient_vapor: Whether the stationary point is the vapour-like phase,
            i.e. whether phase II of the split is the incipient one.
        ln_capital_w: ``ln W`` at the stationary point, taken from the
            stability test rather than rebuilt from ``w`` and ``tpd_min``
            (ADR-0025). It is passed only where the rebuild cannot work - a
            stationary point whose ``w`` has already rounded to an exact
            ``0.0`` because ``ln W`` spans 1450 - and where it is passed the
            rebuild is not executed at all, so every seed that existed before
            ADR-0025 is formed by the same expression it was formed by then.
            The two agree to rounding wherever both are defined; ``ln W`` is
            the accurate one, being what the iteration carried.

    Returns:
        ``u = ln n`` for the components present in the feed; entries for absent
        components are ``-inf`` and are never read by the stage.
    """
    active = z > 0.0
    tiny = float(np.finfo(float).tiny)
    with np.errstate(divide="ignore"):
        ln_w = np.log(np.maximum(np.asarray(w, dtype=float), tiny))
        ln_z = np.log(np.where(active, z, 1.0))
    if ln_capital_w is None:
        ln_capital_w = ln_w - (tpd_min if math.isfinite(tpd_min) else 0.0)
    else:
        ln_capital_w = np.asarray(ln_capital_w, dtype=float)
    ln_k = (ln_capital_w - ln_z) if incipient_vapor else (ln_z - ln_capital_w)
    ln_k = np.clip(ln_k, -_SEED_LN_K_CLAMP, _SEED_LN_K_CLAMP)

    beta = _SEED_PHASE_FRACTION
    ln_t = np.logaddexp(math.log1p(-beta), math.log(beta) + ln_k)
    seed = math.log(beta) + ln_k + ln_z - ln_t
    # ``beta K / ((1 - beta) + beta K) < 1`` exactly, but it rounds to one for
    # a ``K`` past about ``e^36``; the ceiling keeps the *stored* seed strictly
    # inside the box, where the stage needs it.
    ceiling = ln_z + math.log1p(-1e-12)
    return np.where(active, np.minimum(seed, ceiling), -math.inf)


def seed_from_iterate(*, z: np.ndarray, x_ii: np.ndarray, beta: float) -> np.ndarray:
    """Phase-II log mole numbers from a linear iterate ``(x_ii, beta)``.

    Used when the linear stage has run and failed. ``beta`` outside ``(0, 1)``
    is a negative-flash iterate and is pulled back the same way
    :func:`chemthermo.flash._detect._phi_phi_second_order` pulls it back, then
    the mole numbers are taken logarithmically so that a component the linear
    stage had already driven under its floor keeps its magnitude.
    """
    active = z > 0.0
    values = np.asarray(x_ii, dtype=float)
    fraction = float(beta)
    if not 0.0 < fraction < 1.0 or np.any(fraction * values >= z):
        fraction = float(np.min(np.where(values > 0.0, z / np.maximum(values, 1e-300), 1.0))) * 0.5
        fraction = min(max(fraction, 1e-8), 1.0 - 1e-8)
    moles = fraction * values
    # Nothing below the smallest positive double survived the linear stage, so
    # a zero here is a component it could not represent at all; it is restarted
    # a factor 1e-8 below the feed rather than at minus infinity.
    floor = 1e-8 * z
    with np.errstate(divide="ignore"):
        seed = np.log(np.where(moles > 0.0, moles, floor))
    ceiling = np.log(np.where(active, z, 1.0)) + math.log1p(-1e-12)
    return np.where(active, np.minimum(seed, ceiling), -math.inf)


def log_space_split(
    *,
    z: np.ndarray,
    u0: np.ndarray,
    terms_i: Callable[[np.ndarray], np.ndarray],
    terms_ii: Callable[[np.ndarray], np.ndarray],
    settings: FlashSettings,
) -> LogSpaceSplit:
    """Damped Newton on the two-phase split in ``u = ln n`` (ADR-0024).

    The objective, the gradient and the acceptance test are those of
    :func:`chemthermo.flash._second_order._second_order_split`; only the
    variable changes. See the module docstring for why the Newton system is
    written on the residual rather than on the true gradient.

    Args:
        z: Feed mole fractions (normalized).
        u0: Starting log mole numbers of phase II; ``-inf`` for components
            absent from the feed.
        terms_i: Phase I's tangent-plane fugacity terms at a normalized
            composition.
        terms_ii: Phase II's terms. Called with a composition that may contain
            exact zeros.
        settings: ``second_order_tol`` and ``second_order_max_iter``.

    Returns:
        The converged (or best) split. The caller decides whether
        ``residual`` is good enough.

    Raises:
        ModelError: If the starting point itself is inadmissible - the box is
            violated, or neither phase can be evaluated there.
    """
    active = z > 0.0
    index = np.flatnonzero(active)
    size = index.size
    identity = np.eye(size)

    def evaluate(u: np.ndarray) -> _Iterate:
        with np.errstate(over="ignore", under="ignore"):
            moles = np.where(active, np.exp(u), 0.0)
        if not np.all(np.isfinite(moles)):
            raise ValueError("mole numbers left the representable range")
        reference = z - moles
        if np.any(reference[active] <= 0.0):
            raise ValueError("phase I emptied of a component")
        total_i = float(np.sum(reference[active]))
        if total_i <= 0.0:
            raise ValueError("phase I is empty")
        x_i = reference / total_i
        ln_x_ii = np.where(active, u - _log_sum_exp(u[active]), -math.inf)
        with np.errstate(under="ignore"):
            x_ii = np.where(active, np.exp(ln_x_ii), 0.0)
        ln_f_i = terms_i(x_i)
        ln_f_ii = terms_ii(x_ii)
        activity_i = np.log(x_i[active]) + ln_f_i[active]
        activity_ii = ln_x_ii[active] + ln_f_ii[active]
        energy = float(reference[active] @ activity_i + moles[active] @ activity_ii)
        if not math.isfinite(energy):
            raise ValueError("non-finite Gibbs energy")
        return _Iterate(
            energy=energy,
            residual=activity_ii - activity_i,
            x_i=x_i,
            x_ii=x_ii,
            ln_x_ii=ln_x_ii,
            ln_k=activity_ii - activity_i - ln_f_ii[active] + ln_f_i[active],
            moles=moles,
            ln_f_i=ln_f_i,
            ln_f_ii=ln_f_ii,
        )

    try:
        current = evaluate(np.asarray(u0, dtype=float))
    except (ValueError, ModelError) as exc:
        raise ModelError(f"The log-space split cannot start from this seed: {exc}.") from exc

    u = np.asarray(u0, dtype=float).copy()
    residual = float(np.max(np.abs(current.residual)))
    max_delta_k = math.inf
    iterations = 0

    for iteration in range(1, settings.second_order_max_iter + 1):
        if residual < settings.second_order_tol:
            break
        iterations = iteration

        jacobian = np.zeros((size, size), dtype=float)
        usable = True
        for column, component in enumerate(index):
            plus = u.copy()
            minus = u.copy()
            plus[component] += _HESSIAN_STEP
            minus[component] -= _HESSIAN_STEP
            try:
                forward = evaluate(plus).residual
                backward = evaluate(minus).residual
            except (ValueError, ModelError):
                usable = False
                break
            jacobian[:, column] = (forward - backward) / (2.0 * _HESSIAN_STEP)

        # ``-r`` is the log-space successive-substitution step and is the
        # documented fallback whenever the Newton direction is unavailable or
        # is not a descent direction for ``g``.
        direction = -current.residual
        if usable:
            try:
                newton = np.linalg.solve(jacobian, -current.residual)
            except np.linalg.LinAlgError:
                newton = np.linalg.solve(
                    jacobian + 1e-12 * identity, -current.residual
                )  # pragma: no cover - a singular finite-difference Jacobian
            if np.all(np.isfinite(newton)):
                direction = newton

        gradient = current.moles[index] * current.residual
        slope = float(gradient @ direction)
        if slope >= 0.0:
            direction = -current.residual
            slope = float(gradient @ direction)

        scale = 1.0
        accepted = False
        while scale >= _MIN_LINE_SEARCH_SCALE:
            candidate = u.copy()
            candidate[index] = u[index] + scale * direction
            try:
                trial = evaluate(candidate)
            except (ValueError, ModelError):
                scale *= 0.5
                continue
            trial_residual = float(np.max(np.abs(trial.residual)))
            if (
                trial.energy < current.energy + _ARMIJO_C * scale * slope
                or trial_residual < residual
            ):
                with np.errstate(over="ignore", under="ignore"):
                    max_delta_k = float(np.max(np.abs(np.exp(trial.ln_k) - np.exp(current.ln_k))))
                u = candidate
                current = trial
                residual = trial_residual
                accepted = True
                break
            scale *= 0.5

        if not accepted:
            break

    return LogSpaceSplit(
        x_i=current.x_i,
        x_ii=current.x_ii,
        ln_x_ii=current.ln_x_ii,
        beta=float(np.sum(current.moles[active])),
        ln_f_i=current.ln_f_i,
        ln_f_ii=current.ln_f_ii,
        residual=residual,
        iterations=iterations,
        max_delta_k=max_delta_k,
    )
