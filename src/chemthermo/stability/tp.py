"""Isothermal-isobaric phase stability by Michelsen's tangent-plane criterion.

Theory
------
References: M. L. Michelsen, "The isothermal flash problem. Part I. Stability",
Fluid Phase Equilibria 9 (1982) 1-19; M. L. Michelsen and J. M. Mollerup,
"Thermodynamic Models: Fundamentals and Computational Aspects", chapter on
stability analysis. The equations below are restated in our own notation; the
implementation is written from the mathematics, not ported from other code.

A feed of composition ``z`` at fixed ``(T, P)`` is thermodynamically stable when
the molar Gibbs energy surface lies entirely above (or on) the tangent
hyperplane constructed at ``z``. Writing the chemical potential as
``mu_i = g_i^0(T) + R T ln(x_i phi_i(x) P / P^0)``, the vertical distance from
the surface at a trial composition ``w`` to that plane, divided by ``R T``, is
the *reduced tangent plane distance*

    d_i    = ln z_i + ln phi_i(z)                                        (1)
    tpd(w) = sum_i w_i [ ln w_i + ln phi_i(w) - d_i ]                     (2)

with ``sum_i w_i = 1``. The feed is stable if and only if ``tpd(w) >= 0`` for
every admissible ``w``. Because the constrained minimization of (2) is awkward,
Michelsen removes the normalization constraint by introducing unnormalized mole
numbers ``W_i >= 0`` with ``w = W / sum_j W_j`` and minimizing

    tm(W) = 1 + sum_i W_i [ ln W_i + ln phi_i(w) - d_i - 1 ]              (3)

Using the Gibbs-Duhem relation ``sum_i w_i d ln phi_i = 0`` at fixed ``T, P``,
the gradient of (3) collapses to

    d tm / d W_k = ln W_k + ln phi_k(w) - d_k                             (4)

so the stationary points of (3) satisfy

    ln W_i + ln phi_i(w) - d_i = 0   for all i.                           (5)

Two consequences are used throughout this module and checked numerically in
``tests/test_stability_tp.py``. Substituting (5) into (3) gives

    tm* = 1 - sum_i W_i                                                   (6)

and, writing ``S = sum_i W_i`` and ``ln W_i = ln w_i + ln S``, condition (5)
gives ``ln w_i + ln phi_i(w) - d_i = -ln S`` for every ``i``, hence from (2)

    tpd(w) = -ln(S) = -ln(sum_i W_i)      (at a stationary point)         (7)

and therefore ``tm* = 1 - exp(-tpd)``. So ``sum_i W_i > 1`` <=> ``tpd < 0``
<=> ``tm* < 0`` <=> the feed is unstable, and the three criteria always agree in
sign. At the trivial solution ``W = z`` we have ``S = 1``, ``tpd = 0`` and
``tm = 0``.

Activity-coefficient models
---------------------------
For two *liquid* phases described by an activity-coefficient model at fixed
``T`` (and nominal ``P``), the reduced molar Gibbs energy of mixing is
``m(x) = sum_i x_i ln x_i + g^E(x)`` and the chemical potential is
``mu_i = mu_i^0(T, P) + R T ln(x_i gamma_i(x))``. Both phases share the same
pure-liquid reference ``mu_i^0``, so it cancels from the difference that defines
the tangent-plane distance and equations (1)-(7) hold verbatim with

    ln gamma_i(w)   in place of   ln phi_i(w).

Nothing else changes: the same ``tm(W)``, the same stationarity condition (5),
the same relations (6) and (7). This is exactly the tangent-plane distance
``D(x)`` of Tessier, Brennecke & Stadtherr, Chem. Eng. Sci. 55 (2000) 1785,
section 2, whose published stationary points are reproduced in
``tests/validation/test_stability_nrtl_tessier2000.py``.

The one thing that is *not* shared between the two families is the fugacity
term itself and the set of sensible initial estimates. Those two differences
are isolated behind the internal evaluator contract in
``chemthermo.stability._evaluator`` (ADR-0007), so the solver below never
branches on the model family.

Combining an activity model for the liquid with an EOS for the vapor
(gamma-phi stability) is deliberately **not** supported and raises
``ModelError``; see ADR-0007 for the reason.

Root selection (equation-of-state models only)
----------------------------------------------
For a cubic equation of state the fugacity coefficients are multivalued: each
real compressibility root gives a different ``phi``. Michelsen and Mollerup
require the root of *lowest Gibbs energy* at that ``(T, P, composition)`` for
both the feed and every trial evaluation, otherwise the tangent plane itself is
not the physical one. At fixed ``T``, ``P`` and ``w`` the molar Gibbs energy is

    G/RT = sum_i w_i [ g_i^0/RT + ln(w_i P / P^0) ] + sum_i w_i ln phi_i(w)

and only the last term depends on which root was taken (the ideal-mixing part is
root independent). The last term is exactly the residual Gibbs energy
``G^res/(R T)``. Hence the minimum-Gibbs root is the one minimizing
``sum_i w_i ln phi_i(w)``. The EOS evaluator therefore evaluates the model's
``fugacity_coefficients`` with ``phase="vapor"`` and ``phase="liquid"`` and
keeps the branch with the smaller ``sum_i w_i ln phi_i``; when only one real
root exists both calls return the same values and the choice is immaterial. The
selected branch is recorded in the result diagnostics. An activity-coefficient
model has a single branch, so no selection is performed and the branch labels
are None.

Trial sets
----------
For an equation of state the trial set is the two Wilson-based estimates
``W = K_wilson * z`` (vapor-like) and ``W = z / K_wilson`` (liquid-like), plus
one pure-component-dominant trial per component. For an activity model the
Wilson estimates are dropped: they are a vapor-liquid construction built from
``Tc``, ``Pc`` and ``omega`` and carry no information about a liquid-liquid
split. Pure-component-dominant trials are what Michelsen (1982) recommends for
liquid-liquid stability, and later robustness work (for example Li and
Firoozabadi, AIChE J. 58 (2012) 2244-2258) reaches the same conclusion. Every
trial set is deterministic; there are no randomized or adaptive restarts.

Two-stage solution
------------------
Stage 1 is successive substitution on equation (5),

    ln W_i^(k+1) = d_i - ln phi_i(w^(k)),    w^(k) = W^(k) / sum_j W_j^(k) (8)

which is the fixed-point form Michelsen recommends for the first stage of a
stability test. Iteration stops on the residual of (5), or early when the trial
collapses onto the trivial solution ``w -> z`` (detected through
``sum_i (ln(W_i / z_i))**2 < trivial_tol``).

Near a plait point the fixed-point map (8) has a contraction ratio close to one:
it converges linearly and very slowly, or drifts to the trivial solution.
Michelsen and Mollerup therefore recommend switching to a second-order method
after a few substitutions. Stage 2 here is a damped Newton solve of the
stationarity system (5) in the variables ``u_i = ln W_i`` (which keeps
``W_i > 0`` for free), with a central-difference Jacobian of

    g_i(u) = u_i + ln phi_i(w(u)) - d_i,     w(u) = exp(u) / sum_j exp(u_j)

a cap on ``|delta u|`` and a backtracking line search on ``max_i |g_i|``. It is
entered only after ``settings.ssi_iterations`` substitutions have failed to meet
``settings.tol``, so trials that converge quickly are untouched by it.

Limits of the test
------------------
A negative tangent-plane distance is a *proof* of instability. The converse is
not true here: reporting ``stable`` only means that no negative tangent-plane
distance was found from this deterministic trial set. Michelsen's test is a
local stationary-point search; a stationary point that no initial estimate
reaches can hide an instability.
"""

from __future__ import annotations

import math

import numpy as np

from ..core import Mixture
from ..exceptions import CompositionError, ModelError
from ..models import ActivityModel, EquationOfState
from ..validation import validate_pressure, validate_temperature
from ._evaluator import (
    _ActivityTangentPlane,
    _EOSTangentPlane,
    _ln_phi_min_gibbs,  # noqa: F401  (re-exported: tests import root selection from here)
    _TangentPlaneEvaluator,
)
from .results import StabilityResult, StabilityTrial
from .settings import StabilitySettings

_LN_W_MIN = -700.0
_LN_W_MAX = 700.0
_JACOBIAN_STEP = 1e-6
_MIN_LINE_SEARCH_SCALE = 1e-12

__all__ = ["stability_tp"]


def stability_tp(
    mixture: Mixture,
    *,
    temperature_K: float,
    pressure_Pa: float,
    eos: EquationOfState | None = None,
    activity_model: ActivityModel | None = None,
    settings: StabilitySettings | None = None,
) -> StabilityResult:
    """Test a feed for phase stability at fixed temperature, pressure and composition.

    Args:
        mixture: Mixture with a mole-fraction composition (the feed ``z``).
        temperature_K: Temperature in K.
        pressure_Pa: Pressure in Pa. Required for API uniformity and always
            validated; it does not affect the result when ``activity_model`` is
            used (``diagnostics["pressure_dependent"]`` records which case
            applies).
        eos: Equation-of-state model supplying fugacity coefficients. Both the
            ``"vapor"`` and ``"liquid"`` branches are evaluated and the
            minimum-Gibbs branch is used (see the module docstring).
        activity_model: Activity-coefficient model supplying ``gamma`` for a
            liquid-liquid stability test. Exactly one of ``eos`` and
            ``activity_model`` must be given.
        settings: Iteration controls; defaults to ``StabilitySettings()``.

    Returns:
        StabilityResult with the verdict, the minimum reduced tangent-plane
        distance, the minimizing trial composition, the implied incipient-phase
        K-values (``w_i / z_i``) and a per-trial record.

    Raises:
        InputRangeError: If temperature or pressure is non-physical.
        ModelError: If neither or both of ``eos`` and ``activity_model`` are
            given (the combined gamma-phi case is not yet supported), if the
            composition basis is not molar, or if the model returns non-finite
            or non-positive values at the feed composition.
        CompositionError: If the mixture composition is empty or non-positive.

    Notes:
        ``status`` is ``"unstable"`` when some trial converged to a non-trivial
        stationary point with ``tpd < -settings.tpd_tol``; ``"stable"`` when at
        least one trial converged and none did; ``"inconclusive"`` when no trial
        converged at all. See the honesty note on
        :class:`chemthermo.StabilityResult`.

        The analysis is deterministic for fixed inputs and settings.
    """

    temperature = validate_temperature(temperature_K)
    pressure = validate_pressure(pressure_Pa)

    if eos is not None and activity_model is not None:
        raise ModelError(
            "stability_tp accepts exactly one of 'eos' and 'activity_model'. Combined "
            "gamma-phi stability (activity-coefficient liquid against an equation-of-state "
            "vapor) is not yet supported."
        )
    if eos is None and activity_model is None:
        raise ModelError("stability_tp requires a model: pass either 'eos' or 'activity_model'.")
    if mixture.basis != "mole":
        raise ModelError("stability_tp currently requires mole-fraction compositions.")

    settings = settings or StabilitySettings()

    z_raw = np.array(mixture.fractions, dtype=float)
    if z_raw.size == 0:
        raise CompositionError("Mixture composition must be non-empty.")
    if np.any(z_raw < 0.0):
        raise CompositionError("Feed composition fractions must be non-negative.")
    total = float(np.sum(z_raw))
    if total <= 0.0:
        raise CompositionError("Feed composition fractions must sum to a positive value.")
    z = z_raw / total

    evaluator: _TangentPlaneEvaluator
    if eos is not None:
        evaluator = _EOSTangentPlane(
            eos, mixture=mixture, temperature=temperature, pressure=pressure
        )
    else:
        assert activity_model is not None
        evaluator = _ActivityTangentPlane(activity_model, mixture=mixture, temperature=temperature)

    active = z > 0.0
    n_active = int(np.count_nonzero(active))

    ln_f_feed, feed_branch = evaluator.ln_fugacity_terms(z)

    # d_i = ln z_i + ln phi_i(z); inactive components are excluded rather than
    # floored, which keeps W_i = 0 for them at every iteration (equation (8)).
    d = np.full(z.shape, -np.inf, dtype=float)
    d[active] = np.log(z[active]) + ln_f_feed[active]

    trials: list[StabilityTrial] = []
    for label, w0 in evaluator.initial_estimates(z, active):
        trials.append(
            _run_trial(
                label=label,
                w0=w0,
                z=z,
                d=d,
                active=active,
                evaluator=evaluator,
                settings=settings,
            )
        )

    return _summarize(
        z=z,
        temperature=temperature,
        pressure=pressure,
        feed_branch=feed_branch,
        trials=tuple(trials),
        settings=settings,
        n_active=n_active,
        evaluator=evaluator,
    )


def _run_trial(
    *,
    label: str,
    w0: np.ndarray,
    z: np.ndarray,
    d: np.ndarray,
    active: np.ndarray,
    evaluator: _TangentPlaneEvaluator,
    settings: StabilitySettings,
) -> StabilityTrial:
    """Successive substitution on equation (8), then an optional Newton stage."""
    w = w0
    ln_w_capital = np.where(active, np.log(np.where(active, w, 1.0)), -np.inf)
    branch: str | None = None
    residual = math.inf
    sum_w_capital = 1.0
    iterations = 0

    ssi_budget = (
        min(settings.ssi_iterations, settings.max_iter)
        if settings.second_order
        else settings.max_iter
    )

    for iteration in range(1, ssi_budget + 1):
        iterations = iteration
        try:
            ln_f, branch = evaluator.ln_fugacity_terms(w)
        except ModelError as exc:
            return _failed_trial(label, iterations, f"model_error: {exc}")

        residual = float(np.max(np.abs(ln_w_capital[active] + ln_f[active] - d[active])))
        if residual < settings.tol:
            tpd = _tpd(w, ln_f, d, active)
            return StabilityTrial(
                label=label,
                converged=True,
                iterations=iterations,
                tpd=tpd,
                sum_W=sum_w_capital,
                trivial=_is_trivial(ln_w_capital, z, active, settings.trivial_tol),
                residual=residual,
                phase_branch=branch,
                composition=tuple(w.tolist()),
                termination_reason="stationarity_met",
                ssi_iterations=iterations,
                second_order_iterations=0,
                converged_stage="successive-substitution",
            )

        # Equation (8): ln W_i <- d_i - ln phi_i(w).
        ln_w_new = np.full(z.shape, -np.inf, dtype=float)
        ln_w_new[active] = np.clip(d[active] - ln_f[active], _LN_W_MIN, _LN_W_MAX)
        if np.any(~np.isfinite(ln_w_new[active])):
            return _failed_trial(label, iterations, "non_finite_ln_W")

        w_capital = np.where(active, np.exp(np.where(active, ln_w_new, 0.0)), 0.0)
        sum_w_capital = float(np.sum(w_capital))
        if not math.isfinite(sum_w_capital) or sum_w_capital <= 0.0:
            return _failed_trial(label, iterations, "degenerate_sum_W")

        ln_w_capital = ln_w_new
        w = w_capital / sum_w_capital

        if _is_trivial(ln_w_capital, z, active, settings.trivial_tol):
            try:
                ln_f, branch = evaluator.ln_fugacity_terms(w)
            except ModelError as exc:
                return _failed_trial(label, iterations, f"model_error: {exc}")
            return StabilityTrial(
                label=label,
                converged=True,
                iterations=iterations,
                tpd=_tpd(w, ln_f, d, active),
                sum_W=sum_w_capital,
                trivial=True,
                residual=float(np.max(np.abs(ln_w_capital[active] + ln_f[active] - d[active]))),
                phase_branch=branch,
                composition=tuple(w.tolist()),
                termination_reason="trivial_solution",
                ssi_iterations=iterations,
                second_order_iterations=0,
                converged_stage="successive-substitution",
            )

    if not settings.second_order:
        return StabilityTrial(
            label=label,
            converged=False,
            iterations=iterations,
            tpd=float("nan"),
            sum_W=sum_w_capital,
            trivial=False,
            residual=residual,
            phase_branch=branch,
            composition=tuple(w.tolist()),
            termination_reason="max_iter",
            ssi_iterations=iterations,
            second_order_iterations=0,
            converged_stage=None,
        )

    return _second_order_stage(
        label=label,
        ln_w_capital=ln_w_capital,
        z=z,
        d=d,
        active=active,
        evaluator=evaluator,
        settings=settings,
        ssi_iterations=iterations,
    )


def _second_order_stage(
    *,
    label: str,
    ln_w_capital: np.ndarray,
    z: np.ndarray,
    d: np.ndarray,
    active: np.ndarray,
    evaluator: _TangentPlaneEvaluator,
    settings: StabilitySettings,
    ssi_iterations: int,
) -> StabilityTrial:
    """Damped Newton solve of equation (5) in the variables ``u = ln W``.

    The gradient of ``tm`` with respect to ``W_k`` is equation (4); in the
    ``ln W`` variables the residual solved here is

        g_i(u) = u_i + ln phi_i(w(u)) - d_i,   w(u) = exp(u) / sum_j exp(u_j)

    whose Jacobian ``dg/du`` is built by central differences, so it shares no
    analytic derivative code with the model. Steps are capped at
    ``settings.second_order_max_step`` and then backtracked until
    ``max_i |g_i|`` decreases, which keeps ``W > 0`` (automatic in ``ln W``) and
    prevents the iteration from being thrown out of the model's domain.
    """
    index = np.flatnonzero(active)
    d_active = d[index]

    def evaluate(u: np.ndarray) -> tuple[np.ndarray, np.ndarray, str | None]:
        ln_w_full = np.full(z.shape, -np.inf, dtype=float)
        ln_w_full[index] = u
        w_capital = np.zeros(z.shape, dtype=float)
        w_capital[index] = np.exp(np.clip(u, _LN_W_MIN, _LN_W_MAX))
        total = float(np.sum(w_capital))
        if not math.isfinite(total) or total <= 0.0:
            raise ModelError("Second-order stage produced a degenerate sum of mole numbers.")
        w_local = w_capital / total
        terms, branch_local = evaluator.ln_fugacity_terms(w_local)
        return u + terms[index] - d_active, w_local, branch_local

    u = ln_w_capital[index].copy()
    iterations = 0
    reason = "second_order_max_iter"

    try:
        g, w, branch = evaluate(u)
    except ModelError as exc:
        return _failed_trial(label, ssi_iterations, f"model_error: {exc}")
    residual = float(np.max(np.abs(g)))

    for iteration in range(1, settings.second_order_max_iter + 1):
        if residual < settings.tol:
            reason = "stationarity_met"
            break
        iterations = iteration

        jacobian = np.zeros((index.size, index.size), dtype=float)
        try:
            for column in range(index.size):
                plus = u.copy()
                minus = u.copy()
                plus[column] += _JACOBIAN_STEP
                minus[column] -= _JACOBIAN_STEP
                g_plus, _, _ = evaluate(plus)
                g_minus, _, _ = evaluate(minus)
                jacobian[:, column] = (g_plus - g_minus) / (2.0 * _JACOBIAN_STEP)
        except ModelError as exc:
            return _failed_trial(label, ssi_iterations + iterations, f"model_error: {exc}")

        try:
            step = np.linalg.solve(jacobian, -g)
        except np.linalg.LinAlgError:
            reason = "second_order_singular_jacobian"
            break
        if np.any(~np.isfinite(step)):
            reason = "second_order_singular_jacobian"
            break

        largest = float(np.max(np.abs(step)))
        if largest > settings.second_order_max_step:
            step = step * (settings.second_order_max_step / largest)

        scale = 1.0
        accepted = False
        while scale >= _MIN_LINE_SEARCH_SCALE:
            candidate = np.clip(u + scale * step, _LN_W_MIN, _LN_W_MAX)
            try:
                g_trial, w_trial, branch_trial = evaluate(candidate)
            except ModelError:
                scale *= 0.5
                continue
            residual_trial = float(np.max(np.abs(g_trial)))
            if residual_trial < residual:
                u, g, w, branch, residual = (
                    candidate,
                    g_trial,
                    w_trial,
                    branch_trial,
                    residual_trial,
                )
                accepted = True
                break
            scale *= 0.5

        if not accepted:
            reason = "second_order_no_progress"
            break
    else:
        if residual < settings.tol:
            reason = "stationarity_met"

    ln_w_capital_final = np.full(z.shape, -np.inf, dtype=float)
    ln_w_capital_final[index] = u
    sum_w_capital = float(np.sum(np.exp(np.clip(u, _LN_W_MIN, _LN_W_MAX))))
    converged = residual < settings.tol
    trivial = converged and _is_trivial(ln_w_capital_final, z, active, settings.trivial_tol)

    try:
        ln_f, branch = evaluator.ln_fugacity_terms(w)
    except ModelError as exc:
        return _failed_trial(label, ssi_iterations + iterations, f"model_error: {exc}")

    return StabilityTrial(
        label=label,
        converged=converged,
        iterations=ssi_iterations + iterations,
        tpd=_tpd(w, ln_f, d, active) if converged else float("nan"),
        sum_W=sum_w_capital,
        trivial=trivial,
        residual=residual,
        phase_branch=branch,
        composition=tuple(w.tolist()),
        termination_reason=reason,
        ssi_iterations=ssi_iterations,
        second_order_iterations=iterations,
        converged_stage="second-order" if converged else None,
    )


def _failed_trial(label: str, iterations: int, reason: str) -> StabilityTrial:
    return StabilityTrial(
        label=label,
        converged=False,
        iterations=iterations,
        tpd=float("nan"),
        sum_W=float("nan"),
        trivial=False,
        residual=float("nan"),
        phase_branch=None,
        composition=None,
        termination_reason=reason,
        ssi_iterations=iterations,
        second_order_iterations=0,
        converged_stage=None,
    )


def _tpd(w: np.ndarray, ln_f: np.ndarray, d: np.ndarray, active: np.ndarray) -> float:
    """Reduced tangent-plane distance, equation (2), over active components."""
    mask = active & (w > 0.0)
    if not np.any(mask):
        return float("nan")
    return float(np.sum(w[mask] * (np.log(w[mask]) + ln_f[mask] - d[mask])))


def _is_trivial(
    ln_w_capital: np.ndarray, z: np.ndarray, active: np.ndarray, trivial_tol: float
) -> bool:
    """True when ln(W_i / z_i) -> 0 for all active components."""
    ln_k = ln_w_capital[active] - np.log(z[active])
    return bool(float(np.sum(ln_k * ln_k)) < trivial_tol)


def _summarize(
    *,
    z: np.ndarray,
    temperature: float,
    pressure: float,
    feed_branch: str | None,
    trials: tuple[StabilityTrial, ...],
    settings: StabilitySettings,
    n_active: int,
    evaluator: _TangentPlaneEvaluator,
) -> StabilityResult:
    """Reduce per-trial outcomes to a verdict and the minimizing trial."""
    converged = [trial for trial in trials if trial.converged]
    non_trivial = [trial for trial in converged if not trial.trivial and math.isfinite(trial.tpd)]

    best: StabilityTrial | None = None
    for trial in non_trivial:
        if best is None or trial.tpd < best.tpd:
            best = trial

    if best is not None and best.tpd < -settings.tpd_tol:
        status = "unstable"
    elif converged:
        status = "stable"
    else:
        status = "inconclusive"

    tpd_min = best.tpd if best is not None else 0.0

    trial_composition: tuple[float, ...] | None = None
    k_values: tuple[float, ...] | None = None
    phase_branch: str | None = None
    if best is not None and best.composition is not None:
        trial_composition = best.composition
        phase_branch = best.phase_branch
        w = np.array(best.composition, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = np.where(z > 0.0, w / np.where(z > 0.0, z, 1.0), 0.0)
        k_values = tuple(float(value) for value in ratios)

    diagnostics: dict[str, float | int | str | bool] = {
        "trial_count": len(trials),
        "converged_trial_count": len(converged),
        "non_trivial_trial_count": len(non_trivial),
        "trivial_trial_count": sum(1 for trial in converged if trial.trivial),
        "active_component_count": n_active,
        "model_family": evaluator.model_family,
        "pressure_dependent": evaluator.pressure_dependent,
        "status": status,
        "tpd_tol": settings.tpd_tol,
        "tol": settings.tol,
        "trivial_tol": settings.trivial_tol,
        "max_iter": settings.max_iter,
        "second_order_enabled": settings.second_order,
        "ssi_iterations_budget": settings.ssi_iterations,
        "second_order_trial_count": sum(1 for trial in trials if trial.second_order_iterations > 0),
    }
    if feed_branch is not None:
        diagnostics["feed_branch"] = feed_branch
    if best is not None:
        diagnostics["minimizing_trial"] = best.label
        diagnostics["minimizing_trial_iterations"] = best.iterations
        diagnostics["minimizing_trial_ssi_iterations"] = best.ssi_iterations
        diagnostics["minimizing_trial_second_order_iterations"] = best.second_order_iterations
        diagnostics["minimizing_trial_residual"] = best.residual
        if best.converged_stage is not None:
            diagnostics["minimizing_trial_stage"] = best.converged_stage
        diagnostics["sum_W"] = best.sum_W
        if best.sum_W > 0.0 and math.isfinite(best.sum_W):
            # Equation (7); agrees with tpd at a stationary point.
            diagnostics["tpd_from_sum_W"] = -math.log(best.sum_W)
            diagnostics["tm_at_stationary_point"] = 1.0 - best.sum_W

    return StabilityResult(
        temperature_K=temperature,
        pressure_Pa=pressure,
        feed_composition=tuple(z.tolist()),
        stable=status == "stable",
        status=status,
        tpd_min=tpd_min,
        trial_composition=trial_composition,
        k_values=k_values,
        phase_branch=phase_branch,
        feed_branch=feed_branch,
        trials=trials,
        diagnostics=diagnostics,
    )
