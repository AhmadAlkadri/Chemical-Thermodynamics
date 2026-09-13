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

Condition (5) is solved by successive substitution,

    ln W_i^(k+1) = d_i - ln phi_i(w^(k)),    w^(k) = W^(k) / sum_j W_j^(k) (8)

which is the fixed-point form Michelsen recommends for the first stage of a
stability test. Iteration stops on the residual of (5), or early when the trial
collapses onto the trivial solution ``w -> z`` (detected through
``sum_i (ln(W_i / z_i))**2 < trivial_tol``).

Root selection
--------------
For a cubic equation of state the fugacity coefficients are multivalued: each
real compressibility root gives a different ``phi``. Michelsen and Mollerup
require the root of *lowest Gibbs energy* at that ``(T, P, composition)`` for
both the feed and every trial evaluation, otherwise the tangent plane itself is
not the physical one. At fixed ``T``, ``P`` and ``w`` the molar Gibbs energy is

    G/RT = sum_i w_i [ g_i^0/RT + ln(w_i P / P^0) ] + sum_i w_i ln phi_i(w)

and only the last term depends on which root was taken (the ideal-mixing part is
root independent). The last term is exactly the residual Gibbs energy
``G^res/(R T)``. Hence the minimum-Gibbs root is the one minimizing
``sum_i w_i ln phi_i(w)``. This module therefore evaluates the model's
``fugacity_coefficients`` with ``phase="vapor"`` and ``phase="liquid"`` and keeps
the branch with the smaller ``sum_i w_i ln phi_i``; when only one real root
exists both calls return the same values and the choice is immaterial. The
selected branch is recorded in the result diagnostics. This is done generically
over the existing ``EquationOfState`` interface, so no model-specific hook is
required.

Trial set
---------
The trial set is deterministic and small: the two Wilson-based estimates
``W = K_wilson * z`` (vapor-like) and ``W = z / K_wilson`` (liquid-like), plus one
pure-component-dominant trial per component. Pure-component trials are
recommended by Michelsen (1982) and by later robustness work (for example
Li and Firoozabadi, AIChE J. 58 (2012) 2244-2258) because Wilson estimates alone
can miss liquid-liquid stationary points.

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
from ..flash._common import as_float_array, wilson_k
from ..models import EquationOfState
from ..validation import validate_pressure, validate_temperature
from .results import StabilityResult, StabilityTrial
from .settings import StabilitySettings

_LN_W_MIN = -700.0
_LN_W_MAX = 700.0
_PURE_TRIAL_TRACE = 1e-3

__all__ = ["stability_tp"]


def stability_tp(
    mixture: Mixture,
    *,
    temperature_K: float,
    pressure_Pa: float,
    eos: EquationOfState,
    settings: StabilitySettings | None = None,
) -> StabilityResult:
    """Test a feed for phase stability at fixed temperature, pressure and composition.

    Args:
        mixture: Mixture with a mole-fraction composition (the feed ``z``).
        temperature_K: Temperature in K.
        pressure_Pa: Pressure in Pa.
        eos: Equation-of-state model supplying fugacity coefficients. Both the
            ``"vapor"`` and ``"liquid"`` branches are evaluated and the
            minimum-Gibbs branch is used (see the module docstring).
        settings: Iteration controls; defaults to ``StabilitySettings()``.

    Returns:
        StabilityResult with the verdict, the minimum reduced tangent-plane
        distance, the minimizing trial composition, the implied incipient-phase
        K-values (``w_i / z_i``) and a per-trial record.

    Raises:
        InputRangeError: If temperature or pressure is non-physical.
        ModelError: If the EOS is missing, the composition basis is not molar,
            or the model returns non-finite / non-positive fugacity coefficients
            at the feed composition.
        CompositionError: If the mixture composition is empty or non-positive.

    Notes:
        ``status`` is ``"unstable"`` when some trial converged to a non-trivial
        stationary point with ``tpd < -settings.tpd_tol``; ``"stable"`` when at
        least one trial converged and none did; ``"inconclusive"`` when no trial
        converged at all within ``settings.max_iter``. See the honesty note on
        :class:`chemthermo.StabilityResult`.

        The analysis is deterministic for fixed inputs and settings.
    """

    temperature = validate_temperature(temperature_K)
    pressure = validate_pressure(pressure_Pa)

    if eos is None:
        raise ModelError("An equation-of-state model is required for stability_tp.")
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

    active = z > 0.0
    n_active = int(np.count_nonzero(active))

    ln_phi_feed, feed_branch = _ln_phi_min_gibbs(
        eos,
        mixture=mixture,
        temperature=temperature,
        pressure=pressure,
        composition=z,
    )

    # d_i = ln z_i + ln phi_i(z); inactive components are excluded rather than
    # floored, which keeps W_i = 0 for them at every iteration (equation (8)).
    d = np.full(z.shape, -np.inf, dtype=float)
    d[active] = np.log(z[active]) + ln_phi_feed[active]

    trials: list[StabilityTrial] = []
    for label, w0 in _initial_estimates(mixture, z, active, temperature, pressure):
        trials.append(
            _run_trial(
                label=label,
                w0=w0,
                z=z,
                d=d,
                active=active,
                mixture=mixture,
                temperature=temperature,
                pressure=pressure,
                eos=eos,
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
    )


def _ln_phi_min_gibbs(
    eos: EquationOfState,
    *,
    mixture: Mixture,
    temperature: float,
    pressure: float,
    composition: np.ndarray,
) -> tuple[np.ndarray, str]:
    """Return ``ln phi`` on the minimum-Gibbs root branch and the branch name.

    The branch minimizing ``sum_i w_i ln phi_i(w)`` minimizes the molar Gibbs
    energy at fixed ``(T, P, w)`` because that sum is the reduced residual Gibbs
    energy and the ideal-mixing contribution is root independent.
    """
    best_ln_phi: np.ndarray | None = None
    best_branch = ""
    best_g = math.inf
    failures: list[str] = []

    for branch in ("vapor", "liquid"):
        try:
            values = as_float_array(
                eos.fugacity_coefficients(
                    mixture=mixture,
                    temperature_K=temperature,
                    pressure_Pa=pressure,
                    composition=composition.tolist(),
                    phase=branch,
                )
            )
        except Exception as exc:  # noqa: BLE001 - branch may be unavailable
            failures.append(f"{branch}: {exc}")
            continue

        if values.shape != composition.shape:
            failures.append(f"{branch}: inconsistent fugacity coefficient shape")
            continue
        if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
            failures.append(f"{branch}: non-finite or non-positive fugacity coefficients")
            continue

        ln_phi = np.log(values)
        g_res = float(np.sum(composition * ln_phi))
        if not math.isfinite(g_res):
            failures.append(f"{branch}: non-finite reduced residual Gibbs energy")
            continue
        if g_res < best_g:
            best_g = g_res
            best_ln_phi = ln_phi
            best_branch = branch

    if best_ln_phi is None:
        raise ModelError(
            "No usable fugacity-coefficient branch for stability analysis ("
            + "; ".join(failures)
            + ")."
        )
    return best_ln_phi, best_branch


def _initial_estimates(
    mixture: Mixture,
    z: np.ndarray,
    active: np.ndarray,
    temperature: float,
    pressure: float,
) -> list[tuple[str, np.ndarray]]:
    """Deterministic trial-phase initial estimates (normalized compositions)."""
    estimates: list[tuple[str, np.ndarray]] = []

    k = wilson_k(mixture, temperature, pressure)
    estimates.append(("wilson-vapor", _normalized(k * z, active)))
    estimates.append(("wilson-liquid", _normalized(z / k, active)))

    n_active = int(np.count_nonzero(active))
    if n_active > 1:
        trace = _PURE_TRIAL_TRACE / (n_active - 1)
        names = mixture.component_names
        for index in range(z.size):
            if not active[index]:
                continue
            w = np.where(active, trace, 0.0)
            w[index] = 1.0 - _PURE_TRIAL_TRACE
            estimates.append((f"pure-{names[index]}", _normalized(w, active)))

    return estimates


def _normalized(values: np.ndarray, active: np.ndarray) -> np.ndarray:
    """Zero out inactive components and renormalize to sum to one."""
    w = np.where(active, np.asarray(values, dtype=float), 0.0)
    w = np.where(w > 0.0, w, 0.0)
    total = float(np.sum(w))
    if not math.isfinite(total) or total <= 0.0:
        raise ModelError("Stability trial initial estimate has a non-positive total.")
    return w / total


def _run_trial(
    *,
    label: str,
    w0: np.ndarray,
    z: np.ndarray,
    d: np.ndarray,
    active: np.ndarray,
    mixture: Mixture,
    temperature: float,
    pressure: float,
    eos: EquationOfState,
    settings: StabilitySettings,
) -> StabilityTrial:
    """Successive substitution on equation (8) from one initial estimate."""
    w = w0
    ln_w_capital = np.where(active, np.log(np.where(active, w, 1.0)), -np.inf)
    branch: str | None = None
    residual = math.inf
    sum_w_capital = 1.0
    iterations = 0

    for iteration in range(1, settings.max_iter + 1):
        iterations = iteration
        try:
            ln_phi, branch = _ln_phi_min_gibbs(
                eos,
                mixture=mixture,
                temperature=temperature,
                pressure=pressure,
                composition=w,
            )
        except ModelError as exc:
            return _failed_trial(label, iterations, f"model_error: {exc}")

        residual = float(np.max(np.abs(ln_w_capital[active] + ln_phi[active] - d[active])))
        if residual < settings.tol:
            tpd = _tpd(w, ln_phi, d, active)
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
            )

        # Equation (8): ln W_i <- d_i - ln phi_i(w).
        ln_w_new = np.full(z.shape, -np.inf, dtype=float)
        ln_w_new[active] = np.clip(d[active] - ln_phi[active], _LN_W_MIN, _LN_W_MAX)
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
                ln_phi, branch = _ln_phi_min_gibbs(
                    eos,
                    mixture=mixture,
                    temperature=temperature,
                    pressure=pressure,
                    composition=w,
                )
            except ModelError as exc:
                return _failed_trial(label, iterations, f"model_error: {exc}")
            return StabilityTrial(
                label=label,
                converged=True,
                iterations=iterations,
                tpd=_tpd(w, ln_phi, d, active),
                sum_W=sum_w_capital,
                trivial=True,
                residual=float(np.max(np.abs(ln_w_capital[active] + ln_phi[active] - d[active]))),
                phase_branch=branch,
                composition=tuple(w.tolist()),
                termination_reason="trivial_solution",
            )

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
    )


def _tpd(w: np.ndarray, ln_phi: np.ndarray, d: np.ndarray, active: np.ndarray) -> float:
    """Reduced tangent-plane distance, equation (2), over active components."""
    mask = active & (w > 0.0)
    if not np.any(mask):
        return float("nan")
    return float(np.sum(w[mask] * (np.log(w[mask]) + ln_phi[mask] - d[mask])))


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
    feed_branch: str,
    trials: tuple[StabilityTrial, ...],
    settings: StabilitySettings,
    n_active: int,
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
        "feed_branch": feed_branch,
        "status": status,
        "tpd_tol": settings.tpd_tol,
        "tol": settings.tol,
        "trivial_tol": settings.trivial_tol,
        "max_iter": settings.max_iter,
    }
    if best is not None:
        diagnostics["minimizing_trial"] = best.label
        diagnostics["minimizing_trial_iterations"] = best.iterations
        diagnostics["minimizing_trial_residual"] = best.residual
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
