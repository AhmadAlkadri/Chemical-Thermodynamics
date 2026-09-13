"""The shared two-phase split loop: Rachford-Rice plus successive substitution.

This is the first-stage phase split used, unchanged, by all three
``flash_tp`` modes (phi-phi, gamma-phi, gamma-gamma): only the initial ``K``,
the initial ``beta`` and the model call that updates ``K`` differ between
them (see :func:`_solve_k_loop`). Callers seed and interpret the loop
differently (:mod:`chemthermo.flash._detect`, :mod:`chemthermo.flash._legacy`);
this module owns none of that, only the shared numerics.

The phi-phi update also owns one piece of physics rather than pure numerics:
:class:`_PhaseRoot`, which puts each phase on the density/compressibility root
the stability test found **that phase** on, instead of pinning phase I to the
model's liquid branch and phase II to its vapour branch (ADR-0019). That is
what makes a liquid-liquid split expressible from an equation of state at a
state where a vapour root also exists.

Two Rachford-Rice solvers live here, and the difference between them is the
subject of ADR-0016:

- :func:`_rachford_rice` searches ``[0, 1]`` only and reports "no root" when
  the equation does not bracket there. Every path used it before ADR-0016 and
  the legacy ``wilson-heuristic`` path still does, unchanged.
- :func:`_extended_rachford_rice` first *calls* :func:`_rachford_rice` and
  returns its answer bit-for-bit when there is one, and only otherwise solves
  on the wider Leibovici-Neoschil window, where the root may be negative or
  above one (the "negative flash" of Whitson & Michelsen 1989). That is what
  lets successive substitution pass through iterates whose ``K`` has no
  in-window root instead of failing there.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np

from ..core import Mixture
from ..exceptions import ConvergenceError, ModelError
from ..models import ActivityModel, EquationOfState
from ._common import as_float_array, normalize_composition
from .settings import FlashSettings


class _SplitSolution:
    """Converged two-phase split plus the quantities needed to verify it.

    ``ln_f_x`` / ``ln_f_y`` are the tangent-plane fugacity terms of the two
    phases at the returned compositions: ``ln phi`` for phi-phi and gamma-phi,
    ``ln gamma`` for gamma-gamma.
    """

    __slots__ = (
        "K",
        "converged",
        "iterations",
        "ln_f_x",
        "ln_f_y",
        "max_delta",
        "negative_flash_steps",
        "vapor_fraction",
        "x",
        "y",
    )

    def __init__(
        self,
        *,
        x: np.ndarray,
        y: np.ndarray,
        vapor_fraction: float,
        K: np.ndarray,
        ln_f_x: np.ndarray,
        ln_f_y: np.ndarray,
        iterations: int,
        max_delta: float,
        converged: bool = True,
        negative_flash_steps: int = 0,
    ) -> None:
        self.x = x
        self.y = y
        self.vapor_fraction = vapor_fraction
        self.K = K
        self.ln_f_x = ln_f_x
        self.ln_f_y = ln_f_y
        self.iterations = iterations
        self.max_delta = max_delta
        self.converged = converged
        #: Iterations whose Rachford-Rice root fell outside ``[0, 1]``. Zero
        #: unless ``_solve_k_loop`` was given ``extended_rachford_rice=True``.
        self.negative_flash_steps = negative_flash_steps


def _solve_k_loop(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState | None,
    activity_model: ActivityModel | None,
    mode: str,
    settings: FlashSettings,
    z: np.ndarray,
    K: np.ndarray,
    vapor_fraction: float,
    max_iter: int | None = None,
    allow_unconverged: bool = False,
    extended_rachford_rice: bool = False,
    terms_x: Callable[[np.ndarray], np.ndarray] | None = None,
    terms_y: Callable[[np.ndarray], np.ndarray] | None = None,
    roots_x: _PhaseRoot | None = None,
    roots_y: _PhaseRoot | None = None,
) -> _SplitSolution:
    """Successive substitution on K with Rachford-Rice updates of ``beta``.

    This is the first-stage phase split, shared by every mode; only its
    initial ``K``, its ``beta`` and the model call that updates ``K`` differ:

    - phi-phi: ``K = phi^x / phi^y``, each phase's ``phi`` evaluated on that
      phase's **own** density/compressibility root when the caller supplies the
      two :class:`_PhaseRoot` holders (``roots_x`` / ``roots_y``, the
      tangent-plane path, ADR-0019), and on the fixed ``"liquid"`` /
      ``"vapor"`` branches when it does not (the legacy ``wilson-heuristic``
      path and gamma-phi, both unchanged);
    - gamma-phi: ``K = gamma^L phi^L / phi^V``;
    - gamma-gamma: ``K = gamma^I / gamma^II`` (``x`` is phase I, ``y`` phase II
      and ``vapor_fraction`` is the mole fraction of phase II);
    - modified-raoult: each phase carries the tangent-plane term of the phase
      *candidate* the stability test assigned to it (``terms_x`` for ``x``,
      ``terms_y`` for ``y``), and the equal-fugacity condition
      ``ln x_i + t_i^x(x) = ln y_i + t_i^y(y)`` gives
      ``K = exp(t^x(x) - t^y(y))``. That is ``gamma_i Psat_i / P`` when ``x``
      is a liquid and ``y`` an ideal vapor, and ``gamma_i^I / gamma_i^II`` when
      both are liquids - one update rule, two regimes.

    ``allow_unconverged`` returns the last iterate instead of raising when the
    budget runs out, which is how the liquid-liquid, modified-Raoult and (since
    ADR-0016) phi-phi paths hand over to their second-order stage.

    ``extended_rachford_rice`` lets the ``beta`` update leave ``[0, 1]``
    (:func:`_extended_rachford_rice`) instead of failing there. It is
    bit-identical whenever the old solver found a root, so it changes only
    those iterates on which the old loop raised; see ADR-0016.
    """
    budget = settings.max_iter if max_iter is None else max_iter
    max_delta = float("inf")
    negative_flash_steps = 0
    iteration = 0
    x = np.array(z, dtype=float)
    y = np.array(z, dtype=float)
    ln_f_x = np.zeros_like(z)
    ln_f_y = np.zeros_like(z)

    for iteration in range(1, budget + 1):
        x = z / (1.0 + vapor_fraction * (K - 1.0))
        x = normalize_composition(x, label="liquid", error_cls=ConvergenceError)

        y = K * x
        y = normalize_composition(y, label="vapor", error_cls=ConvergenceError)

        if mode == "modified-raoult":
            assert terms_x is not None and terms_y is not None
            ln_f_x = terms_x(x)
            ln_f_y = terms_y(y)
            if ln_f_x.shape != K.shape or ln_f_y.shape != K.shape:
                raise ModelError("Phase candidate returned inconsistent term shapes.")
            K_new = np.exp(ln_f_x - ln_f_y)
            if np.any(~np.isfinite(K_new)) or np.any(K_new <= 0.0):
                raise ModelError("Non-finite or non-positive K-values from the phase candidates.")
        elif mode == "gamma-gamma":
            assert activity_model is not None
            gamma_x = _activity_coefficients(activity_model, mixture, temperature, x)
            gamma_y = _activity_coefficients(activity_model, mixture, temperature, y)
            if gamma_x.shape != K.shape or gamma_y.shape != K.shape:
                raise ModelError("Activity model returned inconsistent coefficient shapes.")
            K_new = gamma_x / gamma_y
            ln_f_x = np.log(gamma_x)
            ln_f_y = np.log(gamma_y)
        else:
            assert eos is not None
            if roots_x is not None and roots_y is not None:
                # ADR-0019: each phase on the root the stability test put it
                # on. The two calls keep the historical order (``y`` first) and
                # the arithmetic below is untouched, so a state whose phases
                # are pinned to the historical ``("liquid", "vapor")`` pair
                # reproduces the pre-slice doubles exactly.
                phi_v = roots_y.fugacity_coefficients(y)
                phi_l = roots_x.fugacity_coefficients(x)
            else:
                phi_v = as_float_array(
                    eos.fugacity_coefficients(
                        mixture=mixture,
                        temperature_K=temperature,
                        pressure_Pa=pressure,
                        composition=y.tolist(),
                        phase="vapor",
                    )
                )
                phi_l = as_float_array(
                    eos.fugacity_coefficients(
                        mixture=mixture,
                        temperature_K=temperature,
                        pressure_Pa=pressure,
                        composition=x.tolist(),
                        phase="liquid",
                    )
                )

            if phi_v.shape != phi_l.shape or phi_v.shape != K.shape:
                raise ModelError("EOS returned inconsistent fugacity coefficient shapes.")
            if np.any(phi_v <= 0.0) or np.any(phi_l <= 0.0):
                raise ModelError("EOS returned non-positive fugacity coefficients.")

            if mode == "gamma-phi":
                assert activity_model is not None
                gamma_l = as_float_array(
                    activity_model.activity_coefficients(
                        mixture=mixture,
                        temperature_K=temperature,
                        composition=x.tolist(),
                    )
                )
                if gamma_l.shape != K.shape:
                    raise ModelError("Activity model returned inconsistent coefficient shapes.")
                if np.any(gamma_l <= 0.0):
                    raise ModelError("Activity model returned non-positive activity coefficients.")
                K_new = gamma_l * phi_l / phi_v
            else:
                K_new = phi_l / phi_v
            ln_f_x = np.log(phi_l)
            ln_f_y = np.log(phi_v)

        max_delta = float(np.max(np.abs(K_new - K)))
        if max_delta < settings.tol:
            return _SplitSolution(
                x=x,
                y=y,
                vapor_fraction=vapor_fraction,
                K=K_new,
                ln_f_x=ln_f_x,
                ln_f_y=ln_f_y,
                iterations=iteration,
                max_delta=max_delta,
                negative_flash_steps=negative_flash_steps,
            )

        if settings.damping is None:
            K = K_new
        else:
            K = K + settings.damping * (K_new - K)

        if np.any(K <= 0.0):
            raise ModelError("Non-positive K-values encountered during iteration.")

        if extended_rachford_rice:
            next_vapor_fraction, rr_status = _extended_rachford_rice(z, K)
            if next_vapor_fraction is None:
                # The updated K has no admissible vapor fraction at all - every
                # K on one side of 1, so the substitution map has walked out of
                # the two-phase region. There is no next iterate to form, so
                # successive substitution stops here; a caller with a
                # second-order stage takes the last admissible iterate from
                # there (ADR-0016), and one without it fails as before.
                if not allow_unconverged:
                    raise ConvergenceError(
                        "Rachford-Rice failed to bracket a vapor fraction even on the "
                        f"extended (negative-flash) window: {rr_status}."
                    )
                break
            if rr_status == "negative-flash":
                negative_flash_steps += 1
        else:
            next_vapor_fraction, _f0, _f1 = _rachford_rice(z, K)
            if next_vapor_fraction is None:
                raise ConvergenceError("Rachford-Rice failed to bracket a vapor fraction.")
        vapor_fraction = next_vapor_fraction

    if allow_unconverged:
        return _SplitSolution(
            x=x,
            y=y,
            vapor_fraction=vapor_fraction,
            K=K,
            ln_f_x=ln_f_x,
            ln_f_y=ln_f_y,
            iterations=iteration,
            max_delta=max_delta,
            converged=False,
            negative_flash_steps=negative_flash_steps,
        )

    raise ConvergenceError(
        f"flash_tp did not converge within the iteration limit; max_delta_k={max_delta:.3e}."
    )


def _activity_coefficients(
    activity_model: ActivityModel,
    mixture: Mixture,
    temperature: float,
    composition: np.ndarray,
) -> np.ndarray:
    """Activity coefficients at ``composition``, validated."""
    values = as_float_array(
        activity_model.activity_coefficients(
            mixture=mixture,
            temperature_K=temperature,
            composition=composition.tolist(),
        )
    )
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ModelError("Activity model returned non-positive activity coefficients.")
    return values


def _ln_gamma_function(
    activity_model: ActivityModel, mixture: Mixture, temperature: float
) -> Callable[[np.ndarray], np.ndarray]:
    """Return ``ln gamma(x)`` as a plain callable on normalized compositions."""

    def ln_gamma(composition: np.ndarray) -> np.ndarray:
        values = np.asarray(composition, dtype=float)
        total = float(np.sum(values))
        if total <= 0.0:
            raise ModelError("Activity model called with a non-positive composition.")
        return np.log(_activity_coefficients(activity_model, mixture, temperature, values / total))

    return ln_gamma


#: The density/compressibility branches an ``EquationOfState`` exposes, in the
#: order :func:`chemthermo.stability._evaluator._cubic_root_candidates` builds
#: them, which is the order the lowest-Gibbs fallback below tries them in.
_ROOT_BRANCHES = ("vapor", "liquid")


class _PhaseRoot:
    """The density/compressibility root **one** phase of a phi-phi split sits on.

    ADR-0019. The pre-slice split hard-wired phase I to the model's
    ``"liquid"`` branch and phase II to its ``"vapor"`` branch, which makes a
    liquid-liquid split inexpressible at any state where a vapour root also
    exists - the water / n-hexane defect of validation Case P-7(iii). What
    replaces it is *per-phase*: each phase is pinned to the branch the
    tangent-plane stability test found **that** phase on (``feed_branch`` for
    the feed-like phase, ``phase_branch`` for the incipient one), so the pair
    may legitimately be ``("liquid", "liquid")``.

    The branch is then **held** for the whole split, which is ADR-0012's rule
    for the stability module's own trials ("a trial pinned to one candidate
    iterates on that candidate's Gibbs surface") applied to the split, and for
    the same reason: re-selecting the lowest-Gibbs root at every *iterate*
    lets successive substitution walk a phase onto its partner's branch and
    collapse to the trivial solution. That is not hypothetical - it is measured
    on the ADR-0016 reference state (PC-SAFT carbon dioxide / n-decane,
    ``z = (0.9, 0.1)``, 240 K, 1.0 MPa), where per-iterate re-selection moves
    the liquid phase onto the vapour root at the fifteenth iterate and the
    split then runs away to ``beta = -3.2e+09``. See ADR-0019 "Alternatives
    considered".

    The lowest-Gibbs rule (Michelsen & Mollerup: each phase on the root of
    lowest Gibbs energy) is the **fallback**, used only where the pinned branch
    is not evaluable at a composition, and it is what the post-split stability
    test independently re-applies to every converged phase - so a split that
    did converge onto a higher-Gibbs root is refused there rather than
    returned.

    Attributes:
        branch: The pinned branch, or None when the caller had none to pin.
        selected: The branch actually used at the last call.
        fallbacks: How often the pinned branch was unusable and the
            lowest-Gibbs rule had to decide instead.
    """

    __slots__ = (
        "_eos",
        "_mixture",
        "_pressure",
        "_temperature",
        "branch",
        "fallbacks",
        "selected",
    )

    def __init__(
        self,
        eos: EquationOfState,
        mixture: Mixture,
        temperature: float,
        pressure: float,
        *,
        branch: str | None = None,
    ) -> None:
        self._eos = eos
        self._mixture = mixture
        self._temperature = temperature
        self._pressure = pressure
        self.branch = branch if branch in _ROOT_BRANCHES else None
        self.selected = self.branch
        self.fallbacks = 0

    def fugacity_coefficients(self, composition: np.ndarray) -> np.ndarray:
        """``phi(w)`` on this phase's root, as the model returned it.

        ``composition`` must already be normalized - the split loop hands over
        exactly the array it built with
        :func:`chemthermo.flash._common.normalize_composition` and it reaches
        the model unchanged, because re-normalizing a vector whose sum is one
        only to the last bit would perturb the model's argument and with it
        every double downstream. :meth:`ln_fugacity_terms` is the entry point
        that does normalize.

        The raw coefficients are returned rather than their logarithm so that
        the caller's ``K = phi^x / phi^y`` update stays the floating-point
        expression it was before ADR-0019.
        """
        w = np.asarray(composition, dtype=float)
        failures: list[str] = []

        if self.branch is not None:
            phi = self._branch_coefficients(self.branch, w, failures)
            if phi is not None:
                self.selected = self.branch
                return phi
            self.fallbacks += 1

        best: np.ndarray | None = None
        best_label = ""
        best_g = math.inf
        for label in _ROOT_BRANCHES:
            phi = self._branch_coefficients(label, w, failures)
            if phi is None:
                continue
            g_res = float(np.sum(w * np.log(phi)))
            if not math.isfinite(g_res):
                failures.append(f"{label}: non-finite reduced residual Gibbs energy")
                continue
            if g_res < best_g:
                best_g = g_res
                best = phi
                best_label = label

        if best is None:
            raise ModelError(
                "No usable fugacity-coefficient root for this phase of the split ("
                + "; ".join(failures)
                + ")."
            )
        self.selected = best_label
        return best

    def ln_fugacity_terms(self, composition: np.ndarray) -> np.ndarray:
        """``ln phi(w)`` on this phase's root, normalizing ``w`` first.

        The tangent-plane fugacity term the second-order stage
        (:func:`chemthermo.flash._second_order._second_order_split`) consumes,
        with the same normalize-then-evaluate arithmetic the pre-ADR-0019
        single-branch callable used.
        """
        values = np.asarray(composition, dtype=float)
        total = float(np.sum(values))
        if total <= 0.0:
            raise ModelError("Equation of state called with a non-positive composition.")
        return np.log(self.fugacity_coefficients(values / total))

    def _branch_coefficients(
        self, branch: str, w: np.ndarray, failures: list[str]
    ) -> np.ndarray | None:
        """``phi`` on one branch, or None when that branch is unusable here."""
        try:
            values = as_float_array(
                self._eos.fugacity_coefficients(
                    mixture=self._mixture,
                    temperature_K=self._temperature,
                    pressure_Pa=self._pressure,
                    composition=w.tolist(),
                    phase=branch,
                )
            )
        except Exception as exc:  # noqa: BLE001 - a root may be absent here
            failures.append(f"{branch}: {exc}")
            return None
        if values.shape != w.shape:
            raise ModelError("EOS returned inconsistent fugacity coefficient shapes.")
        if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
            failures.append(f"{branch}: non-finite or non-positive fugacity coefficients")
            return None
        return values


def _rachford_rice(z: np.ndarray, K: np.ndarray) -> tuple[float | None, float, float]:
    """Solve the Rachford-Rice equation; returns (vapor_fraction, f0, f1)."""

    def f(v: float) -> float:
        denom = 1.0 + v * (K - 1.0)
        if np.any(denom <= 0.0):
            return float("nan")
        return float(np.sum(z * (K - 1.0) / denom))

    f0 = f(0.0)
    f1 = f(1.0)
    if not math.isfinite(f0) or not math.isfinite(f1):
        return None, f0, f1

    if f0 * f1 > 0.0:
        return None, f0, f1

    low, high = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (low + high)
        value = f(mid)
        if not math.isfinite(value):
            return None, f0, f1
        if abs(value) < 1e-12:
            return mid, f0, f1
        if value * f0 > 0.0:
            low = mid
            f0 = value
        else:
            high = mid
    return mid, f0, f1


#: Relative margins tried, smallest first, when stepping inside the open
#: Leibovici-Neoschil window to find a sign change of ``f``. The window's
#: endpoints are poles of ``f`` whenever the component attaining ``K_min`` /
#: ``K_max`` is present in the feed, so the smallest usable margin gives the
#: widest bracket; a larger one is only needed when it is not.
_WINDOW_MARGINS = (1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)
#: Bisection / Newton iterations allowed on the extended window.
_EXTENDED_MAX_ITER = 200


def _rachford_rice_window(K: np.ndarray) -> tuple[float, float] | None:
    """The open Leibovici-Neoschil window of admissible vapor fractions.

    ``x_i = z_i / (1 + beta (K_i - 1))`` and ``y_i = K_i x_i`` are non-negative
    exactly where every denominator ``t_i = 1 + beta (K_i - 1)`` is positive,
    which for ``K_i > 1`` means ``beta > 1 / (1 - K_i)`` (a negative bound) and
    for ``K_i < 1`` means ``beta < 1 / (1 - K_i)`` (a bound above one). The
    binding pair is ``K_max`` and ``K_min``:

        1 / (1 - K_max)  <  beta  <  1 / (1 - K_min)

    (Leibovici & Neoschil 1992; the two-phase case of the ``t_i > 0`` region
    restated in :mod:`chemthermo.flash._multiphase_rr`). The window always
    contains ``[0, 1]`` when ``K_max > 1 > K_min``.

    ``None`` is returned when all ``K`` lie on one side of 1: ``f`` is then of
    one sign everywhere the compositions are admissible and there is no root -
    the single-phase case the callers already handle.

    The bounds are taken over **all** components, not only those present in the
    feed. A component with ``z_i = 0`` contributes nothing to ``f``, so its
    bound can only shrink the window and never introduce a spurious root; the
    payoff is that every composition the loop forms inside the window is
    non-negative componentwise, with no special case for ``z_i = 0``.
    """
    k_min = float(np.min(K))
    k_max = float(np.max(K))
    if not (k_max > 1.0 > k_min):
        return None
    return 1.0 / (1.0 - k_max), 1.0 / (1.0 - k_min)


def _extended_rachford_rice(z: np.ndarray, K: np.ndarray) -> tuple[float | None, str]:
    """Rachford-Rice on the Leibovici-Neoschil window (the "negative flash").

    ``f(beta) = sum_i z_i (K_i - 1) / (1 + beta (K_i - 1))`` is strictly
    decreasing wherever it is defined (``f' = -sum_i z_i (K_i - 1)^2 / t_i^2``),
    so the window of :func:`_rachford_rice_window` holds **at most one** root,
    and exactly one when both endpoints are poles. A root outside ``[0, 1]`` is
    a legitimate iterate, not a failure: it is the negative flash of Whitson &
    Michelsen, *Fluid Phase Equilibria* **53** (1989) 51-71, which is what a
    successive-substitution sequence needs in order to pass through ``K`` sets
    whose implied split is momentarily outside the physical range.

    **Bit-identity.** The in-window solver is called first and its answer is
    returned unchanged whenever it has one, so for every ``(z, K)`` on which
    the pre-ADR-0016 loop made progress this function returns the *same double*
    by construction. The extended window is only ever reached where the old
    code raised.

    Returns:
        ``(beta, status)`` with ``status`` one of ``"bracketed"`` (a root in
        ``[0, 1]``, from :func:`_rachford_rice`), ``"negative-flash"`` (a root
        in the window but outside ``[0, 1]``), ``"single-phase"`` (all ``K`` on
        one side of 1) or ``"no-root"`` (a window with no sign change).
        ``beta`` is ``None`` for the last two.
    """
    in_window, _f0, _f1 = _rachford_rice(z, K)
    if in_window is not None:
        return in_window, "bracketed"

    window = _rachford_rice_window(K)
    if window is None:
        return None, "single-phase"
    lower, upper = window
    span = upper - lower

    slope = K - 1.0

    def f(v: float) -> float:
        denom = 1.0 + v * slope
        if np.any(denom <= 0.0):
            return float("nan")
        return float(np.sum(z * slope / denom))

    def derivative(v: float) -> float:
        denom = 1.0 + v * slope
        return float(-np.sum(z * slope * slope / (denom * denom)))

    low = high = math.nan
    f_low = f_high = math.nan
    for margin in _WINDOW_MARGINS:
        step = margin * span
        candidate_low = lower + step
        candidate_high = upper - step
        if candidate_low >= candidate_high:
            break
        value_low = f(candidate_low)
        value_high = f(candidate_high)
        if math.isfinite(value_low) and math.isfinite(value_high) and value_low > 0.0 > value_high:
            low, high, f_low, f_high = candidate_low, candidate_high, value_low, value_high
            break
    if not math.isfinite(f_low) or not math.isfinite(f_high):
        return None, "no-root"

    # Safeguarded Newton: a Newton step from the current iterate when it stays
    # inside the bracket, a bisection step otherwise. The bracket only ever
    # shrinks, so the iteration cannot leave the window.
    beta = 0.5 * (low + high)
    for _ in range(_EXTENDED_MAX_ITER):
        value = f(beta)
        if not math.isfinite(value):  # pragma: no cover - beta stays inside the bracket
            beta = 0.5 * (low + high)
            continue
        if value > 0.0:
            low, f_low = beta, value
        else:
            high, f_high = beta, value
        # Scale-free stopping test: |f| is a sum of terms that cancel at the
        # root, so its attainable size is set by their magnitude, not by 1.
        denom = 1.0 + beta * slope
        magnitude = float(np.sum(np.abs(z * slope / denom)))
        if abs(value) <= 1e-14 * max(magnitude, 1.0):
            break
        if high - low <= 1e-15 * max(1.0, abs(low), abs(high)):
            break
        gradient = derivative(beta)
        step_taken = beta - value / gradient if gradient < 0.0 else math.nan
        if math.isfinite(step_taken) and low < step_taken < high:
            beta = step_taken
        else:
            beta = 0.5 * (low + high)

    return beta, ("bracketed" if 0.0 <= beta <= 1.0 else "negative-flash")
