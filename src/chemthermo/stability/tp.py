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

Modified Raoult: an activity liquid against an ideal vapor
----------------------------------------------------------
``vapor="ideal"`` (ADR-0010) supplies the pure-liquid reference fugacity that
ADR-0007 recorded as missing, in the one regime where it is honest to write it
down: at low pressure, with ``f_i^0 = Psat_i(T)`` (``phi_i^sat = 1``, Poynting
= 1) and an ideal-gas vapor (``phi_i^V = 1``). The equilibrium condition is then
modified Raoult's law,

    y_i P = x_i gamma_i(x) Psat_i(T)                                      (9)

and both phases can be put on **one** Gibbs surface because both are measured
against the same reference ``ln( f_i / (x_i P) )``:

    liquid:  ln gamma_i(w) + ln( Psat_i(T) / P )
    vapor:   0                                                           (10)

Equations (1)-(7) are then used verbatim with those terms, and at every
composition the candidate of lower Gibbs energy is the one that counts (the
same rule that selects the minimum-Gibbs cubic root; see
:mod:`chemthermo.stability._evaluator`). So a single tangent-plane test detects
a vapor-liquid split, a liquid-liquid split, or neither, and the *label* of the
stationary point says which. Antoine coefficients come from the packaged
databank and are refused outside their stated validity range.

Its limits are exactly the limits of the model: ideal vapor (no ``phi^V``, so
no high pressure), no Poynting correction, no ``phi^sat``, and a temperature
range bounded by the Antoine fits.

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

Trial surfaces (ADR-0012, ADR-0021)
-----------------------------------
An initial estimate may also name the phase *candidate* its trial belongs to,
and the iteration then uses that candidate's terms at every step instead of
re-selecting the lowest-Gibbs one. Re-selecting inside an iteration makes the
successive-substitution map discontinuous where the two candidates exchange
Gibbs energy, and a vapor-like trial started on the far side of that crossing
is dragged onto the liquid candidate and collapses to the trivial solution (or
to the feed's partner liquid) even though a vapor stationary point with
``tpd < 0`` exists. ADR-0012 fixed the surfaces for the modified-Raoult pair
after that failure was measured there; ADR-0021 fixes them for the density
roots of an equation of state after the same failure was measured there
(validation Case P-9 (iv), a water / n-hexane feed above the three-phase
temperature). Every trial of both families is therefore pinned to one surface.

The two families differ in one way, and it is the reason ADR-0012 did not do
both at once. The modified-Raoult candidates are two different models and both
exist at every composition, so a pinned trial always has its surface. A density
root can simply *not exist* at an iterate, and then there is no surface to
walk: both ``phase`` labels name the only root the model has there, so the
trial walks that one. That is not a degradation - where the model has one root,
the pinned surface and the lowest-Gibbs surface are the same surface - but it
is reported rather than left invisible, in
``StabilityTrial.surface_fallback`` / ``surface_fallback_count``, measured
where the solver compares the branches anyway (the trial's stopping point; see
ADR-0021 decision 4 for why not at every iterate, and for how often it fires).

The remaining unpinned trials are an activity-coefficient model on its own
(one candidate, nothing to pin) and the degenerate single-active-component
feed, where there is no composition degree of freedom and pinning would test a
root the feed is not on.

Two things are still reported off the pinned surface:

- the stationarity residual (5) is measured on the trial's own surface, because
  that is the equation the trial is solving, and
- the tangent-plane distance (2) is evaluated with the **lowest-Gibbs**
  candidate at the converged composition, because the distance from the Gibbs
  surface to the tangent plane is by definition the minimum over candidates. A
  trial records both: ``surface`` is what it iterated on, ``phase_branch`` is
  the lowest-Gibbs candidate where it stopped. They normally agree.

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

Where the handover is too early (ADR-0028)
------------------------------------------
"After a few substitutions" is a heuristic, and the robustness map of ADR-0027
found the two states where it costs the analysis its verdict: a 15 wt%
polyethylene / n-pentane feed at 10.5 and 10.8 MPa, just inside the pressure at
which the non-trivial stationary point merges with the trivial one. There every
trial reaches ``second_order_no_progress`` - the Newton line search finds *no*
admissible step, at a residual of 0.685 and 1.746 - from a 50-substitution
iterate that was still travelling, and the whole analysis is ``inconclusive``,
which is the one verdict a flash cannot use.

Neither stage is at fault. Successive substitution on its own never converges
there either (3000 substitutions leave a residual of 53); the same unchanged
Newton stage finishes in seven iterations when it is handed the iterate at
substitution 300 instead of 50. So the analysis is re-run once, from the same
trial compositions, with ``ssi_iterations`` raised to ``max_iter`` - and only
when the first pass ended inconclusive, which makes the retry dormant on every
verdict this module has ever returned. It is reported as
``diagnostics["substitution_budget_retry"]``. See validation Case P-17.

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
from dataclasses import replace
from typing import Literal, Mapping

import numpy as np

from .._eos_memo import scoped as _scoped_eos_memo
from ..core import Mixture
from ..exceptions import CompositionError, ModelError
from ..models import ActivityModel, EquationOfState
from ..validation import validate_pressure, validate_temperature
from ._evaluator import (
    _ActivityTangentPlane,
    _EOSTangentPlane,
    _ln_phi_min_gibbs,  # noqa: F401  (re-exported: tests import root selection from here)
    _ModifiedRaoultTangentPlane,
    _SurfaceTerms,
    _TangentPlaneEvaluator,
)
from .results import StabilityResult, StabilityTrial
from .settings import StabilitySettings

#: The range of ``ln W`` over which ``W = exp(ln W)`` is an ordinary double.
#: Before ADR-0025 these were a **clamp** on the iterate; they are now only the
#: gate that decides which of two algebraically identical normalizations is
#: executed (see :func:`_normalize`), so no iterate is bounded by them.
_LN_W_MIN = -700.0
_LN_W_MAX = 700.0
_JACOBIAN_STEP = 1e-6
_MIN_LINE_SEARCH_SCALE = 1e-12
#: ADR-0035: two trials whose tpd differ by at most this (times max(1, |tpd|))
#: reached the same stationary point. Measured ties: <= 5e-16.
_TIE_TPD = 1e-12
#: ADR-0035: a tied trial replaces the lowest-tpd one only when its
#: stationarity residual is smaller by at least this factor.
_TIE_RESIDUAL_ADVANTAGE = 1e3

#: Admissible values of the ``vapor`` keyword.
VAPOR_CANDIDATES = ("none", "ideal")

__all__ = ["stability_tp"]


@_scoped_eos_memo
def stability_tp(
    mixture: Mixture,
    *,
    temperature_K: float,
    pressure_Pa: float,
    eos: EquationOfState | None = None,
    activity_model: ActivityModel | None = None,
    vapor: Literal["none", "ideal"] = "none",
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
        activity_model: Activity-coefficient model supplying ``gamma``. On its
            own (``vapor="none"``) the test is liquid-liquid. Exactly one of
            ``eos`` and ``activity_model`` must be given.
        vapor: Which vapor phase competes with the activity-model liquid.

            - ``"none"`` (default): no vapor candidate. The test is
              liquid-liquid and the pure-liquid reference cancels.
            - ``"ideal"``: an **ideal-gas** vapor competes with the liquid at
              every trial composition (the modified-Raoult model, ADR-0010).
              The liquid's term becomes
              ``ln gamma_i(w) + ln(Psat_i(T) / P)`` with ``Psat`` from the
              databank Antoine coefficients, the vapor's term is ``0``, and the
              lower-Gibbs candidate is used at each composition, so one
              tangent plane covers vapor-liquid *and* liquid-liquid behavior.
              ``feed_branch`` and ``phase_branch`` then report which candidate
              won (``"liquid"`` / ``"vapor"``).

            Only valid with ``activity_model``; ``vapor="ideal"`` with ``eos``
            is a ``ModelError``.
        settings: Iteration controls; defaults to ``StabilitySettings()``.

    Returns:
        StabilityResult with the verdict, the minimum reduced tangent-plane
        distance, the minimizing trial composition, the implied incipient-phase
        K-values (``w_i / z_i``) and a per-trial record.

    Raises:
        InputRangeError: If temperature or pressure is non-physical, or if
            ``vapor="ideal"`` and the temperature is outside the Antoine
            validity range of a component.
        ModelError: If neither or both of ``eos`` and ``activity_model`` are
            given (the combined *EOS-vapor* gamma-phi case is still not
            supported), if ``vapor`` is not a supported value or is combined
            with an ``eos``, if the composition basis is not molar, or if the
            model returns non-finite or non-positive values at the feed
            composition.
        PropertyNotFoundError: If ``vapor="ideal"`` and a component has no
            Antoine record.
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
            "gamma-phi stability against an equation-of-state vapor is not yet supported; "
            "for a low-pressure ideal vapor use activity_model=... with vapor='ideal'."
        )
    if eos is None and activity_model is None:
        raise ModelError("stability_tp requires a model: pass either 'eos' or 'activity_model'.")
    if vapor not in VAPOR_CANDIDATES:
        raise ModelError(f"vapor must be one of {VAPOR_CANDIDATES}; got {vapor!r}.")
    if vapor != "none" and eos is not None:
        raise ModelError(
            "vapor='ideal' adds an ideal-gas candidate to an activity-coefficient liquid "
            "(the modified-Raoult model) and is only valid with 'activity_model'. An "
            "equation of state already supplies its own vapor branch."
        )
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
    elif vapor == "ideal":
        assert activity_model is not None
        evaluator = _ModifiedRaoultTangentPlane(
            activity_model, mixture=mixture, temperature=temperature, pressure=pressure
        )
    else:
        assert activity_model is not None
        evaluator = _ActivityTangentPlane(activity_model, mixture=mixture, temperature=temperature)

    active = z > 0.0
    n_active = int(np.count_nonzero(active))

    ln_f_feed, feed_branch = evaluator.ln_fugacity_terms(z)
    # ADR-0017: replace the min-Gibbs tie-break label with a compressibility
    # identity where the model supports one. Relabels only; `ln_f_feed` (the
    # numeric terms `d` below is built from) is untouched.
    feed_branch = evaluator.identity_label(z, feed_branch)

    # d_i = ln z_i + ln phi_i(z); inactive components are excluded rather than
    # floored, which keeps W_i = 0 for them at every iteration (equation (8)).
    d = np.full(z.shape, -np.inf, dtype=float)
    d[active] = np.log(z[active]) + ln_f_feed[active]

    def run_all(active_settings: StabilitySettings) -> tuple[StabilityTrial, ...]:
        return tuple(
            _run_trial(
                label=estimate.label,
                w0=estimate.composition,
                surface=estimate.surface,
                z=z,
                d=d,
                active=active,
                evaluator=evaluator,
                settings=active_settings,
            )
            for estimate in evaluator.initial_estimates(z, active)
        )

    result = _summarize(
        z=z,
        temperature=temperature,
        pressure=pressure,
        feed_branch=feed_branch,
        trials=run_all(settings),
        settings=settings,
        n_active=n_active,
        evaluator=evaluator,
    )
    if result.status != "inconclusive":
        return result

    retry = _substitution_budget_retry(settings)
    if retry is None:
        return result
    # ADR-0028. Every trial ran out of iterations without reaching a
    # stationary point, which is the one outcome `flash_tp` cannot use at all:
    # it can neither split the feed nor call it one phase. The whole analysis
    # is therefore re-run once, from the very same trial compositions, with
    # successive substitution given its full `max_iter` budget before the
    # second-order stage takes over instead of `ssi_iterations` of it.
    #
    # Why that and not a change to the Newton stage: the two states this was
    # measured on (validation Case P-17) stall in the Newton stage at a
    # residual of 0.685 and 1.746 with `second_order_no_progress`, i.e. with no
    # admissible step at all, from a substitution iterate that was still
    # travelling. Substitution alone never converges there either - it is the
    # *handover* that is early, not either stage that is broken - and 300
    # substitutions put the iterate where the same unchanged Newton stage
    # finishes in seven. Nothing here is reached unless the first pass already
    # left the caller with an inconclusive verdict, so no verdict, no
    # `tpd_min`, no stationary point and no diagnostic that this function has
    # ever returned can move.
    retried = _summarize(
        z=z,
        temperature=temperature,
        pressure=pressure,
        feed_branch=feed_branch,
        trials=run_all(retry),
        settings=retry,
        n_active=n_active,
        evaluator=evaluator,
        extra_diagnostics={
            "substitution_budget_retry": True,
            "substitution_budget_retry_from": settings.ssi_iterations,
        },
    )
    if retried.status == "inconclusive":
        return result
    return retried


def _substitution_budget_retry(settings: StabilitySettings) -> StabilitySettings | None:
    """The ADR-0028 retry settings, or ``None`` when there is nothing to retry.

    The only field that moves is ``ssi_iterations``, and it moves to
    ``max_iter`` - the budget successive substitution already has when the
    second-order stage is switched off. No tolerance, no step bound and no
    iteration ceiling changes, so the retry can only reach a stationary point
    the first pass would have reached with more substitutions.
    """
    if not settings.second_order or settings.ssi_iterations >= settings.max_iter:
        return None
    return replace(settings, ssi_iterations=settings.max_iter)


def _terms(evaluator: _TangentPlaneEvaluator, w: np.ndarray, surface: str | None) -> _SurfaceTerms:
    """Terms used *inside* a trial iteration: the trial's surface, or min-Gibbs.

    ``surface is None`` reproduces the pre-ADR-0012 call exactly, so evaluators
    that name no surface run bit-identically to before; the terms in hand then
    *are* the lowest-Gibbs ones, which is what ``min_gibbs`` records.
    """
    if surface is None:
        terms, branch = evaluator.ln_fugacity_terms(w)
        return _SurfaceTerms(terms, branch, False, (terms, branch))
    return evaluator.ln_terms_on_surface(w, surface)


def _reported_terms(
    evaluator: _TangentPlaneEvaluator,
    w: np.ndarray,
    surface: str | None,
    evaluated: _SurfaceTerms,
) -> tuple[np.ndarray, str | None, int]:
    """Terms used to *report* a converged trial: always the lowest-Gibbs ones.

    The tangent-plane distance is the distance from the Gibbs surface, which is
    the lower envelope of the candidates, so it is evaluated with the
    minimum-Gibbs candidate at the converged composition even when the trial
    iterated on a pinned surface. For an unpinned trial the terms in hand
    already *are* the lowest-Gibbs ones and are reused, which keeps those
    families bit-identical and spares a model call.

    A pinned trial asks the evaluator (:meth:`ln_report_terms`), because this is
    the one place where it looks at the whole envelope: the third return value
    is 1 when the pinned candidate was not available as a distinct surface at
    ``w`` and 0 otherwise, which is how the ADR-0021 fallback is counted (see
    :class:`chemthermo.stability.results.StabilityTrial`).

    The *label* is then passed through :meth:`_TangentPlaneEvaluator.identity_label`
    (ADR-0017): a no-op for every family except the EOS one, where it may
    replace a min-Gibbs tie-break with a compressibility-measured identity.
    This changes only the reported ``StabilityTrial.phase_branch`` string,
    never ``terms``, so it cannot move a tangent-plane distance or a verdict.
    """
    if surface is None:
        report = evaluated
    else:
        report = evaluator.ln_report_terms(w, surface)
    # `min_gibbs` is populated by every evaluator on both paths; the fallback
    # keeps the contract honest for a future one that cannot supply it.
    if report.min_gibbs is None:
        terms, branch = evaluator.ln_fugacity_terms(w)
    else:
        terms, branch = report.min_gibbs
    return terms, evaluator.identity_label(w, branch), int(surface is not None and report.fell_back)


def _ln_of(total: float) -> float:
    """``ln(total)``, extended to the values ``sum_W`` can take in log space.

    ``0.0`` and ``inf`` are not failures once the sum is allowed to leave the
    exponential's range; they are the two ends of it, and their logarithms are
    the finite numbers the caller actually wants.
    """
    if math.isfinite(total) and total > 0.0:
        return math.log(total)
    if total == 0.0:
        return -math.inf
    return math.inf if total > 0.0 else math.nan


def _logsumexp(values: np.ndarray) -> float:
    """``ln sum_i exp(values_i)`` without ever forming ``exp(values_i)``.

    ``-inf`` entries contribute nothing (an absent component). An all-``-inf``
    input returns ``-inf``; a ``+inf`` or ``nan`` entry is propagated, and both
    are refused by the caller.
    """
    largest = float(np.max(values))
    if not math.isfinite(largest):
        return largest
    with np.errstate(under="ignore"):
        return largest + float(np.log(np.sum(np.exp(values - largest))))


def _normalize(
    ln_w_capital: np.ndarray, active: np.ndarray
) -> tuple[np.ndarray, float, float, bool] | None:
    """``w = W / sum_j W_j`` from ``ln W``, in whichever arithmetic survives.

    Returns ``(w, sum_W, ln_sum_W, log_space)``, or None when the sum is
    degenerate (not positive, or not finite in ``ln``).

    Two algebraically identical routes, and which one runs is the whole of the
    ADR-0025 gate:

    - **every** ``ln W_i`` inside ``[-700, 700]``: form ``W = exp(ln W)``, sum,
      divide. This is the pre-ADR-0025 expression, character for character,
      and it runs wherever the pre-ADR-0025 ``np.clip(..., -700, 700)`` would
      have returned its argument untouched - which is to say wherever the old
      clamp was dormant. Every result that existed before this slice is
      therefore bit-identical by construction.
    - otherwise: ``ln sum_W = logsumexp(ln W)`` and
      ``w_i = exp(ln W_i - ln sum_W)``, which is accurate for any magnitude.
      ``w_i`` may underflow to an exact ``0.0`` (the melt over a
      polymer/solvent feed has ``ln W`` spanning ~1450, so the solvent's
      normalized mole fraction is ``exp(-1449)``) and ``sum_W`` itself may
      overflow to ``inf`` or underflow to ``0.0``; ``ln_sum_W`` is the one that
      is always meaningful, and equation (7) ``tpd = -ln sum_W`` is read from
      it rather than from ``sum_W``.

    The old clamp was not protecting a result: it was keeping ``exp`` in range
    so that the *normalization* stayed finite. Doing the normalization in logs
    removes the need for it, which is why it is gone rather than widened.
    """
    values = ln_w_capital[active]
    if np.any(values < _LN_W_MIN) or np.any(values > _LN_W_MAX):
        ln_total = _logsumexp(values)
        if not math.isfinite(ln_total):
            return None
        w = np.zeros(ln_w_capital.shape, dtype=float)
        with np.errstate(over="ignore", under="ignore"):
            w[active] = np.exp(values - ln_total)
            total = float(np.exp(ln_total))
        return w, total, ln_total, True

    w_capital = np.where(active, np.exp(np.where(active, ln_w_capital, 0.0)), 0.0)
    total = float(np.sum(w_capital))
    if not math.isfinite(total) or total <= 0.0:
        return None
    return w_capital / total, total, math.log(total), False


def _run_trial(
    *,
    label: str,
    w0: np.ndarray,
    surface: str | None,
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
    ln_sum_w_capital = 0.0
    log_space = False
    iterations = 0
    fallbacks = 0

    ssi_budget = (
        min(settings.ssi_iterations, settings.max_iter)
        if settings.second_order
        else settings.max_iter
    )

    for iteration in range(1, ssi_budget + 1):
        iterations = iteration
        try:
            evaluated = _terms(evaluator, w, surface)
        except ModelError as exc:
            return _failed_trial(label, iterations, f"model_error: {exc}", surface)
        ln_f, branch = evaluated.terms, evaluated.label
        fallbacks += int(evaluated.fell_back)

        residual = float(np.max(np.abs(ln_w_capital[active] + ln_f[active] - d[active])))
        if residual < settings.tol:
            try:
                report_f, branch, reported_fallback = _reported_terms(
                    evaluator, w, surface, evaluated
                )
            except ModelError as exc:
                return _failed_trial(label, iterations, f"model_error: {exc}", surface)
            fallbacks += reported_fallback
            tpd = _tpd(w, report_f, d, active)
            return StabilityTrial(
                label=label,
                converged=True,
                iterations=iterations,
                tpd=tpd,
                sum_W=sum_w_capital,
                ln_W=tuple(ln_w_capital.tolist()),
                ln_sum_W=ln_sum_w_capital,
                log_space=log_space,
                trivial=_is_trivial(ln_w_capital, z, active, settings.trivial_tol),
                residual=residual,
                phase_branch=branch,
                composition=tuple(w.tolist()),
                termination_reason="stationarity_met",
                ssi_iterations=iterations,
                second_order_iterations=0,
                converged_stage="successive-substitution",
                surface=surface,
                surface_fallback=fallbacks > 0,
                surface_fallback_count=fallbacks,
            )

        # Equation (8): ln W_i <- d_i - ln phi_i(w). ADR-0025: the proposal is
        # *not* clamped to [-700, 700] any more - a stationary point may
        # legitimately sit at ln W ~ 1450 - and the normalization that used to
        # need the clamp is done by `_normalize` in whichever arithmetic the
        # magnitudes allow.
        proposed = d[active] - ln_f[active]
        if not np.all(np.isfinite(proposed)):
            # The pre-ADR-0025 clamp also absorbed an infinite proposal, and
            # left a NaN to the check below. Neither has been observed - the
            # evaluators refuse a non-finite term - so this keeps that
            # behaviour rather than changing an unreachable branch.
            proposed = np.clip(proposed, _LN_W_MIN, _LN_W_MAX)
        ln_w_new = np.full(z.shape, -np.inf, dtype=float)
        ln_w_new[active] = proposed
        if np.any(~np.isfinite(ln_w_new[active])):
            return _failed_trial(label, iterations, "non_finite_ln_W", surface)

        normalized = _normalize(ln_w_new, active)
        if normalized is None:
            return _failed_trial(label, iterations, "degenerate_sum_W", surface)
        w, sum_w_capital, ln_sum_w_capital, step_in_log_space = normalized
        log_space = log_space or step_in_log_space

        ln_w_capital = ln_w_new

        if _is_trivial(ln_w_capital, z, active, settings.trivial_tol):
            try:
                evaluated = _terms(evaluator, w, surface)
                report_f, report_branch, reported_fallback = _reported_terms(
                    evaluator, w, surface, evaluated
                )
            except ModelError as exc:
                return _failed_trial(label, iterations, f"model_error: {exc}", surface)
            ln_f, branch = evaluated.terms, evaluated.label
            fallbacks += int(evaluated.fell_back) + reported_fallback
            return StabilityTrial(
                label=label,
                converged=True,
                iterations=iterations,
                tpd=_tpd(w, report_f, d, active),
                sum_W=sum_w_capital,
                ln_W=tuple(ln_w_capital.tolist()),
                ln_sum_W=ln_sum_w_capital,
                log_space=log_space,
                trivial=True,
                residual=float(np.max(np.abs(ln_w_capital[active] + ln_f[active] - d[active]))),
                phase_branch=report_branch,
                composition=tuple(w.tolist()),
                termination_reason="trivial_solution",
                ssi_iterations=iterations,
                second_order_iterations=0,
                converged_stage="successive-substitution",
                surface=surface,
                surface_fallback=fallbacks > 0,
                surface_fallback_count=fallbacks,
            )

    if not settings.second_order:
        return StabilityTrial(
            label=label,
            converged=False,
            iterations=iterations,
            tpd=float("nan"),
            sum_W=sum_w_capital,
            ln_W=tuple(ln_w_capital.tolist()),
            ln_sum_W=ln_sum_w_capital,
            log_space=log_space,
            trivial=False,
            residual=residual,
            phase_branch=branch,
            composition=tuple(w.tolist()),
            termination_reason="max_iter",
            ssi_iterations=iterations,
            second_order_iterations=0,
            converged_stage=None,
            surface=surface,
            surface_fallback=fallbacks > 0,
            surface_fallback_count=fallbacks,
        )

    return _second_order_stage(
        label=label,
        surface=surface,
        surface_fallback_count=fallbacks,
        log_space=log_space,
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
    surface: str | None,
    surface_fallback_count: int,
    log_space: bool,
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

    ADR-0025 removed the ``[-700, 700]`` clamp from both places it appeared
    here - the normalization inside ``evaluate`` and the line search's
    candidate - because it was not a safeguard on the *step* but on ``exp``.
    ``u`` is the natural variable of this stage and a stationary point at
    ``u ~ 1450`` is an ordinary double; what could not survive was
    ``exp(1450)``, and :func:`_normalize` no longer forms it. The line search
    keeps a finiteness test instead of a box: a candidate carrying ``inf`` or
    ``nan`` is rejected and the scale halved, which is what the clamp was
    doing for such a candidate anyway.
    """
    index = np.flatnonzero(active)
    d_active = d[index]
    fallbacks = surface_fallback_count
    used_log_space = log_space

    def evaluate(u: np.ndarray) -> tuple[np.ndarray, np.ndarray, str | None]:
        nonlocal fallbacks, used_log_space
        ln_w_full = np.full(z.shape, -np.inf, dtype=float)
        ln_w_full[index] = u
        normalized = _normalize(ln_w_full, active)
        if normalized is None:
            raise ModelError("Second-order stage produced a degenerate sum of mole numbers.")
        w_local, _total, _ln_total, step_in_log_space = normalized
        used_log_space = used_log_space or step_in_log_space
        evaluated = _terms(evaluator, w_local, surface)
        fallbacks += int(evaluated.fell_back)
        return u + evaluated.terms[index] - d_active, w_local, evaluated.label

    u = ln_w_capital[index].copy()
    iterations = 0
    reason = "second_order_max_iter"

    try:
        g, w, branch = evaluate(u)
    except ModelError as exc:
        return _failed_trial(label, ssi_iterations, f"model_error: {exc}", surface)
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
            return _failed_trial(label, ssi_iterations + iterations, f"model_error: {exc}", surface)

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
            candidate = u + scale * step
            if not np.all(np.isfinite(candidate)):
                scale *= 0.5
                continue
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
    if np.any(u < _LN_W_MIN) or np.any(u > _LN_W_MAX):
        used_log_space = True
        ln_sum_w_capital = _logsumexp(u)
        with np.errstate(over="ignore", under="ignore"):
            sum_w_capital = float(np.exp(ln_sum_w_capital))
    else:
        # The pre-ADR-0025 expression, untouched: `np.clip` returned its
        # argument here whenever this branch is the one taken.
        sum_w_capital = float(np.sum(np.exp(u)))
        ln_sum_w_capital = _ln_of(sum_w_capital)
    converged = residual < settings.tol
    trivial = converged and _is_trivial(ln_w_capital_final, z, active, settings.trivial_tol)

    try:
        evaluated = _terms(evaluator, w, surface)
        ln_f, branch, reported_fallback = _reported_terms(evaluator, w, surface, evaluated)
    except ModelError as exc:
        return _failed_trial(label, ssi_iterations + iterations, f"model_error: {exc}", surface)
    fallbacks += int(evaluated.fell_back) + reported_fallback

    return StabilityTrial(
        label=label,
        converged=converged,
        iterations=ssi_iterations + iterations,
        tpd=_tpd(w, ln_f, d, active) if converged else float("nan"),
        sum_W=sum_w_capital,
        ln_W=tuple(ln_w_capital_final.tolist()),
        ln_sum_W=ln_sum_w_capital,
        log_space=used_log_space,
        trivial=trivial,
        residual=residual,
        phase_branch=branch,
        composition=tuple(w.tolist()),
        termination_reason=reason,
        ssi_iterations=ssi_iterations,
        second_order_iterations=iterations,
        converged_stage="second-order" if converged else None,
        surface=surface,
        surface_fallback=fallbacks > 0,
        surface_fallback_count=fallbacks,
    )


def _failed_trial(
    label: str, iterations: int, reason: str, surface: str | None = None
) -> StabilityTrial:
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
        surface=surface,
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
    extra_diagnostics: Mapping[str, float | int | str | bool] | None = None,
) -> StabilityResult:
    """Reduce per-trial outcomes to a verdict and the minimizing trial.

    ``extra_diagnostics`` is written only by the ADR-0028 retry, which is
    itself reached only from an inconclusive first pass, so a result assembled
    without it carries the keys it carried before ADR-0028 and no others.
    """
    converged = [trial for trial in trials if trial.converged]
    non_trivial = [trial for trial in converged if not trial.trivial and math.isfinite(trial.tpd)]

    best: StabilityTrial | None = None
    for trial in non_trivial:
        if best is None or trial.tpd < best.tpd:
            best = trial

    # ADR-0035: several trials often stop at the *same* stationary point, and
    # which of them has the lowest tpd is then last-bit noise. Report a tied
    # trial instead only when it is converged decisively better (residual
    # smaller by _TIE_RESIDUAL_ADVANTAGE); a residual ratio that large is a
    # different stopping point of the iteration, not rounding, so the rule
    # cannot itself flip on noise.
    tie_broken_by_residual = False
    if best is not None:
        window = _TIE_TPD * max(1.0, abs(best.tpd))
        tied = [trial for trial in non_trivial if abs(trial.tpd - best.tpd) <= window]
        sharpest = min(tied, key=lambda trial: trial.residual)
        if sharpest.residual * _TIE_RESIDUAL_ADVANTAGE < best.residual:
            best = sharpest
            tie_broken_by_residual = True

    if best is not None and best.tpd < -settings.tpd_tol:
        status = "unstable"
    elif converged:
        status = "stable"
    else:
        status = "inconclusive"

    tpd_min = best.tpd if best is not None else 0.0

    trial_composition: tuple[float, ...] | None = None
    trial_ln_capital_w: tuple[float, ...] | None = None
    k_values: tuple[float, ...] | None = None
    phase_branch: str | None = None
    if best is not None and best.composition is not None:
        trial_composition = best.composition
        trial_ln_capital_w = best.ln_W
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
        **evaluator.diagnostics,
        **(extra_diagnostics or {}),
    }

    surface_counts: dict[str, int] = {}
    for trial in trials:
        if trial.surface is not None:
            surface_counts[trial.surface] = surface_counts.get(trial.surface, 0) + 1
    log_space_trials = [trial for trial in trials if trial.log_space]
    if log_space_trials:
        # Written only when the ADR-0025 route actually engaged, so its absence
        # is the assertion that nothing left the pre-ADR-0025 arithmetic.
        diagnostics["log_space_trial_count"] = len(log_space_trials)
        diagnostics["log_space_trials"] = ",".join(trial.label for trial in log_space_trials)

    if surface_counts:
        # Deterministic: insertion order is the deterministic trial order.
        diagnostics["trial_surfaces"] = ",".join(
            f"{name}:{count}" for name, count in surface_counts.items()
        )
        diagnostics["surface_fallback_trial_count"] = sum(
            1 for trial in trials if trial.surface_fallback
        )
        diagnostics["surface_fallback_evaluation_count"] = sum(
            trial.surface_fallback_count for trial in trials
        )
    if feed_branch is not None:
        diagnostics["feed_branch"] = feed_branch
    if best is not None:
        diagnostics["minimizing_trial"] = best.label
        if tie_broken_by_residual:
            # Conditional, like log_space_trial_count: its absence says the
            # lowest-tpd trial was reported, as before ADR-0035.
            diagnostics["minimizing_trial_tie_break"] = "residual"
        if best.log_space:
            # Conditional for the same reason as `log_space_trial_count`: its
            # absence is what a consumer reads as "the normalized `w` carries
            # this stationary point, as it always did".
            diagnostics["minimizing_trial_log_space"] = True
        if best.surface is not None:
            diagnostics["minimizing_trial_surface"] = best.surface
        diagnostics["minimizing_trial_iterations"] = best.iterations
        diagnostics["minimizing_trial_ssi_iterations"] = best.ssi_iterations
        diagnostics["minimizing_trial_second_order_iterations"] = best.second_order_iterations
        diagnostics["minimizing_trial_residual"] = best.residual
        if best.converged_stage is not None:
            diagnostics["minimizing_trial_stage"] = best.converged_stage
        diagnostics["sum_W"] = best.sum_W
        diagnostics["ln_sum_W"] = best.ln_sum_W
        if math.isfinite(best.ln_sum_W):
            # Equation (7), evaluated on the surface the trial iterated on. It
            # equals `tpd_min` whenever the trial's surface is also the
            # lowest-Gibbs candidate where it stopped - always so for a cubic
            # or a lone activity liquid, and measured to 3.0e-12 over every
            # *unstable* verdict of the ternary grid of validation Case V-5.
            # Where a pinned surface's stationary point lies above the other
            # candidate (only seen on stable verdicts, whose `tpd_min` is a
            # positive number that decides nothing) the two differ, and
            # `tpd_min` is the smaller, correct one: the tangent-plane distance
            # is measured to the lower envelope of the candidates.
            #
            # ADR-0025: both are read from `ln_sum_W`, which is `math.log` of
            # `sum_W` exactly wherever `sum_W` is a positive, finite double -
            # so this is the same number as before wherever the old guard let
            # it be written - and is finite where `sum_W` has overflowed.
            # `tm*` itself cannot be: `1 - exp(1452)` is `-inf`, which is
            # reported as such rather than omitted, because the sign is the
            # verdict and the magnitude is in `tpd_from_sum_W`.
            diagnostics["tpd_from_sum_W"] = -best.ln_sum_W
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
        trial_ln_W=trial_ln_capital_w,
    )
