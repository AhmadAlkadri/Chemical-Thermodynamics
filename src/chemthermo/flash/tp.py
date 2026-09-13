"""TP flash calculations (T,P) in SI units.

Supports phi-phi (vapor-liquid, equation of state), gamma-phi (activity-model
liquid against an equation-of-state vapor) and gamma-gamma (liquid-liquid, both
phases described by one activity model). Given the same inputs and settings,
results are deterministic.

Phase detection (phi-phi and gamma-gamma)
-----------------------------------------
Both reference paths decide one phase versus two from Michelsen's tangent-plane
stability criterion rather than from Wilson K-value bounds (ADR-0008,
ADR-0009). The flow is

    flash_tp -> stability_tp(feed) -> single phase | seeded split
             -> post-split stability of every converged phase

1. ``stability_tp`` is run on the feed at the same ``(T, P)`` with the same
   model (``eos=`` for phi-phi, ``activity_model=`` for gamma-gamma).
2. ``status == "stable"``: a single-phase ``FlashResult`` is returned with
   ``termination_reason = "feed_stable_tangent_plane"``. For phi-phi its phase
   name is the minimum-Gibbs root branch the stability test selected
   (``feed_branch``); when the cubic has a single real root both branches
   coincide and the label is a *convention*, not a measurement. For
   gamma-gamma the single phase is a liquid and is named ``"liquid"``.
3. ``status == "unstable"``: the converged stationary point seeds the K-values
   (see :func:`_stability_k_seed`) and the successive-substitution /
   Rachford-Rice loop - Michelsen's recommended first-stage phase split - runs
   from there. For phi-phi, if the seeded K-values give no Rachford-Rice root
   the Wilson estimate is tried as a documented fallback, and
   ``diagnostics["k_seed"]`` records which seed was actually used.
4. ``status == "inconclusive"``: a :class:`chemthermo.ConvergenceError` is
   raised. A stability search that could not converge must not silently produce
   a single-phase answer.

The converged split is then verified (material balance, phase fractions, equal
fugacities or equal activities, and a negative Gibbs-energy change against the
single-phase feed) and every residual is reported in ``diagnostics``.

Post-split stability
--------------------
Every converged two-phase result on these two paths is re-tested: each phase is
fed back into ``stability_tp`` with the same model. Two coexisting phases share
one tangent plane, so each of them is *marginally* stable with respect to the
other - a trial that converges onto the partner phase has ``tpd = 0`` up to the
split's own convergence tolerance and is classified ``"marginal"`` here, not as
an instability (validation Case S-3). A phase whose tangent-plane minimum is
genuinely negative somewhere else means the two-phase answer is not a stable
phase set: a third phase is needed, which this release cannot produce, so
``flash_tp`` raises :class:`chemthermo.ConvergenceError` instead of returning
it. ``FlashSettings(post_split_stability=False)`` returns the result anyway
with the failure recorded in ``diagnostics``.

The liquid-liquid split
-----------------------
With two liquid phases described by one activity model at the same pure-liquid
reference state, ``mu_i^0`` cancels and the equilibrium condition is equality of
activities,

    x_i^I gamma_i^I = x_i^II gamma_i^II                                   (1)

so the K-values of the shared split loop are ``K_i = x_i^II / x_i^I =
gamma_i^I / gamma_i^II``. Phase I is the feed-like phase and phase II the
incipient-like one; ``beta`` is the mole fraction of phase II, which is what
Rachford-Rice returns for that K convention. The two phases are named
``"liquid1"`` (feed-like) and ``"liquid2"`` (incipient-like): those labels are
*roles assigned by the seed*, carry no physical identity, and may swap between
two feeds on the same tie-line (see :func:`flash_tp`).

Successive substitution on (1) converges linearly with a ratio close to one
near a plait point - hundreds to thousands of iterations on the Tessier et al.
(2000) Problem 1 feeds - so a second-order stage is required; see
:func:`_second_order_split` for the derivation.

Limits of this slice
--------------------
At most two phases are returned. A feed that needs three is now *reported*
(``ConvergenceError`` from the post-split check) rather than silently returned
as two, but it is still not solved; multiphase flash is the next slice. A
negative ``tpd_min`` proves a feed is not one phase; ``"stable"`` only means no
negative tangent-plane distance was found from the deterministic trial set.
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

import numpy as np

from ..core import Composition, Mixture
from ..exceptions import CompositionError, ConvergenceError, ModelError
from ..models import ActivityModel, EquationOfState
from ..validation import COMPOSITION_SUM_TOL, validate_pressure, validate_temperature
from ._common import as_float_array, normalize_composition, wilson_k
from .results import FlashResult, PhaseResult
from .settings import FlashSettings

#: Seed K-value used for components absent from the feed (``z_i == 0``). Those
#: components have ``x_i = y_i = 0`` at every iteration and are rewritten from
#: the model on the first update, so the seed value cannot affect the result.
_INERT_SEED_K = 1.0

#: Supported ``flash_mode`` values.
FLASH_MODES = ("phi-phi", "gamma-phi", "gamma-gamma")

#: Phase names of a liquid-liquid result. Roles, not identities; see
#: :func:`flash_tp`.
_LIQUID_I = "liquid1"
_LIQUID_II = "liquid2"

#: Central-difference step for the second-order stage's Hessian.
_HESSIAN_STEP = 1e-7
#: Smallest accepted backtracking scale in the second-order line search.
_MIN_LINE_SEARCH_SCALE = 1e-14
#: Armijo constant for the line search on the two-phase Gibbs energy.
_ARMIJO_C = 1e-4


def flash_tp(
    mixture: Mixture,
    *,
    temperature_K: float,
    pressure_Pa: float,
    eos: EquationOfState | None = None,
    activity_model: ActivityModel | None = None,
    flash_mode: str | None = None,
    settings: FlashSettings | None = None,
) -> FlashResult:
    """Perform a TP flash calculation using phi-phi, gamma-phi or gamma-gamma.

    Args:
        mixture: Mixture with mole-fraction composition.
        temperature_K: Temperature in K.
        pressure_Pa: Pressure in Pa.
        eos: Equation-of-state model used for fugacity coefficients. Required
            for ``"phi-phi"`` and ``"gamma-phi"``, and must be omitted for
            ``"gamma-gamma"``.
        activity_model: Activity model. Required for ``"gamma-phi"`` (liquid
            phase) and for ``"gamma-gamma"`` (both liquid phases).
        flash_mode: Case-insensitive mode: ``"phi-phi"``, ``"gamma-phi"`` or
            ``"gamma-gamma"``. ``None`` (the default) infers the mode from the
            models supplied: ``"gamma-gamma"`` when only ``activity_model`` is
            given, ``"phi-phi"`` otherwise.
        settings: Iteration controls (tolerance, damping, max iterations,
            phase-detection mode, stability settings, post-split check,
            second-order stage).

    Returns:
        FlashResult with phase compositions and fractions. Phase names follow:
        VLE -> ``"liquid"``/``"vapor"``, LLE -> ``"liquid1"``/``"liquid2"``,
        single phase -> ``"liquid"`` or ``"vapor"``. ``vapor_fraction`` is
        ``None`` for every gamma-gamma result: neither phase is a vapor, so
        reporting a number there would be fiction.

    Diagnostics:
        Diagnostics keys are implementation details. Current stable keys:

        - Always: ``flash_mode``, ``phase_detection``, ``iterations``,
          ``converged``, ``termination_reason``, ``phase_count``,
          ``phase_state``, ``phase_regime``.
        - Tangent-plane paths (phi-phi default and gamma-gamma):
          ``stability_status``, ``tpd_min``, ``stability_trials``, and
          ``feed_branch`` when the model reports one. A single phase adds
          nothing else and uses
          ``termination_reason = "feed_stable_tangent_plane"``. A two-phase
          result adds ``k_seed``, ``mass_balance_residual``,
          ``delta_g_split_rt``, the equilibrium residual
          (``fugacity_residual`` for phi-phi, ``equilibrium_residual`` for
          gamma-gamma) and the post-split keys ``post_split_checked``,
          ``post_split_stable``, ``post_split_status``,
          ``post_split_tpd_min``, ``phase_stability_<name>`` and
          ``phase_stability_tpd_min_<name>``. Phi-phi additionally reports
          ``incipient_phase``, ``max_delta_k``, ``k_min`` and ``k_max``;
          gamma-gamma additionally reports ``ssi_iterations``,
          ``second_order_iterations`` and ``converged_stage``.
        - Legacy heuristic path: ``k_min``, ``k_max``, ``max_delta_k``,
          ``k_seed`` (``"wilson"``), and, for its single-phase fallbacks,
          ``rr_f0``, ``rr_f1``, ``rr_status``.

    Raises:
        InputRangeError: If temperature or pressure is non-physical.
        ModelError: If required models are missing, if models are combined in an
            unsupported way, or if a model returns invalid values.
        CompositionError: If the mixture composition is invalid.
        ConvergenceError: If iteration fails to converge; or (tangent-plane
            modes only) if the stability analysis is inconclusive, if an
            unstable feed admits no Rachford-Rice root from either seed, or if
            the converged two-phase result fails the post-split stability check
            (a third phase is required).

    Notes:
        With ``phase_detection="tangent-plane"`` (the default for phi-phi and
        the only option for gamma-gamma) a single-phase result means the feed
        was *found stable* by Michelsen's test. With
        ``phase_detection="wilson-heuristic"`` it only means an initial-estimate
        heuristic said so.

        **Vapor/liquid labelling convention.** For a single-phase phi-phi result
        the name is the minimum-Gibbs compressibility root branch of the feed.
        When the cubic has a single real root (dense or supercritical fluids)
        both branches return identical fugacity coefficients, the branch label
        is a tie-break, and the reported ``"vapor"`` / ``"liquid"`` name is
        therefore a naming convention rather than a phase identification. For a
        two-phase phi-phi result the two converged phases are named by
        volatility ordering: the phase enriched (relative to the feed) in the
        component with the largest Wilson K relative to the one with the
        smallest is named ``"vapor"``. ``EquationOfState`` exposes no molar
        volume, so no density-based identification is available; this ordering
        decides the *name* only, never the verdict or the compositions.

        **Liquid-liquid phase names are roles, not identities.** ``"liquid1"``
        is the phase the split was started from as feed-like and ``"liquid2"``
        the one started from the tangent-plane minimizer. Nothing distinguishes
        two liquids the way volatility distinguishes a vapor from a liquid, so
        no attempt is made to name them by composition: two feeds on the same
        tie-line can come back with the same pair of compositions under swapped
        labels. Compare the phase *set*, not ``result.phases["liquid1"]``.
    """

    temperature = validate_temperature(temperature_K)
    pressure = validate_pressure(pressure_Pa)

    mode = _resolve_flash_mode(flash_mode, eos=eos, activity_model=activity_model)
    settings = settings or FlashSettings()

    if mode == "vlle":
        raise ModelError(
            "VLLE support is provided by the optional chemthermo_vlle plugin. "
            "Install chemthermo_vlle to enable VLLE support."
        )
    if mode not in FLASH_MODES:
        raise ModelError(f"Unsupported flash_mode '{flash_mode}'.")

    if mode == "gamma-gamma":
        if activity_model is None:
            raise ModelError("An activity model is required for gamma-gamma (liquid-liquid) flash.")
        if eos is not None:
            raise ModelError(
                "gamma-gamma flash describes both phases with the activity model, so 'eos' "
                "must not be given. Combined gamma-phi equilibrium is flash_mode='gamma-phi'."
            )
    else:
        if eos is None:
            raise ModelError("An equation-of-state model is required for flash_tp.")
        if mode == "gamma-phi" and activity_model is None:
            raise ModelError("An activity model is required for gamma-phi flash.")
        if mode == "phi-phi" and activity_model is not None:
            raise ModelError(
                "activity_model is only used when flash_mode='gamma-phi' or 'gamma-gamma'."
            )

    if mixture.basis != "mole":
        raise ModelError("flash_tp currently requires mole-fraction compositions.")

    z = np.array(mixture.fractions, dtype=float)
    if z.size == 0:
        raise CompositionError("Mixture composition must be non-empty.")

    if mode == "gamma-gamma":
        assert activity_model is not None
        return _flash_tp_liquid_liquid(
            mixture,
            temperature,
            pressure,
            activity_model=activity_model,
            settings=settings,
            z=z,
        )

    assert eos is not None
    if mode == "phi-phi" and settings.phase_detection == "tangent-plane":
        return _flash_tp_tangent_plane(
            mixture, temperature, pressure, eos=eos, mode=mode, settings=settings, z=z
        )

    return _flash_tp_wilson_heuristic(
        mixture,
        temperature,
        pressure,
        eos=eos,
        activity_model=activity_model,
        mode=mode,
        settings=settings,
        z=z,
    )


def _resolve_flash_mode(
    flash_mode: str | None,
    *,
    eos: EquationOfState | None,
    activity_model: ActivityModel | None,
) -> str:
    """Resolve the requested mode, inferring it when it was not given.

    ``flash_mode=None`` means "use the models I passed": an activity model on
    its own is a liquid-liquid problem, anything else keeps the historical
    ``"phi-phi"`` default. Naming a mode explicitly always wins, so
    ``flash_mode="phi-phi"`` without an ``eos`` is still an error rather than a
    silent reinterpretation.
    """
    if flash_mode is None:
        if eos is None and activity_model is not None:
            return "gamma-gamma"
        return "phi-phi"
    return flash_mode.strip().casefold()


def _flash_tp_tangent_plane(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState,
    mode: str,
    settings: FlashSettings,
    z: np.ndarray,
) -> FlashResult:
    """Phi-phi TP flash whose 1-vs-2 phase decision is a stability test."""
    # Imported here, not at module scope: chemthermo.stability._evaluator imports
    # chemthermo.flash._common, so a module-level import would make the two
    # packages' import order significant.
    from ..stability import stability_tp
    from ..stability.tp import _ln_phi_min_gibbs

    stability = stability_tp(
        mixture,
        temperature_K=temperature,
        pressure_Pa=pressure,
        eos=eos,
        settings=settings.stability_settings,
    )

    base: dict[str, float | int | str | bool] = {
        "flash_mode": mode,
        "phase_detection": "tangent-plane",
        "stability_status": stability.status,
        "tpd_min": float(stability.tpd_min),
        "stability_trials": len(stability.trials),
    }
    if stability.feed_branch is not None:
        base["feed_branch"] = stability.feed_branch

    if stability.status == "inconclusive":
        raise ConvergenceError(
            "Tangent-plane stability analysis was inconclusive (no trial converged), so "
            "flash_tp cannot decide whether the feed is one phase or two. Loosen "
            "FlashSettings.stability_settings, or pass "
            "FlashSettings(phase_detection='wilson-heuristic') to use the legacy "
            "K-bound heuristic instead."
        )

    if stability.status == "stable":
        phase_name = stability.feed_branch or "vapor"
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name=phase_name,
            vapor_fraction=1.0 if phase_name == "vapor" else 0.0,
            diagnostics={
                **base,
                "iterations": 0,
                "converged": True,
                "termination_reason": "feed_stable_tangent_plane",
                "phase_count": 1,
                "phase_state": phase_name,
                "phase_regime": "single-phase",
            },
        )

    trial = stability.trial_composition
    if trial is None:
        raise ConvergenceError(
            "Tangent-plane stability reported an unstable feed without a minimizing "
            "trial composition; no phase split can be seeded."
        )

    k_seed, incipient_phase = _stability_k_seed(
        mixture,
        temperature,
        pressure,
        z=z,
        w=np.array(trial, dtype=float),
        tpd_min=float(stability.tpd_min),
    )

    seed_label = "stability"
    vapor_fraction, _f0, _f1 = _rachford_rice(z, k_seed)
    if vapor_fraction is None:
        # Documented fallback: the stationary point is a valid starting phase but
        # its K-values need not bracket a Rachford-Rice root in every geometry.
        seed_label = "wilson"
        k_seed = wilson_k(mixture, temperature, pressure)
        vapor_fraction, _f0, _f1 = _rachford_rice(z, k_seed)
        if vapor_fraction is None:
            raise ConvergenceError(
                "Feed is unstable (tpd_min="
                f"{stability.tpd_min:.6e}) but neither the stability-seeded nor the "
                "Wilson K-values bracket a Rachford-Rice root."
            )

    split = _solve_k_loop(
        mixture,
        temperature,
        pressure,
        eos=eos,
        activity_model=None,
        mode=mode,
        settings=settings,
        z=z,
        K=k_seed,
        vapor_fraction=vapor_fraction,
    )

    ln_phi_feed, _feed_branch = _ln_phi_min_gibbs(
        eos, mixture=mixture, temperature=temperature, pressure=pressure, composition=z
    )
    checks = _verify_split(
        z=z,
        x=split.x,
        y=split.y,
        beta=split.vapor_fraction,
        ln_f_x=split.ln_f_x,
        ln_f_y=split.ln_f_y,
        ln_f_feed=ln_phi_feed,
        residual_key="fugacity_residual",
    )

    post_split = _post_split_stability(
        mixture,
        temperature,
        pressure,
        eos=eos,
        activity_model=None,
        phases=((("liquid"), split.x), ("vapor", split.y)),
        settings=settings,
    )

    return _two_phase_result(
        mixture,
        temperature,
        pressure,
        split.x,
        split.y,
        split.vapor_fraction,
        names=("liquid", "vapor"),
        vapor_fraction=split.vapor_fraction,
        diagnostics={
            **base,
            "iterations": split.iterations,
            "converged": True,
            "termination_reason": "tolerance_met",
            "max_delta_k": split.max_delta,
            "k_min": float(np.min(split.K)),
            "k_max": float(np.max(split.K)),
            "phase_count": 2,
            "phase_state": "two_phase",
            "phase_regime": "VLE",
            "k_seed": seed_label,
            "incipient_phase": incipient_phase,
            **checks,
            **post_split,
        },
    )


def _flash_tp_liquid_liquid(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    activity_model: ActivityModel,
    settings: FlashSettings,
    z: np.ndarray,
) -> FlashResult:
    """Liquid-liquid TP flash driven entirely by an activity model (ADR-0009).

    The phase count is an *output*: the feed is tested with
    ``stability_tp(..., activity_model=...)`` and a split is attempted only when
    that test finds a negative tangent-plane distance. The split itself is the
    shared Rachford-Rice / successive-substitution loop with
    ``K_i = gamma_i^I / gamma_i^II``, followed by the second-order stage of
    :func:`_second_order_split` when the equal-activity residual is still above
    ``settings.second_order_tol``.
    """
    from ..stability import stability_tp

    stability = stability_tp(
        mixture,
        temperature_K=temperature,
        pressure_Pa=pressure,
        activity_model=activity_model,
        settings=settings.stability_settings,
    )

    base: dict[str, float | int | str | bool] = {
        "flash_mode": "gamma-gamma",
        "phase_detection": "tangent-plane",
        "stability_status": stability.status,
        "tpd_min": float(stability.tpd_min),
        "stability_trials": len(stability.trials),
    }

    if stability.status == "inconclusive":
        raise ConvergenceError(
            "Tangent-plane stability analysis was inconclusive (no trial converged), so "
            "flash_tp cannot decide whether the feed is one liquid phase or two. Loosen "
            "FlashSettings.stability_settings and try again."
        )

    if stability.status == "stable":
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name="liquid",
            vapor_fraction=None,
            diagnostics={
                **base,
                "iterations": 0,
                "converged": True,
                "termination_reason": "feed_stable_tangent_plane",
                "phase_count": 1,
                "phase_state": "liquid",
                "phase_regime": "single-phase",
            },
        )

    trial = stability.trial_composition
    if trial is None:
        raise ConvergenceError(
            "Tangent-plane stability reported an unstable feed without a minimizing "
            "trial composition; no phase split can be seeded."
        )

    # Same seed as the phi-phi path: Michelsen's unnormalized mole numbers
    # W = w exp(-tpd), so that f_RR(0) = sum_i W_i - 1 > 0 is bracketable.
    # Phase I is feed-like, phase II incipient-like, K_i = x_i^II / x_i^I.
    w = np.array(trial, dtype=float)
    sum_capital_w = math.exp(-float(stability.tpd_min))
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = (w * sum_capital_w) / z
    usable = np.isfinite(ratio) & (ratio > 0.0)
    k_seed = np.where(usable, ratio, _INERT_SEED_K)

    beta, _f0, _f1 = _rachford_rice(z, k_seed)
    if beta is None:
        raise ConvergenceError(
            f"Feed is unstable (tpd_min={stability.tpd_min:.6e}) but the stability-seeded "
            "K-values do not bracket a Rachford-Rice root, so no liquid-liquid split can "
            "be started."
        )

    ln_gamma = _ln_gamma_function(activity_model, mixture, temperature)

    split = _solve_k_loop(
        mixture,
        temperature,
        pressure,
        eos=None,
        activity_model=activity_model,
        mode="gamma-gamma",
        settings=settings,
        z=z,
        K=k_seed,
        vapor_fraction=beta,
        max_iter=(
            min(settings.ssi_iterations, settings.max_iter)
            if settings.second_order
            else settings.max_iter
        ),
        allow_unconverged=settings.second_order,
    )

    x_i, x_ii, beta = split.x, split.y, split.vapor_fraction
    ln_f_x, ln_f_y = split.ln_f_x, split.ln_f_y
    residual = _equilibrium_residual(x_i, x_ii, ln_f_x, ln_f_y)
    second_order_iterations = 0
    converged_stage = "successive-substitution" if split.converged else None

    if settings.second_order and residual > settings.second_order_tol:
        refined = _second_order_split(
            z=z, x_ii=x_ii, beta=beta, ln_gamma=ln_gamma, settings=settings
        )
        second_order_iterations = refined.iterations
        if refined.residual < residual:
            x_i, x_ii, beta = refined.x_i, refined.x_ii, refined.beta
            ln_f_x, ln_f_y = ln_gamma(x_i), ln_gamma(x_ii)
            residual = _equilibrium_residual(x_i, x_ii, ln_f_x, ln_f_y)
            converged_stage = "second-order"

    if residual > settings.tol and not split.converged:
        raise ConvergenceError(
            "flash_tp did not converge the liquid-liquid split; equal-activity residual="
            f"{residual:.3e} after {split.iterations} successive-substitution and "
            f"{second_order_iterations} second-order iterations."
        )

    checks = _verify_split(
        z=z,
        x=x_i,
        y=x_ii,
        beta=beta,
        ln_f_x=ln_f_x,
        ln_f_y=ln_f_y,
        ln_f_feed=ln_gamma(z / float(np.sum(z))),
        residual_key="equilibrium_residual",
    )

    post_split = _post_split_stability(
        mixture,
        temperature,
        pressure,
        eos=None,
        activity_model=activity_model,
        phases=((_LIQUID_I, x_i), (_LIQUID_II, x_ii)),
        settings=settings,
    )

    return _two_phase_result(
        mixture,
        temperature,
        pressure,
        x_i,
        x_ii,
        beta,
        names=(_LIQUID_I, _LIQUID_II),
        vapor_fraction=None,
        diagnostics={
            **base,
            "iterations": split.iterations + second_order_iterations,
            "converged": True,
            "termination_reason": "tolerance_met",
            "phase_count": 2,
            "phase_state": "two_phase",
            "phase_regime": "LLE",
            "k_seed": "stability",
            "ssi_iterations": split.iterations,
            "second_order_iterations": second_order_iterations,
            "converged_stage": converged_stage or "successive-substitution",
            "k_min": float(np.min(split.K)),
            "k_max": float(np.max(split.K)),
            **checks,
            **post_split,
        },
    )


def _flash_tp_wilson_heuristic(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState,
    activity_model: ActivityModel | None,
    mode: str,
    settings: FlashSettings,
    z: np.ndarray,
) -> FlashResult:
    """Legacy TP VLE solver: Wilson K-bounds and a Rachford-Rice root check.

    Kept unchanged (ADR-0008, ADR-0009) so the pre-``flash-auto-phase-detection``
    behavior stays reachable and testable, and because gamma-phi has no
    stability test to fall back on (ADR-0007). It is deliberately *not* extended
    with the post-split stability check: reproducing the old behavior is the
    whole point of this path.
    """
    K = wilson_k(mixture, temperature, pressure)
    if np.any(K <= 0.0):
        raise ModelError("Non-positive K-values encountered in Wilson estimate.")

    k_min = float(np.min(K))
    k_max = float(np.max(K))
    base: dict[str, float | int | str | bool] = {
        "flash_mode": mode,
        "phase_detection": "wilson-heuristic",
        "k_seed": "wilson",
    }

    if np.all(K <= 1.0) or np.all(K >= 1.0):
        phase_name = "liquid" if np.all(K <= 1.0) else "vapor"
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name=phase_name,
            vapor_fraction=0.0 if phase_name == "liquid" else 1.0,
            diagnostics={
                **base,
                "k_min": k_min,
                "k_max": k_max,
                "iterations": 0,
                "converged": True,
                "termination_reason": "single_phase_k_bounds",
                "max_delta_k": 0.0,
                "phase_count": 1,
                "phase_state": phase_name,
                "phase_regime": "single-phase",
            },
        )

    vapor_fraction, f0, f1 = _rachford_rice(z, K)
    if vapor_fraction is None:
        phase_name = "vapor" if f0 > 0.0 else "liquid"
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name=phase_name,
            vapor_fraction=1.0 if phase_name == "vapor" else 0.0,
            diagnostics={
                **base,
                "k_min": k_min,
                "k_max": k_max,
                "iterations": 0,
                "converged": True,
                "termination_reason": "rr_no_root",
                "max_delta_k": 0.0,
                "phase_count": 1,
                "phase_state": phase_name,
                "phase_regime": "single-phase",
                "rr_f0": float(f0),
                "rr_f1": float(f1),
                "rr_status": "no_root",
            },
        )

    split = _solve_k_loop(
        mixture,
        temperature,
        pressure,
        eos=eos,
        activity_model=activity_model,
        mode=mode,
        settings=settings,
        z=z,
        K=K,
        vapor_fraction=vapor_fraction,
    )
    return _two_phase_result(
        mixture,
        temperature,
        pressure,
        split.x,
        split.y,
        split.vapor_fraction,
        names=("liquid", "vapor"),
        vapor_fraction=split.vapor_fraction,
        diagnostics={
            **base,
            "iterations": split.iterations,
            "converged": True,
            "termination_reason": "tolerance_met",
            "max_delta_k": split.max_delta,
            "k_min": float(np.min(split.K)),
            "k_max": float(np.max(split.K)),
            "phase_count": 2,
            "phase_state": "two_phase",
            "phase_regime": "VLE",
            "post_split_checked": False,
            "post_split_skipped_reason": (
                "gamma_phi_stability_unsupported"
                if mode == "gamma-phi"
                else "legacy_wilson_heuristic_path"
            ),
        },
    )


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
) -> _SplitSolution:
    """Successive substitution on K with Rachford-Rice updates of ``beta``.

    This is the first-stage phase split, shared by all three modes; only its
    initial ``K``, its ``beta`` and the model call that updates ``K`` differ:

    - phi-phi: ``K = phi^L / phi^V`` (``x`` is the liquid, ``y`` the vapor);
    - gamma-phi: ``K = gamma^L phi^L / phi^V``;
    - gamma-gamma: ``K = gamma^I / gamma^II`` (``x`` is phase I, ``y`` phase II
      and ``vapor_fraction`` is the mole fraction of phase II).

    ``allow_unconverged`` returns the last iterate instead of raising when the
    budget runs out, which is how the liquid-liquid path hands over to its
    second-order stage.
    """
    budget = settings.max_iter if max_iter is None else max_iter
    max_delta = float("inf")
    x = np.array(z, dtype=float)
    y = np.array(z, dtype=float)
    ln_f_x = np.zeros_like(z)
    ln_f_y = np.zeros_like(z)

    for iteration in range(1, budget + 1):
        x = z / (1.0 + vapor_fraction * (K - 1.0))
        x = normalize_composition(x, label="liquid", error_cls=ConvergenceError)

        y = K * x
        y = normalize_composition(y, label="vapor", error_cls=ConvergenceError)

        if mode == "gamma-gamma":
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
            )

        if settings.damping is None:
            K = K_new
        else:
            K = K + settings.damping * (K_new - K)

        if np.any(K <= 0.0):
            raise ModelError("Non-positive K-values encountered during iteration.")

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
            iterations=budget,
            max_delta=max_delta,
            converged=False,
        )

    raise ConvergenceError(
        f"flash_tp did not converge within the iteration limit; max_delta_k={max_delta:.3e}."
    )


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


def _equilibrium_residual(
    x: np.ndarray, y: np.ndarray, ln_f_x: np.ndarray, ln_f_y: np.ndarray
) -> float:
    """``max_i |ln(x_i f_i^x) - ln(y_i f_i^y)|`` over components in both phases."""
    both = (x > 0.0) & (y > 0.0)
    if not np.any(both):
        return float("nan")
    return float(
        np.max(np.abs((np.log(x[both]) + ln_f_x[both]) - (np.log(y[both]) + ln_f_y[both])))
    )


def _stability_k_seed(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    z: np.ndarray,
    w: np.ndarray,
    tpd_min: float,
) -> tuple[np.ndarray, str]:
    """Initial K-values from the tangent-plane minimizer.

    The stability test returns the normalized incipient composition ``w`` and
    the reduced tangent-plane distance at that stationary point. Michelsen's
    unnormalized mole numbers follow from the stationary-point identity
    ``tpd = -ln(sum_i W_i)``, so ``W = w * exp(-tpd)``. Using ``W`` rather than
    ``w`` matters: with ``K_i = w_i / z_i`` the Rachford-Rice function at
    ``beta = 0`` is ``sum_i w_i - 1 = 0`` exactly, a degenerate root that cannot
    be bracketed, whereas ``K_i = W_i / z_i`` gives ``sum_i W_i - 1 > 0`` for an
    unstable feed. Both seeds have the same fixed point; only the bracket
    differs.

    The returned K is ``y / x``, so it is ``W / z`` when the incipient phase is
    the vapor-like one and ``z / W`` when it is the liquid-like one. Which of
    the two converged phases is *named* "vapor" is a labelling convention -
    ``EquationOfState`` exposes no molar volume - and is resolved by volatility
    ordering: see :func:`_incipient_is_vapor_like`.

    Returns:
        ``(K, incipient_phase)`` with ``incipient_phase`` in
        ``{"vapor", "liquid"}``.
    """
    sum_capital_w = math.exp(-tpd_min) if math.isfinite(tpd_min) else 1.0
    capital_w = w * sum_capital_w

    incipient_vapor = _incipient_is_vapor_like(mixture, temperature, pressure, z=z, w=w)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = capital_w / z if incipient_vapor else z / capital_w

    usable = np.isfinite(ratio) & (ratio > 0.0)
    K = np.where(usable, ratio, _INERT_SEED_K)
    return K, "vapor" if incipient_vapor else "liquid"


def _incipient_is_vapor_like(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    z: np.ndarray,
    w: np.ndarray,
) -> bool:
    """Name the incipient phase "vapor" or "liquid" by volatility ordering.

    The tangent-plane test finds *a* second phase; it does not say which of the
    two converged phases should be called the vapor. The successive-substitution
    loop is symmetric under swapping the two phases (``x <-> y``,
    ``beta <-> 1 - beta``, ``K <-> 1/K``), so an arbitrary orientation would
    return the mirror-labelled solution roughly half the time.

    The label is decided here by the only volatility ordering the package has:
    the Wilson correlation's ranking of the components (built from ``Tc``,
    ``Pc`` and ``omega``). Writing ``hi`` for the component with the largest
    Wilson K among those present in the feed and ``lo`` for the smallest, the
    incipient phase is called vapor-like when it is enriched in ``hi`` relative
    to ``lo`` compared with the feed, i.e. when
    ``ln(w_hi / z_hi) - ln(w_lo / z_lo) >= 0``.

    Only the *ranking* is used, never the magnitudes, and it decides the name
    only - never the one-versus-two-phase verdict, the compositions, or the
    vapor fraction, all of which come from the tangent-plane test and the
    converged split. The minimum-Gibbs branch label reported by the stability
    test is deliberately *not* used for this: whenever the cubic has a single
    real root (common for dense and near-critical states) both branch calls
    return identical fugacity coefficients and the label is only a tie-break.

    There is no analogue for two liquid phases, which is why the liquid-liquid
    path names its phases by role instead (see :func:`flash_tp`).
    """
    active = z > 0.0
    if int(np.count_nonzero(active)) < 2:
        return True

    k_wilson = wilson_k(mixture, temperature, pressure)
    masked = np.where(active, k_wilson, np.nan)
    hi = int(np.nanargmax(masked))
    lo = int(np.nanargmin(masked))
    if hi == lo:
        return True

    tiny = np.finfo(float).tiny
    score = math.log(max(float(w[hi]), tiny) / float(z[hi])) - math.log(
        max(float(w[lo]), tiny) / float(z[lo])
    )
    return score >= 0.0


def _verify_split(
    *,
    z: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    beta: float,
    ln_f_x: np.ndarray,
    ln_f_y: np.ndarray,
    ln_f_feed: np.ndarray,
    residual_key: str,
) -> dict[str, float | int | str | bool]:
    """Residuals that a converged two-phase solution must satisfy.

    - ``mass_balance_residual``: ``max_i |z_i - (beta y_i + (1 - beta) x_i)|``.
    - ``residual_key``: ``max_i |ln(x_i f_i^x) - ln(y_i f_i^y)|`` with
      ``f = phi`` (phi-phi, key ``fugacity_residual``) or ``f = gamma``
      (gamma-gamma, key ``equilibrium_residual``). The pressure (phi-phi) or the
      pure-liquid reference (gamma-gamma) cancels, so this is the equal-fugacity
      condition in log form.
    - ``delta_g_split_rt``: ``beta sum_i y_i ln(y_i f_i^y)
      + (1 - beta) sum_i x_i ln(x_i f_i^x) - sum_i z_i ln(z_i f_i^feed)``,
      which must be negative for the split to be an improvement on the
      single-phase feed. For phi-phi the feed term uses the *minimum-Gibbs*
      branch (the stability module's convention, i.e. the lowest single-phase
      Gibbs energy available), while each split phase uses the branch the solver
      actually converged on. Since the min-Gibbs branch is never higher in Gibbs
      energy, this combination is a conservative (upper-bound) estimate of the
      true reduction. An activity model has one branch, so for gamma-gamma the
      quantity is exact.
    """
    beta = float(beta)
    balance = float(np.max(np.abs(z - (beta * y + (1.0 - beta) * x))))
    residual = _equilibrium_residual(x, y, ln_f_x, ln_f_y)

    def _reduced_g(fractions: np.ndarray, ln_f: np.ndarray) -> float:
        mask = fractions > 0.0
        return float(np.sum(fractions[mask] * (np.log(fractions[mask]) + ln_f[mask])))

    delta_g = (
        beta * _reduced_g(y, ln_f_y)
        + (1.0 - beta) * _reduced_g(x, ln_f_x)
        - _reduced_g(z, ln_f_feed)
    )

    return {
        "mass_balance_residual": balance,
        residual_key: residual,
        "delta_g_split_rt": delta_g,
    }


def _post_split_stability(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState | None,
    activity_model: ActivityModel | None,
    phases: Sequence[tuple[str, np.ndarray]],
    settings: FlashSettings,
) -> dict[str, float | int | str | bool]:
    """Test every converged phase for stability and adjudicate the phase set.

    Two coexisting phases share one tangent plane, so a stability test run on
    either of them finds its *partner* as a stationary point with ``tpd = 0``
    (validation Case S-3). Up to the split's own convergence tolerance that
    zero is a small negative number - measured worst case -7.0e-09 over the
    in-repo Peng-Robinson grid - so it can dip below ``-tpd_tol`` for a loosely
    converged split. A minimizer that *is* the partner phase is therefore
    reported as ``"marginal"``, not as an instability. "Is the partner" uses the
    same measure as the stability module's trivial-solution test,
    ``sum_i ln(w_i / x_i^partner)^2 < trivial_tol``, rather than a new
    tolerance.

    Anything else negative means the two-phase answer is not a stable phase set
    and a third phase is required, which this release cannot produce.

    The check always runs on the paths that can run it;
    ``settings.post_split_stability`` decides whether a *failure* raises or is
    only reported, so the diagnostics show the failure either way.

    Returns:
        Diagnostics keys ``post_split_checked``, ``post_split_stable``,
        ``post_split_status``, ``post_split_tpd_min``,
        ``phase_stability_<name>`` and ``phase_stability_tpd_min_<name>``.
    """
    from ..stability import StabilitySettings, stability_tp

    stability_settings = settings.stability_settings or StabilitySettings()
    diagnostics: dict[str, float | int | str | bool] = {"post_split_checked": True}
    worst = math.inf
    failures: list[str] = []
    inconclusive: list[str] = []

    for position, (name, composition) in enumerate(phases):
        partner = phases[1 - position][1] if len(phases) == 2 else None
        phase_mixture = Mixture(
            components=mixture.components,
            composition=Composition(
                fractions=tuple(float(value) for value in composition),
                basis=mixture.basis,
                normalize=True,
                tol=COMPOSITION_SUM_TOL,
            ),
        )
        result = stability_tp(
            phase_mixture,
            temperature_K=temperature,
            pressure_Pa=pressure,
            eos=eos,
            activity_model=activity_model,
            settings=stability_settings,
        )

        verdict = result.status
        if verdict == "unstable" and partner is not None and result.trial_composition is not None:
            if _is_same_phase(
                np.array(result.trial_composition, dtype=float),
                partner,
                stability_settings.trivial_tol,
            ):
                verdict = "marginal"

        diagnostics[f"phase_stability_{name}"] = verdict
        diagnostics[f"phase_stability_tpd_min_{name}"] = float(result.tpd_min)
        worst = min(worst, float(result.tpd_min))
        if verdict == "unstable":
            failures.append(name)
        elif verdict == "inconclusive":
            inconclusive.append(name)

    diagnostics["post_split_tpd_min"] = worst if math.isfinite(worst) else float("nan")
    status = "unstable" if failures else ("inconclusive" if inconclusive else "stable")
    diagnostics["post_split_status"] = status
    diagnostics["post_split_stable"] = status == "stable"

    if status != "stable" and settings.post_split_stability:
        detail = ", ".join(failures or inconclusive)
        raise ConvergenceError(
            "The converged two-phase solution is not a stable phase set: the post-split "
            f"stability test reports '{status}' for phase(s) {detail} "
            f"(most negative post-split tpd = {worst:.6e}). A third phase is required, and "
            "flash_tp returns at most two phases in this release (multiphase flash is the "
            "next slice). Pass FlashSettings(post_split_stability=False) to receive the "
            "two-phase result anyway, with this failure recorded in diagnostics."
        )

    return diagnostics


def _is_same_phase(w: np.ndarray, other: np.ndarray, trivial_tol: float) -> bool:
    """True when ``w`` is the composition ``other`` in the trivial-solution metric."""
    mask = (w > 0.0) & (other > 0.0)
    if not np.any(mask):
        return False
    if np.any((w > 0.0) != (other > 0.0)):
        return False
    ln_ratio = np.log(w[mask] / other[mask])
    return bool(float(np.sum(ln_ratio * ln_ratio)) < trivial_tol)


def _single_phase_result(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    *,
    phase_name: str,
    vapor_fraction: float | None,
    diagnostics: dict[str, float | int | str | bool],
) -> FlashResult:
    composition = Composition(
        fractions=mixture.fractions, basis=mixture.basis, normalize=False, tol=COMPOSITION_SUM_TOL
    )
    phase = PhaseResult(name=phase_name, composition=composition)
    phase_fractions = {phase_name: 1.0}
    return FlashResult(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        phases={phase_name: phase},
        vapor_fraction=vapor_fraction,
        phase_fractions=phase_fractions,
        diagnostics=diagnostics,
    )


def _two_phase_result(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    x: np.ndarray,
    y: np.ndarray,
    beta: float,
    *,
    names: tuple[str, str],
    vapor_fraction: float | None,
    diagnostics: dict[str, float | int | str | bool],
) -> FlashResult:
    first_name, second_name = names
    first = PhaseResult(
        name=first_name,
        composition=Composition(
            fractions=tuple(x.tolist()),
            basis=mixture.basis,
            normalize=False,
            tol=COMPOSITION_SUM_TOL,
        ),
    )
    second = PhaseResult(
        name=second_name,
        composition=Composition(
            fractions=tuple(y.tolist()),
            basis=mixture.basis,
            normalize=False,
            tol=COMPOSITION_SUM_TOL,
        ),
    )
    phase_fractions = {first_name: 1.0 - float(beta), second_name: float(beta)}
    return FlashResult(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        phases={first_name: first, second_name: second},
        vapor_fraction=None if vapor_fraction is None else float(vapor_fraction),
        phase_fractions=phase_fractions,
        diagnostics=diagnostics,
    )


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
