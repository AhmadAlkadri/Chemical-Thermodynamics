"""Tangent-plane phase detection: decide 1-vs-2 phases and seed the split.

Both reference paths - phi-phi (:func:`_flash_tp_tangent_plane`) and
gamma-gamma (:func:`_flash_tp_liquid_liquid`) - decide one phase versus two
from Michelsen's tangent-plane stability criterion (ADR-0008, ADR-0009)
rather than from Wilson K-value bounds:

    stability_tp(feed) -> single phase | seeded split
                        -> post-split stability of every converged phase

``stability_tp`` is run on the feed; a stable feed returns a single-phase
result immediately, an unstable feed seeds the shared split loop
(:mod:`chemthermo.flash._split`) from the stationary point
(:func:`_stability_k_seed`), and the converged split is verified and
post-split-checked (:mod:`chemthermo.flash._verify`) before being assembled
into a `FlashResult` (:mod:`chemthermo.flash._assemble`). An inconclusive
stability result raises rather than guessing.
"""

from __future__ import annotations

import math

import numpy as np

from ..core import Mixture
from ..exceptions import ConvergenceError
from ..models import ActivityModel, EquationOfState
from ._assemble import _single_phase_result, _two_phase_result
from ._common import wilson_k
from ._multiphase import _flash_tp_phase_addition, _phase_set_label
from ._second_order import _second_order_split
from ._split import _ln_gamma_function, _rachford_rice, _solve_k_loop
from ._verify import (
    _equilibrium_residual,
    _post_split_report,
    _post_split_stability,
    _reduced_g,
    _verify_split,
)
from .results import FlashResult
from .settings import FlashSettings

#: Seed K-value used for components absent from the feed (``z_i == 0``). Those
#: components have ``x_i = y_i = 0`` at every iteration and are rewritten from
#: the model on the first update, so the seed value cannot affect the result.
_INERT_SEED_K = 1.0

#: Phase names of a liquid-liquid result. Roles, not identities; see
#: :func:`chemthermo.flash.tp.flash_tp`.
_LIQUID_I = "liquid1"
_LIQUID_II = "liquid2"

#: Phase-candidate labels of the modified-Raoult pair, which are also the phase
#: names of a vapor-liquid result on that path.
_LIQUID = "liquid"
_VAPOR = "vapor"
_MODIFIED_RAOULT = "modified-raoult"


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
    :func:`chemthermo.flash._second_order._second_order_split` when the
    equal-activity residual is still above ``settings.second_order_tol``.
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
            z=z, x_ii=x_ii, beta=beta, terms_i=ln_gamma, settings=settings
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


def _flash_tp_modified_raoult(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    activity_model: ActivityModel,
    settings: FlashSettings,
    z: np.ndarray,
) -> FlashResult:
    """Low-pressure gamma-phi TP flash from one tangent plane (ADR-0010).

    The feed is tested against **two** phase candidates at once - an
    activity-coefficient liquid with Antoine pure-liquid reference fugacities
    and an ideal-gas vapor - so the same call returns a vapor-liquid split, a
    liquid-liquid split, or a single phase, and the *candidate label* of the
    stationary point is what decides which:

    ========================  ========================  ====================
    feed candidate            incipient candidate       result
    ========================  ========================  ====================
    liquid                    vapor                     VLE, bubble side
    vapor                     liquid                    VLE, dew side
    liquid                    liquid                    LLE
    ========================  ========================  ====================

    The split loop then evaluates each phase with the candidate it was assigned
    (:func:`chemthermo.flash._split._solve_k_loop`), which makes
    ``K_i = gamma_i Psat_i / P`` for a vapor-liquid pair and
    ``K_i = gamma_i^I / gamma_i^II`` for a liquid-liquid pair without the loop
    knowing the difference. Every converged phase is then re-tested against
    **both** candidates, so a state that needs all three phases raises rather
    than being returned (see :func:`chemthermo.flash._verify._post_split_stability`).
    """
    from ..stability import stability_tp
    from ..stability._evaluator import modified_raoult_candidates

    stability = stability_tp(
        mixture,
        temperature_K=temperature,
        pressure_Pa=pressure,
        activity_model=activity_model,
        vapor="ideal",
        settings=settings.stability_settings,
    )

    base: dict[str, float | int | str | bool] = {
        "flash_mode": _MODIFIED_RAOULT,
        "phase_detection": "tangent-plane",
        "stability_status": stability.status,
        "tpd_min": float(stability.tpd_min),
        "stability_trials": len(stability.trials),
    }
    for key in ("antoine_valid_Tmin_K", "antoine_valid_Tmax_K"):
        if key in stability.diagnostics:
            base[key] = stability.diagnostics[key]
    if stability.feed_branch is not None:
        base["feed_branch"] = stability.feed_branch

    if stability.status == "inconclusive":
        raise ConvergenceError(
            "Tangent-plane stability analysis was inconclusive (no trial converged), so "
            "flash_tp cannot decide whether the feed is one phase or two. Loosen "
            "FlashSettings.stability_settings and try again."
        )

    feed_label = stability.feed_branch or _LIQUID
    if stability.status == "stable":
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name=feed_label,
            vapor_fraction=1.0 if feed_label == _VAPOR else 0.0,
            diagnostics={
                **base,
                "iterations": 0,
                "converged": True,
                "termination_reason": "feed_stable_tangent_plane",
                "phase_count": 1,
                "phase_state": feed_label,
                "phase_regime": "single-phase",
            },
        )

    trial = stability.trial_composition
    if trial is None:
        raise ConvergenceError(
            "Tangent-plane stability reported an unstable feed without a minimizing "
            "trial composition; no phase split can be seeded."
        )
    incipient_label = stability.phase_branch or _LIQUID
    if feed_label == _VAPOR and incipient_label == _VAPOR:
        raise ConvergenceError(
            "The tangent-plane minimizer is a second ideal-gas vapor, which cannot "
            "coexist with the first: an ideal gas has no composition range over which "
            "it demixes. This indicates a model or numerical failure, not a phase split."
        )

    liquid, vapor = modified_raoult_candidates(
        activity_model, mixture=mixture, temperature=temperature, pressure=pressure
    )
    candidates = {_LIQUID: liquid, _VAPOR: vapor}
    terms_x = candidates[feed_label].ln_fugacity_terms
    terms_y = candidates[incipient_label].ln_fugacity_terms

    # Same seed as the other tangent-plane paths: Michelsen's unnormalized mole
    # numbers W = w exp(-tpd), with phase x feed-like and phase y incipient, so
    # K_i = y_i / x_i = W_i / z_i and f_RR(0) = sum_i W_i - 1 > 0 is bracketable.
    w = np.array(trial, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = (w * math.exp(-float(stability.tpd_min))) / z
    usable = np.isfinite(ratio) & (ratio > 0.0)
    k_seed = np.where(usable, ratio, _INERT_SEED_K)

    beta, _f0, _f1 = _rachford_rice(z, k_seed)
    if beta is None:
        raise ConvergenceError(
            f"Feed is unstable (tpd_min={stability.tpd_min:.6e}) but the stability-seeded "
            "K-values do not bracket a Rachford-Rice root, so no split can be started."
        )

    split = _solve_k_loop(
        mixture,
        temperature,
        pressure,
        eos=None,
        activity_model=activity_model,
        mode=_MODIFIED_RAOULT,
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
        terms_x=terms_x,
        terms_y=terms_y,
    )

    x, y, beta = split.x, split.y, split.vapor_fraction
    ln_f_x, ln_f_y = split.ln_f_x, split.ln_f_y
    residual = _equilibrium_residual(x, y, ln_f_x, ln_f_y)
    second_order_iterations = 0
    converged_stage = "successive-substitution" if split.converged else None

    if settings.second_order and residual > settings.second_order_tol:
        refined = _second_order_split(
            z=z,
            x_ii=y,
            beta=beta,
            terms_i=terms_x,
            terms_ii=terms_y,
            settings=settings,
        )
        second_order_iterations = refined.iterations
        if refined.residual < residual:
            x, y, beta = refined.x_i, refined.x_ii, refined.beta
            ln_f_x, ln_f_y = terms_x(x), terms_y(y)
            residual = _equilibrium_residual(x, y, ln_f_x, ln_f_y)
            converged_stage = "second-order"

    if residual > settings.tol and not split.converged:
        raise ConvergenceError(
            "flash_tp did not converge the modified-Raoult split; equilibrium residual="
            f"{residual:.3e} after {split.iterations} successive-substitution and "
            f"{second_order_iterations} second-order iterations."
        )

    checks = _verify_split(
        z=z,
        x=x,
        y=y,
        beta=beta,
        ln_f_x=ln_f_x,
        ln_f_y=ln_f_y,
        ln_f_feed=terms_x(z / float(np.sum(z))),
        residual_key="equilibrium_residual",
    )

    if feed_label == _LIQUID and incipient_label == _LIQUID:
        names = (_LIQUID_I, _LIQUID_II)
        regime = "LLE"
        vapor_fraction: float | None = None
    else:
        names = (feed_label, incipient_label)
        regime = "VLE"
        vapor_fraction = float(beta) if incipient_label == _VAPOR else 1.0 - float(beta)

    report = _post_split_report(
        mixture,
        temperature,
        pressure,
        eos=None,
        activity_model=activity_model,
        phases=((names[0], x), (names[1], y)),
        settings=settings,
        vapor="ideal",
    )

    if report.status != "stable" and settings.post_split_stability:
        if settings.max_phases < 3:
            detail = ", ".join(
                [failure.phase_name for failure in report.instabilities] or report.inconclusive
            )
            raise ConvergenceError(
                "The converged two-phase solution is not a stable phase set: the "
                f"post-split stability test reports '{report.status}' for phase(s) "
                f"{detail} (most negative post-split tpd = {report.tpd_min:.6e}). A third "
                "phase is required, and FlashSettings.max_phases = "
                f"{settings.max_phases} forbids it. Raise max_phases (the default 3 "
                "resolves this state), or pass FlashSettings(post_split_stability=False) "
                "to receive the two-phase result anyway, with this failure recorded in "
                "diagnostics."
            )
        if report.status == "inconclusive":
            raise ConvergenceError(
                "A post-split stability test was inconclusive for phase(s) "
                f"{', '.join(report.inconclusive)}, so flash_tp cannot decide whether the "
                "two-phase set is the answer."
            )
        # The phase set is provably not the answer, and the stability minimizer
        # found on the failing phase is the incipient third phase. Hand over to
        # the phase addition / removal search (ADR-0011).
        failure = report.instabilities[0]
        history = [_phase_set_label((feed_label,)), _phase_set_label(names)]
        history.append(_phase_set_label((*names, failure.branch or _LIQUID)))
        two_phase_g_rt = (1.0 - float(beta)) * _reduced_g(x, ln_f_x) + float(beta) * _reduced_g(
            y, ln_f_y
        )
        return _flash_tp_phase_addition(
            mixture,
            temperature,
            pressure,
            z=z,
            candidates={
                _LIQUID: candidates[_LIQUID].ln_fugacity_terms,
                _VAPOR: candidates[_VAPOR].ln_fugacity_terms,
            },
            labels=(feed_label, incipient_label, failure.branch or _LIQUID),
            compositions=(x, y, failure.composition),
            history=history,
            ln_f_feed=terms_x(z / float(np.sum(z))),
            two_phase_g_rt=two_phase_g_rt,
            activity_model=activity_model,
            vapor="ideal",
            settings=settings,
            base=base,
        )

    post_split = report.diagnostics

    return _two_phase_result(
        mixture,
        temperature,
        pressure,
        x,
        y,
        beta,
        names=names,
        vapor_fraction=vapor_fraction,
        diagnostics={
            **base,
            "iterations": split.iterations + second_order_iterations,
            "converged": True,
            "termination_reason": "tolerance_met",
            "phase_count": 2,
            "phase_state": "two_phase",
            "phase_regime": regime,
            "k_seed": "stability",
            "incipient_phase": incipient_label,
            "ssi_iterations": split.iterations,
            "second_order_iterations": second_order_iterations,
            "converged_stage": converged_stage or "successive-substitution",
            "k_min": float(np.min(split.K)),
            "k_max": float(np.max(split.K)),
            **checks,
            **post_split,
        },
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
    path names its phases by role instead (see
    :func:`chemthermo.flash.tp.flash_tp`).
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
