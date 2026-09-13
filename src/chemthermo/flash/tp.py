"""TP flash calculations (T,P) in SI units.

Supports phi-phi (vapor-liquid, equation of state), modified-raoult
(low-pressure vapor-liquid *and* liquid-liquid from an activity model with
Antoine reference fugacities and an ideal vapor), gamma-gamma (liquid-liquid,
both phases described by one activity model) and the deprecated gamma-phi
(activity-model liquid against an equation-of-state vapor). Given the same
inputs and settings, results are deterministic.

Phase detection (phi-phi, gamma-gamma and modified-raoult)
-----------------------------------------------------------
Every reference path decides one phase versus two from Michelsen's tangent-plane
stability criterion rather than from Wilson K-value bounds (ADR-0008, ADR-0009,
ADR-0010). The flow is

    flash_tp -> stability_tp(feed) -> single phase | seeded split
             -> post-split stability of every converged phase

1. ``stability_tp`` is run on the feed at the same ``(T, P)`` with the same
   model (``eos=`` for phi-phi, ``activity_model=`` for gamma-gamma,
   ``activity_model=`` plus ``vapor="ideal"`` for modified-raoult).
2. ``status == "stable"``: a single-phase ``FlashResult`` is returned with
   ``termination_reason = "feed_stable_tangent_plane"``. For phi-phi its phase
   name is the minimum-Gibbs root branch the stability test selected
   (``feed_branch``); when the cubic has a single real root both branches
   coincide and the label is a *convention*, not a measurement. For
   gamma-gamma the single phase is a liquid and is named ``"liquid"``.
3. ``status == "unstable"``: the converged stationary point seeds the K-values
   (see :func:`chemthermo.flash._detect._stability_k_seed`) and the
   successive-substitution / Rachford-Rice loop - Michelsen's recommended
   first-stage phase split - runs from there. For phi-phi, if the seeded
   K-values give no Rachford-Rice root the Wilson estimate is tried as a
   documented fallback, and ``diagnostics["k_seed"]`` records which seed was
   actually used.
4. ``status == "inconclusive"``: a :class:`chemthermo.ConvergenceError` is
   raised. A stability search that could not converge must not silently produce
   a single-phase answer.

The converged split is then verified (material balance, phase fractions, equal
fugacities or equal activities, and a negative Gibbs-energy change against the
single-phase feed) and every residual is reported in ``diagnostics``.

Post-split stability, and phase addition / removal
--------------------------------------------------
Every converged result on these paths is re-tested: each phase is fed back into
``stability_tp`` with the same model. Coexisting phases share one tangent
plane, so each of them is *marginally* stable with respect to the others - a
trial that converges onto a partner phase has ``tpd = 0`` up to the split's own
convergence tolerance and is classified ``"marginal"`` here, not as an
instability (validation Case S-3). A phase whose tangent-plane minimum is
genuinely negative somewhere else means the phase set is not the answer.

On the ``modified-raoult`` path that failure is now *resolved* rather than
refused (ADR-0011): the minimizer found on the failing phase is the incipient
new phase, it is added, and the enlarged set is re-solved with the multiphase
Rachford-Rice of :mod:`chemthermo.flash._multiphase_rr`; a phase whose fraction
converges to zero or below is removed again. The search is bounded by
``FlashSettings.max_phases`` (default 3) and the sets it visited are reported
in ``diagnostics["phase_set_history"]``. ``max_phases=2`` reproduces the
pre-ADR-0011 behavior, which raises
:class:`chemthermo.ConvergenceError`; so does the phi-phi and gamma-gamma path
at any ``max_phases``, because no state in this repository exercises a third
phase there (ADR-0011 "What remains").
``FlashSettings(post_split_stability=False)`` returns the two-phase result
anyway with the failure recorded in ``diagnostics``, without searching.

Modified Raoult (low-pressure gamma-phi)
----------------------------------------
``flash_mode="modified-raoult"`` puts an activity-coefficient liquid and an
ideal-gas vapor on **one** Gibbs surface, using the pure-liquid reference
fugacity ``f_i^0 = Psat_i(T)`` from the databank Antoine coefficients
(``phi_i^sat = 1``, Poynting = 1, ``phi_i^V = 1``). The equilibrium condition is
modified Raoult's law ``y_i P = x_i gamma_i(x) Psat_i(T)``; in tangent-plane
form the two candidates contribute

    liquid: ln gamma_i(w) + ln( Psat_i(T) / P )      vapor: 0

against the common reference ``ln( f_i / (x_i P) )``. One stability test
therefore detects a vapor-liquid split, a liquid-liquid split or neither, and
the candidate label of the stationary point says which; the split then
evaluates each phase with the candidate assigned to it, so
``K_i = gamma_i Psat_i / P`` (VLE) and ``K_i = gamma_i^I / gamma_i^II`` (LLE)
are the same update rule. See :func:`chemthermo.flash._detect._flash_tp_modified_raoult`.

``flash_mode="gamma-phi"`` is **deprecated** in favour of this mode; see
:func:`flash_tp`.

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
:func:`chemthermo.flash._second_order._second_order_split` for the derivation.

Limits
------
Up to ``FlashSettings.max_phases`` phases are returned on the
``modified-raoult`` path; the phi-phi and gamma-gamma paths still stop at two
and raise when a third is needed. A negative ``tpd_min`` proves a feed is not
one phase; ``"stable"`` only means no negative tangent-plane distance was found
from the deterministic trial set - so a phase count is never more reliable than
the stability test that produced it. Each trial of that set runs on one fixed
phase candidate (ADR-0012), which is what lets a thin three-phase region be
found near a plait point; over the ternary grid of validation Case V-5 the
verdict agrees with an independent lowest-Gibbs classifier at every feed. That
is evidence, not a global proof.

Module layout
-------------
This module is the thin public orchestrator: it validates inputs, resolves the
mode, and dispatches to the internal module that implements it - ``_detect``
(phase detection and split seeding), ``_split`` (the shared K-loop), ``_second_order``
(the liquid-liquid Newton stage), ``_multiphase_rr`` (the multiphase
Rachford-Rice), ``_multiphase`` (the multiphase split and the phase
addition/removal loop), ``_verify`` (residuals and post-split stability),
``_assemble`` (``FlashResult`` construction) and ``_legacy`` (the
``wilson-heuristic`` path). All are internal (ADR-0001): none is re-exported.
"""

from __future__ import annotations

import numpy as np

from ..core import Mixture
from ..exceptions import CompositionError, ModelError
from ..models import ActivityModel, EquationOfState
from ..validation import validate_pressure, validate_temperature
from ._detect import (
    _flash_tp_liquid_liquid,
    _flash_tp_modified_raoult,
    _flash_tp_tangent_plane,
)
from ._legacy import _flash_tp_wilson_heuristic
from .results import FlashResult
from .settings import FlashSettings

#: Supported ``flash_mode`` values.
#:
#: ``"gamma-phi"`` is **deprecated** in favour of ``"modified-raoult"``; see the
#: :func:`flash_tp` docstring and ADR-0010. It is not removed and its numbers
#: are unchanged.
FLASH_MODES = ("phi-phi", "gamma-phi", "gamma-gamma", "modified-raoult")


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
            ``"gamma-gamma"`` and ``"modified-raoult"``.
        activity_model: Activity model. Required for ``"gamma-phi"`` (liquid
            phase), ``"gamma-gamma"`` (both liquid phases) and
            ``"modified-raoult"`` (the liquid candidate).
        flash_mode: Case-insensitive mode: ``"phi-phi"``, ``"modified-raoult"``,
            ``"gamma-gamma"`` or the deprecated ``"gamma-phi"``. ``None`` (the
            default) infers the mode from the models supplied:
            ``"gamma-gamma"`` when only ``activity_model`` is given,
            ``"phi-phi"`` otherwise. ``"modified-raoult"`` and ``"gamma-phi"``
            are never inferred and must be named: which vapor model applies at
            a given pressure is the caller's physical judgement, not
            something the package should guess.
        settings: Iteration controls (tolerance, damping, max iterations,
            phase-detection mode, stability settings, post-split check,
            second-order stage).

    Returns:
        FlashResult with phase compositions and fractions. Phase names follow:
        VLE -> ``"liquid"``/``"vapor"``, LLE -> ``"liquid1"``/``"liquid2"``,
        VLLE -> ``"liquid1"``/``"liquid2"``/``"vapor"``, single phase ->
        ``"liquid"`` or ``"vapor"``. ``vapor_fraction`` is ``None`` for every
        gamma-gamma result and for any result with no ``"vapor"`` phase in it:
        reporting a number there would be fiction. A result that does contain a
        vapor carries that phase's mole fraction.

    Diagnostics:
        Diagnostics keys are implementation details. Current stable keys:

        - Always: ``flash_mode``, ``phase_detection``, ``iterations``,
          ``converged``, ``termination_reason``, ``phase_count``,
          ``phase_state``, ``phase_regime``.
        - Results that went through the phase addition/removal search
          (``modified-raoult`` only, and only when a converged phase set failed
          its post-split test) add ``phase_set_history``, ``phases_added``,
          ``phases_removed``, ``rachford_rice_iterations`` and, when the search
          started from a converged two-phase set, ``delta_g_vs_two_phase_rt``.
          Those keys are **absent** from every other result, deliberately: the
          two-phase numbers of the earlier slices are unchanged down to the
          last bit, diagnostics included.
        - Tangent-plane paths (phi-phi default, gamma-gamma and
          modified-raoult):
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
          ``second_order_iterations`` and ``converged_stage``. Modified-raoult
          reports ``incipient_phase``, ``k_min``, ``k_max``, the three stage
          keys, and ``antoine_valid_Tmin_K`` / ``antoine_valid_Tmax_K`` (the
          intersection of the components' Antoine validity ranges) on every
          result.
        - Legacy heuristic path: ``k_min``, ``k_max``, ``max_delta_k``,
          ``k_seed`` (``"wilson"``), and, for its single-phase fallbacks,
          ``rr_f0``, ``rr_f1``, ``rr_status``.

    Raises:
        InputRangeError: If temperature or pressure is non-physical.
        ModelError: If required models are missing, if models are combined in an
            unsupported way, or if a model returns invalid values.
        InputRangeError: (modified-raoult) If the temperature is outside the
            Antoine validity range of a component.
        PropertyNotFoundError: (modified-raoult) If a component has no Antoine
            record.
        CompositionError: If the mixture composition is invalid.
        ConvergenceError: If iteration fails to converge; or (tangent-plane
            modes only) if the stability analysis is inconclusive, if an
            unstable feed admits no Rachford-Rice root from either seed, or if
            a converged phase set fails the post-split stability check and no
            further phase may be added (phi-phi and gamma-gamma always, or
            ``FlashSettings.max_phases`` reached on the modified-Raoult path).

    Notes:
        **``flash_mode="gamma-phi"`` is DEPRECATED** in favour of
        ``"modified-raoult"``. It is not removed, nothing about it changed, and
        its removal would need its own ADR. It is deprecated because it is not
        a consistent model: it sets ``K_i = gamma_i phi_i^L / phi_i^V`` with
        ``gamma`` from the activity model *and* ``phi^L`` from the equation of
        state evaluated on the liquid mixture, so the liquid's nonideality is
        counted twice, and it carries no pure-liquid reference fugacity at all
        (no ``Psat_i``, no ``phi_i^sat``, no Poynting), so the two phases are
        not on one Gibbs surface. That is also why it has no stability test and
        no post-split check (ADR-0007, ADR-0009). ``"modified-raoult"`` is the
        low-pressure model written down correctly; for high pressure use
        ``"phi-phi"``.

        **Modified-Raoult limits.** Ideal vapor (no ``phi^V``), so low pressure
        only; no Poynting correction and no ``phi^sat``; and the temperature
        must lie inside every component's Antoine validity range, which is
        enforced rather than extrapolated.

        With ``phase_detection="tangent-plane"`` (the default for phi-phi and
        the only option for gamma-gamma and modified-raoult) a single-phase
        result means the feed was *found stable* by Michelsen's test. With
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
        labels. Compare the phase *set*, not ``result.phases["liquid1"]``. The
        same holds for the two liquids of a three-phase result, whose ordering
        follows the order in which the search happened to create them; only the
        ``"vapor"`` name carries a model-level meaning (it is the phase the
        ideal-gas candidate describes).
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
    elif mode == "modified-raoult":
        if activity_model is None:
            raise ModelError("An activity model is required for modified-Raoult flash.")
        if eos is not None:
            raise ModelError(
                "modified-raoult flash describes the vapor as an ideal gas and the liquid "
                "with the activity model plus Antoine reference fugacities, so 'eos' must "
                "not be given."
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

    if mode == "modified-raoult":
        assert activity_model is not None
        return _flash_tp_modified_raoult(
            mixture,
            temperature,
            pressure,
            activity_model=activity_model,
            settings=settings,
            z=z,
        )

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
