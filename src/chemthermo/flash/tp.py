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
:func:`chemthermo.flash._second_order._second_order_split` for the derivation.

Limits of this slice
--------------------
At most two phases are returned. A feed that needs three is now *reported*
(``ConvergenceError`` from the post-split check) rather than silently returned
as two, but it is still not solved; multiphase flash is the next slice. A
negative ``tpd_min`` proves a feed is not one phase; ``"stable"`` only means no
negative tangent-plane distance was found from the deterministic trial set.

Module layout
-------------
This module is the thin public orchestrator: it validates inputs, resolves the
mode, and dispatches to the internal module that implements it - ``_detect``
(phase detection and split seeding), ``_split`` (the shared K-loop), ``_second_order``
(the liquid-liquid Newton stage), ``_verify`` (residuals and post-split
stability), ``_assemble`` (``FlashResult`` construction) and ``_legacy`` (the
``wilson-heuristic`` path). All are internal (ADR-0001): none is re-exported.
"""

from __future__ import annotations

import numpy as np

from ..core import Mixture
from ..exceptions import CompositionError, ModelError
from ..models import ActivityModel, EquationOfState
from ..validation import validate_pressure, validate_temperature
from ._detect import _flash_tp_liquid_liquid, _flash_tp_tangent_plane
from ._legacy import _flash_tp_wilson_heuristic
from .results import FlashResult
from .settings import FlashSettings

#: Supported ``flash_mode`` values.
FLASH_MODES = ("phi-phi", "gamma-phi", "gamma-gamma")


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
