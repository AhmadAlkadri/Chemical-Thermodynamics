"""TP flash calculations (T,P) for VLE in SI units.

Supports phi-phi and gamma-phi. Given the same inputs and settings, results are
deterministic.

Phase detection (phi-phi)
-------------------------
The reference phi-phi path decides one phase versus two from Michelsen's
tangent-plane stability criterion rather than from Wilson K-value bounds
(ADR-0008). The flow is

    flash_tp -> stability_tp(feed) -> single phase | seeded two-phase split

1. ``stability_tp`` is run on the feed at the same ``(T, P)`` with the same
   equation of state.
2. ``status == "stable"``: a single-phase ``FlashResult`` is returned. Its phase
   name is the minimum-Gibbs root branch that the stability test selected for
   the feed (``feed_branch``). When the cubic has a single real root both
   branches coincide and the label is a *convention*, not a measurement; see the
   note on :func:`flash_tp`.
3. ``status == "unstable"``: the converged stationary point seeds the K-values
   (see :func:`_stability_k_seed`) and the existing successive-substitution /
   Rachford-Rice loop - Michelsen's recommended first-stage phase split - runs
   unchanged. If the seeded K-values give no Rachford-Rice root the Wilson
   estimate is tried as a documented fallback, and ``diagnostics["k_seed"]``
   records which seed was actually used.
4. ``status == "inconclusive"``: a :class:`chemthermo.ConvergenceError` is
   raised. A stability search that could not converge must not silently produce
   a single-phase answer.

The converged split is then verified (material balance, phase fractions, equal
fugacities, and a negative Gibbs-energy change against the single-phase feed)
and every residual is reported in ``diagnostics``.

``FlashSettings(phase_detection="wilson-heuristic")`` restores the legacy
K-bound / Rachford-Rice heuristic. Gamma-phi always uses it in this release.

Limits of this slice
--------------------
At most two phases are returned, and the converged phases are **not** themselves
re-tested for stability (phase addition/removal is the next slice). A negative
``tpd_min`` proves the feed is not one phase; ``"stable"`` only means no
negative tangent-plane distance was found from the deterministic trial set.
"""

from __future__ import annotations

import math

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


def flash_tp(
    mixture: Mixture,
    *,
    temperature_K: float,
    pressure_Pa: float,
    eos: EquationOfState | None,
    activity_model: ActivityModel | None = None,
    flash_mode: str = "phi-phi",
    settings: FlashSettings | None = None,
) -> FlashResult:
    """Perform a TP flash calculation using phi-phi or gamma-phi.

    Args:
        mixture: Mixture with mole-fraction composition.
        temperature_K: Temperature in K.
        pressure_Pa: Pressure in Pa.
        eos: Equation-of-state model used for fugacity coefficients.
        activity_model: Activity model (required for gamma-phi).
        flash_mode: Case-insensitive mode: "phi-phi" or "gamma-phi".
        settings: Iteration controls (tolerance, damping, max iterations,
            phase-detection mode, stability settings).

    Returns:
        FlashResult with phase compositions and fractions. Phase names follow:
        VLE -> "liquid"/"vapor", single-phase -> "liquid" or "vapor".

    Diagnostics:
        Diagnostics keys are implementation details. Current stable keys:

        - Always: ``flash_mode``, ``phase_detection``, ``iterations``,
          ``converged``, ``termination_reason``, ``phase_count``,
          ``phase_state``, ``phase_regime``.
        - Tangent-plane path (phi-phi default): ``stability_status``,
          ``tpd_min``, ``stability_trials``, and ``feed_branch`` when the EOS
          reports one. Single phase adds nothing else and uses
          ``termination_reason = "feed_stable_tangent_plane"``. A two-phase
          result adds ``k_seed`` ("stability"/"wilson"), ``incipient_phase``,
          ``max_delta_k``, ``k_min``, ``k_max``, ``mass_balance_residual``,
          ``fugacity_residual`` and ``delta_g_split_rt``.
        - Legacy heuristic path: ``k_min``, ``k_max``, ``max_delta_k``,
          ``k_seed`` ("wilson"), and, for its single-phase fallbacks,
          ``rr_f0``, ``rr_f1``, ``rr_status``.

    Raises:
        InputRangeError: If temperature or pressure is non-physical.
        ModelError: If required models are missing or return invalid values.
        CompositionError: If the mixture composition is invalid.
        ConvergenceError: If iteration fails to converge, or (tangent-plane
            mode only) if the stability analysis is inconclusive, or if an
            unstable feed admits no Rachford-Rice root from either seed.

    Notes:
        With ``phase_detection="tangent-plane"`` (the default for phi-phi) a
        single-phase result means the feed was *found stable* by Michelsen's
        test. With ``phase_detection="wilson-heuristic"`` it only means an
        initial-estimate heuristic said so.

        **Vapor/liquid labelling convention.** For a single-phase result the
        name is the minimum-Gibbs compressibility root branch of the feed. When
        the cubic has a single real root (dense or supercritical fluids) both
        branches return identical fugacity coefficients, the branch label is a
        tie-break, and the reported ``"vapor"`` / ``"liquid"`` name is therefore
        a naming convention rather than a phase identification. For a two-phase
        result the two converged phases are named by volatility ordering: the
        phase enriched (relative to the feed) in the component with the largest
        Wilson K relative to the one with the smallest is named ``"vapor"``.
        ``EquationOfState`` exposes no molar volume, so no density-based
        identification is available; this ordering decides the *name* only,
        never the verdict or the compositions.
    """

    temperature = validate_temperature(temperature_K)
    pressure = validate_pressure(pressure_Pa)

    if eos is None:
        raise ModelError("An equation-of-state model is required for flash_tp.")

    settings = settings or FlashSettings()
    mode = flash_mode.strip().casefold()

    if mode in {"phi-phi", "gamma-phi"}:
        if mode == "gamma-phi" and activity_model is None:
            raise ModelError("An activity model is required for gamma-phi flash.")
        if mode != "gamma-phi" and activity_model is not None:
            raise ModelError("activity_model is only used when flash_mode='gamma-phi'.")
        return _flash_tp_vle(
            mixture,
            temperature,
            pressure,
            eos=eos,
            activity_model=activity_model,
            mode=mode,
            settings=settings,
        )

    if mode == "vlle":
        raise ModelError(
            "VLLE support is provided by the optional chemthermo_vlle plugin. "
            "Install chemthermo_vlle to enable VLLE support."
        )

    raise ModelError(f"Unsupported flash_mode '{flash_mode}'.")


def _flash_tp_vle(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState,
    activity_model: ActivityModel | None,
    mode: str,
    settings: FlashSettings,
) -> FlashResult:
    """Dispatch to the tangent-plane path or the legacy heuristic path."""
    if mixture.basis != "mole":
        raise ModelError("flash_tp currently requires mole-fraction compositions.")

    z = np.array(mixture.fractions, dtype=float)
    if z.size == 0:
        raise CompositionError("Mixture composition must be non-empty.")

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
        vapor_fraction=split.vapor_fraction,
        ln_phi_liquid=np.log(split.phi_liquid),
        ln_phi_vapor=np.log(split.phi_vapor),
        ln_phi_feed=ln_phi_feed,
    )

    return _two_phase_result(
        mixture,
        temperature,
        pressure,
        split.x,
        split.y,
        split.vapor_fraction,
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

    Kept unchanged (ADR-0008) so the pre-``flash-auto-phase-detection`` behavior
    stays reachable and testable, and because gamma-phi has no stability test to
    fall back on (ADR-0007).
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
        },
    )


class _SplitSolution:
    """Converged two-phase split plus the quantities needed to verify it."""

    __slots__ = (
        "K",
        "iterations",
        "max_delta",
        "phi_liquid",
        "phi_vapor",
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
        phi_liquid: np.ndarray,
        phi_vapor: np.ndarray,
        iterations: int,
        max_delta: float,
    ) -> None:
        self.x = x
        self.y = y
        self.vapor_fraction = vapor_fraction
        self.K = K
        self.phi_liquid = phi_liquid
        self.phi_vapor = phi_vapor
        self.iterations = iterations
        self.max_delta = max_delta


def _solve_k_loop(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState,
    activity_model: ActivityModel | None,
    mode: str,
    settings: FlashSettings,
    z: np.ndarray,
    K: np.ndarray,
    vapor_fraction: float,
) -> _SplitSolution:
    """Successive substitution on K with Rachford-Rice updates of ``beta``.

    This is the unchanged first-stage phase split; only its initial ``K`` and
    ``vapor_fraction`` depend on how the caller detected two phases.
    """
    max_delta = float("inf")
    for iteration in range(1, settings.max_iter + 1):
        x = z / (1.0 + vapor_fraction * (K - 1.0))
        x = normalize_composition(x, label="liquid", error_cls=ConvergenceError)

        y = K * x
        y = normalize_composition(y, label="vapor", error_cls=ConvergenceError)

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

        max_delta = float(np.max(np.abs(K_new - K)))
        if max_delta < settings.tol:
            return _SplitSolution(
                x=x,
                y=y,
                vapor_fraction=vapor_fraction,
                K=K_new,
                phi_liquid=phi_l,
                phi_vapor=phi_v,
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

    raise ConvergenceError(
        f"flash_tp did not converge within the iteration limit; max_delta_k={max_delta:.3e}."
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
    vapor_fraction: float,
    ln_phi_liquid: np.ndarray,
    ln_phi_vapor: np.ndarray,
    ln_phi_feed: np.ndarray,
) -> dict[str, float | int | str | bool]:
    """Residuals that a converged two-phase solution must satisfy.

    - ``mass_balance_residual``: ``max_i |z_i - (beta y_i + (1 - beta) x_i)|``.
    - ``fugacity_residual``: ``max_i |ln(x_i phi_i^L) - ln(y_i phi_i^V)|``. The
      pressure cancels, so this is the equal-fugacity condition in log form.
    - ``delta_g_split_rt``: ``beta sum_i y_i ln(y_i phi_i^V)
      + (1 - beta) sum_i x_i ln(x_i phi_i^L) - sum_i z_i ln(z_i phi_i^feed)``,
      which must be negative for the split to be an improvement on the
      single-phase feed. The feed term uses the *minimum-Gibbs* branch (the
      stability module's convention, i.e. the lowest single-phase Gibbs energy
      available), while each split phase uses the branch the solver actually
      converged on. Since the min-Gibbs branch is never higher in Gibbs energy,
      this combination is a conservative (upper-bound) estimate of the true
      reduction.
    """
    beta = float(vapor_fraction)
    balance = float(np.max(np.abs(z - (beta * y + (1.0 - beta) * x))))

    both = (x > 0.0) & (y > 0.0)
    if np.any(both):
        fugacity = float(
            np.max(
                np.abs(
                    (np.log(x[both]) + ln_phi_liquid[both]) - (np.log(y[both]) + ln_phi_vapor[both])
                )
            )
        )
    else:
        fugacity = float("nan")

    def _reduced_g(fractions: np.ndarray, ln_phi: np.ndarray) -> float:
        mask = fractions > 0.0
        return float(np.sum(fractions[mask] * (np.log(fractions[mask]) + ln_phi[mask])))

    delta_g = (
        beta * _reduced_g(y, ln_phi_vapor)
        + (1.0 - beta) * _reduced_g(x, ln_phi_liquid)
        - _reduced_g(z, ln_phi_feed)
    )

    return {
        "mass_balance_residual": balance,
        "fugacity_residual": fugacity,
        "delta_g_split_rt": delta_g,
    }


def _single_phase_result(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    *,
    phase_name: str,
    vapor_fraction: float,
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
    vapor_fraction: float,
    diagnostics: dict[str, float | int | str | bool],
) -> FlashResult:
    liquid = PhaseResult(
        name="liquid",
        composition=Composition(
            fractions=tuple(x.tolist()),
            basis=mixture.basis,
            normalize=False,
            tol=COMPOSITION_SUM_TOL,
        ),
    )
    vapor = PhaseResult(
        name="vapor",
        composition=Composition(
            fractions=tuple(y.tolist()),
            basis=mixture.basis,
            normalize=False,
            tol=COMPOSITION_SUM_TOL,
        ),
    )
    phase_fractions = {"liquid": 1.0 - float(vapor_fraction), "vapor": float(vapor_fraction)}
    return FlashResult(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        phases={"liquid": liquid, "vapor": vapor},
        vapor_fraction=float(vapor_fraction),
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
