"""Result containers for tangent-plane phase-stability analysis."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

from ..exceptions import InputRangeError
from ..validation import validate_pressure, validate_temperature

STABILITY_STATUSES = ("stable", "unstable", "inconclusive")


@dataclass(frozen=True)
class StabilityTrial:
    """Outcome of a single trial-phase stationary-point search.

    Attributes:
        label: Deterministic identifier of the initial estimate, for example
            ``"wilson-vapor"``, ``"wilson-liquid"`` or ``"pure-Methane"``.
        converged: True when the stationarity residual met ``settings.tol`` or
            the iteration collapsed onto the trivial solution.
        iterations: Total iterations performed, ``ssi_iterations +
            second_order_iterations``.
        tpd: Reduced tangent-plane distance at the final trial composition
            (dimensionless, units of RT). ``nan`` when the trial failed.
        sum_W: Sum of the unnormalized trial mole numbers ``sum_i W_i`` at the
            final point. At a stationary point ``tpd = -ln(sum_W)``.
        trivial: True when the trial collapsed onto the feed composition.
        residual: Final stationarity residual ``max_i |ln W_i + ln phi_i(w) - d_i|``
            (``ln gamma_i`` for an activity model).
        phase_branch: Label of the lowest-Gibbs phase candidate selected at the
            final trial composition: the compressibility root for an EOS, or
            ``"liquid"`` / ``"vapor"`` for the modified-Raoult candidate pair
            (``vapor="ideal"``). None for a single-candidate evaluator (an
            activity-coefficient model on its own) or if unavailable.
        surface: Label of the phase candidate this trial *iterated on*, held
            fixed for every successive-substitution and Newton step (ADR-0012).
            None when the trial re-selected the lowest-Gibbs candidate at every
            iterate, which is every trial of an equation of state and of an
            activity-coefficient model on its own. It normally equals
            ``phase_branch``; a difference means the trial converged to a
            stationary point of its own surface at a composition where the
            *other* candidate has the lower Gibbs energy, and ``tpd`` is then
            the (smaller) lowest-Gibbs value, not the surface's own.
        surface_fallback: True when ``surface`` was named but the candidate was
            not evaluable at some iterate and the lowest-Gibbs candidate had to
            be used there instead. Always False for the modified-Raoult pair,
            whose two candidates are both evaluable at every composition.
        composition: Normalized trial composition ``w`` at the final point.
        termination_reason: Short machine-readable reason string.
        ssi_iterations: Successive-substitution iterations performed.
        second_order_iterations: Second-order (Newton) iterations performed;
            0 when the trial finished during successive substitution or when
            the second-order stage is disabled.
        converged_stage: ``"successive-substitution"`` or ``"second-order"``
            for a converged trial, None otherwise.
    """

    label: str
    converged: bool
    iterations: int
    tpd: float
    sum_W: float
    trivial: bool
    residual: float
    phase_branch: str | None
    composition: tuple[float, ...] | None
    termination_reason: str
    ssi_iterations: int = 0
    second_order_iterations: int = 0
    converged_stage: str | None = None
    surface: str | None = None
    surface_fallback: bool = False

    def __post_init__(self) -> None:
        if not self.label.strip():
            raise ValueError("Trial label must be non-empty.")
        if self.iterations < 0:
            raise InputRangeError("Trial iterations must be non-negative.")
        if self.ssi_iterations < 0 or self.second_order_iterations < 0:
            raise InputRangeError("Trial stage iteration counts must be non-negative.")


@dataclass(frozen=True)
class StabilityResult:
    """Outcome of a tangent-plane stability analysis at fixed T, P and z.

    Attributes:
        temperature_K: Temperature in K.
        pressure_Pa: Pressure in Pa.
        feed_composition: Normalized feed mole fractions ``z``.
        stable: True when ``status == "stable"``.
        status: One of ``"stable"``, ``"unstable"`` or ``"inconclusive"``.
        tpd_min: Smallest reduced tangent-plane distance found over the
            converged non-trivial trials (dimensionless, units of RT). ``0.0``
            when every converged trial collapsed onto the trivial solution.
        trial_composition: Normalized minimizing trial composition ``w``, or
            None when no non-trivial stationary point was found.
        k_values: Implied K-values of the incipient phase, ``w_i / z_i``.
            Values > 1 mean the incipient phase is enriched in component ``i``
            relative to the feed; the incipient phase is the *new* phase, so
            ``k_values`` maps feed -> incipient, never the reverse.
        phase_branch: Label of the lowest-Gibbs phase candidate at the
            minimizing trial composition - the *incipient* phase's identity.
            None for a single-candidate evaluator.
        feed_branch: Label of the lowest-Gibbs phase candidate at the feed: the
            compressibility root for an EOS, ``"liquid"`` / ``"vapor"`` for the
            modified-Raoult pair. None for an activity-coefficient model on its
            own (there is a single candidate, so no selection is performed).
        trials: Per-trial records in deterministic order.
        diagnostics: Diagnostic metadata (implementation detail keys). When any
            trial was pinned to a phase-candidate surface (ADR-0012),
            ``trial_surfaces`` holds the per-surface trial counts as a
            deterministic ``"<label>:<count>"`` string in order of first
            appearance, and ``minimizing_trial_surface`` names the surface the
            minimizing trial ran on.

    Honesty note:
        ``stable`` means "no negative tangent-plane distance was found from the
        deterministic trial set used here". It is not a global proof of
        stability: Michelsen's test is a local stationary-point search and a
        missed stationary point can hide an instability.
    """

    temperature_K: float
    pressure_Pa: float
    feed_composition: tuple[float, ...]
    stable: bool
    status: str
    tpd_min: float
    trial_composition: tuple[float, ...] | None = None
    k_values: tuple[float, ...] | None = None
    phase_branch: str | None = None
    feed_branch: str | None = None
    trials: tuple[StabilityTrial, ...] = ()
    diagnostics: Mapping[str, float | int | str | bool] = field(default_factory=dict)

    def __post_init__(self) -> None:
        validate_temperature(self.temperature_K)
        validate_pressure(self.pressure_Pa)

        if self.status not in STABILITY_STATUSES:
            raise ValueError(f"status must be one of {STABILITY_STATUSES}; got {self.status!r}.")
        if self.stable != (self.status == "stable"):
            raise ValueError("stable must be True if and only if status == 'stable'.")
        if not self.feed_composition:
            raise ValueError("feed_composition must be non-empty.")
        if self.trial_composition is not None and len(self.trial_composition) != len(
            self.feed_composition
        ):
            raise ValueError("trial_composition length must match feed_composition length.")
        if self.k_values is not None and len(self.k_values) != len(self.feed_composition):
            raise ValueError("k_values length must match feed_composition length.")
