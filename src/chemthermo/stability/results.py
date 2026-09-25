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
            final point. At a stationary point ``tpd = -ln(sum_W)``. It may be
            ``inf`` (or an exact ``0.0``) when the stationary point's mole
            numbers leave the exponential's range - a polymer melt against a
            solvent-vapour feed has ``ln W`` of order ``1450`` - in which case
            ``ln_sum_W`` is the one to read. See ADR-0025.
        ln_W: ``ln W_i``, the *unnormalized* log mole numbers at the final
            point, ``-inf`` for a component absent from the feed. This is the
            quantity the iteration actually carries, and it stays exact where
            ``W_i`` and the normalized ``w_i`` do not: for the melt above,
            ``ln W = (1452.2, 3.6)`` while ``w`` rounds to ``(1.0, 0.0)``.
            None for a trial that failed before taking a step.
        ln_sum_W: ``ln sum_i W_i`` at the final point. At a stationary point
            this is ``-tpd`` (equation (7)) and it is finite for any magnitude
            of ``ln W``. ``nan`` for a failed trial.
        log_space: True when this trial's normalization had to be done in logs
            at least once, i.e. when some ``ln W_i`` left ``[-700, 700]``.
            False means the trial ran the pre-ADR-0025 arithmetic character
            for character, so its numbers are bit-identical to what it
            returned before that decision record.
        trivial: True when the trial collapsed onto the feed composition.
        residual: Final stationarity residual ``max_i |ln W_i + ln phi_i(w) - d_i|``
            (``ln gamma_i`` for an activity model).
        phase_branch: Label of the lowest-Gibbs phase candidate selected at the
            final trial composition: the compressibility root for an EOS, or
            ``"liquid"`` / ``"vapor"`` for the modified-Raoult candidate pair
            (``vapor="ideal"``). None for a single-candidate evaluator (an
            activity-coefficient model on its own) or if unavailable.
        surface: Label of the phase candidate this trial *iterated on*, held
            fixed for every successive-substitution and Newton step (ADR-0012
            for the modified-Raoult pair, ADR-0021 for the density roots of an
            equation of state). None when the trial re-selected the
            lowest-Gibbs candidate at every iterate: every trial of an
            activity-coefficient model on its own, and the degenerate
            single-active-component trials. It normally equals
            ``phase_branch``; a difference means the trial converged to a
            stationary point of its own surface at a composition where the
            *other* candidate has the lower Gibbs energy, and ``tpd`` is then
            the (smaller) lowest-Gibbs value, not the surface's own.
        surface_fallback: True when the trial could not run on ``surface``
            alone. Two things set it, and both mean "there was no choice of
            surface here": the named candidate raised at an iterate (ADR-0012
            decision 5 - never reached through the modified-Raoult pair, whose
            two candidates are evaluable at every composition), and, for an
            equation of state, the model having a **single** admissible density
            root, which both ``phase`` labels then name (ADR-0021). The second
            is measured where the solver compares the branches anyway - the
            trial's stopping point - and not at every iterate, so this flag
            answers "did this trial *stop* in a one-root region", not "did it
            ever meet one". ADR-0021 decision 4 says why, and records how often
            the condition holds over the validation grids (about nine
            evaluations in ten).
        surface_fallback_count: How many model evaluations inside this trial set
            it (0 when ``surface_fallback`` is False). Evaluations, not
            iterations.
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
    surface_fallback_count: int = 0
    # ADR-0025. Declared last, with defaults, so the positional order of every
    # pre-ADR-0025 field is unchanged.
    ln_W: tuple[float, ...] | None = None
    ln_sum_W: float = float("nan")
    log_space: bool = False

    def __post_init__(self) -> None:
        if not self.label.strip():
            raise ValueError("Trial label must be non-empty.")
        if self.iterations < 0:
            raise InputRangeError("Trial iterations must be non-negative.")
        if self.ssi_iterations < 0 or self.second_order_iterations < 0:
            raise InputRangeError("Trial stage iteration counts must be non-negative.")
        if self.surface_fallback_count < 0:
            raise InputRangeError("Trial surface-fallback count must be non-negative.")
        if self.surface_fallback != (self.surface_fallback_count > 0):
            raise ValueError("surface_fallback must be True if and only if its count is positive.")
        if self.ln_W is not None and self.composition is not None:
            if len(self.ln_W) != len(self.composition):
                raise ValueError("ln_W length must match composition length.")


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
            None when no non-trivial stationary point was found. A component
            whose share of the stationary point is below the smallest positive
            double is an exact ``0.0`` here; ``trial_ln_W`` keeps its
            magnitude (ADR-0025).
        trial_ln_W: The minimizing trial's *unnormalized* ``ln W`` - the
            quantity Michelsen's iteration carries and the one a log-space
            split stage wants as a seed. None when there is no minimizing
            trial, or when the trial did not record one.
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
            appearance, ``surface_fallback_trial_count`` counts the trials that
            had to leave their surface at least once,
            ``surface_fallback_evaluation_count`` sums those evaluations, and
            ``minimizing_trial_surface`` names the surface the minimizing trial
            ran on. ``tpd_from_sum_W`` is then equation (7) on
            *that* surface: it equals ``tpd_min`` whenever the minimizing trial
            stopped where its own candidate is the lowest-Gibbs one, which is
            every unstable verdict measured so far. ``ln_sum_W`` is the same
            quantity read from the logarithm rather than from ``sum_W``, so it
            survives a stationary point outside the exponential's range, and
            ``tm_at_stationary_point`` is then ``-inf`` rather than a number
            (ADR-0025). ``minimizing_trial_tie_break`` (``"residual"``)
            appears only when the lowest-``tpd`` trial was passed over for a
            tied one converged at least 1000x better (ADR-0035).
            ``log_space_trial_count`` appears only when at least
            one trial had to normalize in logs; its absence is the statement
            that every trial ran the pre-ADR-0025 arithmetic.

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
    trial_ln_W: tuple[float, ...] | None = None

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
        if self.trial_ln_W is not None and len(self.trial_ln_W) != len(self.feed_composition):
            raise ValueError("trial_ln_W length must match feed_composition length.")
