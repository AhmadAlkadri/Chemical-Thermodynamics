"""Flash calculation settings."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from ..exceptions import InputRangeError

if TYPE_CHECKING:  # pragma: no cover - import kept out of the runtime cycle
    from ..stability import StabilitySettings

#: Admissible values of :attr:`FlashSettings.phase_detection`.
PHASE_DETECTION_MODES = ("tangent-plane", "wilson-heuristic")


@dataclass(frozen=True)
class FlashSettings:
    """Solver settings for flash calculations.

    Attributes:
        max_iter: Maximum number of iterations.
        tol: Convergence tolerance on K-value updates (dimensionless).
        damping: Optional damping factor for K updates in (0, 1].
        phase_detection: How the phi-phi path decides one phase versus two.

            - ``"tangent-plane"`` (default): run Michelsen's tangent-plane
              stability test on the feed (:func:`chemthermo.stability_tp`),
              return a single phase only when the feed is found stable, and
              seed the two-phase split from the stability minimizer. An
              ``"inconclusive"`` stability result raises
              :class:`chemthermo.ConvergenceError` rather than silently
              producing a single-phase answer.
            - ``"wilson-heuristic"``: the legacy path, which declares a single
              phase when all Wilson K-values fall on one side of 1 or when
              Rachford-Rice finds no root for the Wilson K-values. Both are
              heuristics on an *initial estimate*, not thermodynamic criteria.
              Kept reachable so the old behavior stays testable (ADR-0008).

            Gamma-phi always uses ``"wilson-heuristic"`` in this release; see
            ADR-0007 and ADR-0008 for why gamma-phi stability is not available.
        stability_settings: Settings forwarded to
            :func:`chemthermo.stability_tp` when ``phase_detection`` is
            ``"tangent-plane"``, for the feed test and for the post-split test
            of each converged phase. ``None`` uses ``StabilitySettings()``.
        post_split_stability: Refuse to return a two-phase result whose phases
            are not themselves stable (ADR-0009). Every converged phase is fed
            back into :func:`chemthermo.stability_tp` and the outcome is always
            reported in ``diagnostics``; this flag decides what happens when
            that check *fails*. ``True`` (default) raises
            :class:`chemthermo.ConvergenceError` saying that a third phase is
            required; ``False`` returns the two-phase result anyway, with the
            failure visible in ``diagnostics["post_split_status"]``.

            The check runs on the tangent-plane phi-phi path and on the
            liquid-liquid (``"gamma-gamma"``) path. It cannot run for
            ``"gamma-phi"`` (there is no gamma-phi stability test, ADR-0007),
            and it deliberately does not run on the legacy
            ``phase_detection="wilson-heuristic"`` path, whose purpose is to
            reproduce pre-ADR-0008 behavior unchanged. Those two paths report
            ``diagnostics["post_split_checked"] = False`` and a
            ``post_split_skipped_reason``.
        second_order: Run a second-order stage after successive substitution in
            the **liquid-liquid** (``"gamma-gamma"``) split. The stage is a
            damped Newton minimization of the two-phase Gibbs energy whose
            gradient is the equal-activity residual (ADR-0009). Near a plait
            point successive substitution needs thousands of iterations, so the
            stage is what makes those feeds solvable at all. The phi-phi and
            gamma-phi splits are unchanged by this release and never enter it.
        ssi_iterations: Successive-substitution iterations performed in the
            liquid-liquid split before the second-order stage takes over.
            Capped by ``max_iter``.
        second_order_max_iter: Maximum second-order iterations in the
            liquid-liquid split.
        second_order_tol: Target for the second-order stage, measured on the
            equal-activity residual ``max_i |ln(x_i^I gamma_i^I)
            - ln(x_i^II gamma_i^II)|``. It is tighter than ``tol`` because the
            stage converges quadratically (one extra step is cheap) and because
            this residual is what a caller verifies in ``diagnostics``. The
            stage stops early when it can no longer make progress; a split is
            accepted as converged as soon as it meets ``tol``.

    Notes:
        The solver is deterministic for fixed inputs, models, and settings.
    """

    max_iter: int = 100
    tol: float = 1e-8
    damping: float | None = None
    phase_detection: str = "tangent-plane"
    stability_settings: StabilitySettings | None = None
    post_split_stability: bool = True
    second_order: bool = True
    ssi_iterations: int = 50
    second_order_max_iter: int = 100
    second_order_tol: float = 1e-12

    def __post_init__(self) -> None:
        if self.max_iter <= 0:
            raise InputRangeError("max_iter must be positive.")
        if self.tol <= 0.0:
            raise InputRangeError("tol must be positive.")
        if self.damping is not None:
            if not (0.0 < self.damping <= 1.0):
                raise InputRangeError("damping must be in (0, 1] when specified.")
        if self.phase_detection not in PHASE_DETECTION_MODES:
            raise InputRangeError(
                f"phase_detection must be one of {PHASE_DETECTION_MODES}; "
                f"got {self.phase_detection!r}."
            )
        if self.ssi_iterations <= 0:
            raise InputRangeError("ssi_iterations must be positive.")
        if self.second_order_max_iter <= 0:
            raise InputRangeError("second_order_max_iter must be positive.")
        if self.second_order_tol <= 0.0:
            raise InputRangeError("second_order_tol must be positive.")
