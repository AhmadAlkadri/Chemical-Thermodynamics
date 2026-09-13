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
            ``"tangent-plane"``. ``None`` uses ``StabilitySettings()``.

    Notes:
        The solver is deterministic for fixed inputs, models, and settings.
    """

    max_iter: int = 100
    tol: float = 1e-8
    damping: float | None = None
    phase_detection: str = "tangent-plane"
    stability_settings: StabilitySettings | None = None

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
