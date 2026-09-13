"""Settings for tangent-plane phase-stability analysis."""

from __future__ import annotations

from dataclasses import dataclass

from ..exceptions import InputRangeError


@dataclass(frozen=True)
class StabilitySettings:
    """Solver settings for :func:`chemthermo.stability_tp`.

    Attributes:
        max_iter: Maximum successive-substitution iterations per trial.
        tol: Convergence tolerance on the stationarity residual
            ``max_i |ln W_i + ln phi_i(w) - d_i|`` (dimensionless).
        trivial_tol: Trial is declared to have collapsed onto the trivial
            solution (``w -> z``) when ``sum_i (ln(W_i / z_i))**2`` drops below
            this value (dimensionless).
        tpd_tol: A converged trial is only reported as an instability when its
            reduced tangent-plane distance satisfies ``tpd < -tpd_tol``
            (dimensionless, in units of RT).

    Notes:
        The analysis is deterministic for fixed inputs, models, and settings.
    """

    max_iter: int = 300
    tol: float = 1e-10
    trivial_tol: float = 1e-8
    tpd_tol: float = 1e-8

    def __post_init__(self) -> None:
        if self.max_iter <= 0:
            raise InputRangeError("max_iter must be positive.")
        if self.tol <= 0.0:
            raise InputRangeError("tol must be positive.")
        if self.trivial_tol <= 0.0:
            raise InputRangeError("trivial_tol must be positive.")
        if self.tpd_tol <= 0.0:
            raise InputRangeError("tpd_tol must be positive.")
