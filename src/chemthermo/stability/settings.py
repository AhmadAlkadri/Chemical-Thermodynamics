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
            ``max_i |ln W_i + ln phi_i(w) - d_i|`` (dimensionless). For an
            activity model ``ln gamma_i`` replaces ``ln phi_i``.
        trivial_tol: Trial is declared to have collapsed onto the trivial
            solution (``w -> z``) when ``sum_i (ln(W_i / z_i))**2`` drops below
            this value (dimensionless).
        tpd_tol: A converged trial is only reported as an instability when its
            reduced tangent-plane distance satisfies ``tpd < -tpd_tol``
            (dimensionless, in units of RT). The default 1e-8 is deliberately
            far below the smallest published tangent-plane minimum this package
            is validated against (|D| ~ 1e-5 at the near-plait-point feeds of
            Tessier et al. 2000, Table 2), so those instabilities are not
            rounded away.
        second_order: Run a second-order (Newton) stage on the stationarity
            condition when successive substitution has used its
            ``ssi_iterations`` budget without meeting ``tol``. Near a plait
            point successive substitution converges linearly with a ratio close
            to one, or drifts to the trivial solution, so the second stage is
            what makes those cases usable.
        ssi_iterations: Successive-substitution iterations performed before the
            second-order stage takes over. Ignored when ``second_order`` is
            False, in which case successive substitution runs to ``max_iter``.
            The default 50 is larger than any iteration count reached by the
            validated Peng-Robinson cases (17), so enabling the second stage by
            default leaves those results bit-identical.
        second_order_max_iter: Maximum second-order iterations per trial.
        second_order_max_step: Largest allowed ``|delta ln W_i|`` in one
            second-order step, before the backtracking line search. Keeps the
            iterate in a range where ``W > 0`` and the model stays evaluable.

    Notes:
        The analysis is deterministic for fixed inputs, models, and settings.
    """

    max_iter: int = 300
    tol: float = 1e-10
    trivial_tol: float = 1e-8
    tpd_tol: float = 1e-8
    second_order: bool = True
    ssi_iterations: int = 50
    second_order_max_iter: int = 100
    second_order_max_step: float = 4.0

    def __post_init__(self) -> None:
        if self.max_iter <= 0:
            raise InputRangeError("max_iter must be positive.")
        if self.tol <= 0.0:
            raise InputRangeError("tol must be positive.")
        if self.trivial_tol <= 0.0:
            raise InputRangeError("trivial_tol must be positive.")
        if self.tpd_tol <= 0.0:
            raise InputRangeError("tpd_tol must be positive.")
        if self.ssi_iterations <= 0:
            raise InputRangeError("ssi_iterations must be positive.")
        if self.second_order_max_iter <= 0:
            raise InputRangeError("second_order_max_iter must be positive.")
        if self.second_order_max_step <= 0.0:
            raise InputRangeError("second_order_max_step must be positive.")
