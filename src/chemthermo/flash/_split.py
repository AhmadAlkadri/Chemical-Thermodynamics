"""The shared two-phase split loop: Rachford-Rice plus successive substitution.

This is the first-stage phase split used, unchanged, by all three
``flash_tp`` modes (phi-phi, gamma-phi, gamma-gamma): only the initial ``K``,
the initial ``beta`` and the model call that updates ``K`` differ between
them (see :func:`_solve_k_loop`). Callers seed and interpret the loop
differently (:mod:`chemthermo.flash._detect`, :mod:`chemthermo.flash._legacy`);
this module owns none of that, only the shared numerics.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np

from ..core import Mixture
from ..exceptions import ConvergenceError, ModelError
from ..models import ActivityModel, EquationOfState
from ._common import as_float_array, normalize_composition
from .settings import FlashSettings


class _SplitSolution:
    """Converged two-phase split plus the quantities needed to verify it.

    ``ln_f_x`` / ``ln_f_y`` are the tangent-plane fugacity terms of the two
    phases at the returned compositions: ``ln phi`` for phi-phi and gamma-phi,
    ``ln gamma`` for gamma-gamma.
    """

    __slots__ = (
        "K",
        "converged",
        "iterations",
        "ln_f_x",
        "ln_f_y",
        "max_delta",
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
        ln_f_x: np.ndarray,
        ln_f_y: np.ndarray,
        iterations: int,
        max_delta: float,
        converged: bool = True,
    ) -> None:
        self.x = x
        self.y = y
        self.vapor_fraction = vapor_fraction
        self.K = K
        self.ln_f_x = ln_f_x
        self.ln_f_y = ln_f_y
        self.iterations = iterations
        self.max_delta = max_delta
        self.converged = converged


def _solve_k_loop(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState | None,
    activity_model: ActivityModel | None,
    mode: str,
    settings: FlashSettings,
    z: np.ndarray,
    K: np.ndarray,
    vapor_fraction: float,
    max_iter: int | None = None,
    allow_unconverged: bool = False,
) -> _SplitSolution:
    """Successive substitution on K with Rachford-Rice updates of ``beta``.

    This is the first-stage phase split, shared by all three modes; only its
    initial ``K``, its ``beta`` and the model call that updates ``K`` differ:

    - phi-phi: ``K = phi^L / phi^V`` (``x`` is the liquid, ``y`` the vapor);
    - gamma-phi: ``K = gamma^L phi^L / phi^V``;
    - gamma-gamma: ``K = gamma^I / gamma^II`` (``x`` is phase I, ``y`` phase II
      and ``vapor_fraction`` is the mole fraction of phase II).

    ``allow_unconverged`` returns the last iterate instead of raising when the
    budget runs out, which is how the liquid-liquid path hands over to its
    second-order stage.
    """
    budget = settings.max_iter if max_iter is None else max_iter
    max_delta = float("inf")
    x = np.array(z, dtype=float)
    y = np.array(z, dtype=float)
    ln_f_x = np.zeros_like(z)
    ln_f_y = np.zeros_like(z)

    for iteration in range(1, budget + 1):
        x = z / (1.0 + vapor_fraction * (K - 1.0))
        x = normalize_composition(x, label="liquid", error_cls=ConvergenceError)

        y = K * x
        y = normalize_composition(y, label="vapor", error_cls=ConvergenceError)

        if mode == "gamma-gamma":
            assert activity_model is not None
            gamma_x = _activity_coefficients(activity_model, mixture, temperature, x)
            gamma_y = _activity_coefficients(activity_model, mixture, temperature, y)
            if gamma_x.shape != K.shape or gamma_y.shape != K.shape:
                raise ModelError("Activity model returned inconsistent coefficient shapes.")
            K_new = gamma_x / gamma_y
            ln_f_x = np.log(gamma_x)
            ln_f_y = np.log(gamma_y)
        else:
            assert eos is not None
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
            ln_f_x = np.log(phi_l)
            ln_f_y = np.log(phi_v)

        max_delta = float(np.max(np.abs(K_new - K)))
        if max_delta < settings.tol:
            return _SplitSolution(
                x=x,
                y=y,
                vapor_fraction=vapor_fraction,
                K=K_new,
                ln_f_x=ln_f_x,
                ln_f_y=ln_f_y,
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

    if allow_unconverged:
        return _SplitSolution(
            x=x,
            y=y,
            vapor_fraction=vapor_fraction,
            K=K,
            ln_f_x=ln_f_x,
            ln_f_y=ln_f_y,
            iterations=budget,
            max_delta=max_delta,
            converged=False,
        )

    raise ConvergenceError(
        f"flash_tp did not converge within the iteration limit; max_delta_k={max_delta:.3e}."
    )


def _activity_coefficients(
    activity_model: ActivityModel,
    mixture: Mixture,
    temperature: float,
    composition: np.ndarray,
) -> np.ndarray:
    """Activity coefficients at ``composition``, validated."""
    values = as_float_array(
        activity_model.activity_coefficients(
            mixture=mixture,
            temperature_K=temperature,
            composition=composition.tolist(),
        )
    )
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ModelError("Activity model returned non-positive activity coefficients.")
    return values


def _ln_gamma_function(
    activity_model: ActivityModel, mixture: Mixture, temperature: float
) -> Callable[[np.ndarray], np.ndarray]:
    """Return ``ln gamma(x)`` as a plain callable on normalized compositions."""

    def ln_gamma(composition: np.ndarray) -> np.ndarray:
        values = np.asarray(composition, dtype=float)
        total = float(np.sum(values))
        if total <= 0.0:
            raise ModelError("Activity model called with a non-positive composition.")
        return np.log(_activity_coefficients(activity_model, mixture, temperature, values / total))

    return ln_gamma


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
