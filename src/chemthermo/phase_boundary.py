"""VLE calculations (Bubble Point, Dew Point) for diagrams."""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import TypedDict

import numpy as np

from .core import Mixture
from .exceptions import ConvergenceError, ModelError
from .flash._common import as_float_array, normalize_composition, wilson_k
from .models import ActivityModel, EquationOfState
from .validation import COMPOSITION_SUM_TOL, validate_pressure, validate_temperature

__all__ = [
    "bubble_temperature",
    "dew_temperature",
    "bubble_pressure",
    "dew_pressure",
]


class _VLESolverSettings(TypedDict):
    tol: float
    max_iter: int
    inner_max_iter: int
    inner_tol: float


def bubble_temperature(
    mixture: Mixture,
    pressure_Pa: float,
    eos: EquationOfState,
    activity_model: ActivityModel | None = None,
    settings: dict | None = None,
) -> tuple[float, list[float]]:
    """Calculate bubble-point temperature and equilibrium vapor composition."""
    pressure = validate_pressure(pressure_Pa)
    x = _require_mole_mixture(mixture)
    cfg = _parse_settings(settings)

    temperature_guess = 0.7 * float(np.dot(x, [component.tc_k for component in mixture.components]))
    temperature_guess = max(temperature_guess, 200.0)

    def evaluate(temperature: float) -> tuple[float, np.ndarray]:
        temp = validate_temperature(temperature)
        _, y, residual = _solve_k_bubble(
            mixture,
            temp,
            pressure,
            x,
            eos,
            activity_model,
            inner_tol=cfg["inner_tol"],
            inner_max_iter=cfg["inner_max_iter"],
        )
        return residual, y

    temperature, y = _solve_scalar_residual(
        temperature_guess,
        evaluate,
        tol=cfg["tol"],
        max_iter=cfg["max_iter"],
        min_value=1.0,
        label="bubble_temperature",
    )
    return float(temperature), y.tolist()


def dew_temperature(
    mixture: Mixture,
    pressure_Pa: float,
    eos: EquationOfState,
    activity_model: ActivityModel | None = None,
    settings: dict | None = None,
) -> tuple[float, list[float]]:
    """Calculate dew-point temperature and equilibrium liquid composition."""
    pressure = validate_pressure(pressure_Pa)
    y = _require_mole_mixture(mixture)
    cfg = _parse_settings(settings)

    temperature_guess = 0.7 * float(np.dot(y, [component.tc_k for component in mixture.components]))
    temperature_guess = max(temperature_guess, 200.0)

    def evaluate(temperature: float) -> tuple[float, np.ndarray]:
        temp = validate_temperature(temperature)
        _, x, residual = _solve_k_dew(
            mixture,
            temp,
            pressure,
            y,
            eos,
            activity_model,
            inner_tol=cfg["inner_tol"],
            inner_max_iter=cfg["inner_max_iter"],
        )
        return residual, x

    temperature, x = _solve_scalar_residual(
        temperature_guess,
        evaluate,
        tol=cfg["tol"],
        max_iter=cfg["max_iter"],
        min_value=1.0,
        label="dew_temperature",
    )
    return float(temperature), x.tolist()


def bubble_pressure(
    mixture: Mixture,
    temperature_K: float,
    eos: EquationOfState,
    activity_model: ActivityModel | None = None,
    settings: dict | None = None,
) -> tuple[float, list[float]]:
    """Calculate bubble-point pressure and equilibrium vapor composition."""
    temperature = validate_temperature(temperature_K)
    x = _require_mole_mixture(mixture)
    cfg = _parse_settings(settings)

    pressure_guess = 0.0
    for component, x_i in zip(mixture.components, x):
        pressure_guess += (
            x_i
            * component.pc_pa
            * math.exp(5.373 * (1.0 + component.omega) * (1.0 - component.tc_k / temperature))
        )
    pressure_guess = max(pressure_guess, 1.0e5)

    def evaluate(pressure: float) -> tuple[float, np.ndarray]:
        p = validate_pressure(pressure)
        _, y, residual = _solve_k_bubble(
            mixture,
            temperature,
            p,
            x,
            eos,
            activity_model,
            inner_tol=cfg["inner_tol"],
            inner_max_iter=cfg["inner_max_iter"],
        )
        return residual, y

    pressure, y = _solve_scalar_residual(
        pressure_guess,
        evaluate,
        tol=cfg["tol"],
        max_iter=cfg["max_iter"],
        min_value=1.0,
        label="bubble_pressure",
    )
    return float(pressure), y.tolist()


def dew_pressure(
    mixture: Mixture,
    temperature_K: float,
    eos: EquationOfState,
    activity_model: ActivityModel | None = None,
    settings: dict | None = None,
) -> tuple[float, list[float]]:
    """Calculate dew-point pressure and equilibrium liquid composition."""
    temperature = validate_temperature(temperature_K)
    y = _require_mole_mixture(mixture)
    cfg = _parse_settings(settings)

    inverse_pressure_guess = 0.0
    for component, y_i in zip(mixture.components, y):
        pseudo_psat = component.pc_pa * math.exp(
            5.373 * (1.0 + component.omega) * (1.0 - component.tc_k / temperature)
        )
        inverse_pressure_guess += y_i / pseudo_psat
    pressure_guess = 1.0 / inverse_pressure_guess if inverse_pressure_guess > 0.0 else 1.0e5
    pressure_guess = max(pressure_guess, 1.0e5)

    def evaluate(pressure: float) -> tuple[float, np.ndarray]:
        p = validate_pressure(pressure)
        _, x, residual = _solve_k_dew(
            mixture,
            temperature,
            p,
            y,
            eos,
            activity_model,
            inner_tol=cfg["inner_tol"],
            inner_max_iter=cfg["inner_max_iter"],
        )
        return residual, x

    pressure, x = _solve_scalar_residual(
        pressure_guess,
        evaluate,
        tol=cfg["tol"],
        max_iter=cfg["max_iter"],
        min_value=1.0,
        label="dew_pressure",
    )
    return float(pressure), x.tolist()


def _parse_settings(settings: dict | None) -> _VLESolverSettings:
    user_settings = settings or {}

    tol = float(user_settings.get("tol", 1e-8))
    max_iter = int(user_settings.get("max_iter", 100))
    inner_max_iter = int(user_settings.get("inner_max_iter", 60))
    inner_tol = float(user_settings.get("inner_tol", min(1e-8, tol)))

    if tol <= 0.0:
        raise ModelError("VLE tolerance must be positive.")
    if max_iter <= 0:
        raise ModelError("VLE max_iter must be positive.")
    if inner_max_iter <= 0:
        raise ModelError("VLE inner_max_iter must be positive.")
    if inner_tol <= 0.0:
        raise ModelError("VLE inner_tol must be positive.")

    return {
        "tol": tol,
        "max_iter": max_iter,
        "inner_max_iter": inner_max_iter,
        "inner_tol": inner_tol,
    }


def _require_mole_mixture(mixture: Mixture) -> np.ndarray:
    if mixture.basis != "mole":
        raise ModelError("VLE boundary calculations require mole-fraction compositions.")

    fractions = np.array(mixture.fractions, dtype=float)
    if fractions.size == 0:
        raise ModelError("Mixture composition must be non-empty.")
    if np.any(~np.isfinite(fractions)):
        raise ModelError("Mixture composition contains non-finite values.")
    if np.any(fractions < 0.0):
        raise ModelError("Mixture composition contains negative fractions.")

    total = float(np.sum(fractions))
    if abs(total - 1.0) > COMPOSITION_SUM_TOL:
        raise ModelError(
            f"Mixture composition must sum to 1 within {COMPOSITION_SUM_TOL:.1e}; got {total:.8f}."
        )

    return fractions


def _solve_scalar_residual(
    initial_value: float,
    evaluate: Callable[[float], tuple[float, np.ndarray]],
    *,
    tol: float,
    max_iter: int,
    min_value: float,
    label: str,
) -> tuple[float, np.ndarray]:
    value = max(float(initial_value), min_value)
    residual, composition = evaluate(value)

    if abs(residual) <= tol and _is_valid_composition(composition):
        return value, composition

    previous_value = max(value * 1.01, min_value * 1.01)
    previous_residual, _ = _evaluate_with_backtracking(
        previous_value,
        value,
        evaluate,
        min_value=min_value,
        label=label,
    )

    for _ in range(2, max_iter + 1):
        denominator = residual - previous_residual
        if not math.isfinite(denominator) or abs(denominator) < 1e-14:
            # Fallback step opposite the residual sign.
            step = -math.copysign(0.1 * max(value, min_value), residual)
            candidate = value + step
        else:
            candidate = value - residual * (value - previous_value) / denominator

        # Clamp scalar step to keep updates local and positive.
        max_step = 0.2 * max(value, min_value)
        lower = max(min_value, value - max_step)
        upper = value + max_step

        if not math.isfinite(candidate):
            candidate = 0.5 * (lower + upper)
        candidate = min(max(candidate, lower), upper)

        candidate_residual, candidate_composition = _evaluate_with_backtracking(
            candidate,
            value,
            evaluate,
            min_value=min_value,
            label=label,
        )

        previous_value, previous_residual = value, residual
        value, residual = candidate, candidate_residual
        composition = candidate_composition

        if abs(residual) <= tol and _is_valid_composition(composition):
            return value, composition

    raise ConvergenceError(
        f"{label} did not converge within {max_iter} iterations; residual={residual:.3e}."
    )


def _evaluate_with_backtracking(
    candidate: float,
    current: float,
    evaluate: Callable[[float], tuple[float, np.ndarray]],
    *,
    min_value: float,
    label: str,
) -> tuple[float, np.ndarray]:
    probe = max(candidate, min_value)
    for _ in range(8):
        try:
            return evaluate(probe)
        except ModelError:
            probe = max(min_value, 0.5 * (probe + current))

    raise ConvergenceError(f"{label} failed to evaluate stable model states while iterating.")


def _solve_k_bubble(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    x: np.ndarray,
    eos: EquationOfState,
    activity_model: ActivityModel | None,
    *,
    inner_tol: float,
    inner_max_iter: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    K = wilson_k(mixture, temperature_K, pressure_Pa)

    for _ in range(inner_max_iter):
        y = normalize_composition(K * x, label="vapor", error_cls=ConvergenceError)
        K_new = _model_k(
            mixture,
            temperature_K,
            pressure_Pa,
            liquid_composition=x,
            vapor_composition=y,
            eos=eos,
            activity_model=activity_model,
        )

        max_delta = float(np.max(np.abs(K_new - K)))
        K = K_new
        if max_delta <= inner_tol:
            y = normalize_composition(K * x, label="vapor", error_cls=ConvergenceError)
            residual = float(np.sum(K * x) - 1.0)
            return K, y, residual

    raise ConvergenceError(
        f"Bubble-point K-iteration failed to converge within {inner_max_iter} iterations."
    )


def _solve_k_dew(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    y: np.ndarray,
    eos: EquationOfState,
    activity_model: ActivityModel | None,
    *,
    inner_tol: float,
    inner_max_iter: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    K = wilson_k(mixture, temperature_K, pressure_Pa)

    for _ in range(inner_max_iter):
        x = normalize_composition(y / K, label="liquid", error_cls=ConvergenceError)
        K_new = _model_k(
            mixture,
            temperature_K,
            pressure_Pa,
            liquid_composition=x,
            vapor_composition=y,
            eos=eos,
            activity_model=activity_model,
        )

        max_delta = float(np.max(np.abs(K_new - K)))
        K = K_new
        if max_delta <= inner_tol:
            x = normalize_composition(y / K, label="liquid", error_cls=ConvergenceError)
            residual = float(np.sum(y / K) - 1.0)
            return K, x, residual

    raise ConvergenceError(
        f"Dew-point K-iteration failed to converge within {inner_max_iter} iterations."
    )


def _model_k(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    *,
    liquid_composition: np.ndarray,
    vapor_composition: np.ndarray,
    eos: EquationOfState,
    activity_model: ActivityModel | None,
) -> np.ndarray:
    try:
        phi_v = as_float_array(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                composition=vapor_composition.tolist(),
                phase="vapor",
            )
        )
        phi_l = as_float_array(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                composition=liquid_composition.tolist(),
                phase="liquid",
            )
        )
    except Exception as exc:
        raise ModelError(
            "EOS fugacity coefficient evaluation failed "
            f"at T={temperature_K:.6g} K, P={pressure_Pa:.6g} Pa."
        ) from exc

    if phi_v.shape != phi_l.shape:
        raise ModelError("EOS returned mismatched fugacity coefficient shapes.")
    if np.any(~np.isfinite(phi_v)) or np.any(~np.isfinite(phi_l)):
        raise ModelError("EOS returned non-finite fugacity coefficients.")
    if np.any(phi_v <= 0.0) or np.any(phi_l <= 0.0):
        raise ModelError("EOS returned non-positive fugacity coefficients.")

    if activity_model is None:
        gamma_l = np.ones_like(phi_l)
    else:
        try:
            gamma_l = as_float_array(
                activity_model.activity_coefficients(
                    mixture=mixture,
                    temperature_K=temperature_K,
                    composition=liquid_composition.tolist(),
                )
            )
        except Exception as exc:
            raise ModelError(
                f"Activity-model coefficient evaluation failed at T={temperature_K:.6g} K."
            ) from exc

        if gamma_l.shape != phi_l.shape:
            raise ModelError("Activity model returned inconsistent coefficient shapes.")
        if np.any(~np.isfinite(gamma_l)) or np.any(gamma_l <= 0.0):
            raise ModelError("Activity model returned non-positive activity coefficients.")

    K = gamma_l * phi_l / phi_v
    if np.any(~np.isfinite(K)) or np.any(K <= 0.0):
        raise ModelError("K-value update produced non-positive or non-finite values.")

    return K


def _is_valid_composition(values: np.ndarray) -> bool:
    if np.any(~np.isfinite(values)):
        return False
    if np.any(values < -COMPOSITION_SUM_TOL):
        return False
    return abs(float(np.sum(values)) - 1.0) <= 10.0 * COMPOSITION_SUM_TOL
