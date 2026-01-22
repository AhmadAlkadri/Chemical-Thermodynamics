"""TP flash solver using phi-phi approach."""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from ..core import Composition, Mixture
from ..exceptions import CompositionError, ConvergenceError, ModelError
from ..models import EquationOfState
from ..validation import COMPOSITION_SUM_TOL, validate_pressure, validate_temperature
from .results import FlashResult, PhaseResult
from .settings import FlashSettings


def flash_tp(
    mixture: Mixture,
    *,
    temperature_K: float,
    pressure_Pa: float,
    eos: EquationOfState | None,
    settings: FlashSettings | None = None,
) -> FlashResult:
    """Perform a TP flash calculation using a cubic EOS and phi-phi updates."""

    temperature = validate_temperature(temperature_K)
    pressure = validate_pressure(pressure_Pa)

    if eos is None:
        raise ModelError("An equation-of-state model is required for flash_tp.")

    if mixture.basis != "mole":
        raise ModelError("flash_tp currently requires mole-fraction compositions.")

    settings = settings or FlashSettings()

    z = np.array(mixture.fractions, dtype=float)
    if z.size == 0:
        raise CompositionError("Mixture composition must be non-empty.")

    K = _wilson_k(mixture, temperature, pressure)
    if np.any(K <= 0.0):
        raise ModelError("Non-positive K-values encountered in Wilson estimate.")

    k_min = float(np.min(K))
    k_max = float(np.max(K))

    if np.all(K <= 1.0):
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name="liquid",
            vapor_fraction=0.0,
            diagnostics={"k_min": k_min, "k_max": k_max, "iterations": 0, "converged": True},
        )

    if np.all(K >= 1.0):
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name="vapor",
            vapor_fraction=1.0,
            diagnostics={"k_min": k_min, "k_max": k_max, "iterations": 0, "converged": True},
        )

    vapor_fraction, f0, f1 = _rachford_rice(z, K)
    if vapor_fraction is None:
        phase_name = "vapor" if f0 > 0.0 else "liquid"
        vapor_fraction_value = 1.0 if phase_name == "vapor" else 0.0
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name=phase_name,
            vapor_fraction=vapor_fraction_value,
            diagnostics={
                "k_min": k_min,
                "k_max": k_max,
                "iterations": 0,
                "converged": True,
                "rr_f0": float(f0),
                "rr_f1": float(f1),
                "rr_status": "no_root",
            },
        )

    max_delta = float("inf")
    for iteration in range(1, settings.max_iter + 1):
        x = z / (1.0 + vapor_fraction * (K - 1.0))
        if np.any(x < 0.0):
            raise ConvergenceError("Negative liquid composition encountered in flash_tp.")
        x /= x.sum()

        y = K * x
        if np.any(y < 0.0):
            raise ConvergenceError("Negative vapor composition encountered in flash_tp.")
        y /= y.sum()

        phi_v = _as_array(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature,
                pressure_Pa=pressure,
                composition=y.tolist(),
                phase="vapor",
            )
        )
        phi_l = _as_array(
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

        K_new = phi_l / phi_v
        max_delta = float(np.max(np.abs(K_new - K)))
        if max_delta < settings.tol:
            return _two_phase_result(
                mixture,
                temperature,
                pressure,
                x,
                y,
                vapor_fraction,
                diagnostics={
                    "iterations": iteration,
                    "converged": True,
                    "max_delta_k": max_delta,
                    "k_min": float(np.min(K_new)),
                    "k_max": float(np.max(K_new)),
                },
            )

        if settings.damping is None:
            K = K_new
        else:
            K = K + settings.damping * (K_new - K)

        if np.any(K <= 0.0):
            raise ModelError("Non-positive K-values encountered during iteration.")

        vapor_fraction, f0, f1 = _rachford_rice(z, K)
        if vapor_fraction is None:
            raise ConvergenceError("Rachford-Rice failed to bracket a vapor fraction.")

    raise ConvergenceError(
        "flash_tp did not converge within the iteration limit; "
        f"max_delta_k={max_delta:.3e}."
    )


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
    return FlashResult(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        phases={phase_name: phase},
        vapor_fraction=vapor_fraction,
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
    return FlashResult(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        phases={"liquid": liquid, "vapor": vapor},
        vapor_fraction=float(vapor_fraction),
        diagnostics=diagnostics,
    )


def _wilson_k(mixture: Mixture, temperature: float, pressure: float) -> np.ndarray:
    values = []
    for component in mixture.components:
        tc = component.tc_k
        pc = component.pc_pa
        omega = component.omega

        ln_k = math.log(pc / pressure) + 5.373 * (1.0 + omega) * (1.0 - tc / temperature)
        values.append(math.exp(ln_k))

    return np.array(values, dtype=float)


def _rachford_rice(z: np.ndarray, K: np.ndarray) -> tuple[float | None, float, float]:
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


def _as_array(values: Sequence[float]) -> np.ndarray:
    return np.array(list(values), dtype=float)
