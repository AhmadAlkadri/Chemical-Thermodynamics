"""TP flash calculations (T,P) for VLE/VLLE in SI units.

Supports phi-phi, gamma-phi, and a staged VLLE scaffold. The solver uses
Wilson K-value initialization, Rachford-Rice vapor fraction updates, and
fixed-point K updates. Given the same inputs and settings, results are
deterministic.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from ..core import Composition, Mixture
from ..exceptions import CompositionError, ConvergenceError, ModelError
from ..models import ActivityModel, EquationOfState
from ..validation import COMPOSITION_SUM_TOL, validate_pressure, validate_temperature
from .results import FlashResult, PhaseResult
from .settings import FlashSettings


@dataclass(frozen=True)
class _LleResult:
    x1: np.ndarray
    x2: np.ndarray
    beta_l1: float
    beta_l2: float
    iterations: int
    converged: bool
    max_delta_x: float
    tol: float
    composition_residual: float
    termination_reason: str
    separation: float


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
    """Perform a TP flash calculation using phi-phi, gamma-phi, or VLLE scaffolding.

    Args:
        mixture: Mixture with mole-fraction composition.
        temperature_K: Temperature in K.
        pressure_Pa: Pressure in Pa.
        eos: Equation-of-state model used for fugacity coefficients.
        activity_model: Activity model (required for gamma-phi and VLLE).
        flash_mode: Case-insensitive mode: "phi-phi", "gamma-phi", or "vlle".
        settings: Iteration controls (tolerance, damping, max iterations).

    Returns:
        FlashResult with phase compositions and fractions. Phase names follow:
        VLE -> "liquid"/"vapor", single-phase -> "liquid" or "vapor",
        VLLE -> "vapor"/"liquid1"/"liquid2".

    Diagnostics:
        Diagnostics keys are implementation details. Current stable keys include:
        - VLE: iterations, converged, termination_reason, max_delta_k, k_min,
          k_max, phase_count, phase_state, phase_regime, flash_mode.
        - VLE single-phase fallbacks may also include rr_f0, rr_f1, rr_status.
        - VLLE: phase_regime, phase_state, phase_count, flash_mode, vle_iterations,
          lle_iterations, lle_max_delta_x, lle_composition_residual, lle_tol,
          termination_reason, converged, and optional vlle_gamma_max.

    Raises:
        InputRangeError: If temperature or pressure is non-physical.
        ModelError: If required models are missing or return invalid values.
        CompositionError: If the mixture composition is invalid.
        ConvergenceError: If iteration fails to converge.
        ValueError: If VLLE is requested without an activity model.

    Notes:
        Single-phase fallbacks are returned (not raised) when K-bounds or
        Rachford-Rice root checks indicate only one phase is possible.
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
        if activity_model is None:
            raise ValueError("An activity model is required for VLLE flash.")
        return _flash_tp_vlle(
            mixture,
            temperature,
            pressure,
            eos=eos,
            activity_model=activity_model,
            settings=settings,
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
    """Internal TP VLE solver (phi-phi or gamma-phi) with K-value iteration."""
    if mixture.basis != "mole":
        raise ModelError("flash_tp currently requires mole-fraction compositions.")

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
            diagnostics={
                "k_min": k_min,
                "k_max": k_max,
                "iterations": 0,
                "converged": True,
                "termination_reason": "single_phase_k_bounds",
                "max_delta_k": 0.0,
                "phase_count": 1,
                "phase_state": "liquid",
                "phase_regime": "single-phase",
                "flash_mode": mode,
            },
        )

    if np.all(K >= 1.0):
        return _single_phase_result(
            mixture,
            temperature,
            pressure,
            phase_name="vapor",
            vapor_fraction=1.0,
            diagnostics={
                "k_min": k_min,
                "k_max": k_max,
                "iterations": 0,
                "converged": True,
                "termination_reason": "single_phase_k_bounds",
                "max_delta_k": 0.0,
                "phase_count": 1,
                "phase_state": "vapor",
                "phase_regime": "single-phase",
                "flash_mode": mode,
            },
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
                "termination_reason": "rr_no_root",
                "max_delta_k": 0.0,
                "phase_count": 1,
                "phase_state": phase_name,
                "phase_regime": "single-phase",
                "flash_mode": mode,
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

        if mode == "gamma-phi":
            assert activity_model is not None
            gamma_l = _as_array(
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
                    "termination_reason": "tolerance_met",
                    "max_delta_k": max_delta,
                    "k_min": float(np.min(K_new)),
                    "k_max": float(np.max(K_new)),
                    "phase_count": 2,
                    "phase_state": "two_phase",
                    "phase_regime": "VLE",
                    "flash_mode": mode,
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


def _flash_tp_vlle(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    *,
    eos: EquationOfState,
    activity_model: ActivityModel,
    settings: FlashSettings,
) -> FlashResult:
    """Stage gamma-phi VLE, then split the liquid into two LLE phases."""
    z = np.array(mixture.fractions, dtype=float)
    if z.size == 0:
        raise CompositionError("Mixture composition must be non-empty.")

    gamma_max = (
        _gamma_max(
            activity_model=activity_model,
            mixture=mixture,
            temperature_K=temperature,
            composition=z,
        )
        if activity_model is not None
        else None
    )

    try:
        vle_result = _flash_tp_vle(
            mixture,
            temperature,
            pressure,
            eos=eos,
            activity_model=activity_model,
            mode="gamma-phi",
            settings=settings,
        )
    except ConvergenceError as exc:
        raise ConvergenceError("VLLE staging failed during VLE solve.") from exc

    if "liquid" not in vle_result.phases:
        diagnostics = _merge_diagnostics(
            vle_result.diagnostics,
            {
                "phase_regime": "single-phase",
                "termination_reason": "vlle_no_liquid",
                "flash_mode": "vlle",
            },
        )
        return FlashResult(
            temperature_K=vle_result.temperature_K,
            pressure_Pa=vle_result.pressure_Pa,
            phases=vle_result.phases,
            vapor_fraction=vle_result.vapor_fraction,
            phase_fractions=vle_result.phase_fractions,
            diagnostics=diagnostics,
        )

    beta_v = float(vle_result.vapor_fraction or 0.0)
    if beta_v >= 1.0 - COMPOSITION_SUM_TOL:
        diagnostics = _merge_diagnostics(
            vle_result.diagnostics,
            {
                "phase_regime": "single-phase",
                "termination_reason": "vlle_no_liquid",
                "flash_mode": "vlle",
            },
        )
        return FlashResult(
            temperature_K=vle_result.temperature_K,
            pressure_Pa=vle_result.pressure_Pa,
            phases=vle_result.phases,
            vapor_fraction=vle_result.vapor_fraction,
            phase_fractions=vle_result.phase_fractions,
            diagnostics=diagnostics,
        )

    liquid_phase = vle_result.phases["liquid"]
    z_liq = np.array(liquid_phase.composition.fractions, dtype=float)

    lle_result = _solve_lle(
        mixture,
        temperature,
        pressure,
        z_liq,
        eos=eos,
        activity_model=activity_model,
        settings=settings,
    )

    gamma_diag = {"vlle_gamma_max": gamma_max} if gamma_max is not None else {}
    if not lle_result.converged or lle_result.separation < 1e-3:
        diagnostics = _merge_diagnostics(
            vle_result.diagnostics,
            {
                "phase_regime": "VLE",
                "termination_reason": "vlle_lle_failed",
                "lle_iterations": lle_result.iterations,
                "lle_max_delta_x": lle_result.max_delta_x,
                "lle_tol": lle_result.tol,
                "flash_mode": "vlle",
            },
        )
        diagnostics.update(gamma_diag)
        return FlashResult(
            temperature_K=vle_result.temperature_K,
            pressure_Pa=vle_result.pressure_Pa,
            phases=vle_result.phases,
            vapor_fraction=vle_result.vapor_fraction,
            phase_fractions=vle_result.phase_fractions,
            diagnostics=diagnostics,
        )

    beta_l1_total = (1.0 - beta_v) * lle_result.beta_l1
    beta_l2_total = (1.0 - beta_v) * lle_result.beta_l2

    diagnostics = _merge_diagnostics(
        vle_result.diagnostics,
        {
            "phase_regime": "VLLE",
            "phase_count": 3,
            "phase_state": "three_phase",
            "flash_mode": "vlle",
            "vle_iterations": vle_result.diagnostics.get("iterations", 0),
            "lle_iterations": lle_result.iterations,
            "lle_max_delta_x": lle_result.max_delta_x,
            "lle_composition_residual": lle_result.composition_residual,
            "lle_tol": lle_result.tol,
            "termination_reason": lle_result.termination_reason,
            "converged": vle_result.diagnostics.get("converged", False)
            and lle_result.converged,
        },
    )
    diagnostics.update(gamma_diag)

    vapor = vle_result.phases.get("vapor")
    if vapor is None:
        diagnostics = _merge_diagnostics(
            diagnostics,
            {
                "phase_regime": "single-phase",
                "termination_reason": "vlle_no_vapor",
            },
        )
        return FlashResult(
            temperature_K=vle_result.temperature_K,
            pressure_Pa=vle_result.pressure_Pa,
            phases=vle_result.phases,
            vapor_fraction=vle_result.vapor_fraction,
            phase_fractions=vle_result.phase_fractions,
            diagnostics=diagnostics,
        )

    return _three_phase_result(
        mixture,
        temperature,
        pressure,
        lle_result.x1,
        lle_result.x2,
        np.array(vapor.composition.fractions, dtype=float),
        beta_v,
        beta_l1_total,
        beta_l2_total,
        diagnostics=diagnostics,
    )


def _solve_lle(
    mixture: Mixture,
    temperature: float,
    pressure: float,
    z_liq: np.ndarray,
    *,
    eos: EquationOfState,
    activity_model: ActivityModel,
    settings: FlashSettings,
) -> _LleResult:
    """Fixed-point LLE split for a liquid feed composition."""
    x1, x2 = _initialize_liquid_split(z_liq)
    anchor_x1 = x1.copy()
    anchor_x2 = x2.copy()
    max_delta = float("inf")
    iterations = 0
    damping = 0.5 if settings.damping is None else settings.damping
    anchor_weight = 0.1
    lle_tol = max(settings.tol * 10.0, 1e-6)
    for iteration in range(1, settings.max_iter + 1):
        iterations = iteration
        gamma1 = _as_array(
            activity_model.activity_coefficients(
                mixture=mixture,
                temperature_K=temperature,
                composition=x1.tolist(),
            )
        )
        gamma2 = _as_array(
            activity_model.activity_coefficients(
                mixture=mixture,
                temperature_K=temperature,
                composition=x2.tolist(),
            )
        )
        phi1 = _as_array(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature,
                pressure_Pa=pressure,
                composition=x1.tolist(),
                phase="liquid",
            )
        )
        phi2 = _as_array(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature,
                pressure_Pa=pressure,
                composition=x2.tolist(),
                phase="liquid",
            )
        )

        if np.any(gamma1 <= 0.0) or np.any(gamma2 <= 0.0):
            raise ModelError("Activity model returned non-positive activity coefficients.")
        if np.any(phi1 <= 0.0) or np.any(phi2 <= 0.0):
            raise ModelError("EOS returned non-positive fugacity coefficients.")

        ratio = gamma1 * phi1 / (gamma2 * phi2)
        x2_target = _normalize(x1 * ratio)
        x1_target = _normalize(x2 / ratio)

        new_x1 = _normalize((1.0 - damping) * x1 + damping * x1_target)
        new_x2 = _normalize((1.0 - damping) * x2 + damping * x2_target)
        new_x1 = _normalize((1.0 - anchor_weight) * new_x1 + anchor_weight * anchor_x1)
        new_x2 = _normalize((1.0 - anchor_weight) * new_x2 + anchor_weight * anchor_x2)

        max_delta = float(max(np.max(np.abs(new_x1 - x1)), np.max(np.abs(new_x2 - x2))))
        x1 = new_x1
        x2 = new_x2

        if max_delta < lle_tol:
            break

    beta_l1 = _compute_beta_l1(z_liq, x1, x2)
    separation = float(np.max(np.abs(x1 - x2)))

    if beta_l1 is None:
        return _LleResult(
            x1=x1,
            x2=x2,
            beta_l1=0.5,
            beta_l2=0.5,
            iterations=iterations,
            converged=False,
            max_delta_x=max_delta,
            tol=lle_tol,
            composition_residual=float("inf"),
            termination_reason="lle_invalid_split",
            separation=separation,
        )

    beta_l2 = 1.0 - beta_l1
    z_recon = beta_l1 * x1 + beta_l2 * x2
    composition_residual = float(np.max(np.abs(z_recon - z_liq)))

    converged = max_delta < lle_tol
    termination_reason = "lle_converged" if converged else "lle_max_iter"

    return _LleResult(
        x1=x1,
        x2=x2,
        beta_l1=beta_l1,
        beta_l2=beta_l2,
        iterations=iterations,
        converged=converged,
        max_delta_x=max_delta,
        tol=lle_tol,
        composition_residual=composition_residual,
        termination_reason=termination_reason,
        separation=separation,
    )


def _initialize_liquid_split(z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    count = z.size
    if count == 1:
        return z.copy(), z.copy()

    z = _normalize(z)
    delta = 0.15
    if count == 2:
        x1 = np.array([min(1.0 - 1e-6, max(1e-6, z[0] + delta)), 0.0], dtype=float)
        x1[1] = 1.0 - x1[0]
        x2 = np.array([max(1e-6, min(1.0 - 1e-6, z[0] - delta)), 0.0], dtype=float)
        x2[1] = 1.0 - x2[0]
        return x1, x2

    idx = int(np.argmax(z))
    x1 = z.copy()
    x2 = z.copy()
    x1[idx] = min(1.0 - 1e-6, z[idx] + delta)
    x2[idx] = max(1e-6, z[idx] - delta)
    x1 = _normalize(x1)
    x2 = _normalize(x2)
    return x1, x2


def _compute_beta_l1(z: np.ndarray, x1: np.ndarray, x2: np.ndarray) -> float | None:
    values = []
    for zi, x1i, x2i in zip(z, x1, x2):
        denom = x1i - x2i
        if abs(denom) > 1e-8:
            values.append((zi - x2i) / denom)

    if not values:
        return None

    beta = float(np.mean(values))
    if not (0.0 <= beta <= 1.0):
        return None
    return beta


def _gamma_max(
    *,
    activity_model: ActivityModel,
    mixture: Mixture,
    temperature_K: float,
    composition: np.ndarray,
) -> float | None:
    gamma = activity_model.activity_coefficients(
        mixture=mixture,
        temperature_K=temperature_K,
        composition=composition.tolist(),
    )
    return float(max(gamma)) if gamma else None


def _merge_diagnostics(
    base: Mapping[str, float | str | int | bool],
    extra: Mapping[str, float | str | int | bool],
) -> dict[str, float | str | int | bool]:
    merged = dict(base)
    merged.update(extra)
    return merged


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


def _three_phase_result(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    x1: np.ndarray,
    x2: np.ndarray,
    y: np.ndarray,
    beta_v: float,
    beta_l1: float,
    beta_l2: float,
    *,
    diagnostics: dict[str, float | int | str | bool],
) -> FlashResult:
    liquid1 = PhaseResult(
        name="liquid1",
        composition=Composition(
            fractions=tuple(x1.tolist()),
            basis=mixture.basis,
            normalize=False,
            tol=COMPOSITION_SUM_TOL,
        ),
    )
    liquid2 = PhaseResult(
        name="liquid2",
        composition=Composition(
            fractions=tuple(x2.tolist()),
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
    phase_fractions = {
        "vapor": float(beta_v),
        "liquid1": float(beta_l1),
        "liquid2": float(beta_l2),
    }
    return FlashResult(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        phases={"vapor": vapor, "liquid1": liquid1, "liquid2": liquid2},
        vapor_fraction=float(beta_v),
        phase_fractions=phase_fractions,
        diagnostics=diagnostics,
    )


def _wilson_k(mixture: Mixture, temperature: float, pressure: float) -> np.ndarray:
    """Wilson K-value estimate for each component (dimensionless)."""
    values = []
    for component in mixture.components:
        tc = component.tc_k
        pc = component.pc_pa
        omega = component.omega

        ln_k = math.log(pc / pressure) + 5.373 * (1.0 + omega) * (1.0 - tc / temperature)
        values.append(math.exp(ln_k))

    return np.array(values, dtype=float)


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


def _normalize(values: np.ndarray) -> np.ndarray:
    """Clip to positive values and normalize to sum to 1."""
    values = np.clip(values, 1e-12, None)
    total = float(np.sum(values))
    if total <= 0.0:
        raise CompositionError("Composition fractions must sum to a positive value.")
    return values / total


def _as_array(values: Sequence[float]) -> np.ndarray:
    return np.array(list(values), dtype=float)
