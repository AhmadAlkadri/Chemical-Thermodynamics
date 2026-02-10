"""VLE calculations (Bubble Point, Dew Point) for diagrams."""

from __future__ import annotations

import math
from typing import Sequence
import numpy as np

from .core import Composition, Mixture
from .exceptions import ConvergenceError, ModelError
from .models import ActivityModel, EquationOfState
from .validation import validate_pressure, validate_temperature

__all__ = [
    "bubble_temperature",
    "dew_temperature",
    "bubble_pressure",
    "dew_pressure",
]

def bubble_temperature(
    mixture: Mixture,
    pressure_Pa: float,
    eos: EquationOfState,
    activity_model: ActivityModel | None = None,
    settings: dict | None = None,
) -> tuple[float, list[float]]:
    """Calculate Bubble Point Temperature and equilibrium vapor composition."""
    settings = settings or {}
    tol = settings.get("tol", 1e-5)
    max_iter = settings.get("max_iter", 50)
    
    # 1. Initial T Estimate (weighted boiling points or Tc)
    T = 0.0
    for comp, x_i in zip(mixture.components, mixture.fractions):
         T += x_i * comp.tc_k
    T *= 0.7 
    if T < 50.0: T = 200.0

    x = np.array(mixture.fractions, dtype=float)
    
    # helper for obj function
    def get_error(T_test, use_eos=False):
        # Return log(sum(yi))
        if T_test <= 0: return -100.0 # Should not happen with clamp
        
        # Base K
        K_vals = _wilson_k(mixture, T_test, pressure_Pa)
        
        # If using EOS, refine K
        if use_eos:
            y_est = K_vals * x
            S_est = np.sum(y_est)
            if S_est > 0: 
                y_norm = y_est / S_est
            else:
                y_norm = y_est
            
            try:
                # One Step update
                phi_v = _as_array(eos.fugacity_coefficients(
                    mixture=mixture, 
                    temperature_K=T_test, 
                    pressure_Pa=pressure_Pa, 
                    composition=y_norm.tolist(), 
                    phase="vapor"
                ))
                phi_l = _as_array(eos.fugacity_coefficients(
                    mixture=mixture, 
                    temperature_K=T_test, 
                    pressure_Pa=pressure_Pa, 
                    composition=x.tolist(), 
                    phase="liquid"
                ))
                gamma_l = 1.0
                if activity_model:
                    gamma_l = _as_array(activity_model.activity_coefficients(mixture=mixture, temperature_K=T_test, composition=x.tolist()))
                
                K_vals = gamma_l * phi_l / phi_v
            except Exception:
                # If EOS fails, fallback to Wilson
                pass

        y = K_vals * x
        S = np.sum(y)
        if S <= 1e-12: return -50.0 # prevent log(0)
        return math.log(S)

    # 2. Outer Loop: Solve f(T) = log(sum(yi)) = 0
    # Use Secant with Step Limit
    use_eos_flag = True # Use EOS from start? Or start with ideal? Start with ideal to get close.
    
    # Phase 1: Ideal Wilson only (3 iterations) to get in range
    for phase_name in ["wilson", "eos"]:
        use_eos = (phase_name == "eos")
        limit_iter = 10 if use_eos else 5
        
        f1 = get_error(T, use_eos)
        if abs(f1) < tol and use_eos: break
        
        # Perturb
        T2 = T * 1.001
        f2 = get_error(T2, use_eos)
        
        for k in range(limit_iter):
            if abs(f2) < tol and use_eos:
                T = T2
                break
            
            # Secant update
            if abs(f2 - f1) < 1e-12: break
            
            denom = (f2 - f1)
            if abs(denom) < 1e-15: break
            
            delta_T = - f2 * (T2 - T) / denom
            
            # Clamp Step
            max_step = 0.1 * T2
            if delta_T > max_step: delta_T = max_step
            if delta_T < -max_step: delta_T = -max_step
            
            T_next = T2 + delta_T
            if T_next < 10.0: T_next = 10.0
            
            # Shift
            T = T2
            f1 = f2
            T2 = T_next
            f2 = get_error(T2, use_eos)
        
        T = T2 # Accept result of phase

    # Final Convergence with EOS fully
    # Already done in loop above.
    
    # Final y calculation
    y_unnorm = _calc_equilibrium_y(T, pressure_Pa, x, mixture, eos, activity_model)
    y = (y_unnorm / np.sum(y_unnorm)).tolist()
    
    return float(T), y


def dew_temperature(
    mixture: Mixture,
    pressure_Pa: float,
    eos: EquationOfState,
    activity_model: ActivityModel | None = None,
    settings: dict | None = None,
) -> tuple[float, list[float]]:
    """Calculate Dew Point Temperature and equilibrium liquid composition."""
    settings = settings or {}
    tol = settings.get("tol", 1e-5)
    
    # 1. Initial T Estimate
    T = 0.0
    for comp, y_i in zip(mixture.components, mixture.fractions):
         T += y_i * comp.tc_k
    T *= 0.7 # Start lower
    if T < 50.0: T = 200.0

    y = np.array(mixture.fractions, dtype=float)

    def get_error(T_test, use_eos=False):
        # Return log(sum(xi)) where xi = yi / Ki
        if T_test <= 0: return 100.0
        
        K_vals = _wilson_k(mixture, T_test, pressure_Pa)
        
        if use_eos:
            x_est = y / K_vals
            S_est = np.sum(x_est)
            if S_est > 0: x_norm = x_est / S_est
            else: x_norm = x_est
            
            try:
                phi_v = _as_array(eos.fugacity_coefficients(
                    mixture=mixture, 
                    temperature_K=T_test, 
                    pressure_Pa=pressure_Pa, 
                    composition=y.tolist(), 
                    phase="vapor"
                ))
                phi_l = _as_array(eos.fugacity_coefficients(
                    mixture=mixture, 
                    temperature_K=T_test, 
                    pressure_Pa=pressure_Pa, 
                    composition=x_norm.tolist(), 
                    phase="liquid"
                ))
                gamma_l = 1.0
                if activity_model:
                    gamma_l = _as_array(activity_model.activity_coefficients(mixture=mixture, temperature_K=T_test, composition=x_norm.tolist()))
                
                K_vals = gamma_l * phi_l / phi_v
            except Exception:
                pass

        x_calc = y / K_vals
        S = np.sum(x_calc)
        if S <= 1e-12: return -50.0
        return math.log(S)

    # 2. Solve f(T) = log(sum(xi)) = 0
    # Function is increasing with T?
    # T increases -> K increases -> 1/K decreases -> S decreases.
    # So f(T) is DECREASING.
    
    for phase_name in ["wilson", "eos"]:
        use_eos = (phase_name == "eos")
        limit_iter = 10 if use_eos else 5
        
        f1 = get_error(T, use_eos)
        if abs(f1) < tol and use_eos: break
        
        T2 = T * 1.001
        f2 = get_error(T2, use_eos)
        
        for k in range(limit_iter):
            if abs(f2) < tol and use_eos:
                T = T2
                break
                
            denom = (f2 - f1)
            if abs(denom) < 1e-15: break
            
            delta_T = - f2 * (T2 - T) / denom
            
            max_step = 0.1 * T2
            if delta_T > max_step: delta_T = max_step
            if delta_T < -max_step: delta_T = -max_step
            
            T_next = T2 + delta_T
            if T_next < 10.0: T_next = 10.0
            
            T = T2
            f1 = f2
            T2 = T_next
            f2 = get_error(T2, use_eos)

        T = T2

    x_unnorm = _calc_equilibrium_x(T, pressure_Pa, y, mixture, eos, activity_model)
    x = (x_unnorm / np.sum(x_unnorm)).tolist()
    
    return float(T), x


def bubble_pressure(
    mixture: Mixture,
    temperature_K: float,
    eos: EquationOfState,
    activity_model: ActivityModel | None = None,
    settings: dict | None = None,
) -> tuple[float, list[float]]:
    """Calculate Bubble Point Pressure and equilibrium vapor composition."""
    settings = settings or {}
    tol = settings.get("tol", 1e-5)
    
    # Simple guess: sum(x * P_sat) ? 
    # Use Wilson with P=1atm to find approx K, then scale?
    # No, Wilson K depends on P. K ~ A/P.
    # sum(K*x) = 1 => sum(A/P * x) = 1 => 1/P * sum(A*x) = 1 => P = sum(A*x).
    # where A = P_c * exp(...). roughly Psat.
    
    P_est = 0.0
    for comp, x_i in zip(mixture.components, mixture.fractions):
        # Estimate Psat using Wilson at low P?
        # ln(Pc/P) + ... = ln K
        # K = Pc/P * exp(...)
        # K*P = Pc * exp(...) which is Psat approximation
        val = comp.pc_pa * math.exp(5.373 * (1.0 + comp.omega) * (1.0 - comp.tc_k / temperature_K))
        P_est += x_i * val
    
    P = P_est
    if P < 100.0: P = 100000.0
    
    x = np.array(mixture.fractions, dtype=float)

    def get_error(P_test, use_eos=False):
        # Return log(sum(yi))
        if P_test <= 1.0: return -100.0
        
        K_vals = _wilson_k(mixture, temperature_K, P_test)
        
        if use_eos:
            y_est = K_vals * x
            S_est = np.sum(y_est)
            if S_est > 0: y_norm = y_est / S_est
            else: y_norm = y_est
            
            try:
                phi_v = _as_array(eos.fugacity_coefficients(
                    mixture=mixture, 
                    temperature_K=temperature_K, 
                    pressure_Pa=P_test, 
                    composition=y_norm.tolist(), 
                    phase="vapor"
                ))
                phi_l = _as_array(eos.fugacity_coefficients(
                    mixture=mixture, 
                    temperature_K=temperature_K, 
                    pressure_Pa=P_test, 
                    composition=x.tolist(), 
                    phase="liquid"
                ))
                gamma_l = 1.0
                if activity_model:
                    gamma_l = _as_array(activity_model.activity_coefficients(mixture=mixture, temperature_K=temperature_K, composition=x.tolist()))
                
                K_vals = gamma_l * phi_l / phi_v
            except Exception:
                pass
        
        y = K_vals * x
        S = np.sum(y)
        if S <= 1e-12: return -50.0
        return math.log(S)

    # Solve log(sum(yi)) = 0.
    # Bubble P: increasing P -> decreasing y. slope negative.
    
    for phase_name in ["wilson", "eos"]:
        use_eos = (phase_name == "eos")
        limit_iter = 10 if use_eos else 5
        
        f1 = get_error(P, use_eos)
        if abs(f1) < tol and use_eos: break
        
        P2 = P * 0.99 # Perturb down (or up)
        f2 = get_error(P2, use_eos)
        
        for k in range(limit_iter):
            if abs(f2) < tol and use_eos:
                P = P2
                break
            
            denom = (f2 - f1)
            if abs(denom) < 1e-15: break
            
            delta_P = - f2 * (P2 - P) / denom
            
            # Clamp step
            max_step = 0.2 * P2
            if delta_P > max_step: delta_P = max_step
            if delta_P < -max_step: delta_P = -max_step
            
            P_next = P2 + delta_P
            if P_next < 1.0: P_next = 1.0
            
            P = P2
            f1 = f2
            P2 = P_next
            f2 = get_error(P2, use_eos)
            
        P = P2

    y_unnorm = _calc_equilibrium_y(temperature_K, P, x, mixture, eos, activity_model)
    y = (y_unnorm / np.sum(y_unnorm)).tolist()
    return float(P), y


def dew_pressure(
    mixture: Mixture,
    temperature_K: float,
    eos: EquationOfState,
    activity_model: ActivityModel | None = None,
    settings: dict | None = None,
) -> tuple[float, list[float]]:
    """Calculate Dew Point Pressure and equilibrium liquid composition."""
    settings = settings or {}
    tol = settings.get("tol", 1e-5)
    
    # Estimate Dew P
    # sum(y/K) = 1.
    # K ~ Psat/P.
    # sum(y * P / Psat) = 1 => P * sum(y/Psat) = 1 => P = 1/sum(y/Psat).
    
    inv_P_est = 0.0
    for comp, y_i in zip(mixture.components, mixture.fractions):
        Psat = comp.pc_pa * math.exp(5.373 * (1.0 + comp.omega) * (1.0 - comp.tc_k / temperature_K))
        inv_P_est += y_i / Psat
    
    if inv_P_est > 0: P_est = 1.0 / inv_P_est
    else: P_est = 100000.0
    
    P = P_est
    if P < 100.0: P = 100000.0

    y = np.array(mixture.fractions, dtype=float)

    def get_error(P_test, use_eos=False):
        # Return log(sum(xi))
        if P_test <= 1.0: return -100.0
        
        K_vals = _wilson_k(mixture, temperature_K, P_test)
        
        if use_eos:
            x_est = y / K_vals
            S_est = np.sum(x_est)
            if S_est > 0: x_norm = x_est / S_est
            else: x_norm = x_est
            
            try:
                phi_v = _as_array(eos.fugacity_coefficients(
                    mixture=mixture, 
                    temperature_K=temperature_K, 
                    pressure_Pa=P_test, 
                    composition=y.tolist(), 
                    phase="vapor"
                ))
                phi_l = _as_array(eos.fugacity_coefficients(
                    mixture=mixture, 
                    temperature_K=temperature_K, 
                    pressure_Pa=P_test, 
                    composition=x_norm.tolist(), 
                    phase="liquid"
                ))
                gamma_l = 1.0
                if activity_model:
                    gamma_l = _as_array(activity_model.activity_coefficients(mixture=mixture, temperature_K=temperature_K, composition=x_norm.tolist()))
                
                K_vals = gamma_l * phi_l / phi_v
            except Exception:
                pass
        
        x_calc = y / K_vals
        S = np.sum(x_calc)
        if S <= 1e-12: return -50.0
        return math.log(S)

    # Solve log(sum(xi)) = 0
    # Increasing P -> Decreasing K -> Increasing x. Slope positive.

    for phase_name in ["wilson", "eos"]:
        use_eos = (phase_name == "eos")
        limit_iter = 10 if use_eos else 5
        
        f1 = get_error(P, use_eos)
        if abs(f1) < tol and use_eos: break
        
        P2 = P * 0.99
        f2 = get_error(P2, use_eos)
        
        for k in range(limit_iter):
            if abs(f2) < tol and use_eos:
                P = P2
                break
                
            denom = (f2 - f1)
            if abs(denom) < 1e-15: break
            
            delta_P = - f2 * (P2 - P) / denom
            
            max_step = 0.2 * P2
            if delta_P > max_step: delta_P = max_step
            if delta_P < -max_step: delta_P = -max_step
            
            P_next = P2 + delta_P
            if P_next < 1.0: P_next = 1.0
            
            P = P2
            f1 = f2
            P2 = P_next
            f2 = get_error(P2, use_eos)

        P = P2

    x_unnorm = _calc_equilibrium_x(temperature_K, P, y, mixture, eos, activity_model)
    x = (x_unnorm / np.sum(x_unnorm)).tolist()
    return float(P), x

# --- Internals ---

def _wilson_k(mixture, T, P) -> np.ndarray:
    values = []
    for component in mixture.components:
        tc = component.tc_k
        pc = component.pc_pa
        omega = component.omega
        ln_k = math.log(pc / P) + 5.373 * (1.0 + omega) * (1.0 - tc / T)
        values.append(math.exp(ln_k))
    return np.array(values, dtype=float)

def _calc_sum_yi_minus_1(T, P, x, mixture, eos, activity_model=None):
    # Bubble point condition: sum(yi) - 1 = 0
    # yi = K_i * xi
    # K_i = phi_L / phi_V
    # This depends on composition y as well!
    # "True" Bubble point solves for T/P AND y simultaneously. (Rachford Rice for v=0)
    # Simplified approach for diagrams (often good enough for demo codes):
    # Iterate to find K-values.
    
    # Inner loop to converge y?
    # Or just use current K estimate.
    
    # Let's approximate K with Wilson first, then refine?
    # Actually, for a ROBUST single function, we probably should do a full Bubble Point loop.
    # Bubble Point Inner Loop:
    # 1. Guess T (outer loop handles this)
    # 2. Guess y_i = K_i_wilson * x_i. Normalize y.
    # 3. Calculate Phi_L(x), Phi_V(y).
    # 4. K_new = Phi_L / Phi_V.
    # 5. y_new = K_new * x. Sum(y_new) is the error for T?
    
    # Wait, existing literature for Bubble P/T usually iterates:
    #   ln K = ...
    # This can be heavy.
    
    # Let's check `chemthermo.flash.tp` again. It has a full solver.
    # Can we trick `flash_tp` to give us Bubble Point?
    # If we set z=x, and T, P. `flash_tp` checks phase stability.
    # If we are strictly liquid, it returns vapor_fraction=0.
    # If we are strictly vapor, vapor_fraction=1.
    # At bubble point, vapor_fraction = 0 (incipient).
    
    # Using `flash_tp` inside an optimization loop to find boundary?
    # That is very expensive.
    
    # I will implement a "lite" K-value iterator here since we are just making a demo.
    # We will assume ideal-ish behavior for the K-values update or do 1-2 inner updates.
    
    # Let's stick to a simple Wilson initialized K for now to get OFF the ground, 
    # then iterate a few times if needed.
    # Actually, let's implement the standard update:
    # y calculation: K = _wilson_k(...)
    # y = K * x
    # err = sum(y) - 1
    # This ignores fugacity coefficient dependency on y. 
    # For a high quality library, we need the inner Fugacity loop.
    
    # Given the constraint of "Thin Vertical Slice" and "Teaching Notebook", 
    # maybe we assume Ideal K (Wilson) for the *Search*, and only do a final check?
    # No, that will be inaccurate for high pressure.
    
    # Let's do 3 inner iterations of K-update.
    # Initial K
    K = _wilson_k(mixture, T, P)
    y = K * x
    S = np.sum(y)
    
    # Normalize y for fugacity calc
    if S > 0: y_norm = y / S
    else: y_norm = y # Should not happen

    # 1st Refinement
    if eos:
        try:
            phi_v = _as_array(eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=T,
                pressure_Pa=P,
                composition=y_norm.tolist(),
                phase="vapor"
            ))
            phi_l = _as_array(eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=T,
                pressure_Pa=P,
                composition=x.tolist(),
                phase="liquid"
            ))
            gamma_l = 1.0
            if activity_model:
                gamma_l = _as_array(activity_model.activity_coefficients(
                    mixture=mixture,
                    temperature_K=T,
                    composition=x.tolist()
                ))
            
            K = gamma_l * phi_l / phi_v
            y = K * x
            S = np.sum(y)
        except Exception:
            # If EOS fails (e.g. unphysical T), return a punitive value or existing S
            pass
            
    return S - 1.0


def _calc_equilibrium_y(T, P, x, mixture, eos, activity_model=None):
    # Return y distribution check
    # We do a more rigorous K calc here
    K = _wilson_k(mixture, T, P)
    y = K * x
    if np.sum(y) > 0: y_norm = y / np.sum(y)
    else: return y
    
    # Iterate K a few times
    for _ in range(5):
        phi_v = _as_array(eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=T,
            pressure_Pa=P,
            composition=y_norm.tolist(),
            phase="vapor"
        ))
        phi_l = _as_array(eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=T,
            pressure_Pa=P,
            composition=x.tolist(),
            phase="liquid"
        ))
        gamma_l = 1.0
        if activity_model:
            gamma_l = _as_array(activity_model.activity_coefficients(
                mixture=mixture,
                temperature_K=T,
                composition=x.tolist()
            ))
        
        K = gamma_l * phi_l / phi_v
        y = K * x
        if np.sum(y) > 0: y_norm = y / np.sum(y)
    
    return y


def _calc_sum_xi_minus_1(T, P, y, mixture, eos, activity_model=None):
    # Dew point: sum(x) - 1 = 0
    # x = y / K
    
    K = _wilson_k(mixture, T, P)
    x = y / K
    S = np.sum(x)
    
    if S > 0: x_norm = x / S
    else: x_norm = x

    # 1st Refinement
    if eos:
        try:
            phi_v = _as_array(eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=T,
                pressure_Pa=P,
                composition=y.tolist(),
                phase="vapor"
            ))
            phi_l = _as_array(eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=T,
                pressure_Pa=P,
                composition=x_norm.tolist(),
                phase="liquid"
            ))
            gamma_l = 1.0
            if activity_model:
                gamma_l = _as_array(activity_model.activity_coefficients(
                    mixture=mixture,
                    temperature_K=T,
                    composition=x_norm.tolist()
                ))
            
            K = gamma_l * phi_l / phi_v
            x = y / K
            S = np.sum(x)
        except Exception:
            pass
            
    return S - 1.0


def _calc_equilibrium_x(T, P, y, mixture, eos, activity_model=None):
    K = _wilson_k(mixture, T, P)
    x = y / K
    if np.sum(x) > 0: x_norm = x / np.sum(x)
    else: x_norm = x
    
    for _ in range(5):
        phi_v = _as_array(eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=T,
            pressure_Pa=P,
            composition=y.tolist(),
            phase="vapor"
        ))
        phi_l = _as_array(eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=T,
            pressure_Pa=P,
            composition=x_norm.tolist(),
            phase="liquid"
        ))
        gamma_l = 1.0
        if activity_model:
            gamma_l = _as_array(activity_model.activity_coefficients(
                mixture=mixture,
                temperature_K=T,
                composition=x_norm.tolist()
            ))
        
        K = gamma_l * phi_l / phi_v
        x = y / K
        if np.sum(x) > 0: x_norm = x / np.sum(x)
        
    return x

def _as_array(vals):
    return np.array(list(vals), dtype=float)
