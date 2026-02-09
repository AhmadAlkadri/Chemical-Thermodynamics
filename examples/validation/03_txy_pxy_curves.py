"""
Slice 3: Txy or Pxy curve extraction
------------------------------------
Compares chemthermo vs thermo for binary VLE curves.
Choice: Txy diagram at fixed P.

Scope:
- System: Methane (1) + Ethane (2)
- State: Fixed P = 20 bar
- Sweep: x1/y1 in [0.0, 0.1, ..., 1.0]

Compute:
- Bubble Temperature T_bub(x)
- Dew Temperature T_dew(y)

Metrics:
- RMS difference in T_bub and T_dew
"""

import sys
import csv
import warnings
from pathlib import Path

import numpy as np
import chemthermo as ct
from chemthermo.models.peng_robinson import PengRobinsonEOS

try:
    import thermo
    from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL
except ImportError:
    print("Error: 'thermo' package not found.")
    sys.exit(1)

def solve_bubble_point(mixture, P, x, eos, T_guess=150.0, tol=1e-5):
    """Simple bubble point solver for T."""
    # This is a naive implementation for validation.
    # At bubble point: sum(K*x) = 1.
    # We iterate T. T affects K.
    # New T = Old T * (sum(Kx))^-0.1 (simple dumping) or Newton.
    
    T = T_guess
    for _ in range(50):
        try:
            # We need K-values. We can get them from fugacities at P, T, x (liquid) and assume y=Kx.
            # But y is unknown.
            # Standard algorithm involves inner loop for y.
            # Let's try probing via flash_tp?
            # If we flash at T, P, z=x, we get phase split.
            # If V/F close to 0, we are at bubble point?
            # No, flash is robust.
            # Let's iterate T until flash gives V/F = 0 (saturated liquid).
            # Actually, standard bubble T algo:
            # 1. Guess T.
            # 2. Calc K values (assume yi = Ki*xi).
            # 3. Calc sum_yi = sum(Ki*xi).
            # 4. Check convergence: sum_yi = 1.
            # ... this requires ideal solution estimate first.
            
            # Use chemthermo loop if we had one. We don't have bubble_t exposed yet.
            # Let's brute force Search using flash_tp?
            # Start low T (liquid), increase until 2-phase?
            # Better: Use thermo as reference helper? No, that defeats validation.
            
            # Use a robust bisection/secant on sum(Kx) - 1?
            # Need K values.
            # K_i = phi_l / phi_v. Phi_v depends on y. y depends on K.
            # This is complex to implement from scratch in a script without risk.
            
            # Alternative: Use "Flash at fixed P, z=x".
            # If T < T_bubble, VF=0. If T > T_dew, VF=1. If T_bubble < T < T_dew, 0 < VF < 1.
            # Bubble point is where VF -> 0+.
            # So find T such that flash_tp(T, P, z=x).vapor_fraction = 1e-6 approx.
            
            # Simple Bisection:
            T_min, T_max = 100.0, 300.0
            for _ in range(20):
                T_mid = 0.5 * (T_min + T_max)
                res = None
                try:
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore", category=RuntimeWarning, module="chemicals")
                        res = ct.flash_tp(mixture, temperature_K=T_mid, pressure_Pa=P, eos=eos)
                    vf = res.vapor_fraction
                except ct.exceptions.ConvergenceError:
                    vf = None

                # If single phase liquid (vf=0 or None/liquid), we are too cold. T_mid is too low.
                if vf == 0.0 or (vf is None and res and "liquid" in getattr(res, "phases", {})):
                    T_min = T_mid
                # If single phase vapor (vf=1 or None/vapor), we are too hot.
                elif vf == 1.0 or (vf is None and res and "vapor" in getattr(res, "phases", {})):
                    T_max = T_mid
                # If 2-phase, check VF.
                elif vf is not None:
                    # We want VF ~ 0.
                    if vf > 1e-4:
                        T_max = T_mid # Too hot, too much vapor
                    else:
                        return T_mid # Close enough to bubble point
                else:
                    return None
            
            return 0.5*(T_min + T_max)

        except Exception:
            return None
            
    return None

def solve_dew_point(mixture, P, y, eos, T_guess=150.0):
    # Dew point is where VF -> 1-.
    # Find T such that flash_tp(T, P, z=y).vapor_fraction = 1 - 1e-6.
            
    T_min, T_max = 100.0, 300.0
    for _ in range(20):
        T_mid = 0.5 * (T_min + T_max)
        res = None
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=RuntimeWarning, module="chemicals")
                res = ct.flash_tp(mixture, temperature_K=T_mid, pressure_Pa=P, eos=eos)
            vf = res.vapor_fraction
        except ct.exceptions.ConvergenceError:
            vf = None

        # If single phase liquid, too cold.
        if vf == 0.0 or (vf is None and res and "liquid" in getattr(res, "phases", {})):
            T_min = T_mid
        # If single phase vapor, too hot.
        elif vf == 1.0 or (vf is None and res and "vapor" in getattr(res, "phases", {})):
            T_max = T_mid
        # If 2-phase
        elif vf is not None:
            # We want VF ~ 1.
            if vf < 1.0 - 1e-4:
                T_min = T_mid # Too cold, too much liquid
            else:
                return T_mid # Close enough
        else:
            return None
    
    return 0.5*(T_min + T_max)

def main():
    comps = ["Methane", "Ethane"]
    P = 2.0e6 # 20 bar
    
    print(f"Running Validation Slice 3: {comps[0]} + {comps[1]}")
    print(f"Txy Diagram at P={P/1e5} bar")
    
    # Setup Thermo
    ids = [c.lower() for c in comps]
    constants, properties = ChemicalConstantsPackage.from_IDs(ids)
    kijs = [[0.0]*2 for _ in range(2)]
    eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
    gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
    
    # Setup Chemthermo
    eos = PengRobinsonEOS()
    
    fractions = np.linspace(0.0, 1.0, 11)
    results = []
    
    print(f"{'z1':<5} {'T_bub_ct':<10} {'T_bub_th':<10} {'T_dew_ct':<10} {'T_dew_th':<10} {'Diff_Bub':<10} {'Diff_Dew':<10}")
    print("-" * 80)
    
    for z1 in fractions:
        z2 = 1.0 - z1
        zs = [z1, z2]
        mix = ct.Mixture.from_database(comps, zs)
        
        # --- Bubble Point (x=z) ---
        # Thermo
        try:
            # T_bubble at P, zs
            # thermo doesn't expose strict BubbleT easily via FlashVL without iteration or using
            # specific methods. But flasher.flash(P=P, VF=0, zs=zs) works.
            bub_res = flasher.flash(P=P, VF=0.0, zs=zs)
            ref_T_bub = bub_res.T
        except Exception:
            ref_T_bub = None
            
        # Chemthermo
        ct_T_bub = solve_bubble_point(mix, P, zs, eos)
        
        # --- Dew Point (y=z) ---
        # Thermo
        try:
            dew_res = flasher.flash(P=P, VF=1.0, zs=zs)
            ref_T_dew = dew_res.T
        except Exception:
            ref_T_dew = None
            
        # Chemthermo
        ct_T_dew = solve_dew_point(mix, P, zs, eos)
        
        # Metrics
        d_bub = abs(ct_T_bub - ref_T_bub) if (ct_T_bub and ref_T_bub) else None
        d_dew = abs(ct_T_dew - ref_T_dew) if (ct_T_dew and ref_T_dew) else None
        
        results.append({
            "z1": z1,
            "ct_T_bub": ct_T_bub, "ref_T_bub": ref_T_bub,
            "ct_T_dew": ct_T_dew, "ref_T_dew": ref_T_dew,
            "d_bub": d_bub, "d_dew": d_dew
        })
        
        def fmt(x): return f"{x:.2f}" if x else "-"
        def fmtd(x): return f"{x:.2e}" if x is not None else "-"
        print(f"{z1:<5.1f} {fmt(ct_T_bub):<10} {fmt(ref_T_bub):<10} {fmt(ct_T_dew):<10} {fmt(ref_T_dew):<10} {fmtd(d_bub):<10} {fmtd(d_dew):<10}")

    # Artifact
    out_file = Path(__file__).parent / "03_txy_pxy_curves.csv"
    with open(out_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nArtifact saved to: {out_file}")

if __name__ == "__main__":
    main()
