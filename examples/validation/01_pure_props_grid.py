"""
Slice 1: Pure-component property grid
-------------------------------------
Compares chemthermo vs thermo for:
1. Vapor Pressure Psat(T)
2. Compressibility Factor Z(T, P) for saturated vapor/liquid

Scope:
- Component: Methane
- Temperature: 110 K - 180 K (Critical Temp ~190.56 K)
"""

import sys
import math
import csv
from pathlib import Path

import numpy as np
import chemthermo as ct
from chemthermo.models.peng_robinson import PengRobinsonEOS

# Try to import thermo, or exit with a helpful message
try:
    import thermo
    from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashPureVLS
except ImportError:
    print("Error: 'thermo' package not found. Please install it with 'pip install thermo'.")
    sys.exit(1)


def solve_psat_chemthermo(mixture, T, eos, P_guess, tol=1e-5, max_iter=50):
    """
    Find Psat such that fugacity coefficients of liquid and vapor match.
    Uses a simple Newton-Raphson or Secant on ln(phi_L) - ln(phi_V).
    For robustness here, we'll use a simple iterative update based on ratio.
    """
    P = P_guess
    for _ in range(max_iter):
        # We need both phases to exist as potentially stable or metastable roots
        try:
            phi_l = eos.fugacity_coefficients(
                mixture=mixture, temperature_K=T, pressure_Pa=P, 
                composition=mixture.fractions, phase="liquid"
            )[0]
            phi_v = eos.fugacity_coefficients(
                mixture=mixture, temperature_K=T, pressure_Pa=P, 
                composition=mixture.fractions, phase="vapor"
            )[0]
        except Exception:
            # If we fail to find roots (e.g. strict single phase region), we can't find Psat easily
            return None
        
        # At equilibrium: phi_L * P * x = phi_V * P * y. Pure => x=y=1 => phi_L = phi_V.
        # K = phi_L / phi_V. We want K = 1.
        ratio = phi_l / phi_v
        
        if abs(ratio - 1.0) < tol:
            return P
        
        # Simple update: P_new = P_old * ratio
        # Rationale: If phi_L > phi_V, liquid is more stable (lower energetic cost?), 
        # normally roughly P_sat ~= P_current * (phi_l/phi_v) for some EOS.
        # Actually standard successive substitution is P_new = P_old * (phi_l / phi_v)
        P = P * ratio
        
    return None

def main():
    # Setup Component
    name = "Methane"
    # Methane critical properties (approx): Tc=190.56 K, Pc=4.599 MPa
    mixture = ct.Mixture.from_database([name], [1.0])
    eos = PengRobinsonEOS()
    
    # Thermo setup
    # We use PRMIX to strictly match Peng-Robinson behavior if possible,
    # or just use Chemical for a broad comparison. 
    # To be "apples to apples" with chemthermo PR, we should use thermo's PR.
    constants, properties = ChemicalConstantsPackage.from_IDs([name])
    eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
    gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    flasher = FlashPureVLS(constants, properties, gas=gas, liquids=[liquid], solids=[])
    
    # Grid
    temperatures = np.linspace(110.0, 180.0, 15)
    
    # Storage
    results = []
    
    print(f"Running Validation Slice 1: {name}")
    print(f"{'T (K)':<10} {'Type':<10} {'Thermo':<15} {'ChemThermo':<15} {'RelDiff':<15}")
    print("-" * 70)
    
    for T in temperatures:
        # --- 1. Vapor Pressure Psat ---
        # Thermo
        try:
            ref_flash = flasher.flash(T=T, VF=0.5)
            ref_Psat = ref_flash.P
        except Exception:
            ref_Psat = None

        # Chemthermo
        # Initial guess from Antoine or just loose correlation? 
        # For Methane @ 110K ~ 0.8 bar. @ 180K ~ 40 bar.
        # Let's use ref_Psat as guess if available to speed up, or a simple guess
        guess_P = ref_Psat if ref_Psat else 1e5
        calc_Psat = solve_psat_chemthermo(mixture, T, eos, guess_P)
        
        if ref_Psat and calc_Psat:
            err = abs(calc_Psat - ref_Psat) / ref_Psat
            results.append({
                "T": T, "Property": "Psat", 
                "Thermo": ref_Psat, "ChemThermo": calc_Psat, "RelError": err
            })
            print(f"{T:<10.1f} {'Psat':<10} {ref_Psat:<15.4e} {calc_Psat:<15.4e} {err:<15.2e}")
        else:
            print(f"{T:<10.1f} {'Psat':<10} {'FAIL':<15} {'FAIL':<15} {'-':<15}")

            # Z factor
            # Use the Z from the flash result directly
            if ref_Psat and calc_Psat and ref_flash:
                P_check = ref_Psat 
                
                # Liquid Z
                try:
                    ref_Zl = ref_flash.liquid0.Z
                    calc_Zl = eos.compressibility_factor(
                        mixture=mixture, temperature_K=T, pressure_Pa=P_check, 
                        composition=[1.0], phase="liquid"
                    )
                    err_l = abs(calc_Zl - ref_Zl) / abs(ref_Zl)
                    results.append({
                        "T": T, "Property": "Z_liq", 
                        "Thermo": ref_Zl, "ChemThermo": calc_Zl, "RelError": err_l
                    })
                except Exception:
                    pass

                # Vapor Z
                try:
                    ref_Zv = ref_flash.gas.Z
                    calc_Zv = eos.compressibility_factor(
                        mixture=mixture, temperature_K=T, pressure_Pa=P_check, 
                        composition=[1.0], phase="vapor"
                    )
                    err_v = abs(calc_Zv - ref_Zv) / abs(ref_Zv)
                    results.append({
                        "T": T, "Property": "Z_vap", 
                        "Thermo": ref_Zv, "ChemThermo": calc_Zv, "RelError": err_v
                    })
                except Exception:
                    pass

    # Save artifact
    out_file = Path(__file__).parent / "01_pure_props_grid.csv"
    with open(out_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["T", "Property", "Thermo", "ChemThermo", "RelError"])
        writer.writeheader()
        writer.writerows(results)
    
    # Summary
    if not results:
        print("No successful comparisons.")
        sys.exit(1)
        
    errs = np.array([r["RelError"] for r in results])
    print("\nSummary Statistics (Relative Error):")
    print(f"Count:  {len(errs)}")
    print(f"Max:    {np.max(errs):.2e}")
    print(f"Mean:   {np.mean(errs):.2e}")
    print(f"Median: {np.median(errs):.2e}")
    print(f"RMS:    {np.sqrt(np.mean(errs**2)):.2e}")
    
    print(f"\nArtifact saved to: {out_file}")

if __name__ == "__main__":
    main()
