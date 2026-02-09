"""
Slice 2: Binary K-value surface (TP flash)
------------------------------------------
Compares chemthermo vs thermo for a binary mixture TP flash.
Sweeps composition z1 from 0.05 to 0.95.

Scope:
- System: Methane (1) + Ethane (2)
- State: Fixed T=160 K, P=10 bar (Chosen to likely be in VLE region)
- Sweep: z1 in [0.05, ..., 0.95]

Metrics:
- Vapor fraction beta
- K-values K1, K2
"""

import sys
import csv
from pathlib import Path

import numpy as np
import chemthermo as ct
from chemthermo.models.peng_robinson import PengRobinsonEOS

# Limit imports to what we need
try:
    import thermo
    from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, PRMIX, FlashVL
except ImportError:
    print("Error: 'thermo' package not found.")
    sys.exit(1)

def main():
    # Setup
    comps = ["Methane", "Ethane"]
    T = 160.0
    P = 1.0e6  # 10 bar
    
    # Chemthermo setup
    eos = PengRobinsonEOS()
    
    # Thermo setup
    ids = [c.lower() for c in comps]
    constants, properties = ChemicalConstantsPackage.from_IDs(ids)
    # Ensure no kij interaction parameters for base comparison (chemthermo defaults to 0)
    kijs = [[0.0]*2 for _ in range(2)]
    eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
    
    gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
    
    # Sweep
    z1_values = np.linspace(0.05, 0.95, 19)
    results = []
    
    print(f"Running Validation Slice 2: {comps[0]} + {comps[1]}")
    print(f"T={T} K, P={P/1e5} bar")
    print(f"{'z1':<8} {'beta_ct':<10} {'beta_th':<10} {'K1_ct':<10} {'K1_th':<10} {'K2_ct':<10} {'K2_th':<10} {'MaxDiff':<10}")
    print("-" * 90)
    
    for z1 in z1_values:
        z2 = 1.0 - z1
        zs = [z1, z2]
        
        # --- Thermo ---
        try:
            ref_res = flasher.flash(T=T, P=P, zs=zs)
            ref_beta = ref_res.VF
            ref_x = ref_res.liquid0.zs
            ref_y = ref_res.gas.zs
            ref_K = [y/x for y, x in zip(ref_y, ref_x)]
        except Exception as e:
            ref_beta = None
            ref_K = [None, None]
            
        # --- Chemthermo ---
        try:
            mix = ct.Mixture.from_database(comps, zs)
            ct_res = ct.flash_tp(mix, temperature_K=T, pressure_Pa=P, eos=eos)
            
            ct_beta = ct_res.vapor_fraction
            if ct_beta is None:
                # Single phase
                if "liquid" in ct_res.phases: ct_beta = 0.0
                elif "vapor" in ct_res.phases: ct_beta = 1.0
                
            # Extract K-values if VLE
            if ct_res.vapor_fraction is not None and 0 < ct_res.vapor_fraction < 1:
                ct_x = ct_res.phases["liquid"].composition.fractions
                ct_y = ct_res.phases["vapor"].composition.fractions
                ct_K = [y/x for y, x in zip(ct_y, ct_x)]
            else:
                ct_K = [None, None]
                
        except Exception as e:
            ct_beta = None
            ct_K = [None, None]
            
        # Compare
        row = {
            "z1": z1,
            "ct_beta": ct_beta, "ref_beta": ref_beta,
            "ct_K1": ct_K[0], "ref_K1": ref_K[0],
            "ct_K2": ct_K[1], "ref_K2": ref_K[1]
        }
        
        # Metrics
        diffs = []
        if ct_beta is not None and ref_beta is not None:
            diffs.append(abs(ct_beta - ref_beta))
        if ct_K[0] and ref_K[0]:
            diffs.append(abs(ct_K[0] - ref_K[0]))
            diffs.append(abs(ct_K[1] - ref_K[1]))
        
        max_diff = max(diffs) if diffs else 0.0
        row["max_diff"] = max_diff
        results.append(row)
        
        # Print
        def fmt(x): return f"{x:.4f}" if x is not None else "-"
        print(f"{z1:<8.2f} {fmt(ct_beta):<10} {fmt(ref_beta):<10} {fmt(ct_K[0]):<10} {fmt(ref_K[0]):<10} {fmt(ct_K[1]):<10} {fmt(ref_K[1]):<10} {max_diff:<10.2e}")

    # Save artifact
    out_file = Path(__file__).parent / "02_binary_k_surface.csv"
    with open(out_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
        
    print(f"\nArtifact saved to: {out_file}")
    
    # Summary assertions could be here, but we just want the artifact for now.
    worst_case = max(r["max_diff"] for r in results)
    print(f"Worst-case disagreement: {worst_case:.2e}")

if __name__ == "__main__":
    main()
