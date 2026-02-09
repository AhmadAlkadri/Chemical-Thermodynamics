"""
Slice 4: Multicomponent TP-flash sweep
--------------------------------------
Compares chemthermo vs thermo for a 3-component mixture at random compositions.

Scope:
- System: Methane (1) + Ethane (2) + Propane (3)
- State: Fixed T=200 K, P=20 bar
- Sweep: 20 random compositions (Dirichlet distribution)

Metrics:
- Phase count agreement (1 vs 2 phases)
- Vapor fraction difference
- Composition differences (L2 norm)
"""

import sys
import csv
import random
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

def main():
    # Setup
    comps = ["Methane", "Ethane", "Propane"]
    T = 200.0
    P = 2.0e6
    num_samples = 20
    seed = 42
    
    print(f"Running Validation Slice 4: {' + '.join(comps)}")
    print(f"T={T} K, P={P/1e5} bar, {num_samples} samples")
    
    # Chemthermo
    eos = PengRobinsonEOS()
    
    # Thermo
    ids = [c.lower() for c in comps]
    constants, properties = ChemicalConstantsPackage.from_IDs(ids)
    kijs = [[0.0]*3 for _ in range(3)]
    eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
    gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
    
    # Generate samples
    rng = np.random.default_rng(seed)
    samples = rng.dirichlet(alpha=[1.0, 1.0, 1.0], size=num_samples)
    
    results = []
    
    print(f"{'Sample':<6} {'Phase_CT':<10} {'Phase_TH':<10} {'beta_diff':<10} {'dx_norm':<10} {'dy_norm':<10}")
    print("-" * 70)
    
    for i, zs in enumerate(samples):
        # Thermo
        try:
            ref_res = flasher.flash(T=T, P=P, zs=zs.tolist())
            ref_vf = ref_res.VF
            ref_phases = 2 if (0 < ref_vf < 1) else 1
            if ref_phases == 2:
                ref_x = np.array(ref_res.liquid0.zs)
                ref_y = np.array(ref_res.gas.zs)
            else:
                ref_x = None
                ref_y = None
        except Exception:
            ref_vf = None
            ref_phases = -1
            ref_x, ref_y = None, None

        # Chemthermo
        try:
            mix = ct.Mixture.from_database(comps, zs)
            ct_res = ct.flash_tp(mix, temperature_K=T, pressure_Pa=P, eos=eos)
            ct_vf = ct_res.vapor_fraction
            
            ct_phases = 0
            if ct_vf is not None:
                if 0 < ct_vf < 1: ct_phases = 2
                else: ct_phases = 1
            
            ct_x, ct_y = None, None
            if ct_phases == 2:
                ct_x = np.array(ct_res.phases["liquid"].composition.fractions)
                ct_y = np.array(ct_res.phases["vapor"].composition.fractions)
                
        except Exception:
            ct_vf = None
            ct_phases = -1
        
        # Metrics
        beta_diff = abs(ct_vf - ref_vf) if (ct_vf is not None and ref_vf is not None) else None
        
        dx_norm = None
        if ct_x is not None and ref_x is not None:
            dx_norm = np.linalg.norm(ct_x - ref_x)
            
        dy_norm = None
        if ct_y is not None and ref_y is not None:
            dy_norm = np.linalg.norm(ct_y - ref_y)
            
        results.append({
            "sample_id": i,
            "z1": zs[0], "z2": zs[1], "z3": zs[2],
            "phases_ct": ct_phases, "phases_th": ref_phases,
            "beta_ct": ct_vf, "beta_th": ref_vf,
            "beta_diff": beta_diff,
            "dx_norm": dx_norm, "dy_norm": dy_norm
        })
        
        def fmtd(x): return f"{x:.2e}" if x is not None else "-"
        print(f"{i:<6} {ct_phases:<10} {ref_phases:<10} {fmtd(beta_diff):<10} {fmtd(dx_norm):<10} {fmtd(dy_norm):<10}")

    # Artifact
    out_file = Path(__file__).parent / "04_multicomponent_flash_sweep.csv"
    with open(out_file, "w", newline="") as f:
        # csv.DictWriter helper?
        # keys need to be consistent. DictWriter requires all keys present?
        # We'll just write what we have, handling missing optional keys if we were being strict.
        # But here results have uniform keys.
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
        
    print(f"\nArtifact saved to: {out_file}")
    
    # Summary
    mismatches = sum(1 for r in results if r["phases_ct"] != r["phases_th"])
    print(f"Phase classification mismatches: {mismatches} / {num_samples}")

if __name__ == "__main__":
    main()
