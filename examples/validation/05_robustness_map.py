"""
Slice 5: Robustness / failure-mode map
--------------------------------------
Maps convergence and phase stability across a T-P grid for a binary mixture.

Scope:
- System: Methane (1) + Ethane (2)
- Composition: Fixed z = [0.5, 0.5]
- Grid: T=[100, 300] (10 pts), P=[1e5, 1e7] (10 pts) -> 100 points

Metrics:
- Convergence status (Success/Fail) for both
- Phase count agreement
- Classification buckets: MATCH, PHASE_MISMATCH, CT_FAIL, TH_FAIL, BOTH_FAIL
"""

import sys
import csv
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
    comps = ["Methane", "Ethane"]
    zs = [0.5, 0.5]
    
    T_range = np.linspace(100.0, 300.0, 10)
    P_range = np.geomspace(1.0e5, 1.0e7, 10) # 1 bar to 100 bar
    
    print(f"Running Validation Slice 5: Robustness Map")
    print(f"System: {' + '.join(comps)} @ z={zs}")
    print(f"Grid: 10x10 points (T: 100-300 K, P: 1-100 bar)")
    
    # Chemthermo
    eos = PengRobinsonEOS()
    mix = ct.Mixture.from_database(comps, zs)
    
    # Thermo
    ids = [c.lower() for c in comps]
    constants, properties = ChemicalConstantsPackage.from_IDs(ids)
    kijs = [[0.0]*2 for _ in range(2)]
    eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas, 'kijs': kijs}
    gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
    
    results = []
    
    counts = {
        "MATCH": 0,
        "PHASE_MISMATCH": 0,
        "CT_FAIL": 0,
        "TH_FAIL": 0,
        "BOTH_FAIL": 0
    }
    
    print(f"{'T(K)':<8} {'P(bar)':<8} {'Status':<15} {'Phase_CT':<10} {'Phase_TH':<10}")
    print("-" * 60)
    
    for T in T_range:
        for P in P_range:
            # Thermo
            th_status = "OK"
            th_phases = -1
            try:
                ref_res = flasher.flash(T=T, P=P, zs=zs)
                th_phases = 2 if (0 < ref_res.VF < 1) else 1
            except Exception:
                th_status = "FAIL"
                
            # Chemthermo
            ct_status = "OK"
            ct_phases = -1
            try:
                ct_res = ct.flash_tp(mix, temperature_K=T, pressure_Pa=P, eos=eos)
                if ct_res.vapor_fraction is not None:
                     ct_phases = 2 if (0 < ct_res.vapor_fraction < 1) else 1
                else: 
                     # Should not happen for successful flash_tp return type unless I misinterpreted
                     # actually flash_tp always sets vapor_fraction (0 or 1 for single phase)
                     # Wait, single phase fallback sets vf=0 or 1.
                     pass 
            except Exception:
                ct_status = "FAIL"
            
            # Classification
            outcome = "UNKNOWN"
            if ct_status == "FAIL" and th_status == "FAIL":
                outcome = "BOTH_FAIL"
            elif ct_status == "FAIL":
                outcome = "CT_FAIL"
            elif th_status == "FAIL":
                outcome = "TH_FAIL"
            else:
                if ct_phases == th_phases:
                    outcome = "MATCH"
                else:
                    outcome = "PHASE_MISMATCH"
            
            counts[outcome] += 1
            
            results.append({
                "T": T, "P": P,
                "Outcome": outcome,
                "CT_Status": ct_status, "TH_Status": th_status,
                "CT_Phases": ct_phases, "TH_Phases": th_phases
            })
            
            if outcome != "MATCH":
                print(f"{T:<8.1f} {P/1e5:<8.1f} {outcome:<15} {ct_phases:<10} {th_phases:<10}")

    # Artifact
    out_file = Path(__file__).parent / "05_robustness_map.csv"
    with open(out_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    
    print("-" * 60)
    print(f"Artifact saved to: {out_file}")
    print("\nSummary Counts:")
    for k, v in counts.items():
        print(f"  {k:<15}: {v}")

if __name__ == "__main__":
    main()
