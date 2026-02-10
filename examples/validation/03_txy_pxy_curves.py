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

import csv
import argparse
import os
import sys
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

def _parse_outdir() -> Path | None:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--outdir", type=Path, default=None)
    args, _ = parser.parse_known_args()

    if args.outdir is not None:
        return args.outdir

    env_outdir = os.environ.get("CHEMTHERMO_OUTDIR")
    if env_outdir:
        return Path(env_outdir)

    return None

def main():
    outdir = _parse_outdir()
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
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=RuntimeWarning, module="chemicals")
                bub_res = flasher.flash(P=P, VF=0.0, zs=zs)
            ref_T_bub = bub_res.T
        except Exception:
            ref_T_bub = None
            
        # Chemthermo
        ct_T_bub, _ = ct.bubble_temperature(mix, P, eos)
        
        # --- Dew Point (y=z) ---
        # Thermo
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=RuntimeWarning, module="chemicals")
                dew_res = flasher.flash(P=P, VF=1.0, zs=zs)
            ref_T_dew = dew_res.T
        except Exception:
            ref_T_dew = None
            
        # Chemthermo
        ct_T_dew, _ = ct.dew_temperature(mix, P, eos)
        
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

    # Save artifact (optional).
    if outdir is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        out_file = outdir / "03_txy_pxy_curves.csv"
        with out_file.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nArtifact saved to: {out_file}")
    else:
        print("\nNo --outdir/CHEMTHERMO_OUTDIR provided; skipping artifact write.")

if __name__ == "__main__":
    main()
