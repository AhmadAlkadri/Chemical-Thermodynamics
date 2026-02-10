"""Minimal deterministic external reference check against thermo."""

from __future__ import annotations

import numpy as np

import chemthermo as ct

try:
    import thermo
except ImportError:
    print("Error: 'thermo' package is required. Install with: pip install -e '.[validation]'")
    raise SystemExit(1)


CEOSGas = thermo.CEOSGas
CEOSLiquid = thermo.CEOSLiquid
ChemicalConstantsPackage = thermo.ChemicalConstantsPackage
FlashVL = thermo.FlashVL
PRMIX = thermo.PRMIX


def _get_vf(state: object) -> float:
    vf = getattr(state, "VF", None)
    if vf is None:
        raise RuntimeError("Reference state does not expose VF.")
    value = vf() if callable(vf) else vf
    return float(value)


def _reference_flasher(ids: list[str]) -> object:
    constants, properties = ChemicalConstantsPackage.from_IDs(ids)
    kijs = [[0.0 for _ in ids] for _ in ids]
    eos_kwargs = {
        "Pcs": constants.Pcs,
        "Tcs": constants.Tcs,
        "omegas": constants.omegas,
        "kijs": kijs,
    }
    gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    return FlashVL(constants, properties, liquid=liquid, gas=gas)


def main() -> int:
    names = ["Methane", "Ethane", "Propane"]
    zs = [0.50, 0.30, 0.20]
    temperature_K = 240.0
    pressure_Pa = 3.0e6

    mixture = ct.Mixture.from_database(names, zs, normalize=True)
    result = ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=ct.PengRobinsonEOS(),
    )

    if "liquid" not in result.phases or "vapor" not in result.phases:
        print("FAIL: chemthermo did not produce a two-phase result for the reference case.")
        return 1

    ref = _reference_flasher([name.lower() for name in names]).flash(
        T=temperature_K,
        P=pressure_Pa,
        zs=zs,
    )
    ref_vf = _get_vf(ref)
    if not (0.0 < ref_vf < 1.0):
        print("FAIL: thermo did not produce a two-phase result for the reference case.")
        return 1

    beta = float(result.vapor_fraction)
    x = np.array(result.phases["liquid"].composition.fractions)
    y = np.array(result.phases["vapor"].composition.fractions)

    ref_x = np.array(ref.liquid0.zs)
    ref_y = np.array(ref.gas.zs)

    beta_diff = abs(beta - ref_vf)
    x_diff = np.abs(x - ref_x)
    y_diff = np.abs(y - ref_y)

    beta_tol = 5e-3
    comp_atol = 1e-2
    comp_rtol = 5e-2

    print("Reference validation (chemthermo vs thermo)")
    print(f"Case: names={names}, T={temperature_K:.1f} K, P={pressure_Pa:.3e} Pa")
    print(f"beta chemthermo={beta:.8f}, beta thermo={ref_vf:.8f}, |delta|={beta_diff:.3e}")
    print(f"max |x - x_ref| = {float(np.max(x_diff)):.3e}")
    print(f"max |y - y_ref| = {float(np.max(y_diff)):.3e}")

    ok_beta = beta_diff <= beta_tol
    ok_x = np.allclose(x, ref_x, rtol=comp_rtol, atol=comp_atol)
    ok_y = np.allclose(y, ref_y, rtol=comp_rtol, atol=comp_atol)

    if ok_beta and ok_x and ok_y:
        print("PASS: reference case is within tolerances.")
        return 0

    print("FAIL: reference case exceeded tolerance thresholds.")
    print(f"Thresholds: beta_tol={beta_tol}, comp_rtol={comp_rtol}, comp_atol={comp_atol}")
    return 1


if __name__ == "__main__":
    exit_code = main()
    if exit_code != 0:
        raise SystemExit(exit_code)
