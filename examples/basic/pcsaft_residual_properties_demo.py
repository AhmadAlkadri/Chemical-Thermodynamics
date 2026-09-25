"""PC-SAFT residual enthalpy and entropy at saturation (ADR-0034).

Prints the reduced residual properties of the liquid and vapour roots at the
saturation pressure, and their difference - the enthalpy of vaporization,
because the ideal-gas parts cancel at equal T and P - for n-hexane at 300 K
(non-associating) and water at 373.15 K (2B association, the C2 amendment).
No optional dependency.
"""

from __future__ import annotations

from chemthermo.eos import PCSAFTEOS
from chemthermo.eos.pcsaft import R_J_PER_MOL_K

TEMPERATURE_K = 300.0
#: n-hexane's PC-SAFT saturation pressure at 300 K (ledger Case P-2).
P_SAT_PA = 21858.084278856164
WATER_T_K = 373.15


def _saturation_pressure(eos: PCSAFTEOS, temperature: float, guess: float) -> float:
    """Secant on P for ln phi(vapour root) = ln phi(liquid root)."""

    def mismatch(pressure: float) -> float:
        roots = eos.density_roots(
            temperature_K=temperature, pressure_Pa=pressure, composition=[1.0]
        )
        ln_phi = [
            eos.ln_fugacity_coefficients(
                temperature_K=temperature, density_mol_m3=rho, composition=[1.0]
            )[0]
            for rho in (roots[0], roots[-1])
        ]
        return float(ln_phi[0] - ln_phi[1])

    p0, p1 = 0.9 * guess, 1.1 * guess
    f0, f1 = mismatch(p0), mismatch(p1)
    for _ in range(50):
        if abs(f1) < 1e-14 or f1 == f0:
            break
        p0, p1, f0 = p1, p1 - f1 * (p1 - p0) / (f1 - f0), f1
        f1 = mismatch(p1)
    return p1


def main() -> None:
    report("n-Hexane", TEMPERATURE_K, P_SAT_PA)
    water = PCSAFTEOS(components=("Water",))
    p_sat = _saturation_pressure(water, WATER_T_K, 1.0e5)
    print(f"\nPC-SAFT water saturation pressure at {WATER_T_K} K: {p_sat:.3f} Pa")
    report("Water", WATER_T_K, p_sat)
    print(
        "Remark, not a check: water's measured enthalpy of vaporization at its normal "
        "boiling point is about 40.65 kJ/mol."
    )


def report(name: str, temperature: float, pressure: float) -> None:
    eos = PCSAFTEOS(components=(name,))
    roots = eos.density_roots(temperature_K=temperature, pressure_Pa=pressure, composition=[1.0])
    vapor_rho, liquid_rho = roots[0], roots[-1]
    rt = R_J_PER_MOL_K * temperature
    rows = {}
    for label, rho in (("vapor", vapor_rho), ("liquid", liquid_rho)):
        props = eos.residual_properties(
            temperature_K=temperature, density_mol_m3=rho, composition=[1.0]
        )
        rows[label] = props
        print(
            f"{label:<6} rho={rho:12.4f} mol/m^3  Z={props['z']:.6f}  "
            f"H^res={props['h_res'] * rt / 1000.0:9.4f} kJ/mol  "
            f"S^res(T,P)={props['s_res_tp'] * R_J_PER_MOL_K:9.4f} J/(mol K)"
        )
    dh_vap = (rows["vapor"]["h_res"] - rows["liquid"]["h_res"]) * rt / 1000.0
    ds_vap = (rows["vapor"]["s_res_tp"] - rows["liquid"]["s_res_tp"]) * R_J_PER_MOL_K
    print(f"{name}: PC-SAFT enthalpy of vaporization at {temperature} K: {dh_vap:.4f} kJ/mol")
    # At saturation dG = 0, so dH_vap = T dS_vap: a consistency check of the pair.
    print(f"T * dS_vap: {temperature * ds_vap / 1000.0:.4f} kJ/mol (equal at saturation)")


if __name__ == "__main__":
    main()
