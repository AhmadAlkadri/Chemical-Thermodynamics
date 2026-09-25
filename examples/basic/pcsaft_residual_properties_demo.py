"""PC-SAFT residual enthalpy and entropy along a saturated-liquid-like state (ADR-0034).

Prints the reduced residual properties of liquid and vapour n-hexane at 300 K
and the difference of the residual enthalpies between the two roots at the
saturation pressure - which is the enthalpy of vaporization, because the
ideal-gas parts cancel at equal T and P. No optional dependency.
"""

from __future__ import annotations

from chemthermo.eos import PCSAFTEOS
from chemthermo.eos.pcsaft import R_J_PER_MOL_K

TEMPERATURE_K = 300.0
#: n-hexane's PC-SAFT saturation pressure at 300 K (ledger Case P-2).
P_SAT_PA = 21858.084278856164


def main() -> None:
    eos = PCSAFTEOS(components=("n-Hexane",))
    vapor_rho, liquid_rho = eos.density_roots(
        temperature_K=TEMPERATURE_K, pressure_Pa=P_SAT_PA, composition=[1.0]
    )
    rt = R_J_PER_MOL_K * TEMPERATURE_K
    rows = {}
    for label, rho in (("vapor", vapor_rho), ("liquid", liquid_rho)):
        props = eos.residual_properties(
            temperature_K=TEMPERATURE_K, density_mol_m3=rho, composition=[1.0]
        )
        rows[label] = props
        print(
            f"{label:<6} rho={rho:12.4f} mol/m^3  Z={props['z']:.6f}  "
            f"H^res={props['h_res'] * rt / 1000.0:9.4f} kJ/mol  "
            f"S^res(T,P)={props['s_res_tp'] * R_J_PER_MOL_K:9.4f} J/(mol K)"
        )
    dh_vap = (rows["vapor"]["h_res"] - rows["liquid"]["h_res"]) * rt / 1000.0
    ds_vap = (rows["vapor"]["s_res_tp"] - rows["liquid"]["s_res_tp"]) * R_J_PER_MOL_K
    print(f"PC-SAFT enthalpy of vaporization at {TEMPERATURE_K} K: {dh_vap:.4f} kJ/mol")
    # At saturation dG = 0, so dH_vap = T dS_vap: a consistency check of the pair.
    print(f"T * dS_vap: {TEMPERATURE_K * ds_vap / 1000.0:.4f} kJ/mol (equal at saturation)")


if __name__ == "__main__":
    main()
