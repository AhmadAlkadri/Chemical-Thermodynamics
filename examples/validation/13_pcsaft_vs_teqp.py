"""PC-SAFT against teqp (NIST), with PASS/FAIL on every state.

What is being validated
-----------------------
``chemthermo.eos.PCSAFTEOS`` and ``teqp`` implement the same published model -
Gross & Sadowski, Ind. Eng. Chem. Res. 40 (2001) 1244, non-associating - from
the same 42 universal constants, but they do not share a line of derivative
code:

  * teqp writes one ``alphar`` and obtains **every** derivative by automatic
    differentiation (generalized complex step / autodiff types);
  * chemthermo writes ``A^res/RT`` and then **analytic** density and
    composition derivatives, assembled by chain rule over the ``zeta`` moments,
    ``mbar``, ``m2es3`` and ``m2e2s3``.

So agreement on ``A^res/RT`` tests the model, and agreement on ``Z`` and
``ln phi_i`` tests the hand-written derivatives, which is where mistakes
actually live.

What is checked here
--------------------
1. ``A^res/RT``, ``Z``, ``P`` and ``ln phi_i`` at fourteen states: pure
   n-hexane at four ``(T, rho)``; binaries at gas-like and liquid-like
   densities; a nonzero per-pair ``k_ij`` (methane / n-decane, 0.03); two
   ternaries; and one state inside the spinodal where ``Z < 0``, where
   chemthermo refuses to produce fugacity coefficients rather than returning a
   ``nan``.
2. The internal identity ``sum_i x_i ln phi_i = A^res/RT + Z - 1 - ln Z``,
   which follows from Euler's theorem and must hold to round-off.
3. Pure-component **saturation** at 300 K and 400 K: a small bisection density
   root finder written in this script (chemthermo ships no density solver yet)
   against ``teqp``'s own ``pure_VLE_T`` Newton solve.
4. A negative control: perturbing one ``sigma`` by 1 % must break the
   agreement, so the comparison cannot be passing vacuously.

Model versus experiment
-----------------------
The 300 K saturation pressure is also printed next to the commonly tabulated
n-hexane vapour pressure of about 21.7 kPa. That number was **not** verified
against a primary source here and is a remark, not an assertion: PC-SAFT with
the published parameters is being compared to another implementation of
PC-SAFT, not to experiment.

Requires the optional ``teqp`` dependency::

    pip install -e ".[validation]"

The script prints a message and exits 0 when it is missing.
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

import numpy as np

try:
    import teqp
except ImportError:  # pragma: no cover - exercised only without the extra
    teqp = None  # type: ignore[assignment]

from chemthermo import ModelError
from chemthermo.eos import PCSAFTEOS
from chemthermo.eos.pcsaft import R_J_PER_MOL_K

# Gross & Sadowski (2001) Table 1; written out here so the reference model is
# built from this file rather than from the package under test.
PARAMETERS: dict[str, tuple[float, float, float]] = {
    "Methane": (1.0000, 3.7039, 150.03),
    "Nitrogen": (1.2053, 3.3130, 90.96),
    "Carbon dioxide": (2.0729, 2.7852, 169.21),
    "n-Hexane": (3.0576, 3.7983, 236.77),
    "n-Decane": (4.6627, 3.8384, 243.87),
}

TOLERANCE = 1e-10

STATES: list[tuple[str, tuple[str, ...], object, list[float], float, float]] = [
    ("C6 pure, 300 K, gas", ("n-Hexane",), 0.0, [1.0], 300.0, 100.0),
    ("C6 pure, 300 K, liquid", ("n-Hexane",), 0.0, [1.0], 300.0, 7700.0),
    ("C6 pure, 400 K, liquid", ("n-Hexane",), 0.0, [1.0], 400.0, 6800.0),
    ("C6 pure, 500 K, dense", ("n-Hexane",), 0.0, [1.0], 500.0, 3000.0),
    ("C1/C6 0.5/0.5, 300 K, gas", ("Methane", "n-Hexane"), 0.0, [0.5, 0.5], 300.0, 200.0),
    (
        "C1/C6 0.5/0.5, 300 K, liquid",
        ("Methane", "n-Hexane"),
        0.0,
        [0.5, 0.5],
        300.0,
        11000.0,
    ),
    ("C1/C6 0.2/0.8, 450 K, liquid", ("Methane", "n-Hexane"), 0.0, [0.2, 0.8], 450.0, 8000.0),
    (
        "C1/C6 0.2/0.8, 300 K, spinodal",
        ("Methane", "n-Hexane"),
        0.0,
        [0.2, 0.8],
        300.0,
        8000.0,
    ),
    (
        "C1/C10 kij=0.03, 350 K, gas",
        ("Methane", "n-Decane"),
        {("Methane", "n-Decane"): 0.03},
        [0.3, 0.7],
        350.0,
        100.0,
    ),
    (
        "C1/C10 kij=0.03, 350 K, liquid",
        ("Methane", "n-Decane"),
        {("Methane", "n-Decane"): 0.03},
        [0.3, 0.7],
        350.0,
        6500.0,
    ),
    ("N2/C1 0.4/0.6, 150 K, liquid", ("Nitrogen", "Methane"), 0.0, [0.4, 0.6], 150.0, 20000.0),
    (
        "CO2/C10 0.4/0.6, 320 K, liquid",
        ("Carbon dioxide", "n-Decane"),
        0.0,
        [0.4, 0.6],
        320.0,
        8000.0,
    ),
    (
        "C1/C6/N2 0.3/0.4/0.3, 250 K, gas",
        ("Methane", "n-Hexane", "Nitrogen"),
        0.0,
        [0.3, 0.4, 0.3],
        250.0,
        500.0,
    ),
    (
        "C1/C6/N2 0.1/0.7/0.2, 250 K, liquid",
        ("Methane", "n-Hexane", "Nitrogen"),
        0.0,
        [0.1, 0.7, 0.2],
        250.0,
        10000.0,
    ),
]


def build_reference(components: Sequence[str], kij_matrix: np.ndarray):
    coefficients = [
        {
            "name": name,
            "m": PARAMETERS[name][0],
            "sigma_Angstrom": PARAMETERS[name][1],
            "epsilon_over_k": PARAMETERS[name][2],
            "BibTeXKey": "Gross-IECR-2001",
        }
        for name in components
    ]
    return teqp.make_model(
        {"kind": "PCSAFT", "model": {"coeffs": coefficients, "kmat": kij_matrix.tolist()}}
    )


def bisect(f: Callable[[float], float], lo: float, hi: float, iterations: int = 200) -> float:
    f_lo, f_hi = f(lo), f(hi)
    if f_lo == 0.0:
        return lo
    if f_hi == 0.0:
        return hi
    if f_lo * f_hi > 0.0:
        raise ValueError(f"bracket does not straddle a root: {f_lo!r}, {f_hi!r}")
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if f_mid == 0.0:
            return mid
        if f_lo * f_mid < 0.0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)


def pure_saturation(
    eos: PCSAFTEOS, temperature: float, density_max: float = 9500.0
) -> tuple[float, float, float]:
    """``(Psat, rho_liquid, rho_vapour)`` from equal fugacity on the two roots.

    Written in this script on purpose: chemthermo does not ship a density
    solver yet (that is the next slice). Scan the isotherm for the spinodal
    extrema, bisect ``P(rho) = P`` on each branch, then bisect on
    ``ln phi(liquid) = ln phi(vapour)``.
    """

    def pressure(rho: float) -> float:
        return eos.pressure_Pa(temperature_K=temperature, density_mol_m3=rho, composition=[1.0])

    def ln_phi(rho: float) -> float:
        return eos.ln_fugacity_coefficients(
            temperature_K=temperature, density_mol_m3=rho, composition=[1.0]
        )[0]

    grid = np.geomspace(1e-2, density_max, 6000)
    curve = np.array([pressure(float(rho)) for rho in grid])
    top = int(np.argmax(np.where(grid < 0.6 * density_max, curve, -np.inf)))
    bottom = top + int(np.argmin(curve[top:]))

    def roots(p: float) -> tuple[float, float]:
        rho_v = bisect(lambda rho: pressure(rho) - p, 1e-8, float(grid[top]))
        rho_l = bisect(lambda rho: pressure(rho) - p, float(grid[bottom]), density_max)
        return rho_l, rho_v

    def residual(p: float) -> float:
        rho_l, rho_v = roots(p)
        return ln_phi(rho_l) - ln_phi(rho_v)

    p_sat = bisect(residual, max(float(curve[bottom]), 1e-3), 0.999999 * float(curve[top]))
    rho_l, rho_v = roots(p_sat)
    return p_sat, rho_l, rho_v


def main() -> None:
    if teqp is None:
        print("teqp is not installed; skipping the PC-SAFT cross-check.")
        print('Install it with:  pip install -e ".[validation]"')
        return

    print(f"PC-SAFT (Gross & Sadowski 2001, non-associating) vs teqp {teqp.__version__}")
    print("chemthermo: analytic derivatives. teqp: automatic differentiation.")
    print(f"tolerance on A^res/RT, Z, P and ln phi: {TOLERANCE:.0e} (relative and absolute)\n")

    failures: list[str] = []
    worst_a_res = 0.0
    worst_z = 0.0
    worst_ln_phi = 0.0
    worst_identity = 0.0

    def record(label: str, ok: bool) -> None:
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
        if not ok:
            failures.append(label)

    print("=" * 78)
    print("1) Residual properties at 14 states")
    print("=" * 78)
    header = f"  {'state':<34}{'|dA^res|':>12}{'|dZ|':>12}{'max|dlnphi|':>14}"
    print(header)
    for label, components, kij, x, temperature, density in STATES:
        eos = PCSAFTEOS(components=components, kij=kij)  # type: ignore[arg-type]
        model = build_reference(components, eos.kij_matrix())
        z_array = np.array(x, dtype=float)

        reference_a_res = float(model.get_Ar00(temperature, density, z_array))
        reference_z = 1.0 + float(model.get_Ar01(temperature, density, z_array))

        a_res = eos.residual_helmholtz(
            temperature_K=temperature, volume_m3=1.0 / density, composition=x
        )
        z_factor = eos.compressibility_factor(
            temperature_K=temperature, density_mol_m3=density, composition=x
        )
        pressure = eos.pressure_Pa(temperature_K=temperature, density_mol_m3=density, composition=x)

        d_a_res = abs(a_res - reference_a_res)
        d_z = abs(z_factor - reference_z)
        d_p = abs(pressure - reference_z * density * R_J_PER_MOL_K * temperature)
        worst_a_res = max(worst_a_res, d_a_res)
        worst_z = max(worst_z, d_z)

        ok = d_a_res <= TOLERANCE * max(1.0, abs(reference_a_res))
        ok = ok and d_z <= TOLERANCE * max(1.0, abs(reference_z))
        ok = ok and d_p <= TOLERANCE * max(1.0, abs(pressure))

        if reference_z > 0.0:
            reference_ln_phi = np.log(
                np.asarray(model.get_fugacity_coefficients(temperature, density * z_array))
            )
            ln_phi = np.array(
                eos.ln_fugacity_coefficients(
                    temperature_K=temperature, density_mol_m3=density, composition=x
                )
            )
            d_ln_phi = float(np.max(np.abs(ln_phi - reference_ln_phi)))
            worst_ln_phi = max(worst_ln_phi, d_ln_phi)
            ok = ok and d_ln_phi <= TOLERANCE * max(1.0, float(np.max(np.abs(reference_ln_phi))))

            identity = float(z_array @ ln_phi)
            expected = a_res + z_factor - 1.0 - math.log(z_factor)
            worst_identity = max(worst_identity, abs(identity - expected))
            ok = ok and abs(identity - expected) < 1e-12
            shown = f"{d_ln_phi:>14.2e}"
        else:
            try:
                eos.ln_fugacity_coefficients(
                    temperature_K=temperature, density_mol_m3=density, composition=x
                )
            except ModelError:
                shown = f"{'refused (Z<0)':>14}"
            else:
                shown = f"{'NOT REFUSED':>14}"
                ok = False

        print(f"  {label:<34}{d_a_res:>12.2e}{d_z:>12.2e}{shown}")
        if not ok:
            failures.append(label)

    print(
        f"\n  worst |dA^res/RT| = {worst_a_res:.2e}, worst |dZ| = {worst_z:.2e}, "
        f"worst max|d ln phi| = {worst_ln_phi:.2e}"
    )
    print(f"  worst |sum_i x_i ln phi_i - (A^res/RT + Z - 1 - ln Z)| = {worst_identity:.2e}")
    record("every state agrees with teqp inside the tolerance", not failures)

    print("\n" + "=" * 78)
    print("2) Negative control: a 1 % change in sigma must break the agreement")
    print("=" * 78)
    perturbed = teqp.make_model(
        {
            "kind": "PCSAFT",
            "model": {
                "coeffs": [
                    {
                        "name": "n-Hexane",
                        "m": PARAMETERS["n-Hexane"][0],
                        "sigma_Angstrom": PARAMETERS["n-Hexane"][1] * 1.01,
                        "epsilon_over_k": PARAMETERS["n-Hexane"][2],
                        "BibTeXKey": "perturbed",
                    }
                ],
                "kmat": [[0.0]],
            },
        }
    )
    hexane = PCSAFTEOS(components=("n-Hexane",))
    gap = abs(
        hexane.residual_helmholtz(temperature_K=300.0, volume_m3=1.0 / 7700.0, composition=[1.0])
        - float(perturbed.get_Ar00(300.0, 7700.0, np.array([1.0])))
    )
    print(f"  |dA^res/RT| with sigma * 1.01 = {gap:.3e}  (tolerance {TOLERANCE:.0e})")
    record("the comparison is not vacuous", gap > 1e-3)

    print("\n" + "=" * 78)
    print("3) Pure n-hexane saturation: bisection here vs teqp's pure_VLE_T")
    print("=" * 78)
    for temperature, guesses in ((300.0, (7700.0, 10.0)), (400.0, (6800.0, 150.0))):
        model = build_reference(("n-Hexane",), hexane.kij_matrix())
        rho_l_ref, rho_v_ref = model.pure_VLE_T(temperature, guesses[0], guesses[1], 200)
        p_ref = (
            rho_l_ref
            * R_J_PER_MOL_K
            * temperature
            * (1.0 + model.get_Ar01(temperature, rho_l_ref, np.array([1.0])))
        )
        p_sat, rho_l, rho_v = pure_saturation(hexane, temperature)
        print(f"  T = {temperature:.1f} K")
        print(
            f"    Psat  : {p_sat:>16.6f} Pa   teqp {p_ref:>16.6f}   "
            f"rel {abs(p_sat - p_ref) / p_ref:.2e}"
        )
        print(
            f"    rho_L : {rho_l:>16.6f} mol/m^3  teqp {rho_l_ref:>16.6f}   "
            f"rel {abs(rho_l - rho_l_ref) / rho_l_ref:.2e}"
        )
        print(
            f"    rho_V : {rho_v:>16.6f} mol/m^3  teqp {rho_v_ref:>16.6f}   "
            f"rel {abs(rho_v - rho_v_ref) / rho_v_ref:.2e}"
        )
        record(
            f"saturation at {temperature:.0f} K agrees with teqp to 1e-6 relative",
            abs(p_sat - p_ref) / p_ref < 1e-6
            and abs(rho_l - rho_l_ref) / rho_l_ref < 1e-6
            and abs(rho_v - rho_v_ref) / rho_v_ref < 1e-6,
        )
        if temperature == 300.0:
            print(
                f"    model vs experiment: PC-SAFT gives {p_sat / 1e3:.3f} kPa here; the "
                "n-hexane vapour\n"
                "    pressure at 300 K is commonly tabulated near 21.7 kPa (unverified "
                "against a\n"
                "    primary source). Reported as a remark; nothing above asserts it."
            )

    print("\n" + "=" * 78)
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
