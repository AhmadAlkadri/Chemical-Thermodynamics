"""PC-SAFT association against FeOs, with PASS/FAIL on every check.

What is being validated
-----------------------
``chemthermo.eos.PCSAFTEOS`` and FeOs (feos-org/feos, MIT OR Apache-2.0)
implement the same published model - Gross & Sadowski, Ind. Eng. Chem. Res. 41
(2002) 5510, hard chain + dispersion + **association** - and share no
derivative code:

  * FeOs writes one Helmholtz energy in Rust and obtains every derivative by
    automatic differentiation (``num-dual``);
  * chemthermo writes ``A^res/RT`` and then analytic derivatives, with the
    association term's first derivatives coming from Michelsen & Hendriks's
    stationary ``Q`` function and its second density derivative from the small
    linear solve that goes with it.

teqp, the reference used for the non-associating slices, does **not** implement
association, which is why FeOs is here.

What is checked
---------------
1. **Term by term** at eighteen states - pure water and pure ethanol at four
   each (gas-like, two liquid-like, high temperature) and ten mixture states
   including water/ethanol (cross association), water/n-hexane (one
   associating component) and a ternary: the hard-chain, dispersion and
   association contributions separately, then ``A^res/RT``, ``Z`` and
   ``ln phi_i``.
2. **The ``sigma^3`` versus ``d^3`` question** in the association strength
   ``Delta``. Both spellings are in circulation and the 2002 paper is
   paywalled (HTTP 403 from here), so the convention was settled numerically:
   this script recomputes the pure-water association term both ways and shows
   which one reproduces FeOs.
3. **Non-associating n-hexane is unchanged**, bit for bit, against the numbers
   pinned in validation Case P-1.
4. **Equilibrium**: pure-water saturation at 373.15 K against FeOs's own
   ``PhaseEquilibrium.pure``; a water/ethanol vapour-liquid flash at 351 K; and
   a water/n-hexane liquid-liquid split, with FeOs's fugacities evaluated at
   chemthermo's phases.
5. **A negative control**: perturbing one association parameter by 1 % must
   break the agreement.

One shared input is deliberately not shared
-------------------------------------------
chemthermo packages the 42 universal constants of the 2001 paper **as
printed**, to ten figures, which is also what teqp uses. FeOs hard-codes them
to fourteen figures. The two tables differ by up to 4.8e-09, which floors any
comparison of the *dispersion* term at about 1e-9 - an input difference, not an
implementation difference. Every dispersion-dependent quantity is therefore
reported twice: as shipped, and with FeOs's constants substituted in. The
association term does not depend on those constants at all.

Model versus experiment
-----------------------
Two remarks are printed and **nothing asserts them**: water's saturation
pressure at 373.15 K (101.325 kPa by the definition of the normal boiling
point) and the water/n-hexane mutual solubilities (commonly tabulated near
5e-4 mole fraction water in hexane and 2e-6 hexane in water at 298 K, not
verified against a primary source here). With ``k_ij = 0`` PC-SAFT is known to
be poor for the latter; this script compares one implementation of the model
to another, not the model to measurement.

Requires the optional ``feos`` dependency::

    pip install -e ".[validation]"

The script prints a message and exits 0 when it is missing.
"""

from __future__ import annotations

import json
import math

import numpy as np

try:
    import si_units as si
    from feos import (
        Contributions,
        EquationOfState,
        Parameters,
        PhaseEquilibrium,
        PureRecord,
        State,
    )
except ImportError:  # pragma: no cover - exercised only without the extra
    si = None  # type: ignore[assignment]

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.eos import pcsaft as pcsaft_module
from chemthermo.eos.pcsaft import R_J_PER_MOL_K

# Gross & Sadowski (2002) Table 1 for the associating four, (2001) Table 1 for
# n-hexane; written out here so the reference model is built from this file.
# ``name -> (MW, m, sigma/A, eps/k, kappa^AB or None, eps^AB/k)``.
PARAMETERS: dict[str, tuple[float, float, float, float, float | None, float]] = {
    "Water": (18.015, 1.0656, 3.0007, 366.51, 0.034868, 2500.7),
    "Ethanol": (46.069, 2.3827, 3.1771, 198.24, 0.032384, 2653.4),
    "Methanol": (32.042, 1.5255, 3.2300, 188.90, 0.035176, 2899.5),
    "n-Hexane": (86.177, 3.0576, 3.7983, 236.77, None, 0.0),
}

#: FeOs's own universal constants (``crates/feos/src/pcsaft/eos/dispersion.rs``).
FEOS_A = np.array(
    [
        [
            0.91056314451539,
            0.63612814494991,
            2.68613478913903,
            -26.5473624914884,
            97.7592087835073,
            -159.591540865600,
            91.2977740839123,
        ],
        [
            -0.30840169182720,
            0.18605311591713,
            -2.50300472586548,
            21.4197936296668,
            -65.2558853303492,
            83.3186804808856,
            -33.7469229297323,
        ],
        [
            -0.09061483509767,
            0.45278428063920,
            0.59627007280101,
            -1.72418291311787,
            -4.13021125311661,
            13.7766318697211,
            -8.67284703679646,
        ],
    ]
)
FEOS_B = np.array(
    [
        [
            0.72409469413165,
            2.23827918609380,
            -4.00258494846342,
            -21.00357681484648,
            26.8556413626615,
            206.5513384066188,
            -355.60235612207947,
        ],
        [
            -0.57554980753450,
            0.69950955214436,
            3.89256733895307,
            -17.21547164777212,
            192.6722644652495,
            -161.8264616487648,
            -165.2076934555607,
        ],
        [
            0.09768831158356,
            -0.25575749816100,
            -9.15585615297321,
            20.64207597439724,
            -38.80443005206285,
            93.6267740770146,
            -29.66690558514725,
        ],
    ]
)

#: ``(label, components, x, T/K, rho/(mol/m^3))``.
STATES: list[tuple[str, tuple[str, ...], list[float], float, float]] = [
    ("water 300 K / 55000", ("Water",), [1.0], 300.0, 55000.0),
    ("water 350 K / 50000", ("Water",), [1.0], 350.0, 50000.0),
    ("water 373.15 K / 100", ("Water",), [1.0], 373.15, 100.0),
    ("water 550 K / 41241", ("Water",), [1.0], 550.0, 41241.1907),
    ("water 550 K / 1250", ("Water",), [1.0], 550.0, 1250.012),
    ("ethanol 300 K / 17000", ("Ethanol",), [1.0], 300.0, 17000.0),
    ("ethanol 350 K / 15929", ("Ethanol",), [1.0], 350.0, 15928.6188),
    ("ethanol 450 K / 691", ("Ethanol",), [1.0], 450.0, 691.1054),
    ("ethanol 500 K / 9929", ("Ethanol",), [1.0], 500.0, 9929.2725),
    ("H2O/EtOH 0.2/0.8 320 K / 20000", ("Water", "Ethanol"), [0.2, 0.8], 320.0, 20000.0),
    ("H2O/EtOH 0.5/0.5 320 K / 25401", ("Water", "Ethanol"), [0.5, 0.5], 320.0, 25401.3411),
    ("H2O/EtOH 0.8/0.2 320 K / 36508", ("Water", "Ethanol"), [0.8, 0.2], 320.0, 36507.9501),
    ("H2O/EtOH 0.5/0.5 400 K / 94.6", ("Water", "Ethanol"), [0.5, 0.5], 400.0, 94.6239),
    ("H2O/EtOH 0.2/0.8 351 K / 18575", ("Water", "Ethanol"), [0.2, 0.8], 351.0, 18574.9573),
    ("H2O/C6 0.3/0.7 298 K / 10229", ("Water", "n-Hexane"), [0.3, 0.7], 298.15, 10229.0998),
    ("H2O/C6 0.9/0.1 298 K / 33473", ("Water", "n-Hexane"), [0.9, 0.1], 298.15, 33472.7843),
    ("H2O/C6 0.5/0.5 400 K / 93.9", ("Water", "n-Hexane"), [0.5, 0.5], 400.0, 93.9327),
    (
        "MeOH/H2O/C6 0.3/0.4/0.3 320 K / 16685",
        ("Methanol", "Water", "n-Hexane"),
        [0.3, 0.4, 0.3],
        320.0,
        16684.835,
    ),
]

#: Case P-1's pinned non-associating n-hexane numbers.
PINNED_HEXANE = (-5.783742760059239, 0.661534529144653, -5.709015132378622)

TIGHT = 1e-10
DISPERSION_LIMITED = 1e-8
LN_PHI_LIMITED = 1e-5

failures: list[str] = []


def record(label: str, ok: bool) -> None:
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        failures.append(label)


def pure_record(name: str) -> "PureRecord":
    mw, m, sigma, epsilon, kappa, epsilon_ab = PARAMETERS[name]
    payload: dict[str, object] = {
        "identifier": {"name": name},
        "molarweight": mw,
        "m": m,
        "sigma": sigma,
        "epsilon_k": epsilon,
    }
    if kappa is not None:
        payload["association_sites"] = [
            {"kappa_ab": kappa, "epsilon_k_ab": epsilon_ab, "na": 1.0, "nb": 1.0}
        ]
    return PureRecord.from_json_str(json.dumps(payload))


def feos_eos(names: tuple[str, ...]) -> "EquationOfState":
    """FeOs's PC-SAFT with ``k_ij = 0`` (an empty binary-record list)."""
    records = [pure_record(name) for name in names]
    if len(records) == 1:
        return EquationOfState.pcsaft(Parameters.new_pure(records[0]))
    return EquationOfState.pcsaft(Parameters.from_records(records))


def feos_report(
    names: tuple[str, ...], temperature: float, density: float, x: list[float]
) -> tuple[dict[str, float], float, np.ndarray]:
    eos = feos_eos(names)
    unit = si.MOL / si.METER**3
    if len(names) == 1:
        state = State(eos, temperature=temperature * si.KELVIN, density=density * unit)
    else:
        state = State(
            eos,
            temperature=temperature * si.KELVIN,
            density=density * unit,
            composition=np.asarray(x),
        )
    factor = R_J_PER_MOL_K * temperature
    contributions = {
        label: (value / si.JOULE * si.MOL) / factor
        for label, value in state.residual_molar_helmholtz_energy_contributions()
    }
    z_factor = (state.pressure() / si.PASCAL) / (density * factor)
    mu_res = np.atleast_1d(state.chemical_potential(Contributions.Residual) / si.JOULE * si.MOL)
    return contributions, float(z_factor), mu_res / factor - math.log(z_factor)


def term_by_term() -> None:
    print("\n1) Term by term against FeOs, at 18 states")
    print(
        "   Two runs: 'shipped' uses the published 10-figure universal constants,\n"
        "   'matched' substitutes FeOs's 14-figure ones into chemthermo."
    )
    shipped = {key: 0.0 for key in ("hc", "disp", "assoc", "a_res", "Z", "ln_phi")}
    matched = dict(shipped)
    published_a = pcsaft_module.A_UNIVERSAL
    published_b = pcsaft_module.B_UNIVERSAL
    try:
        for mode, worst in (("shipped", shipped), ("matched", matched)):
            if mode == "matched":
                pcsaft_module.A_UNIVERSAL = FEOS_A
                pcsaft_module.B_UNIVERSAL = FEOS_B
            for label, names, x, temperature, density in STATES:
                eos = PCSAFTEOS(components=names)
                contributions, z_reference, ln_phi_reference = feos_report(
                    names, temperature, density, x
                )
                terms = eos.residual_helmholtz_terms(
                    temperature_K=temperature, volume_m3=1.0 / density, composition=x
                )
                z_factor = eos.compressibility_factor(
                    temperature_K=temperature, density_mol_m3=density, composition=x
                )
                ln_phi = np.asarray(
                    eos.ln_fugacity_coefficients(
                        temperature_K=temperature, density_mol_m3=density, composition=x
                    )
                )
                deltas = {
                    "hc": abs(
                        terms["hard-chain"]
                        - (contributions["Hard Sphere"] + contributions["Hard Chain"])
                    ),
                    "disp": abs(terms["dispersion"] - contributions["Dispersion"]),
                    "assoc": abs(terms["association"] - contributions["Association"]),
                    "a_res": abs(terms["total"] - sum(contributions.values())),
                    "Z": abs(z_factor - z_reference),
                    "ln_phi": float(np.max(np.abs(ln_phi - ln_phi_reference))),
                }
                for key, value in deltas.items():
                    worst[key] = max(worst[key], value)
                if mode == "matched":
                    print(
                        f"    {label:<38} a_assoc {terms['association']:>12.6f}  "
                        f"|d| {deltas['assoc']:.1e}  |dlnphi| {deltas['ln_phi']:.1e}"
                    )
    finally:
        pcsaft_module.A_UNIVERSAL = published_a
        pcsaft_module.B_UNIVERSAL = published_b

    print("\n    worst |difference| over the 18 states:")
    header = f"    {'quantity':<10}{'shipped constants':>22}{'FeOs constants':>20}"
    print(header)
    for key in ("hc", "disp", "assoc", "a_res", "Z", "ln_phi"):
        print(f"    {key:<10}{shipped[key]:>22.3e}{matched[key]:>20.3e}")

    record(
        f"association term agrees to {TIGHT:g} with either table "
        f"(worst {max(shipped['assoc'], matched['assoc']):.2e})",
        max(shipped["assoc"], matched["assoc"]) < TIGHT,
    )
    record(
        f"hard-chain term agrees to {TIGHT:g} with either table (worst "
        f"{max(shipped['hc'], matched['hc']):.2e})",
        max(shipped["hc"], matched["hc"]) < TIGHT,
    )
    record(
        f"with FeOs's constants every quantity agrees to {TIGHT:g} "
        f"(worst {max(matched.values()):.2e})",
        max(matched.values()) < TIGHT,
    )
    record(
        f"as shipped, A^res/RT and Z agree to {DISPERSION_LIMITED:g} "
        f"(worst {max(shipped['a_res'], shipped['Z']):.2e})",
        max(shipped["a_res"], shipped["Z"]) < DISPERSION_LIMITED,
    )
    record(
        f"as shipped, ln phi agrees to {LN_PHI_LIMITED:g} (worst {shipped['ln_phi']:.2e})",
        shipped["ln_phi"] < LN_PHI_LIMITED,
    )
    difference = max(
        float(np.max(np.abs(FEOS_A - published_a))), float(np.max(np.abs(FEOS_B - published_b)))
    )
    print(
        f"    The two tables of universal constants differ by up to {difference:.1e}, which is\n"
        "    where the 'shipped' column's floor comes from. Nothing above is a"
        " disagreement\n    about the model."
    )


def sigma_versus_d() -> None:
    print("\n2) The association strength's prefactor: sigma^3 or d^3?")
    temperature, density = 300.0, 55000.0
    _, m, sigma, epsilon_k, kappa, epsilon_ab = PARAMETERS["Water"]
    assert kappa is not None
    rho = density * 6.02214076e23 * 1e-30
    d = sigma * (1.0 - 0.12 * math.exp(-3.0 * epsilon_k / temperature))
    zeta_2 = math.pi / 6.0 * rho * m * d**2
    eta = math.pi / 6.0 * rho * m * d**3
    u = 1.0 - eta
    contact = 0.5 * d
    g = 1.0 / u + 3.0 * contact * zeta_2 / u**2 + 2.0 * contact**2 * zeta_2**2 / u**3

    def association(prefactor: float) -> float:
        strength = rho * prefactor * g * kappa * math.expm1(epsilon_ab / temperature)
        site = (-1.0 + math.sqrt(1.0 + 4.0 * strength)) / (2.0 * strength)
        return 2.0 * math.log(site) - site + 1.0

    contributions, _, _ = feos_report(("Water",), temperature, density, [1.0])
    reference = contributions["Association"]
    with_sigma = association(sigma**3)
    with_d = association(d**3)
    print(f"    FeOs                      a_assoc = {reference:.12f}")
    print(
        f"    sigma^3 (sigma = {sigma:.4f} A)  a_assoc = {with_sigma:.12f}  "
        f"|d| {abs(with_sigma - reference):.1e}"
    )
    print(
        f"    d^3     (d     = {d:.4f} A)  a_assoc = {with_d:.12f}  "
        f"|d| {abs(with_d - reference):.1e}"
    )
    record("sigma^3 reproduces FeOs to 1e-12", abs(with_sigma - reference) < 1e-12)
    record("d^3 does not (the choice is unambiguous)", abs(with_d - reference) > 1e-3)
    shipped = PCSAFTEOS(components=("Water",)).residual_helmholtz_terms(
        temperature_K=temperature, volume_m3=1.0 / density, composition=[1.0]
    )["association"]
    record(
        "the shipped model uses the sigma^3 form",
        abs(shipped - with_sigma) < 1e-13,
    )
    print(
        "    The 2002 paper was NOT read (pubs.acs.org returns HTTP 403 from here), so\n"
        "    this numerical test is the evidence for the convention, not a citation."
    )


def non_associating_unchanged() -> None:
    print("\n3) Non-associating n-hexane is untouched (validation Case P-1's pins)")
    eos = PCSAFTEOS(components=("n-Hexane",))
    a_res = eos.residual_helmholtz(temperature_K=300.0, volume_m3=1.0 / 7700.0, composition=[1.0])
    z_factor = eos.compressibility_factor(
        temperature_K=300.0, density_mol_m3=7700.0, composition=[1.0]
    )
    ln_phi = eos.ln_fugacity_coefficients(
        temperature_K=300.0, density_mol_m3=7700.0, composition=[1.0]
    )[0]
    print(f"    A^res/RT {a_res!r}\n    Z        {z_factor!r}\n    ln phi   {ln_phi!r}")
    record(
        "n-hexane at 300 K / 7700 mol/m^3 is bit-identical to the pinned values",
        (a_res, z_factor, ln_phi) == PINNED_HEXANE,
    )
    record("and the model reports itself non-associating", not eos.associates())


def saturation() -> None:
    print("\n4) Pure water saturation at 373.15 K")
    temperature = 373.15
    eos = PCSAFTEOS(components=("Water",))

    def gap(pressure: float) -> float:
        roots = eos.density_roots(
            temperature_K=temperature, pressure_Pa=pressure, composition=[1.0]
        )
        return (
            eos.ln_fugacity_coefficients(
                temperature_K=temperature, density_mol_m3=roots[-1], composition=[1.0]
            )[0]
            - eos.ln_fugacity_coefficients(
                temperature_K=temperature, density_mol_m3=roots[0], composition=[1.0]
            )[0]
        )

    low, high = 5.0e4, 2.0e5
    low_value = gap(low)
    for _ in range(90):
        middle = 0.5 * (low + high)
        value = gap(middle)
        if low_value * value <= 0.0:
            high = middle
        else:
            low, low_value = middle, value
    pressure = 0.5 * (low + high)
    roots = eos.density_roots(temperature_K=temperature, pressure_Pa=pressure, composition=[1.0])

    equilibrium = PhaseEquilibrium.pure(feos_eos(("Water",)), temperature * si.KELVIN)
    unit = si.MOL / si.METER**3
    reference = (
        equilibrium.liquid.pressure() / si.PASCAL,
        equilibrium.liquid.density / unit,
        equilibrium.vapor.density / unit,
    )
    mine = (pressure, roots[-1], roots[0])
    for name, value, expected in zip(("Psat", "rho_L", "rho_V"), mine, reference):
        print(
            f"    {name:<6} {value:>18,.6f}   FeOs {expected:>18,.6f}   "
            f"rel {abs(value - expected) / expected:.2e}"
        )
    record(
        "saturation matches FeOs's own Newton solve to 1e-6 relative",
        all(abs(a - b) / b < 1e-6 for a, b in zip(mine, reference)),
    )
    print(
        f"    Model versus experiment (a remark): 373.15 K is water's normal boiling\n"
        f"    point, so the experimental saturation pressure there is 101,325 Pa by\n"
        f"    definition; this model is {(101325.0 - pressure) / 101325.0 * 100:.2f} %"
        " low. Nothing asserts that."
    )


def equal_fugacity(
    names: tuple[str, ...],
    temperature: float,
    phases: list[list[float]],
    densities: list[float],
) -> float:
    logs = []
    for x, density in zip(phases, densities):
        _, _, ln_phi = feos_report(names, temperature, density, x)
        logs.append(np.log(np.asarray(x)) + np.atleast_1d(ln_phi))
    return float(np.max(np.abs(logs[0] - logs[1])))


def vapor_liquid_flash() -> None:
    print("\n5) Water / ethanol vapour-liquid flash at 351 K and 80 kPa, z = 0.7 / 0.3")
    names = ("Water", "Ethanol")
    temperature, pressure = 351.0, 80.0e3
    mixture = ct.Mixture.from_database(list(names), [0.7, 0.3])
    published_a, published_b = pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL
    try:
        # FeOs evaluates the fugacities at chemthermo's own densities, so the
        # universal constants have to match or the residual is floored at the
        # 1e-6 level by the table difference, not by the solver.
        pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = FEOS_A, FEOS_B
        eos = PCSAFTEOS()
        stability = ct.stability_tp(
            mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos
        )
        print(f"    stability_tp -> {stability.status!r}, tpd_min = {stability.tpd_min:.6e}")
        result = ct.flash_tp(mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos)
        phases: list[list[float]] = []
        densities: list[float] = []
        for name in ("liquid", "vapor"):
            fractions = list(result.phases[name].composition.fractions)
            roots = eos.density_roots(
                temperature_K=temperature,
                pressure_Pa=pressure,
                composition=fractions,
                mixture=mixture,
            )
            density = roots[-1] if name == "liquid" else roots[0]
            phases.append(fractions)
            densities.append(density)
            print(
                f"    {name:<7} x = ({fractions[0]:.8f}, {fractions[1]:.8f})  "
                f"rho = {density:>12,.4f} mol/m^3"
            )
        print(
            f"    vapor_fraction = {result.vapor_fraction:.8f}, "
            f"dG/RT = {result.diagnostics['delta_g_split_rt']:.3e}, "
            f"fugacity residual = {result.diagnostics['fugacity_residual']:.2e}"
        )
        residual = equal_fugacity(names, temperature, phases, densities)
        print(f"    FeOs equal-fugacity residual at chemthermo's phases: {residual:.2e}")
        record("the split is a Gibbs decrease", result.diagnostics["delta_g_split_rt"] < 0.0)
        record(
            "both phases pass the post-split stability test",
            result.diagnostics["post_split_status"] == "stable",
        )
        record("FeOs agrees the phases are in equilibrium to 1e-8", residual < 1e-8)
    finally:
        pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = published_a, published_b


def liquid_liquid_split() -> None:
    print("\n6) Water / n-hexane at 298.15 K, z = 0.5 / 0.5")
    names = ("Water", "n-Hexane")
    temperature = 298.15
    mixture = ct.Mixture.from_database(list(names), [0.5, 0.5])

    eos = PCSAFTEOS()
    stability = ct.stability_tp(mixture, temperature_K=temperature, pressure_Pa=101325.0, eos=eos)
    print(f"    at 1 atm: stability_tp -> {stability.status!r}, tpd_min = {stability.tpd_min:.6e}")
    raised = False
    try:
        ct.flash_tp(mixture, temperature_K=temperature, pressure_Pa=101325.0, eos=eos)
    except ct.ConvergenceError:
        raised = True
    record("at 1 atm the feed is unstable", stability.status == "unstable")
    record(
        "at 1 atm flash_tp refuses rather than returning a spurious vapour-liquid pair",
        raised,
    )
    print(
        "    That refusal is a known limitation, not a bug in the model: the phi-phi\n"
        "    split can only pair a vapour-root phase with a liquid-root one, and at\n"
        "    1 atm the isotherm still has a vapour root. FeOs's own two-phase flash at\n"
        "    the same state returns two liquids. See ADR-0018."
    )

    published_a, published_b = pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL
    try:
        pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = FEOS_A, FEOS_B
        eos = PCSAFTEOS()
        pressure = 1.0e6
        roots = eos.density_roots(
            temperature_K=temperature,
            pressure_Pa=pressure,
            composition=[0.5, 0.5],
            mixture=mixture,
        )
        print(
            f"\n    at 1 MPa the feed has {len(roots)} density root(s): the vapour branch is gone"
        )
        result = ct.flash_tp(mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos)
        phases: list[list[float]] = []
        densities: list[float] = []
        identities: list[str] = []
        for name in ("liquid", "vapor"):
            fractions = list(result.phases[name].composition.fractions)
            density = eos.density_roots(
                temperature_K=temperature,
                pressure_Pa=pressure,
                composition=fractions,
                mixture=mixture,
            )[0]
            identity = eos.phase_identity(
                mixture=mixture,
                temperature_K=temperature,
                pressure_Pa=pressure,
                composition=fractions,
                phase="liquid",
            )
            phases.append(fractions)
            densities.append(density)
            identities.append(identity)
            print(
                f"    labelled {name!r:<8} x = ({fractions[0]:.8f}, {fractions[1]:.8f})  "
                f"rho = {density:>12,.2f}  measured identity {identity!r}"
            )
        print(
            f"    phase_label_method = {result.diagnostics['phase_label_method']!r}, "
            f"vapor_fraction = {result.vapor_fraction:.6f}"
        )
        residual = equal_fugacity(names, temperature, phases, densities)
        print(f"    FeOs equal-fugacity residual at chemthermo's phases: {residual:.2e}")
        record("both converged phases measure as liquids", identities == ["liquid", "liquid"])
        record("the split is a Gibbs decrease", result.diagnostics["delta_g_split_rt"] < 0.0)
        record("FeOs agrees the two liquids are in equilibrium to 1e-8", residual < 1e-8)
        print(
            "    Labels: ADR-0017 measures both phases as liquids, but the phi-phi path\n"
            "    has no 'liquid1' / 'liquid2' naming, so it falls back to the Wilson\n"
            "    ranking and 'vapor_fraction' is the hexane-rich LIQUID's fraction.\n"
            "    Recorded as a limitation; fixing it is the next slice's job."
        )
        print(
            "    Mutual solubilities from this model, against the commonly tabulated\n"
            "    experimental figures (which were NOT verified against a primary source\n"
            "    here, and which nothing above asserts):\n"
            f"      hexane in the water-rich phase  {phases[0][1]:.3e}   experiment ~2e-6\n"
            f"      water  in the hexane-rich phase {phases[1][0]:.3e}   experiment ~5e-4\n"
            "    With k_ij = 0 PC-SAFT is known to be poor for water / hydrocarbon mutual\n"
            "    solubilities, and it is an order of magnitude out on both. This section\n"
            "    is a check of the code against FeOs, not of the model against measurement."
        )
    finally:
        pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = published_a, published_b


def negative_control() -> None:
    print("\n7) Negative control")
    temperature, density = 300.0, 55000.0
    contributions, _, _ = feos_report(("Water",), temperature, density, [1.0])
    perturbed = ct.PCSAFTParameters.from_records(
        [
            {
                "name": "Water",
                "m": 1.0656,
                "sigma_A": 3.0007,
                "epsilon_k_K": 366.51,
                "association": {"kappa_ab": 0.034868 * 1.01, "epsilon_ab_k_K": 2500.7},
            }
        ]
    )
    shifted = PCSAFTEOS(components=("Water",), parameters=perturbed).residual_helmholtz_terms(
        temperature_K=temperature, volume_m3=1.0 / density, composition=[1.0]
    )["association"]
    difference = abs(shifted - contributions["Association"])
    print(f"    kappa^AB + 1 % moves a_assoc by {difference:.3e}")
    record("a 1 % change in kappa^AB breaks the agreement", difference > 1e-3)


def main() -> None:
    if si is None:
        print("feos is not installed; skipping. Install with: pip install -e '.[validation]'")
        return
    print("PC-SAFT association versus FeOs (validation Cases P-6 and P-7, ADR-0018)")
    print("=" * 78)
    term_by_term()
    sigma_versus_d()
    non_associating_unchanged()
    saturation()
    vapor_liquid_flash()
    liquid_liquid_split()
    negative_control()

    print("\n" + "=" * 78)
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
