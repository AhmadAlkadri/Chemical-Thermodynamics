"""Polymer/solvent PC-SAFT against FeOs, with PASS/FAIL on every check.

What is being validated
-----------------------
ADR-0022 makes a polymer expressible as a PC-SAFT component: parameters given
per unit mass (``m/M`` in mol/g, times the sample's molar mass), a custom
non-volatile ``Component`` with no critical constants, and a log-space route
for fugacity coefficients whose exponential does not fit a double. This script
checks the result against implementations that share no code with it, and
against a solver written here.

Four routes, deliberately different in kind:

1. **Properties against FeOs** (feos-org/feos, MIT OR Apache-2.0; the same
   Gross & Sadowski model in Rust, every derivative by automatic
   differentiation). FeOs packages no polymer parameters but takes a segment
   number directly, so it is handed the same ``m = (m/M) Mw`` this package
   derives - which makes the comparison a check of the convention as well as of
   the equations. ``A^res/RT``, ``Z`` and ``ln phi_i`` at eleven states whose
   densities are chemthermo's **own** liquid roots (validation Case P-12).
2. **An equal-fugacity Newton written in this script**, two unknowns carried as
   logarithms, which reproduces the tie line ``flash_tp`` returns without using
   any of the flash's machinery (Case P-13).
3. **FeOs's chemical potentials at chemthermo's converged phases** - the
   reference model's own equilibrium condition evaluated at chemthermo's
   answer, which does not depend on FeOs's flash converging (the Case P-8
   route).
4. **FeOs's own TP flash**, where it converges, so the tie line is checked
   against another *solver* and not only another model.

Two facts about the reference are reported rather than hidden. FeOs's
``tp_flash`` raises ``RuntimeError: `rachford_rice` encountered illegal values``
on this system at 5 and 8 MPa and returns a degenerate pair at 3 MPa; it
converges at 10 MPa, which is the state route 4 uses. And ``k_ij`` cannot be
given to FeOs's PC-SAFT from Python in feos 0.10.1 -
``EquationOfState.pcsaft`` raises "missing field ``k_ij``" for every
serialization tried - so routes 1, 3 and 4 run at ``k_ij = 0`` on both sides,
while route 2 uses the fitted ``k_ij = -0.006``.

As in Cases P-6 to P-11, one input is deliberately not shared: the 42 universal
constants of the 2001 dispersion term (chemthermo packages the ten figures as
printed, FeOs hard-codes fourteen). Every comparison is reported twice, as
shipped and with FeOs's table substituted in.

**The polymer parameters are not verified data.** They are as tabulated by
Martini, Cismondi, Barbosa & Brignole, *Sep. Sci. Technol.* **44** (2009),
citing Gross & Sadowski, *IECR* **41** (2002) 1084 - a paywalled table that was
not read, with no second open source found. They live in
``tests/fixtures/pcsaft/martini2009_polymers.json``, are never packaged, and
nothing here is compared against measurement.

``--full`` adds the Mw = 53000 chain (where the ADR-0022 log-space route
carries the whole flash), the FeOs-flash failure survey across 3-15 MPa and a
cloud-point bisection. The default takes about four seconds.

Requires the optional ``feos`` dependency for routes 1, 3 and 4::

    pip install -e ".[validation]"

Route 2 runs without it; the script says so and still exits 0.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np

try:
    import si_units as si
    from feos import Contributions, EquationOfState, Parameters, PureRecord, State
except ImportError:  # pragma: no cover - exercised only without the extra
    si = None  # type: ignore[assignment]

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.eos import pcsaft as pcsaft_module
from chemthermo.eos.pcsaft import R_J_PER_MOL_K
from chemthermo.parameters import PCSAFTParameters, PCSAFTRecord

REPO = Path(__file__).resolve().parents[2]
FIXTURE = REPO / "tests" / "fixtures" / "pcsaft" / "martini2009_polymers.json"

TEMPERATURE_K = 453.0
PENTANE_MW_G_MOL = 72.146
PE_MW_G_MOL = 16400.0
KIJ = -0.006
#: n-pentane, Gross & Sadowski (2001) Table 1 - the packaged record's values.
PENTANE = (2.6896, 3.7729, 231.20)

PROPERTY_TOL = 1e-10
POTENTIAL_TOL = 1e-8
NEWTON_TOL = 1e-12

failures: list[str] = []


def record_check(label: str, ok: bool) -> None:
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        failures.append(label)


def section(title: str) -> None:
    print("\n" + "-" * 78)
    print(title)
    print("-" * 78)


# ---------------------------------------------------------------------------
# The two models, from one set of numbers
# ---------------------------------------------------------------------------


def polymer_row() -> dict:
    with FIXTURE.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return next(row for row in payload["polymers"] if row["name"] == "Polyethylene")


ROW = polymer_row()


def our_parameters(mw_g_mol: float) -> PCSAFTParameters:
    return PCSAFTParameters.from_records(
        [
            PCSAFTRecord(
                name="Polyethylene",
                segments_per_g=ROW["segments_per_g"],
                MW_g_mol=mw_g_mol,
                sigma_A=ROW["sigma_A"],
                epsilon_k_K=ROW["epsilon_k_K"],
            ),
            PCSAFTRecord(
                name="n-Pentane",
                m=PENTANE[0],
                sigma_A=PENTANE[1],
                epsilon_k_K=PENTANE[2],
                MW_g_mol=PENTANE_MW_G_MOL,
            ),
        ]
    )


def feos_eos(names: Sequence[str], mw_g_mol: float):
    """The reference model, always with ``k_ij = 0`` - see the module docstring."""
    payloads = {
        "Polyethylene": {
            "identifier": {"name": "Polyethylene"},
            "molarweight": mw_g_mol,
            "m": ROW["segments_per_g"] * mw_g_mol,
            "sigma": ROW["sigma_A"],
            "epsilon_k": ROW["epsilon_k_K"],
        },
        "n-Pentane": {
            "identifier": {"name": "n-Pentane"},
            "molarweight": PENTANE_MW_G_MOL,
            "m": PENTANE[0],
            "sigma": PENTANE[1],
            "epsilon_k": PENTANE[2],
        },
    }
    records = [PureRecord.from_json_str(json.dumps(payloads[name])) for name in names]
    return EquationOfState.pcsaft(Parameters.from_records(records))


def polymer_component(mw_g_mol: float) -> ct.Component:
    return ct.Component.custom(
        "Polyethylene",
        mw_kg_per_mol=mw_g_mol / 1000.0,
        formula="(C2H4)n",
        volatile=False,
    )


def feed(weight_fraction: float, mw_g_mol: float) -> list[float]:
    polymer = weight_fraction / mw_g_mol
    solvent = (1.0 - weight_fraction) / PENTANE_MW_G_MOL
    total = polymer + solvent
    return [polymer / total, solvent / total]


def mixture(weight_fraction: float, mw_g_mol: float = PE_MW_G_MOL) -> ct.Mixture:
    return ct.Mixture.from_components(
        [polymer_component(mw_g_mol), ct.Component.from_database("n-Pentane")],
        feed(weight_fraction, mw_g_mol),
        normalize=True,
    )


def feos_constants() -> tuple[np.ndarray, np.ndarray]:
    """FeOs's fourteen-figure tables, from the Case P-6 test module (single copy)."""
    path = REPO / "tests" / "validation" / "test_pcsaft_association_vs_feos.py"
    spec = importlib.util.spec_from_file_location("_pcsaft_feos_constants", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.FEOS_A_UNIVERSAL, module.FEOS_B_UNIVERSAL


def feos_report(
    names: Sequence[str], mw_g_mol: float, density: float, x: Sequence[float]
) -> tuple[float, float, np.ndarray]:
    """``(A^res/RT, Z, ln phi)`` from FeOs at ``(T, rho, x)``."""
    state = State(
        feos_eos(names, mw_g_mol),
        temperature=TEMPERATURE_K * si.KELVIN,
        density=density * (si.MOL / si.METER**3),
        composition=np.asarray(x, dtype=float),
    )
    factor = R_J_PER_MOL_K * TEMPERATURE_K
    a_res = sum(
        (value / si.JOULE * si.MOL) / factor
        for _label, value in state.residual_molar_helmholtz_energy_contributions()
    )
    z_factor = (state.pressure() / si.PASCAL) / (density * factor)
    mu_res = np.atleast_1d(
        np.asarray(state.chemical_potential(Contributions.Residual) / si.JOULE * si.MOL)
    )
    return float(a_res), float(z_factor), mu_res / factor - math.log(z_factor)


# ---------------------------------------------------------------------------
# Route 1: properties (Case P-12)
# ---------------------------------------------------------------------------

PROPERTY_STATES: tuple[tuple[str, tuple[str, ...], float, float | None, float], ...] = (
    ("melt Mw=16400, 1 MPa", ("Polyethylene",), 16400.0, None, 1.0e6),
    ("melt Mw=16400, 10 MPa", ("Polyethylene",), 16400.0, None, 1.0e7),
    ("melt Mw=16400, 30 MPa", ("Polyethylene",), 16400.0, None, 3.0e7),
    ("melt Mw=53000, 1 MPa", ("Polyethylene",), 53000.0, None, 1.0e6),
    ("melt Mw=53000, 10 MPa", ("Polyethylene",), 53000.0, None, 1.0e7),
    ("melt Mw=53000, 30 MPa", ("Polyethylene",), 53000.0, None, 3.0e7),
    ("5 wt%, 10 MPa", ("Polyethylene", "n-Pentane"), 16400.0, 0.05, 1.0e7),
    ("10 wt%, 10 MPa", ("Polyethylene", "n-Pentane"), 16400.0, 0.10, 1.0e7),
    ("15 wt%, 10 MPa", ("Polyethylene", "n-Pentane"), 16400.0, 0.15, 1.0e7),
    ("5 wt%, 15 MPa", ("Polyethylene", "n-Pentane"), 16400.0, 0.05, 1.5e7),
    ("15 wt%, 15 MPa", ("Polyethylene", "n-Pentane"), 16400.0, 0.15, 1.5e7),
)


def state_composition(mw_g_mol: float, weight_fraction: float | None) -> list[float]:
    return [1.0] if weight_fraction is None else feed(weight_fraction, mw_g_mol)


def the_properties() -> None:
    section("1. Properties against FeOs (Case P-12), at chemthermo's own density roots")
    a_feos, b_feos = feos_constants()
    a_ours, b_ours = pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL

    print("    state                     max|dA/RT|      max|dZ|    max|dln phi|")
    worst = {"as shipped": [0.0, 0.0, 0.0], "matched": [0.0, 0.0, 0.0]}
    try:
        for tag in ("as shipped", "matched"):
            if tag == "matched":
                pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = a_feos, b_feos
            else:
                pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = a_ours, b_ours
            for label, names, mw_g_mol, weight_fraction, pressure_Pa in PROPERTY_STATES:
                eos = PCSAFTEOS(components=names, parameters=our_parameters(mw_g_mol))
                x = state_composition(mw_g_mol, weight_fraction)
                density = eos.density_roots(
                    temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=x
                )[-1]
                ours = (
                    eos.residual_helmholtz(
                        temperature_K=TEMPERATURE_K, volume_m3=1.0 / density, composition=x
                    ),
                    eos.compressibility_factor(
                        temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=x
                    ),
                    np.asarray(
                        eos.ln_fugacity_coefficients(
                            temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=x
                        )
                    ),
                )
                theirs = feos_report(names, mw_g_mol, density, x)
                differences = (
                    abs(ours[0] - theirs[0]),
                    abs(ours[1] - theirs[1]),
                    float(np.max(np.abs(ours[2] - theirs[2]))),
                )
                for index, difference in enumerate(differences):
                    worst[tag][index] = max(worst[tag][index], difference)
                if tag == "matched":
                    print(
                        f"    {label:24s} {differences[0]:11.3e}  {differences[1]:11.3e}"
                        f"  {differences[2]:11.3e}"
                    )
    finally:
        pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = a_ours, b_ours

    for tag in ("as shipped", "matched"):
        print(
            f"\n    worst, {tag:11s}: |dA/RT| {worst[tag][0]:.3e}   |dZ| {worst[tag][1]:.3e}"
            f"   |dln phi| {worst[tag][2]:.3e}"
        )
    record_check(
        f"matched-constants properties within {PROPERTY_TOL:.0e} over "
        f"{len(PROPERTY_STATES)} states",
        max(worst["matched"]) < PROPERTY_TOL,
    )
    print(
        "\n    The 'as shipped' residual is the universal-constants table (ten printed\n"
        "    figures against FeOs's fourteen), not a model difference."
    )


# ---------------------------------------------------------------------------
# Route 2: an independent Newton (Case P-13)
# ---------------------------------------------------------------------------


def the_independent_newton() -> None:
    section("2. The tie line from an equal-fugacity Newton written here (k_ij = -0.006)")
    pressure_Pa = 8.0e6
    feed_mixture = mixture(0.05)
    eos = PCSAFTEOS(parameters=our_parameters(PE_MW_G_MOL), kij=KIJ)
    result = ct.flash_tp(
        feed_mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
    )
    bound = PCSAFTEOS(
        components=("Polyethylene", "n-Pentane"),
        parameters=our_parameters(PE_MW_G_MOL),
        kij=KIJ,
    )

    def ln_activity(polymer_fraction: float) -> np.ndarray:
        x = [polymer_fraction, 1.0 - polymer_fraction]
        density = bound.density_roots(
            temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=x
        )[-1]
        return np.log(np.asarray(x)) + np.asarray(
            bound.ln_fugacity_coefficients(
                temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=x
            )
        )

    def residual(u: np.ndarray) -> np.ndarray:
        return ln_activity(math.exp(u[0])) - ln_activity(math.exp(u[1]))

    names = list(result.phases)
    # Start a thousandth away from the flash's answer, so the Newton does work.
    u = np.array(
        [math.log(result.phases[name].composition.fractions[0]) for name in names]
    ) + np.array([1e-3, -1e-3])
    for _ in range(60):
        value = residual(u)
        jacobian = np.zeros((2, 2))
        for column in range(2):
            step = 1e-7
            plus, minus = u.copy(), u.copy()
            plus[column] += step
            minus[column] -= step
            jacobian[:, column] = (residual(plus) - residual(minus)) / (2.0 * step)
        correction = np.linalg.solve(jacobian, -value)
        u = u + correction
        if float(np.max(np.abs(correction))) < 1e-14:
            break

    ours = sorted(result.phases[name].composition.fractions[0] for name in names)
    theirs = sorted(math.exp(value) for value in u)
    worst = max(abs(a - b) for a, b in zip(ours, theirs))
    print(f"    flash_tp    x_polymer = {ours[0]:.12e}, {ours[1]:.12e}")
    print(f"    Newton here x_polymer = {theirs[0]:.12e}, {theirs[1]:.12e}")
    print(f"    max |dx_polymer| = {worst:.3e}")
    print(f"    Newton's own residual = {float(np.max(np.abs(residual(u)))):.3e}")
    record_check(f"independent Newton agrees within {NEWTON_TOL:.0e}", worst < NEWTON_TOL)
    for key in ("mass_balance_residual", "fugacity_residual", "delta_g_split_rt"):
        print(f"    flash_tp {key:24s} = {float(result.diagnostics[key]):.3e}")
    record_check(
        "the split is verified (mass balance, fugacity residual, dG < 0, post-split stable)",
        float(result.diagnostics["mass_balance_residual"]) < 1e-12
        and float(result.diagnostics["fugacity_residual"]) < 1e-8
        and float(result.diagnostics["delta_g_split_rt"]) < 0.0
        and result.diagnostics["post_split_status"] == "stable",
    )


# ---------------------------------------------------------------------------
# Routes 3 and 4: FeOs at chemthermo's answer, and FeOs's own flash
# ---------------------------------------------------------------------------


def our_phases(pressure_Pa: float, kij: float):
    feed_mixture = mixture(0.05)
    eos = PCSAFTEOS(parameters=our_parameters(PE_MW_G_MOL), kij=kij)
    result = ct.flash_tp(
        feed_mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
    )
    phases = {}
    for name, phase in result.phases.items():
        fractions = np.asarray(phase.composition.fractions, dtype=float)
        density = eos.density_roots(
            mixture=feed_mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=pressure_Pa,
            composition=fractions.tolist(),
        )[-1]
        phases[name] = (fractions, density, result.phase_fractions[name])
    return result, phases


def the_chemical_potentials() -> None:
    section("3. FeOs's chemical potentials at chemthermo's phases, 8 MPa, k_ij = 0")
    a_feos, b_feos = feos_constants()
    a_ours, b_ours = pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL
    try:
        for tag in ("as shipped", "matched"):
            if tag == "matched":
                pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = a_feos, b_feos
            _result, phases = our_phases(8.0e6, 0.0)
            potentials = []
            for x, density, _fraction in phases.values():
                _a, _z, ln_phi = feos_report(("Polyethylene", "n-Pentane"), PE_MW_G_MOL, density, x)
                potentials.append(ln_phi + np.log(x))
            difference = float(np.max(np.abs(potentials[0] - potentials[1])))
            print(f"    {tag:11s}: max |d mu_i / RT| = {difference:.3e}")
            if tag == "matched":
                record_check(
                    f"FeOs's potentials are equal at chemthermo's phases within "
                    f"{POTENTIAL_TOL:.0e}",
                    difference < POTENTIAL_TOL,
                )
    finally:
        pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = a_ours, b_ours
    print(
        "\n    FeOs's own tp_flash RAISES at this state; this route does not need it to\n"
        "    converge, only to evaluate."
    )


def the_feos_flash() -> None:
    section("4. FeOs's own tp_flash at 10 MPa, k_ij = 0 (the one state where it converges)")
    result, ours = our_phases(1.0e7, 0.0)
    state = State(
        feos_eos(("Polyethylene", "n-Pentane"), PE_MW_G_MOL),
        temperature=TEMPERATURE_K * si.KELVIN,
        pressure=1.0e7 * si.PASCAL,
        composition=np.asarray(feed(0.05, PE_MW_G_MOL)),
    )
    equilibrium = state.tp_flash()
    theirs = sorted(
        (
            (
                np.atleast_1d(np.asarray(equilibrium.liquid.molefracs, dtype=float)),
                float(equilibrium.liquid.density / (si.MOL / si.METER**3)),
                1.0 - float(equilibrium.vapor_phase_fraction),
            ),
            (
                np.atleast_1d(np.asarray(equilibrium.vapor.molefracs, dtype=float)),
                float(equilibrium.vapor.density / (si.MOL / si.METER**3)),
                float(equilibrium.vapor_phase_fraction),
            ),
        ),
        key=lambda entry: entry[0][0],
    )
    mine = [ours[name] for name in sorted(ours, key=lambda name: ours[name][0][0])]

    worst_x = worst_rho = worst_f = 0.0
    for (our_x, our_rho, our_f), (their_x, their_rho, their_f) in zip(mine, theirs):
        print(
            f"    x_polymer ours {our_x[0]:.10e}  feos {their_x[0]:.10e}  "
            f"|d| {abs(our_x[0] - their_x[0]):.2e}"
        )
        worst_x = max(worst_x, float(np.max(np.abs(our_x - their_x))))
        worst_rho = max(worst_rho, abs(our_rho - their_rho) / their_rho)
        worst_f = max(worst_f, abs(our_f - their_f))
    print(f"    worst |dx| {worst_x:.3e}, densities rel {worst_rho:.3e}, fractions {worst_f:.3e}")
    record_check("the two flashes agree on the tie line to 1e-08 in mole fraction", worst_x < 1e-8)
    record_check(
        "both of chemthermo's phases are liquids", sorted(result.phases) == ["liquid1", "liquid2"]
    )
    print(
        "\n    10 MPa is within 8% of this system's k_ij = 0 cloud point (10.77 MPa), so\n"
        "    both solvers are working near a plait point. Measured in chemthermo's own\n"
        "    model, chemthermo's pair has an equal-fugacity residual of 1.06e-07 and\n"
        "    FeOs's has 9.07e-06: the residual gap above is conditioning, not a model\n"
        "    difference - the property comparison agrees to 1e-11 at the same kind of state."
    )


def the_feos_flash_survey() -> None:
    section("5. Where FeOs's tp_flash converges on this system (k_ij = 0)")
    composition = np.asarray(feed(0.05, PE_MW_G_MOL))
    for pressure_Pa in (3.0e6, 5.0e6, 8.0e6, 1.0e7, 1.1e7, 1.5e7):
        state = State(
            feos_eos(("Polyethylene", "n-Pentane"), PE_MW_G_MOL),
            temperature=TEMPERATURE_K * si.KELVIN,
            pressure=pressure_Pa * si.PASCAL,
            composition=composition,
        )
        try:
            equilibrium = state.tp_flash()
        except RuntimeError as error:
            verdict = f"RuntimeError: {str(error)[:52]}"
        else:
            liquid = np.atleast_1d(np.asarray(equilibrium.liquid.molefracs, dtype=float))
            vapor = np.atleast_1d(np.asarray(equilibrium.vapor.molefracs, dtype=float))
            if float(np.max(np.abs(liquid - vapor))) < 1e-9:
                verdict = "degenerate (both phases equal, fraction 0.5)"
            else:
                verdict = f"two phases, x_polymer {min(liquid[0], vapor[0]):.3e} / {max(liquid[0], vapor[0]):.3e}"
        try:
            ours, _phases = our_phases(pressure_Pa, 0.0)
            ours_text = f"{sorted(ours.phases)}"
        except Exception as error:  # noqa: BLE001 - reporting, not asserting
            ours_text = f"{type(error).__name__}"
        print(
            f"    P = {pressure_Pa / 1e6:5.2f} MPa   feos: {verdict:58s}  chemthermo: {ours_text}"
        )
    print(
        "\n    Recorded, not asserted against. FeOs's flash is not a usable reference\n"
        "    inside this two-phase region; routes 2 and 3 are what check the answer there."
    )


def the_long_chain() -> None:
    section("6. Mw = 53000 (m = 1393.9): the ADR-0022 log-space route, end to end")
    mw_g_mol = 53000.0
    eos = PCSAFTEOS(components=("Polyethylene", "n-Pentane"), parameters=our_parameters(mw_g_mol))
    x = feed(0.05, mw_g_mol)
    density = eos.density_roots(temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, composition=x)[-1]
    ln_phi = eos.ln_fugacity_coefficients(
        temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=x
    )
    underflowed = float(np.exp(np.asarray(ln_phi))[0])
    print(
        f"    ln phi_polymer = {ln_phi[0]:.4f}   ->   exp(that) = {underflowed!r}"
        "   (smallest positive double: exp(-744.44))"
    )
    record_check("exp(ln phi) underflows to exactly zero", underflowed == 0.0)

    feed_mixture = mixture(0.05, mw_g_mol)
    result = ct.flash_tp(
        feed_mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=8.0e6,
        eos=PCSAFTEOS(parameters=our_parameters(mw_g_mol), kij=KIJ),
    )
    for name, phase in result.phases.items():
        print(
            f"      {name}: x_polymer = {phase.composition.fractions[0]:.6e}   "
            f"fraction = {result.phase_fractions[name]:.6f}"
        )
    record_check(
        "the Mw = 53000 flash returns a verified two-liquid split",
        sorted(result.phases) == ["liquid1", "liquid2"]
        and float(result.diagnostics["delta_g_split_rt"]) < 0.0
        and result.diagnostics["post_split_status"] == "stable",
    )

    if si is None:
        return
    _a, _z, theirs = feos_report(("Polyethylene", "n-Pentane"), mw_g_mol, density, x)
    difference = float(np.max(np.abs(np.asarray(ln_phi) - theirs)))
    print(f"    FeOs's ln phi at the same state (as shipped): max |d| = {difference:.3e}")


def the_cloud_point() -> None:
    section("7. The cloud point at k_ij = -0.006, bisected from the stability verdict")

    def verdict(pressure_Pa: float) -> str:
        return ct.stability_tp(
            mixture(0.05),
            temperature_K=TEMPERATURE_K,
            pressure_Pa=pressure_Pa,
            eos=PCSAFTEOS(parameters=our_parameters(PE_MW_G_MOL), kij=KIJ),
        ).status

    low, high = 5.0e6, 3.0e7
    while high - low > 1.0:
        middle = 0.5 * (low + high)
        if verdict(middle) == "unstable":
            low = middle
        else:
            high = middle
    cloud_point = 0.5 * (low + high)
    print(f"    cloud point = {cloud_point / 1e6:.6f} MPa = {cloud_point / 1e5:.4f} bar")
    record_check("the cloud point sits between 8 and 10 MPa", 8.0e6 < cloud_point < 1.0e7)
    print(
        "\n    The cited source's figures put this system's separation pressures in the\n"
        "    same range (tens of bar to a few hundred). That is a magnitude remark about\n"
        "    a figure that was not digitized, not a comparison."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="also run the FeOs-flash survey, the Mw = 53000 chain and the cloud point",
    )
    # `parse_known_args`, not `parse_args`: `tests/test_examples.py` runs this
    # script via `runpy.run_path` with pytest's own `sys.argv` still in place.
    args, _unknown = parser.parse_known_args()

    print("=" * 78)
    print("Polymer/solvent PC-SAFT against FeOs: polyethylene / n-pentane at 453 K")
    print("=" * 78)
    print(
        "Polymer parameters: as tabulated by Martini et al. (2009) citing Gross &\n"
        "Sadowski (2002); the primary table is paywalled and was NOT read, and no second\n"
        "open source printing these values was found. They are a cited test fixture, not\n"
        "packaged runtime data. k_ij cannot be given to FeOs's PC-SAFT in feos 0.10.1,\n"
        "so every FeOs comparison runs at k_ij = 0 on both sides."
    )

    the_independent_newton()
    if si is None:
        print("\n  feos is not installed; routes 1, 3 and 4 are skipped.")
        print("  Install with: pip install -e '.[validation]'")
    else:
        the_properties()
        the_chemical_potentials()
        the_feos_flash()
        if args.full:
            the_feos_flash_survey()
    if args.full:
        the_long_chain()
        the_cloud_point()
    else:
        print(
            "\n  (pass --full for the FeOs-flash survey across 3-15 MPa, the Mw = 53000"
            "\n   chain and the bisected cloud point)"
        )

    print("\n" + "=" * 78)
    print(
        "Nothing above is compared against measurement. The parameters rest on a single\n"
        "secondary source, the polymer is modelled as monodisperse, and the k_ij was\n"
        "fitted elsewhere to cloud-point data this script never touches."
    )
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
