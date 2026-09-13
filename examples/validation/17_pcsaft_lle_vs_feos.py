"""The EOS liquid-liquid tie line against FeOs, with PASS/FAIL on every check.

What is being validated
-----------------------
ADR-0019 taught ``flash_tp``'s phi-phi split to put each phase on the density
branch the tangent-plane stability test found *that phase* on, instead of
pinning one phase to the model's liquid root and the other to its vapour root.
That is what makes a liquid-liquid split from an equation of state expressible
at a state where a vapour root still exists - water / n-hexane at 298.15 K and
1 atm, validation Case P-7(iii)'s documented failure, now Case P-8.

The reference is **FeOs** (feos-org/feos, MIT OR Apache-2.0): the same Gross &
Sadowski PC-SAFT model implemented independently in Rust, with every derivative
from automatic differentiation (``num-dual``) and with its **own** two-phase
flash. Three separate things are compared, so a coincidence in one does not
carry the others:

1. **Tie line, densities and phase amounts** against FeOs's own
   ``State.tp_flash`` - two solvers, not just two models.
2. **FeOs's chemical potentials evaluated at chemthermo's phases** - this one
   does not depend on FeOs's flash converging at all; it asks whether the state
   chemthermo returned is an equilibrium state *of the reference model*.
3. **The lever rule across three feeds**, and the fact that the tie line does
   not move with the feed.

A negative control is included: perturbing water's association energy by 1 %
must move the tie line, so "the two agree" is not a statement about numbers
that no longer depend on the model.

One shared input is deliberately not shared
-------------------------------------------
chemthermo packages the 42 universal constants of the 2001 dispersion term
**as printed** (ten figures); FeOs hard-codes them to fourteen. The tables
differ by up to 4.8e-09, which is an input difference, not an implementation
difference, and it floors any comparison that evaluates FeOs at chemthermo's
own densities at about 1e-6. Every such quantity is therefore reported twice:
as shipped, and with FeOs's constants substituted in.

What the reference will not do
------------------------------
FeOs's own ``tp_flash`` raises ``"stability analysis did not converge"`` on
this binary at ``z = 0.2/0.8`` (measured with feos 0.10.1). That is printed,
not hidden, and every feed is checked against the tie line from
``z = 0.5/0.5``, which FeOs does return.

Model versus experiment
-----------------------
Printed and **not** asserted: with ``k_ij = 0`` this model gives about 1.7e-05
mole fraction hexane in the water-rich phase and 6.3e-03 water in the
hexane-rich phase, against commonly tabulated experimental figures near 2e-06
and 5e-04 at 298 K (tabulated values, not verified against a primary source
here). This script compares one implementation of the model against another,
not the model against measurement.

Requires the optional ``feos`` dependency::

    pip install -e ".[validation]"

The script prints a message and exits 0 when it is missing.
"""

from __future__ import annotations

import json
import math
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
from chemthermo.models.base import KAPPA_LIQUID_THRESHOLD

NAMES = ("Water", "n-Hexane")
TEMPERATURE_K = 298.15
ATMOSPHERE_PA = 101325.0
HIGH_PRESSURE_PA = 1.0e6

#: Gross & Sadowski (2002) Table 1 for water, (2001) Table 1 for n-hexane,
#: written out here so the reference model is built from this file.
PARAMETERS: dict[str, tuple[float, float, float, float, float | None, float]] = {
    "Water": (18.015, 1.0656, 3.0007, 366.51, 0.034868, 2500.7),
    "n-Hexane": (86.177, 3.0576, 3.7983, 236.77, None, 0.0),
}

failures: list[str] = []


def record(label: str, ok: bool) -> None:
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        failures.append(label)


# ---------------------------------------------------------------------------
# The FeOs side
# ---------------------------------------------------------------------------


def _pure_record(name: str) -> "PureRecord":
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


def feos_eos() -> "EquationOfState":
    """The reference model, always with ``k_ij = 0`` (FeOs's own default)."""
    return EquationOfState.pcsaft(Parameters.from_records([_pure_record(n) for n in NAMES]))


def feos_tp_flash(pressure_Pa: float, z: Sequence[float]):
    """FeOs's own flash: ``((x, rho), (x, rho), fraction of the second)``."""
    state = State(
        feos_eos(),
        temperature=TEMPERATURE_K * si.KELVIN,
        pressure=pressure_Pa * si.PASCAL,
        composition=np.asarray(z, dtype=float),
    )
    equilibrium = state.tp_flash()
    mol_per_m3 = si.MOL / si.METER**3
    dense = (
        np.atleast_1d(np.asarray(equilibrium.liquid.molefracs, dtype=float)),
        float(equilibrium.liquid.density / mol_per_m3),
    )
    light = (
        np.atleast_1d(np.asarray(equilibrium.vapor.molefracs, dtype=float)),
        float(equilibrium.vapor.density / mol_per_m3),
    )
    return dense, light, float(equilibrium.vapor_phase_fraction)


def feos_reduced_potentials(density: float, x: Sequence[float]) -> np.ndarray:
    """``mu_i / RT`` from FeOs at chemthermo's ``(T, rho, x)``, up to a constant."""
    values = np.asarray(x, dtype=float)
    state = State(
        feos_eos(),
        temperature=TEMPERATURE_K * si.KELVIN,
        density=density * (si.MOL / si.METER**3),
        composition=values,
    )
    factor = R_J_PER_MOL_K * TEMPERATURE_K
    z_factor = (state.pressure() / si.PASCAL) / (density * factor)
    mu_residual = np.atleast_1d(
        np.asarray(state.chemical_potential(Contributions.Residual) / si.JOULE * si.MOL)
    )
    return mu_residual / factor - math.log(z_factor) + np.log(values)


# ---------------------------------------------------------------------------
# The chemthermo side
# ---------------------------------------------------------------------------


def our_phases(pressure_Pa: float, z: Sequence[float]):
    """``flash_tp`` plus, per phase, ``(composition, density, phase fraction)``."""
    mixture = ct.Mixture.from_database(list(NAMES), list(z), normalize=True)
    eos = PCSAFTEOS()
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
    phases = {}
    for name, phase in result.phases.items():
        fractions = np.asarray(phase.composition.fractions, dtype=float)
        density = eos.density_roots(
            mixture=mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=pressure_Pa,
            composition=fractions.tolist(),
        )[-1]
        phases[name] = (fractions, density, result.phase_fractions[name])
    return result, phases


def kappa(pressure_Pa: float, composition: Sequence[float]) -> float:
    """``kappa = P / (rho dP/drho)``, from the public pressure routine."""
    mixture = ct.Mixture.from_database(list(NAMES), [0.5, 0.5], normalize=True)
    eos = PCSAFTEOS()
    density = eos.density_roots(
        mixture=mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure_Pa,
        composition=list(composition),
    )[-1]
    bound = PCSAFTEOS(components=NAMES)
    step = 1e-4 * density
    slope = (
        bound.pressure_Pa(
            temperature_K=TEMPERATURE_K,
            density_mol_m3=density + step,
            composition=list(composition),
        )
        - bound.pressure_Pa(
            temperature_K=TEMPERATURE_K,
            density_mol_m3=density - step,
            composition=list(composition),
        )
    ) / (2.0 * step)
    return pressure_Pa / (density * slope)


#: FeOs's own universal constants, loaded from the Case P-6/P-7 script rather
#: than pasted a third time (two copies of a 42-number table already drift).
def feos_universal_constants() -> tuple[np.ndarray, np.ndarray]:
    import importlib.util
    from pathlib import Path

    path = Path(__file__).with_name("16_pcsaft_association_vs_feos.py")
    spec = importlib.util.spec_from_file_location("_pcsaft_feos_constants", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.FEOS_A, module.FEOS_B


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------


def the_state() -> None:
    print("\n1) The state: water / n-hexane, 298.15 K, 1 atm, z = 0.5 / 0.5")
    eos = PCSAFTEOS()
    mixture = ct.Mixture.from_database(list(NAMES), [0.5, 0.5], normalize=True)
    stability = ct.stability_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=ATMOSPHERE_PA, eos=eos
    )
    roots = eos.density_roots(
        mixture=mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=ATMOSPHERE_PA,
        composition=[0.5, 0.5],
    )
    print(
        f"    stability_tp -> {stability.status!r}, tpd_min = {stability.tpd_min:.6e}, "
        f"branches {stability.feed_branch!r} / {stability.phase_branch!r}"
    )
    print(
        f"    density roots at the feed: {len(roots)} -> "
        + ", ".join(f"{value:,.2f}" for value in roots)
        + " mol/m^3"
    )
    record("the feed is unstable", stability.status == "unstable")
    record("a vapour root still exists here (the premise of Case P-7(iii))", len(roots) == 2)


def the_tie_line(section: int, pressure_Pa: float, label: str) -> None:
    print(f"\n{section}) The tie line at {label}, against FeOs's own tp_flash")
    result, ours = our_phases(pressure_Pa, (0.5, 0.5))
    record(
        "flash_tp returns two liquids named liquid1 / liquid2",
        sorted(result.phases) == ["liquid1", "liquid2"] and result.vapor_fraction is None,
    )
    record("phase_regime is LLE", result.diagnostics["phase_regime"] == "LLE")
    record(
        "the names came from the compressibility criterion, not a fallback",
        result.diagnostics["phase_label_method"] == "compressibility",
    )
    record("the split is a Gibbs decrease", float(result.diagnostics["delta_g_split_rt"]) < 0.0)
    record("equal fugacities to 1e-9", float(result.diagnostics["fugacity_residual"]) < 1e-9)
    record("mass balance to 1e-12", float(result.diagnostics["mass_balance_residual"]) < 1e-12)
    record(
        "every converged phase passes its own stability test",
        result.diagnostics["post_split_status"] == "stable",
    )

    dense, light, light_fraction = feos_tp_flash(pressure_Pa, (0.5, 0.5))
    water_rich, hexane_rich = ours["liquid1"], ours["liquid2"]
    print(
        f"    chemthermo liquid1 x = ({water_rich[0][0]:.10f}, {water_rich[0][1]:.10e}), "
        f"rho = {water_rich[1]:,.4f}, fraction {water_rich[2]:.10f}\n"
        f"    FeOs       'liquid' x = ({dense[0][0]:.10f}, {dense[0][1]:.10e}), "
        f"rho = {dense[1]:,.4f}\n"
        f"    chemthermo liquid2 x = ({hexane_rich[0][0]:.10e}, {hexane_rich[0][1]:.10f}), "
        f"rho = {hexane_rich[1]:,.4f}, fraction {hexane_rich[2]:.10f}\n"
        f"    FeOs        'vapor' x = ({light[0][0]:.10e}, {light[0][1]:.10f}), "
        f"rho = {light[1]:,.4f}, fraction {light_fraction:.10f}"
    )
    print("    (FeOs's container names its two phases 'liquid' and 'vapor'; both are liquids)")
    record("FeOs's two phases are both liquid densities", dense[1] > 5000.0 and light[1] > 5000.0)
    record(
        "compositions agree with FeOs's flash to 1e-8",
        float(np.max(np.abs(water_rich[0] - dense[0]))) < 1e-8
        and float(np.max(np.abs(hexane_rich[0] - light[0]))) < 1e-8,
    )
    record(
        "densities agree with FeOs's flash to 1e-6 relative",
        abs(water_rich[1] - dense[1]) / dense[1] < 1e-6
        and abs(hexane_rich[1] - light[1]) / light[1] < 1e-6,
    )
    record(
        "the phase amounts agree with FeOs's flash to 1e-8",
        abs(hexane_rich[2] - light_fraction) < 1e-8,
    )

    for name, (composition, _density, _fraction) in ours.items():
        value = kappa(pressure_Pa, composition)
        print(f"    kappa({name}) = {value:.3e}  (< {KAPPA_LIQUID_THRESHOLD} is a liquid)")
        record(
            f"{name} is a liquid by an independently computed kappa",
            0.0 < value < KAPPA_LIQUID_THRESHOLD,
        )


def the_chemical_potentials(section: int, pressure_Pa: float, label: str) -> None:
    print(f"\n{section}) FeOs's chemical potentials at chemthermo's phases, {label}")
    published_a, published_b = pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL
    feos_a, feos_b = feos_universal_constants()
    residuals: dict[str, float] = {}
    try:
        for tag, tables in (("as shipped", None), ("matched", (feos_a, feos_b))):
            if tables is not None:
                pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = tables
            else:
                pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = published_a, published_b
            _result, ours = our_phases(pressure_Pa, (0.5, 0.5))
            water_rich, hexane_rich = ours["liquid1"], ours["liquid2"]
            residuals[tag] = float(
                np.max(
                    np.abs(
                        feos_reduced_potentials(water_rich[1], water_rich[0])
                        - feos_reduced_potentials(hexane_rich[1], hexane_rich[0])
                    )
                )
            )
    finally:
        pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = published_a, published_b

    for tag, value in residuals.items():
        print(f"    max_i |mu_i^I - mu_i^II| / RT, {tag:<10} constants: {value:.3e}")
    print(
        "    The 'as shipped' number is floored by the universal-constants table (ten\n"
        "    printed figures against FeOs's fourteen), not by either solver."
    )
    record(
        "FeOs's chemical potentials agree to 1e-8 with matched constants",
        residuals["matched"] < 1e-8,
    )


def three_feeds() -> None:
    print("\n6) Three feeds, one tie line")
    dense, light, _fraction = feos_tp_flash(ATMOSPHERE_PA, (0.5, 0.5))
    for feed in ((0.5, 0.5), (0.2, 0.8), (0.8, 0.2)):
        result, ours = our_phases(ATMOSPHERE_PA, feed)
        z = np.asarray(feed, dtype=float) / float(sum(feed))
        beta = ours["liquid2"][2]
        lever = float(
            np.max(np.abs(z - ((1.0 - beta) * ours["liquid1"][0] + beta * ours["liquid2"][0])))
        )
        same = float(np.max(np.abs(ours["liquid1"][0] - dense[0]))) < 1e-8 and (
            float(np.max(np.abs(ours["liquid2"][0] - light[0]))) < 1e-8
        )
        print(
            f"    z = ({feed[0]:.1f}, {feed[1]:.1f}): liquid2 fraction {beta:.8f}, "
            f"lever-rule residual {lever:.1e}"
        )
        record(f"z = {feed} returns FeOs's tie line to 1e-8", same)
        record(f"z = {feed} satisfies the lever rule to 1e-12", lever < 1e-12)

    for feed in ((0.2, 0.8), (0.8, 0.2)):
        try:
            feos_tp_flash(ATMOSPHERE_PA, feed)
        except Exception as error:  # noqa: BLE001 - reporting the reference's own limit
            print(f"    FeOs's own tp_flash at z = {feed} raises: {error}")


def negative_control() -> None:
    print("\n7) Negative control")
    _result, ours = our_phases(ATMOSPHERE_PA, (0.5, 0.5))
    perturbed = ct.PCSAFTParameters.from_records(
        [
            {
                "name": "Water",
                "m": 1.0656,
                "sigma_A": 3.0007,
                "epsilon_k_K": 366.51,
                "association": {"kappa_ab": 0.034868, "epsilon_ab_k_K": 2500.7 * 1.01},
            },
            {"name": "n-Hexane", "m": 3.0576, "sigma_A": 3.7983, "epsilon_k_K": 236.77},
        ]
    )
    shifted = ct.flash_tp(
        ct.Mixture.from_database(list(NAMES), [0.5, 0.5], normalize=True),
        temperature_K=TEMPERATURE_K,
        pressure_Pa=ATMOSPHERE_PA,
        eos=PCSAFTEOS(components=NAMES, parameters=perturbed),
    )
    moved = abs(
        float(shifted.phases["liquid2"].composition.fractions[0]) - float(ours["liquid2"][0][0])
    )
    print(f"    eps^AB + 1 % moves x_water in the hexane-rich phase by {moved:.3e}")
    record("a 1 % change in eps^AB moves the tie line", moved > 1e-4)


def main() -> None:
    if si is None:
        print("feos is not installed; skipping. Install with: pip install -e '.[validation]'")
        return
    print("PC-SAFT liquid-liquid equilibrium versus FeOs (validation Case P-8, ADR-0019)")
    print("=" * 78)
    the_state()
    the_tie_line(2, ATMOSPHERE_PA, "1 atm")
    the_chemical_potentials(3, ATMOSPHERE_PA, "1 atm")
    the_tie_line(4, HIGH_PRESSURE_PA, "1 MPa")
    the_chemical_potentials(5, HIGH_PRESSURE_PA, "1 MPa")
    three_feeds()
    negative_control()

    print("\n" + "=" * 78)
    print(
        "k_ij = 0 throughout. With no fitted binary parameter PC-SAFT is known to be poor\n"
        "for water / hydrocarbon mutual solubilities - about 1.7e-05 mole fraction hexane\n"
        "in water and 6.3e-03 water in hexane here, against commonly tabulated values near\n"
        "2e-06 and 5e-04 at 298 K (not verified against a primary source, and nothing above\n"
        "asserts them). This script compares two implementations of one model."
    )
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
