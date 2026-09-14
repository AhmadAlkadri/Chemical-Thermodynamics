"""A polymer as a PC-SAFT component: polyethylene / n-pentane at 453 K.

What this shows
---------------
Three things that did not exist before ADR-0022, used end to end:

1. **A polymer parameter record.** Polymer PC-SAFT parameters are published
   per unit *mass* - ``m/M`` in mol/g - because the chain length depends on the
   molar mass of the particular sample. ``PCSAFTRecord(segments_per_g=...,
   MW_g_mol=...)`` says exactly that, and derives ``m = (m/M) * Mw`` once, in
   the record. For the ``Mw = 16400 g/mol`` sample below that is
   ``m = 431.32`` segments, against ``2.6896`` for n-pentane: a size ratio of
   160 to 1.
2. **A component the databank does not carry.**
   ``Component.custom("Polyethylene", mw_kg_per_mol=16.4, volatile=False)``
   builds one with no critical constants at all - a polymer has none - and
   ``volatile=False`` gives it the fixed Wilson K-value *estimate* ``1e-10``
   ("essentially absent from the vapour-like trial") so the deterministic
   stability trial set stays complete.
3. **``stability_tp`` and ``flash_tp`` on the result**, which return a verified
   liquid-liquid split: the polymer-lean solvent phase and the polymer-rich
   one, with a mass balance that closes despite a polymer mole fraction of
   1e-05, a negative Gibbs change, and both phases stable when re-tested.

The script prints, in order: the derived segment number, the pure melt density
over 1-30 MPa, a pressure scan of the stability verdict that brackets the
model's cloud point, the split itself with its three verification residuals,
the same scan in *temperature* (the LCST-type direction: this system demixes on
heating), and the effect of ``k_ij``.

**Read this before reading any number.** The polyethylene parameters come from
one open secondary source (Martini, Cismondi, Barbosa & Brignole,
*Sep. Sci. Technol.* **44** (2009), author manuscript in the CONICET
repository) that cites Gross & Sadowski, *IECR* **41** (2002) 1084. That
primary table was **not** read - it is paywalled - and no second open source
printing the same three numbers was found. chemthermo therefore packages **no**
polymer parameters: the values live in
``tests/fixtures/pcsaft/martini2009_polymers.json`` with their provenance
caveat, this script reads them from there, and everything below is a statement
about *this model with these inputs*, never about polyethylene. The polymer is
also treated as **monodisperse** - one chain length, one component - while the
sample this ``k_ij`` was fitted to has a polydispersity of 1.16.

Needs no optional dependency. All inputs are SI; compositions are mole
fractions. ``--full`` adds a second molar mass (``Mw = 53000``, ``m = 1393.9``,
where ``exp(ln phi)`` underflows and ADR-0022's log-space route carries the
flash) and a finer cloud-point bisection; the default takes about five seconds.

Cross-checked against FeOs in
``examples/validation/20_pcsaft_polymer_vs_feos.py`` (validation Cases P-12 and
P-13).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.parameters import PCSAFTParameters, PCSAFTRecord

FIXTURE = (
    Path(__file__).resolve().parents[2]
    / "tests"
    / "fixtures"
    / "pcsaft"
    / "martini2009_polymers.json"
)

TEMPERATURE_K = 453.0
PENTANE_MW_G_MOL = 72.146
#: Table 3 of the cited source, for the Mw = 16400 sample.
KIJ = -0.006


def polymer_row() -> dict:
    with FIXTURE.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return next(row for row in payload["polymers"] if row["name"] == "Polyethylene")


def parameters(mw_g_mol: float) -> PCSAFTParameters:
    """The two-component parameter set, polymer given per unit mass."""
    row = polymer_row()
    return PCSAFTParameters.from_records(
        [
            PCSAFTRecord(
                name="Polyethylene",
                segments_per_g=row["segments_per_g"],
                MW_g_mol=mw_g_mol,
                sigma_A=row["sigma_A"],
                epsilon_k_K=row["epsilon_k_K"],
                source="Martini et al. 2009 Table 1 citing Gross & Sadowski 2002",
            ),
            PCSAFTRecord(
                name="n-Pentane",
                m=2.6896,
                sigma_A=3.7729,
                epsilon_k_K=231.20,
                MW_g_mol=PENTANE_MW_G_MOL,
                source="Gross & Sadowski 2001 Table 1 (the packaged record)",
            ),
        ]
    )


def polymer(mw_g_mol: float) -> ct.Component:
    return ct.Component.custom(
        "Polyethylene",
        mw_kg_per_mol=mw_g_mol / 1000.0,
        formula="(C2H4)n",
        volatile=False,
        source="see tests/fixtures/pcsaft/martini2009_polymers.json",
    )


def feed(weight_fraction: float, mw_g_mol: float) -> list[float]:
    """Mole fractions from a polymer **mass** fraction - the usual polymer basis."""
    moles_polymer = weight_fraction / mw_g_mol
    moles_solvent = (1.0 - weight_fraction) / PENTANE_MW_G_MOL
    total = moles_polymer + moles_solvent
    return [moles_polymer / total, moles_solvent / total]


def mixture(weight_fraction: float, mw_g_mol: float) -> ct.Mixture:
    return ct.Mixture.from_components(
        [polymer(mw_g_mol), ct.Component.from_database("n-Pentane")],
        feed(weight_fraction, mw_g_mol),
        normalize=True,
    )


def section(title: str) -> None:
    print("\n" + "-" * 78)
    print(title)
    print("-" * 78)


def the_record(mw_g_mol: float) -> None:
    section(f"1. The segments-per-mass record, Mw = {mw_g_mol:g} g/mol")
    row = polymer_row()
    record = parameters(mw_g_mol).record("Polyethylene")
    print(f"    m/M   = {row['segments_per_g']} mol/g   (as published)")
    print(f"    Mw    = {mw_g_mol:g} g/mol       (this sample)")
    print(f"    m     = {record.m:.4f}          (derived: m/M * Mw)")
    print(f"    sigma = {record.sigma_A} A,  eps/k = {record.epsilon_k_K} K")
    print(f"    n-pentane for comparison: m = 2.6896  ->  size ratio {record.m / 2.6896:.0f} : 1")


def the_melt(mw_g_mol: float) -> None:
    section(f"2. The pure melt at {TEMPERATURE_K:g} K")
    eos = PCSAFTEOS(components=("Polyethylene",), parameters=parameters(mw_g_mol))
    print("    P / MPa    roots    rho / (mol/m^3)      rho / (g/cm^3)")
    for pressure_Pa in (1.0e6, 1.0e7, 2.0e7, 3.0e7):
        roots = eos.density_roots(
            temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=[1.0]
        )
        density = roots[-1]
        print(
            f"    {pressure_Pa / 1e6:7.1f}   {len(roots):5d}    {density:15.6f}"
            f"      {density * mw_g_mol / 1e6:.6f}"
        )
    print(
        "\n    Commonly tabulated polyethylene melt densities near 450 K are around\n"
        "    0.77-0.80 g/cm^3. That remark is unverified here and is printed only to\n"
        "    say the magnitude is not absurd; nothing is asserted against it."
    )


def verdict(pressure_Pa: float, mw_g_mol: float, kij: float, temperature_K: float) -> str:
    return ct.stability_tp(
        mixture(0.05, mw_g_mol),
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=PCSAFTEOS(parameters=parameters(mw_g_mol), kij=kij),
    ).status


def the_pressure_scan(mw_g_mol: float, full: bool) -> None:
    section(f"3. Stability verdict against pressure, 5 wt% polymer, k_ij = {KIJ}")
    pressures = (5.0e6, 8.0e6, 1.0e7, 1.5e7, 3.0e7)
    for pressure_Pa in pressures:
        print(
            f"    P = {pressure_Pa / 1e6:5.1f} MPa  ->  {verdict(pressure_Pa, mw_g_mol, KIJ, TEMPERATURE_K)}"
        )
    print(
        "\n    The switch between 8 and 10 MPa is the model's cloud point for this feed:\n"
        "    compress the solution and it becomes one phase."
    )
    if not full:
        print("    (pass --full to bisect it)")
        return
    low, high = 8.0e6, 1.0e7
    while high - low > 1.0:
        middle = 0.5 * (low + high)
        if verdict(middle, mw_g_mol, KIJ, TEMPERATURE_K) == "unstable":
            low = middle
        else:
            high = middle
    print(f"    bisected cloud point: {0.5 * (low + high) / 1e6:.6f} MPa")


def the_split(mw_g_mol: float) -> None:
    section(f"4. flash_tp at 8 MPa, 5 wt% polymer, Mw = {mw_g_mol:g}")
    feed_mixture = mixture(0.05, mw_g_mol)
    eos = PCSAFTEOS(parameters=parameters(mw_g_mol), kij=KIJ)
    result = ct.flash_tp(feed_mixture, temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=eos)
    diagnostics = result.diagnostics
    print(f"    feed x_polymer          = {feed_mixture.fractions[0]:.6e}")
    print(f"    phase regime            = {diagnostics['phase_regime']}")
    print(f"    named by                = {diagnostics['phase_label_method']}")
    print(f"    vapor_fraction          = {result.vapor_fraction}")
    for name, phase in result.phases.items():
        fractions = phase.composition.fractions
        density = eos.density_roots(
            mixture=feed_mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=8.0e6,
            composition=list(fractions),
        )[-1]
        print(
            f"      {name}: x_polymer = {fractions[0]:.6e}   fraction = "
            f"{result.phase_fractions[name]:.6f}   rho = {density:.2f} mol/m^3"
        )
    print("    verification:")
    for key in ("mass_balance_residual", "fugacity_residual", "delta_g_split_rt"):
        print(f"      {key:24s} = {float(diagnostics[key]):.3e}")
    print(f"      post_split_status        = {diagnostics['post_split_status']}")
    print(
        "\n    The mass balance closes to 1e-20 even though one phase holds the polymer\n"
        "    at a mole fraction of 1e-05: the extreme size ratio costs nothing here."
    )


def the_temperature_scan(mw_g_mol: float) -> None:
    section("5. Stability verdict against temperature at 10 MPa (LCST-type direction)")
    for temperature_K in (400.0, 425.0, 450.0, 460.0):
        print(f"    T = {temperature_K:6.1f} K  ->  {verdict(1.0e7, mw_g_mol, KIJ, temperature_K)}")
    print(
        "\n    One phase when cold, two when hot: the polymer comes out of solution on\n"
        "    heating. The cited manuscript's Figure 1 describes the same direction (the\n"
        "    region above its cloud-point curve is single phase). That comparison is\n"
        "    qualitative - the figure was not digitized and no number here is fitted to it."
    )


def the_kij_effect(mw_g_mol: float) -> None:
    section("6. The effect of k_ij, at 15 MPa and 453 K")
    for kij in (-0.006, 0.0, 0.02):
        print(f"    k_ij = {kij:+.3f}  ->  {verdict(1.5e7, mw_g_mol, kij, TEMPERATURE_K)}")
    print(
        "\n    A larger k_ij means a worse solvent, so a higher pressure is needed to keep\n"
        "    the polymer dissolved. The source states the same trend."
    )


def the_long_chain() -> None:
    section("7. Mw = 53000 (m = 1393.9): where exp(ln phi) stops existing")
    mw_g_mol = 53000.0
    eos = PCSAFTEOS(components=("Polyethylene", "n-Pentane"), parameters=parameters(mw_g_mol))
    x = feed(0.05, mw_g_mol)
    density = eos.density_roots(temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, composition=x)[-1]
    ln_phi = eos.ln_fugacity_coefficients(
        temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=x
    )
    print(f"    ln phi_polymer = {ln_phi[0]:.4f}")
    print("    exp(that)      = 0.0 exactly - the smallest positive double is exp(-744.44)")
    print(
        "    ADR-0022 keeps the flash in log space where that happens, so the split\n"
        "    below runs instead of the model being refused:"
    )
    feed_mixture = mixture(0.05, mw_g_mol)
    result = ct.flash_tp(
        feed_mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=8.0e6,
        eos=PCSAFTEOS(parameters=parameters(mw_g_mol), kij=KIJ),
    )
    for name, phase in result.phases.items():
        print(
            f"      {name}: x_polymer = {phase.composition.fractions[0]:.6e}   "
            f"fraction = {result.phase_fractions[name]:.6f}"
        )
    print(
        f"      dG/RT = {float(result.diagnostics['delta_g_split_rt']):.3e}, "
        f"post-split {result.diagnostics['post_split_status']}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="also bisect the cloud point and run the Mw = 53000 chain",
    )
    # `parse_known_args`, not `parse_args`: `tests/test_examples.py` runs this
    # script via `runpy.run_path` with pytest's own `sys.argv` still in place.
    args, _unknown = parser.parse_known_args()

    print("=" * 78)
    print("PC-SAFT for a polymer/solvent mixture: polyethylene / n-pentane, 453 K")
    print("=" * 78)
    print(
        "Polymer parameters: as tabulated by Martini et al. (2009) citing Gross &\n"
        "Sadowski (2002); NOT verified against that primary table, which is paywalled,\n"
        "and no second open source printing them was found. They are a cited test\n"
        "fixture, not packaged runtime data - chemthermo packages no polymer\n"
        "parameters. The polymer is modelled as monodisperse. Nothing below is\n"
        "compared against measurement."
    )

    mw_g_mol = 16400.0
    the_record(mw_g_mol)
    the_melt(mw_g_mol)
    the_pressure_scan(mw_g_mol, args.full)
    the_split(mw_g_mol)
    the_temperature_scan(mw_g_mol)
    the_kij_effect(mw_g_mol)
    if args.full:
        the_long_chain()
    else:
        print("\n  (pass --full for the cloud-point bisection and the Mw = 53000 chain)")

    print("\n" + "=" * 78)
    print(
        "Every number above is this model with these inputs. The polymer parameters\n"
        "rest on a single secondary source; the k_ij was fitted by that source to\n"
        "cloud-point data this script never touches; and a real sample is polydisperse\n"
        "while this component is not."
    )


if __name__ == "__main__":
    main()
