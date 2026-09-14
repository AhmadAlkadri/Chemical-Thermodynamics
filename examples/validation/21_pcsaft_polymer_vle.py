"""Polymer/solvent **vapour-liquid** equilibrium, with PASS/FAIL on every check.

What is being validated
-----------------------
ADR-0024 makes the phi-phi split solvable in **log mole numbers**, which is
what a polymer/solvent VLE needs: the vapour above a polyethylene melt holds a
polymer mole fraction of ``exp(-450)``, a number no linear parametrization can
represent, and the tangent-plane minimizer that would seed the split is an
essentially pure melt whose K-values span ``1e+180`` and bracket no
Rachford-Rice root at all. Before this slice ``flash_tp`` raised
``ConvergenceError`` on every state below n-pentane's saturation pressure
(validation Case P-13 (vi), pinned as a defect); it now returns a solvent
vapour in equilibrium with a solvent-swollen melt (Case P-14).

Five routes, deliberately different in kind:

1. **The split itself** at 0.5, 1 and 2 MPa: both phases named from a measured
   compressibility (ADR-0017), the melt's solvent content, the vapour's
   polymer content as ``ln y``, mass balance, equal-fugacity residual, the
   Gibbs change against the single-phase feed, and the post-split stability
   test of both phases.
2. **A one-dimensional equal-fugacity solve written in this script**, with the
   vapour taken to be *exactly* pure solvent: one unknown, ``ln x_solvent`` in
   the melt, solved against the pure-solvent vapour fugacity by a
   finite-difference Newton. It uses no part of the flash - no stability test,
   no Rachford-Rice, no Newton stage, no phase-count logic - only
   ``PCSAFTEOS``'s ``(T, rho, x)`` interface and ``density_roots``.
3. **The polymer's own equal-fugacity condition**, checked in logarithms, which
   is the equation the linear split could not even write down.
4. **FeOs's chemical potentials at chemthermo's converged phases and
   densities** (feos-org/feos, MIT OR Apache-2.0; the same Gross & Sadowski
   model in Rust with every derivative by automatic differentiation). This does
   not depend on FeOs's own flash converging - it does not, on this system.
5. **The ternary** polyethylene / n-pentane / n-hexane at 3 MPa, the second
   state Case P-13 pinned as a runaway.

``k_ij`` cannot be given to FeOs's PC-SAFT from Python in feos 0.10.1
(``EquationOfState.pcsaft`` raises "missing field ``k_ij``" for every
serialization tried), so route 4 runs at ``k_ij = 0`` on **both** sides, and
the flash it checks is re-run at ``k_ij = 0`` for the purpose. Routes 1, 2, 3
and 5 use the fitted ``k_ij = -0.006``. As in Cases P-6 to P-12, the 42
universal constants of the 2001 dispersion term are not shared (chemthermo
packages the ten figures as printed, FeOs hard-codes fourteen), so route 4 is
reported twice, as shipped and with FeOs's table substituted in.

**The polymer parameters are not verified data.** They are as tabulated by
Martini, Cismondi, Barbosa & Brignole, *Sep. Sci. Technol.* **44** (2009),
citing Gross & Sadowski, *IECR* **41** (2002) 1084 - a paywalled table that was
not read. They live in ``tests/fixtures/pcsaft/martini2009_polymers.json``, are
never packaged, and nothing here is compared against measurement.

``--full`` adds the 25-point pressure scan from 0.3 to 12 MPa (the
vapour-liquid / liquid-liquid / single-liquid verdict sequence) and the
bisected vapour-liquid to liquid-liquid boundary. The default takes about ten
seconds.

Requires the optional ``feos`` dependency for route 4 only::

    pip install -e ".[validation]"

Routes 1, 2, 3 and 5 run without it; the script says so and still exits 0.
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
from chemthermo.models.base import KAPPA_LIQUID_THRESHOLD
from chemthermo.parameters import PCSAFTParameters, PCSAFTRecord

REPO = Path(__file__).resolve().parents[2]
FIXTURE = REPO / "tests" / "fixtures" / "pcsaft" / "martini2009_polymers.json"

TEMPERATURE_K = 453.0
PE_MW_G_MOL = 16400.0
PENTANE_MW_G_MOL = 72.146
HEXANE_MW_G_MOL = 86.177
KIJ = -0.006
#: n-pentane and n-hexane, Gross & Sadowski (2001) Table 1.
PENTANE = (2.6896, 3.7729, 231.20)
HEXANE = (3.0576, 3.7983, 236.77)

#: The three states of route 1. All are below the solvent's saturation
#: pressure at 453 K, where the mixture has a vapour density root.
PRESSURES_PA = (5.0e5, 1.0e6, 2.0e6)

NEWTON_TOL = 1e-12
#: The melt's solvent mole fraction from the flash and from the 1-D solve.
COMPOSITION_TOL = 1e-10
#: Equal fugacity in logarithms, including the polymer's.
LOG_FUGACITY_TOL = 1e-8
POTENTIAL_TOL = 1e-8
MASS_BALANCE_TOL = 1e-12

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
# The model, from one set of numbers
# ---------------------------------------------------------------------------


def polymer_row() -> dict:
    with FIXTURE.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return next(row for row in payload["polymers"] if row["name"] == "Polyethylene")


ROW = polymer_row()


def polymer_record() -> PCSAFTRecord:
    return PCSAFTRecord(
        name="Polyethylene",
        segments_per_g=ROW["segments_per_g"],
        MW_g_mol=PE_MW_G_MOL,
        sigma_A=ROW["sigma_A"],
        epsilon_k_K=ROW["epsilon_k_K"],
    )


def binary_parameters() -> PCSAFTParameters:
    return PCSAFTParameters.from_records(
        [
            polymer_record(),
            PCSAFTRecord(
                name="n-Pentane",
                m=PENTANE[0],
                sigma_A=PENTANE[1],
                epsilon_k_K=PENTANE[2],
                MW_g_mol=PENTANE_MW_G_MOL,
            ),
        ]
    )


def ternary_parameters() -> PCSAFTParameters:
    return PCSAFTParameters.from_records(
        [
            polymer_record(),
            PCSAFTRecord(
                name="n-Pentane", m=PENTANE[0], sigma_A=PENTANE[1], epsilon_k_K=PENTANE[2]
            ),
            PCSAFTRecord(name="n-Hexane", m=HEXANE[0], sigma_A=HEXANE[1], epsilon_k_K=HEXANE[2]),
        ]
    )


def polymer_component() -> ct.Component:
    return ct.Component.custom(
        "Polyethylene",
        mw_kg_per_mol=PE_MW_G_MOL / 1000.0,
        formula="(C2H4)n",
        volatile=False,
        source="see tests/fixtures/pcsaft/martini2009_polymers.json",
    )


def binary_mixture(weight_fraction: float = 0.05) -> ct.Mixture:
    polymer = weight_fraction / PE_MW_G_MOL
    solvent = (1.0 - weight_fraction) / PENTANE_MW_G_MOL
    total = polymer + solvent
    return ct.Mixture.from_components(
        [polymer_component(), ct.Component.from_database("n-Pentane")],
        [polymer / total, solvent / total],
        normalize=True,
    )


def ternary_mixture() -> ct.Mixture:
    amounts = [0.05 / PE_MW_G_MOL, 0.475 / PENTANE_MW_G_MOL, 0.475 / HEXANE_MW_G_MOL]
    total = sum(amounts)
    return ct.Mixture.from_components(
        [
            polymer_component(),
            ct.Component.from_database("n-Pentane"),
            ct.Component.from_database("n-Hexane"),
        ],
        [value / total for value in amounts],
        normalize=True,
    )


def binary_eos(kij: float = KIJ) -> PCSAFTEOS:
    return PCSAFTEOS(parameters=binary_parameters(), kij=kij)


def ln_y_polymer(result: ct.FlashResult, phase_name: str) -> float:
    """``ln y_polymer``, from the mole fraction or from the log-space diagnostic.

    ADR-0024 decision 3: a mole fraction the model puts outside the
    exponential's range is reported as an exact ``0.0`` in the composition, and
    its logarithm is in ``diagnostics["log_space_ln_x_min"]``.
    """
    value = float(result.phases[phase_name].composition.fractions[0])
    if value > 0.0:
        return math.log(value)
    return float(result.diagnostics["log_space_ln_x_min"])


def kappa(
    kij: float, pressure_Pa: float, x: Sequence[float]
) -> tuple[tuple[float, ...], list[float]]:
    """``P / (rho dP/drho)`` from the public pressure routine (ADR-0017)."""
    bound = PCSAFTEOS(
        components=("Polyethylene", "n-Pentane"), parameters=binary_parameters(), kij=kij
    )
    roots = bound.density_roots(
        temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=list(x)
    )
    values = []
    for density in roots:
        step = 1e-4 * density
        slope = (
            bound.pressure_Pa(
                temperature_K=TEMPERATURE_K, density_mol_m3=density + step, composition=list(x)
            )
            - bound.pressure_Pa(
                temperature_K=TEMPERATURE_K, density_mol_m3=density - step, composition=list(x)
            )
        ) / (2.0 * step)
        values.append(pressure_Pa / (density * slope))
    return roots, values


# ---------------------------------------------------------------------------
# Route 1: the split
# ---------------------------------------------------------------------------


def the_splits() -> dict[float, ct.FlashResult]:
    section("1. flash_tp on polyethylene(16400) / n-pentane, 5 wt%, 453 K, k_ij = -0.006")
    mixture = binary_mixture()
    eos = binary_eos()
    print(f"  feed z = {tuple(mixture.fractions)}")
    print(
        f"\n  {'P/MPa':>7} {'phases':>15} {'beta_vapor':>13} {'x_C5 (melt)':>14}"
        f" {'ln y_PE':>10} {'resid':>9} {'dG/RT':>11}"
    )
    results: dict[float, ct.FlashResult] = {}
    for pressure_Pa in PRESSURES_PA:
        result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
        results[pressure_Pa] = result
        melt = "liquid" if "liquid" in result.phases else sorted(result.phases)[0]
        vapor = "vapor" if "vapor" in result.phases else sorted(result.phases)[-1]
        print(
            f"  {pressure_Pa / 1e6:7.2f} {'+'.join(sorted(result.phases)):>15}"
            f" {result.vapor_fraction if result.vapor_fraction is not None else float('nan'):13.10f}"
            f" {result.phases[melt].composition.fractions[1]:14.10f}"
            f" {ln_y_polymer(result, vapor):10.2f}"
            f" {float(result.diagnostics['fugacity_residual']):9.1e}"
            f" {float(result.diagnostics['delta_g_split_rt']):11.3e}"
        )

    print()
    for pressure_Pa, result in results.items():
        label = f"{pressure_Pa / 1e6:.1f} MPa"
        record_check(
            f"{label}: a vapour and a liquid", sorted(result.phases) == ["liquid", "vapor"]
        )
        record_check(
            f"{label}: named from a measured compressibility",
            result.diagnostics["phase_label_method"] == "compressibility",
        )
        record_check(
            f"{label}: solved in log mole numbers",
            result.diagnostics["converged_stage"] == "second-order-log",
        )
        record_check(
            f"{label}: mass balance < {MASS_BALANCE_TOL:g}",
            float(result.diagnostics["mass_balance_residual"]) < MASS_BALANCE_TOL,
        )
        record_check(
            f"{label}: equal fugacity < {LOG_FUGACITY_TOL:g}",
            float(result.diagnostics["fugacity_residual"]) < LOG_FUGACITY_TOL,
        )
        record_check(
            f"{label}: dG_split/RT < 0", float(result.diagnostics["delta_g_split_rt"]) < 0.0
        )
        record_check(
            f"{label}: both phases post-split stable",
            result.diagnostics["post_split_status"] == "stable",
        )
        # The compressibility identity, computed here from the public pressure
        # routine rather than from `phase_identity`'s own derivative.
        _roots, vapor_kappa = kappa(KIJ, pressure_Pa, result.phases["vapor"].composition.fractions)
        _melt_roots, melt_kappa = kappa(
            KIJ, pressure_Pa, result.phases["liquid"].composition.fractions
        )
        record_check(
            f"{label}: vapour kappa ~ 1 (measured {vapor_kappa[0]:.4f})",
            vapor_kappa[0] > KAPPA_LIQUID_THRESHOLD,
        )
        record_check(
            f"{label}: melt kappa << 0.5 (measured {melt_kappa[-1]:.4f})",
            0.0 < melt_kappa[-1] < KAPPA_LIQUID_THRESHOLD,
        )
    return results


# ---------------------------------------------------------------------------
# Route 2 and 3: an independent one-dimensional solve
# ---------------------------------------------------------------------------


def the_one_dimensional_solve(results: dict[float, ct.FlashResult]) -> None:
    section("2. An independent 1-D equal-fugacity solve, vapour taken as exactly pure solvent")
    print(
        "  One unknown, ln x_solvent in the melt, against\n"
        "      ln phi_s^V(pure n-pentane vapour) = ln x_s + ln phi_s^L(x)\n"
        "  by a finite-difference Newton. Nothing of flash_tp is used."
    )
    bound = PCSAFTEOS(
        components=("Polyethylene", "n-Pentane"), parameters=binary_parameters(), kij=KIJ
    )

    def ln_phi(pressure_Pa: float, x: Sequence[float], *, vapor: bool) -> np.ndarray:
        roots = bound.density_roots(
            temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=list(x)
        )
        density = roots[0] if vapor else roots[-1]
        return np.asarray(
            bound.ln_fugacity_coefficients(
                temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=list(x)
            )
        )

    print(f"\n  {'P/MPa':>7} {'x_C5 (1-D)':>16} {'x_C5 (flash)':>16} {'|difference|':>13}")
    for pressure_Pa, result in results.items():
        target = float(ln_phi(pressure_Pa, [0.0, 1.0], vapor=True)[1])

        def residual(
            ln_x: float, pressure_Pa: float = pressure_Pa, target: float = target
        ) -> float:
            solvent = math.exp(ln_x)
            x = [1.0 - solvent, solvent]
            return ln_x + float(ln_phi(pressure_Pa, x, vapor=False)[1]) - target

        # Started a thousandth away from the flash's answer so the Newton has
        # real work to do, and bracketed nowhere near it by construction.
        ln_x = math.log(float(result.phases["liquid"].composition.fractions[1])) - 1e-3
        for _ in range(60):
            value = residual(ln_x)
            step = 1e-7
            slope = (residual(ln_x + step) - residual(ln_x - step)) / (2.0 * step)
            correction = -value / slope
            ln_x += correction
            if abs(correction) < NEWTON_TOL:
                break

        ours = float(result.phases["liquid"].composition.fractions[1])
        theirs = math.exp(ln_x)
        difference = abs(ours - theirs)
        print(f"  {pressure_Pa / 1e6:7.2f} {theirs:16.12f} {ours:16.12f} {difference:13.2e}")
        record_check(
            f"{pressure_Pa / 1e6:.1f} MPa: melt x_solvent matches the 1-D solve < {COMPOSITION_TOL:g}",
            difference < COMPOSITION_TOL,
        )
        record_check(
            f"{pressure_Pa / 1e6:.1f} MPa: the 1-D solve is converged",
            abs(residual(ln_x)) < LOG_FUGACITY_TOL,
        )

        # The flash's own vapour-phase solvent fugacity against the pure-solvent
        # one: the vapour is pure to 1e-196, so these must agree to round-off.
        y = list(result.phases["vapor"].composition.fractions)
        flash_vapor = math.log(y[1]) + float(ln_phi(pressure_Pa, y, vapor=True)[1])
        record_check(
            f"{pressure_Pa / 1e6:.1f} MPa: flash vapour solvent fugacity = pure-solvent value",
            abs(flash_vapor - target) < LOG_FUGACITY_TOL,
        )

    section("3. The polymer's own equal-fugacity condition, in logarithms")
    print(
        "  ln x_PE + ln phi_PE^L(x)  =  ln y_PE + ln phi_PE^V(y), with ln y_PE ~ -450.\n"
        "  This is the equation the linear split cannot write down at all."
    )
    print(f"\n  {'P/MPa':>7} {'ln f_PE (melt)':>17} {'ln f_PE (vapour)':>18} {'|difference|':>13}")
    for pressure_Pa, result in results.items():
        x = list(result.phases["liquid"].composition.fractions)
        y = list(result.phases["vapor"].composition.fractions)
        melt = math.log(x[0]) + float(ln_phi(pressure_Pa, x, vapor=False)[0])
        vapor = ln_y_polymer(result, "vapor") + float(ln_phi(pressure_Pa, y, vapor=True)[0])
        print(f"  {pressure_Pa / 1e6:7.2f} {melt:17.10f} {vapor:18.10f} {abs(melt - vapor):13.2e}")
        record_check(
            f"{pressure_Pa / 1e6:.1f} MPa: polymer log fugacities agree < {LOG_FUGACITY_TOL:g}",
            abs(melt - vapor) < LOG_FUGACITY_TOL,
        )


# ---------------------------------------------------------------------------
# Route 4: FeOs chemical potentials
# ---------------------------------------------------------------------------


def feos_constants() -> tuple[np.ndarray, np.ndarray]:
    """FeOs's fourteen-figure tables, from the Case P-6 test module (single copy)."""
    path = REPO / "tests" / "validation" / "test_pcsaft_association_vs_feos.py"
    spec = importlib.util.spec_from_file_location("_pcsaft_feos_constants", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.FEOS_A_UNIVERSAL, module.FEOS_B_UNIVERSAL


def feos_eos():
    """The reference model, always with ``k_ij = 0`` - see the module docstring."""
    payloads = [
        {
            "identifier": {"name": "Polyethylene"},
            "molarweight": PE_MW_G_MOL,
            "m": ROW["segments_per_g"] * PE_MW_G_MOL,
            "sigma": ROW["sigma_A"],
            "epsilon_k": ROW["epsilon_k_K"],
        },
        {
            "identifier": {"name": "n-Pentane"},
            "molarweight": PENTANE_MW_G_MOL,
            "m": PENTANE[0],
            "sigma": PENTANE[1],
            "epsilon_k": PENTANE[2],
        },
    ]
    records = [PureRecord.from_json_str(json.dumps(payload)) for payload in payloads]
    return EquationOfState.pcsaft(Parameters.from_records(records))


def feos_reduced_potentials(density: float, x: Sequence[float]) -> np.ndarray:
    """``mu_i / RT`` from FeOs at chemthermo's ``(T, rho, x)``, up to one constant.

    The omitted constant is the same for both phases of one flash, so it
    cancels in the difference and the comparison is the reference
    implementation's own equilibrium condition at chemthermo's answer.
    """
    values = np.asarray(x, dtype=float)
    state = State(
        feos_eos(),
        temperature=TEMPERATURE_K * si.KELVIN,
        density=density * (si.MOL / si.METER**3),
        composition=values,
    )
    factor = R_J_PER_MOL_K * TEMPERATURE_K
    z_factor = (state.pressure() / si.PASCAL) / (density * factor)
    mu_res = np.atleast_1d(
        np.asarray(state.chemical_potential(Contributions.Residual) / si.JOULE * si.MOL)
    )
    with np.errstate(divide="ignore"):
        return mu_res / factor - math.log(z_factor) + np.log(values)


def the_chemical_potentials() -> None:
    section("4. FeOs chemical potentials at chemthermo's converged phases (k_ij = 0 on both sides)")
    print(
        "  FeOs's own tp_flash does not converge on this system, so this route uses only\n"
        "  chemthermo's two compositions and two densities and lets FeOs supply the\n"
        "  potentials. A polymer mole fraction of ~1e-196 is an ordinary double and FeOs\n"
        "  takes it as given."
    )
    mixture = binary_mixture()
    eos = binary_eos(0.0)
    shipped = (pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL)
    matched = feos_constants()

    for label, constants in (("as shipped", shipped), ("matched constants", matched)):
        pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = constants
        worst = 0.0
        try:
            for pressure_Pa in PRESSURES_PA:
                result = ct.flash_tp(
                    mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
                )
                potentials = []
                for name in ("liquid", "vapor"):
                    x = list(result.phases[name].composition.fractions)
                    roots = eos.density_roots(
                        mixture=mixture,
                        temperature_K=TEMPERATURE_K,
                        pressure_Pa=pressure_Pa,
                        composition=x,
                    )
                    density = roots[-1] if name == "liquid" else roots[0]
                    potentials.append(feos_reduced_potentials(density, x))
                worst = max(worst, float(np.max(np.abs(potentials[0] - potentials[1]))))
        finally:
            pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = shipped
        print(f"  worst |d mu_i / RT| over the three states, {label}: {worst:.3e}")
        if label == "matched constants":
            record_check(
                f"FeOs potentials equal at both phases < {POTENTIAL_TOL:g}", worst < POTENTIAL_TOL
            )


# ---------------------------------------------------------------------------
# Route 5: the ternary
# ---------------------------------------------------------------------------


def the_ternary() -> None:
    section("5. The ternary polyethylene / n-pentane / n-hexane at 3 MPa (Case P-13 (vi) runaway)")
    mixture = ternary_mixture()
    eos = PCSAFTEOS(parameters=ternary_parameters(), kij=KIJ)
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=3.0e6, eos=eos)
    print(f"  feed z      = {tuple(mixture.fractions)}")
    for name in sorted(result.phases):
        print(
            f"  {name:8s}    = {tuple(result.phases[name].composition.fractions)}"
            f"  (phase fraction {result.phase_fractions[name]:.10f})"
        )
    print(f"  converged stage: {result.diagnostics['converged_stage']}")
    record_check("ternary: two phases", len(result.phases) == 2)
    record_check(
        "ternary: mass balance < 1e-12",
        float(result.diagnostics["mass_balance_residual"]) < MASS_BALANCE_TOL,
    )
    record_check(
        "ternary: equal fugacity < 1e-08",
        float(result.diagnostics["fugacity_residual"]) < LOG_FUGACITY_TOL,
    )
    record_check("ternary: dG_split/RT < 0", float(result.diagnostics["delta_g_split_rt"]) < 0.0)
    record_check("ternary: post-split stable", result.diagnostics["post_split_status"] == "stable")


# ---------------------------------------------------------------------------
# --full: the verdict sequence with pressure
# ---------------------------------------------------------------------------


def the_pressure_scan() -> None:
    section("6. (--full) 25 pressures from 0.3 to 12 MPa: the verdict sequence")
    mixture = binary_mixture()
    eos = binary_eos()
    print(f"  {'P/MPa':>7} {'verdict':>17} {'stage':>19} {'x_C5 (PE-rich)':>16} {'dG/RT':>11}")
    verdicts: list[str] = []
    for pressure_Pa in np.linspace(0.3e6, 12.0e6, 25):
        result = ct.flash_tp(
            mixture, temperature_K=TEMPERATURE_K, pressure_Pa=float(pressure_Pa), eos=eos
        )
        regime = str(result.diagnostics["phase_regime"])
        verdicts.append(regime)
        if len(result.phases) == 1:
            print(f"  {pressure_Pa / 1e6:7.3f} {regime:>17} {'-':>19}")
            continue
        rich = max(result.phases, key=lambda n: result.phases[n].composition.fractions[0])
        print(
            f"  {pressure_Pa / 1e6:7.3f} {regime:>17}"
            f" {str(result.diagnostics.get('converged_stage', 'successive-substitution')):>19}"
            f" {result.phases[rich].composition.fractions[1]:16.10f}"
            f" {float(result.diagnostics['delta_g_split_rt']):11.3e}"
        )
    ordered = [
        regime
        for index, regime in enumerate(verdicts)
        if index == 0 or regime != verdicts[index - 1]
    ]
    print(f"\n  verdict sequence: {' -> '.join(ordered)}")
    record_check(
        "the verdict sequence is VLE -> LLE -> single-phase, with no failures",
        ordered == ["VLE", "LLE", "single-phase"],
    )


def the_boundary() -> None:
    section("7. (--full) The vapour-liquid to liquid-liquid boundary, bisected")
    mixture = binary_mixture()
    eos = binary_eos()

    def regime(pressure_Pa: float) -> str:
        return str(
            ct.flash_tp(
                mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
            ).diagnostics["phase_regime"]
        )

    low, high = 2.25e6, 2.8e6
    while high - low > 1.0e3:
        middle = 0.5 * (low + high)
        if regime(middle) == "VLE":
            low = middle
        else:
            high = middle
    print(f"  the split changes character between {low / 1e6:.4f} and {high / 1e6:.4f} MPa")
    print("  (the incipient phase's density root stops being a vapour there)")
    for pressure_Pa, expected in ((low, "VLE"), (high, "LLE")):
        result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
        print(
            f"  {pressure_Pa / 1e6:.4f} MPa: {expected}, dG/RT ="
            f" {float(result.diagnostics['delta_g_split_rt']):.6e},"
            f" post-split {result.diagnostics['post_split_status']}"
        )
        record_check(
            f"{pressure_Pa / 1e6:.4f} MPa: the returned {expected} answer lowers the Gibbs energy",
            float(result.diagnostics["delta_g_split_rt"]) < 0.0
            and result.diagnostics["post_split_status"] == "stable",
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="also run the 25-point pressure scan and the bisected VL/LL boundary",
    )
    # `parse_known_args`, not `parse_args`: `tests/test_examples.py` runs this
    # script via `runpy.run_path` with pytest's own `sys.argv` still in place.
    args, _unknown = parser.parse_known_args()

    print("=" * 78)
    print("Polymer/solvent vapour-liquid equilibrium: polyethylene / n-pentane at 453 K")
    print("=" * 78)
    print(
        "Polymer parameters: as tabulated by Martini et al. (2009) citing Gross &\n"
        "Sadowski (2002); the primary table is paywalled and was NOT read. They are a\n"
        "cited test fixture, not packaged runtime data, and nothing below is compared\n"
        "against measurement."
    )

    results = the_splits()
    the_one_dimensional_solve(results)
    if si is None:
        print("\n  feos is not installed; route 4 is skipped.")
        print("  Install with: pip install -e '.[validation]'")
    else:
        the_chemical_potentials()
    the_ternary()
    if args.full:
        the_pressure_scan()
        the_boundary()
    else:
        print("\n  (pass --full for the 25-point pressure scan and the bisected VL/LL boundary)")

    print("\n" + "=" * 78)
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
