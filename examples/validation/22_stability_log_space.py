"""Tangent-plane stability in log mole numbers, with PASS/FAIL on every check.

What is being validated
-----------------------
ADR-0025 removes a clamp. Michelsen's stability iteration carries the
*unnormalized* mole numbers ``W``, and the only place they were ever needed as
doubles is the normalization ``w = W / sum_j W_j``; before this slice that
normalization was guarded by ``ln W <- clip(ln W, -700, 700)``, which is where
``exp`` stops existing. A stationary point outside that window was therefore
not merely inaccurate - it was unreachable. The iteration parked on the
boundary and the trial was reported ``second_order_no_progress``.

A polymer melt is such a point. A polyethylene of ``Mw = 53 000``
(``m = 1393.9`` segments) has ``ln phi_polymer = -1550`` as a melt at 453 K and
0.5 MPa, and the feed - 5 wt% of it in n-pentane, a vapour at that pressure -
has ``d_polymer = -97.5``, so the stationary point sits at
``ln W_polymer = 1452``. Validation Case P-14 pinned this as a **stability
miss**: the deepest stationary point the trial set could reach was a shallow
vapour-side one (``tpd`` of order 1e-04), the log-space split stage of ADR-0024
was seeded 1300 orders of magnitude from the answer, and ``flash_tp`` raised.

The normalization is now done in logs where, and only where, the old clamp
would have engaged, and the melt is found in three successive substitutions.
This script checks the result five ways, deliberately different in kind:

1. **The stationary point itself**, against its own defining equations (5) and
   (7) re-derived here from ``PCSAFTEOS`` - including that ``sum_i W_i`` has
   overflowed to ``inf`` while ``ln sum_i W_i`` has not, which is the whole
   reason the decision record exists.
2. **The split it seeds**: ``flash_tp`` returns a solvent vapour against a
   solvent-swollen melt, with both phases named from a measured
   compressibility (ADR-0017), mass balance, equal fugacity, the Gibbs change
   against the single-phase feed, and the post-split stability of both phases.
3. **A one-dimensional equal-fugacity solve written in this script**, with the
   vapour taken to be *exactly* pure solvent: one unknown, ``ln x_solvent`` in
   the melt, by a finite-difference Newton. It uses no part of the flash - no
   stability test, no Rachford-Rice, no Newton stage, no phase-count logic -
   only ``PCSAFTEOS``'s ``(T, rho, x)`` interface and ``density_roots``.
4. **FeOs's chemical potentials at chemthermo's converged phases and
   densities** (feos-org/feos, MIT OR Apache-2.0; the same Gross & Sadowski
   model in Rust with every derivative by automatic differentiation).
5. **Dormancy.** The gate is the claim that nothing else moved, so it is
   measured: the shorter ``Mw = 16 400`` chain of Case P-14, a Peng-Robinson
   state and an activity-model state all run the pre-ADR-0025 arithmetic, and
   ``--full`` adds the whole 144-state Peng-Robinson stability grid.

``k_ij`` cannot be given to FeOs's PC-SAFT from Python in feos 0.10.1
(``EquationOfState.pcsaft`` raises "missing field ``k_ij``" for every
serialization tried), so route 4 runs at ``k_ij = 0`` on **both** sides and the
flash it checks is re-run at ``k_ij = 0`` for the purpose. Routes 1, 2, 3 and 5
use the fitted ``k_ij = -0.006``. As in Cases P-6 to P-12 the 42 universal
constants of the 2001 dispersion term are not shared (chemthermo packages the
ten figures as printed, FeOs hard-codes fourteen), so route 4 is reported
twice, as shipped and with FeOs's table substituted in.

**The polymer parameters are not verified data.** They are as tabulated by
Martini, Cismondi, Barbosa & Brignole, *Sep. Sci. Technol.* **44** (2009),
citing Gross & Sadowski, *IECR* **41** (2002) 1084 - a paywalled table that was
not read. They live in ``tests/fixtures/pcsaft/martini2009_polymers.json``, are
never packaged, and nothing here is compared against measurement.

``--full`` adds the 144-state Peng-Robinson stability grid (route 5), a
pressure scan from 0.4 to 2 MPa - the whole region where the feed is a vapour,
every point of which either raised or was seeded from the wrong stationary
point before this slice - and **route 6, the six states ADR-0026 retires**
(validation Case P-16): 0.3 MPa and 2.8 to 3.2 MPa, the last states of the
0.3-3.6 MPa sweep that ``flash_tp`` refused. Each is verified by a solve that
shares no code with the flash - the 1-D equal-fugacity solve of route 3 for the
vapour-liquid state, a two-equation Newton for the five liquid-liquid ones -
and the 0.3 MPa iteration is printed residual by residual, with and without the
curvature safeguard, so the repair is visible rather than asserted. The default
takes about five seconds; ``--full`` about two minutes.

Requires the optional ``feos`` dependency for route 4 only::

    pip install -e ".[validation]"

Routes 1, 2, 3 and 5 run without it; the script says so and still exits 0.
"""

from __future__ import annotations

import argparse
import dataclasses
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
from chemthermo.flash import _detect
from chemthermo.flash._log_space import log_space_seed, log_space_split
from chemthermo.models.base import KAPPA_LIQUID_THRESHOLD
from chemthermo.parameters import PCSAFTParameters, PCSAFTRecord

REPO = Path(__file__).resolve().parents[2]
FIXTURE = REPO / "tests" / "fixtures" / "pcsaft" / "martini2009_polymers.json"

TEMPERATURE_K = 453.0
#: The long chain: `m = 1393.9`, the one whose melt is outside `exp`'s range.
PE_MW_G_MOL = 53000.0
#: The short chain of ADR-0024, which is inside it - the dormancy control.
PE_SHORT_MW_G_MOL = 16400.0
PENTANE_MW_G_MOL = 72.146
KIJ = -0.006
#: n-pentane, Gross & Sadowski (2001) Table 1.
PENTANE = (2.6896, 3.7729, 231.20)

#: The two states validation Case P-14 pinned as a stability miss.
PRESSURES_PA = (5.0e5, 1.0e6)
#: `--full`: the same question across the vapour-liquid region. It stops at
#: 2 MPa because that is where the *feed* stops being a vapour: at 2.2 MPa and
#: above the feed's own lowest-Gibbs branch is the liquid root, there is no
#: melt stationary point distinct from it, and `stability_tp` returns the
#: shallow near-critical point it returned before this slice (measured:
#: `tpd_min = -8.4e-02` at 2.2 MPa, unchanged).
SCAN_PA = (4.0e5, 5.0e5, 7.5e5, 1.0e6, 1.5e6, 2.0e6)
#: The first two pressures past that edge, printed rather than asserted on.
SCAN_BEYOND_PA = (2.2e6, 2.5e6)

#: Where `exp` stops existing, and therefore where the old clamp sat.
LN_W_WINDOW = 700.0

#: `--full` route 6 (ADR-0026, validation Case P-16): the six states of the
#: 0.3-3.6 MPa sweep that raised `ConvergenceError` at HEAD 584c508. The first
#: is vapour-liquid, the other five liquid-liquid.
P16_VLE_PA = 3.0e5
P16_LLE_PA = (2.8e6, 2.9e6, 3.0e6, 3.1e6, 3.2e6)
#: The sweep those six came out of: both molar masses, 0.3-3.6 MPa at 0.1 MPa.
SWEEP_PA = tuple(3.0e5 + 1.0e5 * step for step in range(34))
#: Agreement between the flash and the independent two-equation Newton, on the
#: polymer-rich composition. It is a relative tolerance because the quantity is
#: 9e-04 and both solves are converged to ~1e-13 on their own residuals.
TIE_LINE_REL_TOL = 1e-9

NEWTON_TOL = 1e-12
#: The melt's solvent mole fraction from the flash and from the 1-D solve.
COMPOSITION_TOL = 1e-10
#: Equal fugacity in logarithms, including the polymer's.
LOG_FUGACITY_TOL = 1e-8
#: Equations (5) and (7) at the stationary point.
STATIONARITY_TOL = 1e-8
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


def polymer_record(mw_g_mol: float) -> PCSAFTRecord:
    return PCSAFTRecord(
        name="Polyethylene",
        segments_per_g=ROW["segments_per_g"],
        MW_g_mol=mw_g_mol,
        sigma_A=ROW["sigma_A"],
        epsilon_k_K=ROW["epsilon_k_K"],
    )


def binary_parameters(mw_g_mol: float = PE_MW_G_MOL) -> PCSAFTParameters:
    return PCSAFTParameters.from_records(
        [
            polymer_record(mw_g_mol),
            PCSAFTRecord(
                name="n-Pentane",
                m=PENTANE[0],
                sigma_A=PENTANE[1],
                epsilon_k_K=PENTANE[2],
                MW_g_mol=PENTANE_MW_G_MOL,
            ),
        ]
    )


def binary_mixture(mw_g_mol: float = PE_MW_G_MOL, weight_fraction: float = 0.05) -> ct.Mixture:
    polymer = weight_fraction / mw_g_mol
    solvent = (1.0 - weight_fraction) / PENTANE_MW_G_MOL
    total = polymer + solvent
    return ct.Mixture.from_components(
        [
            ct.Component.custom(
                "Polyethylene",
                mw_kg_per_mol=mw_g_mol / 1000.0,
                formula="(C2H4)n",
                volatile=False,
                source="see tests/fixtures/pcsaft/martini2009_polymers.json",
            ),
            ct.Component.from_database("n-Pentane"),
        ],
        [polymer / total, solvent / total],
        normalize=True,
    )


def binary_eos(kij: float = KIJ, mw_g_mol: float = PE_MW_G_MOL) -> PCSAFTEOS:
    return PCSAFTEOS(parameters=binary_parameters(mw_g_mol), kij=kij)


def kappa(pressure_Pa: float, x: Sequence[float]) -> list[float]:
    """``P / (rho dP/drho)`` from the public pressure routine (ADR-0017)."""
    bound = PCSAFTEOS(
        components=("Polyethylene", "n-Pentane"), parameters=binary_parameters(), kij=KIJ
    )
    values = []
    for density in bound.density_roots(
        temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=list(x)
    ):
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
    return values


# ---------------------------------------------------------------------------
# Route 1: the stationary point
# ---------------------------------------------------------------------------


def the_stationary_points() -> dict[float, ct.StabilityResult]:
    section("1. stability_tp on polyethylene(53000) / n-pentane, 5 wt%, 453 K, k_ij = -0.006")
    mixture = binary_mixture()
    eos = binary_eos()
    z = np.asarray(mixture.fractions, dtype=float)
    print(
        f"  feed z = {tuple(mixture.fractions)}  (m_polymer = {ROW['segments_per_g'] * PE_MW_G_MOL:.1f})"
    )
    print(
        "\n  Before ADR-0025 the three liquid-surface trials ended `second_order_no_progress`\n"
        "  after 51 iterations, parked on the clamp at ln W = 700, and the verdict came from\n"
        "  a shallow vapour-side point with tpd ~ -1e-04."
    )
    print(
        f"\n  {'P/MPa':>7} {'tpd_min':>15} {'ln W_PE':>12} {'ln W_C5':>10}"
        f" {'iters':>6} {'minimizing trial':>20}"
    )

    results: dict[float, ct.StabilityResult] = {}
    for pressure_Pa in PRESSURES_PA:
        result = ct.stability_tp(
            mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
        )
        results[pressure_Pa] = result
        assert result.trial_ln_W is not None
        best = next(
            trial
            for trial in result.trials
            if trial.label == str(result.diagnostics["minimizing_trial"])
        )
        print(
            f"  {pressure_Pa / 1e6:7.2f} {result.tpd_min:15.6f}"
            f" {result.trial_ln_W[0]:12.4f} {result.trial_ln_W[1]:10.4f}"
            f" {best.iterations:6d} {str(result.diagnostics['minimizing_trial']):>20}"
        )

    print()
    for pressure_Pa, result in results.items():
        label = f"{pressure_Pa / 1e6:.1f} MPa"
        assert result.trial_ln_W is not None
        ln_capital_w = np.asarray(result.trial_ln_W, dtype=float)
        w = np.asarray(result.trial_composition, dtype=float)

        record_check(f"{label}: the feed is unstable", result.status == "unstable")
        record_check(f"{label}: the melt is the minimizer", result.phase_branch == "liquid")
        record_check(
            f"{label}: the feed is a vapour",
            result.feed_branch == "vapor",
        )
        record_check(
            f"{label}: ln W_polymer is outside the old clamp ({ln_capital_w[0]:.1f} > {LN_W_WINDOW:.0f})",
            ln_capital_w[0] > LN_W_WINDOW,
        )
        record_check(
            f"{label}: the trial says so ({result.diagnostics.get('log_space_trial_count')} of "
            f"{result.diagnostics['trial_count']} trials in log space)",
            bool(result.diagnostics.get("minimizing_trial_log_space", False)),
        )

        # Equation (5), re-derived here: `ln W_i + ln phi_i(w) - d_i = 0` on the
        # surface the trial iterated on. Nothing of `stability_tp`'s own
        # residual is reused.
        def ln_phi(composition: np.ndarray, phase: str) -> np.ndarray:
            return np.asarray(
                eos.log_fugacity_coefficients(
                    mixture=mixture,
                    temperature_K=TEMPERATURE_K,
                    pressure_Pa=pressure_Pa,
                    composition=composition.tolist(),
                    phase=phase,
                )
            )

        d = np.log(z) + ln_phi(z, "vapor")
        stationarity = float(np.max(np.abs(ln_capital_w + ln_phi(w, "liquid") - d)))
        record_check(
            f"{label}: equation (5) holds at the stationary point ({stationarity:.1e})",
            stationarity < STATIONARITY_TOL,
        )

        # Equation (7): `tpd = -ln sum_i W_i`, with the sum taken in logs.
        largest = float(np.max(ln_capital_w))
        ln_sum = largest + math.log(float(np.sum(np.exp(ln_capital_w - largest))))
        record_check(
            f"{label}: equation (7), tpd = -ln sum_W ({abs(-ln_sum - result.tpd_min):.1e})",
            abs(-ln_sum - result.tpd_min) < STATIONARITY_TOL,
        )
        record_check(
            f"{label}: sum_W itself has overflowed, and is reported as inf",
            result.diagnostics["sum_W"] == math.inf,
        )
        record_check(
            f"{label}: ln sum_W has not ({float(result.diagnostics['ln_sum_W']):.4f})",
            math.isfinite(float(result.diagnostics["ln_sum_W"])),
        )
        # The normalized composition cannot carry this stationary point: the
        # solvent's share of it is exp(-1449).
        record_check(
            f"{label}: the normalized w rounds to (1.0, 0.0), and ln W does not",
            w.tolist() == [1.0, 0.0] and math.isfinite(ln_capital_w[1]),
        )
    return results


# ---------------------------------------------------------------------------
# Route 2: the split it seeds
# ---------------------------------------------------------------------------


def the_splits() -> dict[float, ct.FlashResult]:
    section("2. flash_tp: the split the melt seeds (Case P-14's ConvergenceError, retired)")
    mixture = binary_mixture()
    eos = binary_eos()
    print(
        f"\n  {'P/MPa':>7} {'phases':>15} {'beta_vapor':>13} {'x_C5 (melt)':>14}"
        f" {'wt% C5':>8} {'ln y_PE':>10} {'resid':>9} {'dG/RT':>11}"
    )
    results: dict[float, ct.FlashResult] = {}
    for pressure_Pa in PRESSURES_PA:
        result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
        results[pressure_Pa] = result
        melt = result.phases["liquid"].composition.fractions
        solvent_mass = melt[1] * PENTANE_MW_G_MOL
        polymer_mass = melt[0] * PE_MW_G_MOL
        print(
            f"  {pressure_Pa / 1e6:7.2f} {'+'.join(sorted(result.phases)):>15}"
            f" {result.vapor_fraction if result.vapor_fraction is not None else float('nan'):13.10f}"
            f" {melt[1]:14.10f}"
            f" {100.0 * solvent_mass / (solvent_mass + polymer_mass):8.3f}"
            f" {float(result.diagnostics['log_space_ln_x_min']):10.2f}"
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
            f"{label}: seeded from the tangent-plane melt",
            result.diagnostics["k_seed"] == "stability-log",
        )
        record_check(
            f"{label}: solved in log mole numbers",
            result.diagnostics["converged_stage"] == "second-order-log",
        )
        record_check(
            f"{label}: named from a measured compressibility",
            result.diagnostics["phase_label_method"] == "compressibility",
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
        vapor_kappa = kappa(pressure_Pa, result.phases["vapor"].composition.fractions)
        melt_kappa = kappa(pressure_Pa, result.phases["liquid"].composition.fractions)
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
# Route 3: an independent one-dimensional solve
# ---------------------------------------------------------------------------


def the_one_dimensional_solve(results: dict[float, ct.FlashResult]) -> None:
    section("3. An independent 1-D equal-fugacity solve, vapour taken as exactly pure solvent")
    print(
        "  One unknown, ln x_solvent in the melt, against\n"
        "      ln phi_s^V(pure n-pentane vapour) = ln x_s + ln phi_s^L(x)\n"
        "  by a finite-difference Newton. Nothing of flash_tp or stability_tp is used."
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
            f"{pressure_Pa / 1e6:.1f} MPa: melt x_solvent matches the 1-D solve "
            f"< {COMPOSITION_TOL:g}",
            difference < COMPOSITION_TOL,
        )
        record_check(
            f"{pressure_Pa / 1e6:.1f} MPa: the 1-D solve is converged",
            abs(residual(ln_x)) < LOG_FUGACITY_TOL,
        )

    section("3b. The polymer's own equal-fugacity condition, in logarithms")
    print(
        "  ln x_PE + ln phi_PE^L(x)  =  ln y_PE + ln phi_PE^V(y), with ln y_PE ~ -1500.\n"
        "  This is the equation neither the linear split nor the clamped stability\n"
        "  iteration could write down at all."
    )
    print(f"\n  {'P/MPa':>7} {'ln f_PE (melt)':>17} {'ln f_PE (vapour)':>18} {'|difference|':>13}")
    for pressure_Pa, result in results.items():
        x = list(result.phases["liquid"].composition.fractions)
        y = list(result.phases["vapor"].composition.fractions)
        melt = math.log(x[0]) + float(ln_phi(pressure_Pa, x, vapor=False)[0])
        # ADR-0024 decision 3: `y_polymer` is an exact 0.0 and its logarithm is
        # in the diagnostics.
        vapor = float(result.diagnostics["log_space_ln_x_min"]) + float(
            ln_phi(pressure_Pa, y, vapor=True)[0]
        )
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
        "  potentials. The vapour's polymer mole fraction is an exact 0.0 as a double, so\n"
        "  the ideal part of its potential is -inf on both sides and drops out; what is\n"
        "  compared is the solvent's potential, which is the equation that sets the answer."
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
                difference = np.abs(potentials[0] - potentials[1])
                worst = max(worst, float(np.max(difference[np.isfinite(difference)])))
        finally:
            pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = shipped
        print(f"  worst |d mu_i / RT| over the two states, {label}: {worst:.3e}")
        if label == "matched constants":
            record_check(
                f"FeOs potentials equal at both phases < {POTENTIAL_TOL:g}", worst < POTENTIAL_TOL
            )


# ---------------------------------------------------------------------------
# Route 5: dormancy
# ---------------------------------------------------------------------------


def the_dormancy(full: bool) -> None:
    section("5. Dormancy: the log-space route runs only where the old clamp engaged")
    print(
        "  The gate is the whole bit-identity argument, so it is measured rather than\n"
        "  asserted: a trial records whether it had to normalize in logs, and nothing\n"
        "  that had an answer before ADR-0025 does."
    )

    short = binary_mixture(PE_SHORT_MW_G_MOL)
    short_eos = binary_eos(KIJ, PE_SHORT_MW_G_MOL)
    engaged = 0
    trials = 0
    worst_ln_capital_w = 0.0
    for pressure_Pa in (5.0e5, 1.0e6, 2.0e6):
        result = ct.stability_tp(
            short, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=short_eos
        )
        for trial in result.trials:
            trials += 1
            engaged += int(trial.log_space)
            if trial.ln_W is not None:
                finite = [abs(value) for value in trial.ln_W if math.isfinite(value)]
                worst_ln_capital_w = max(worst_ln_capital_w, max(finite))
    print(
        f"\n  Mw = {PE_SHORT_MW_G_MOL:.0f} (validation Case P-14), 3 states:"
        f" {engaged} of {trials} trials in log space, worst |ln W| = {worst_ln_capital_w:.1f}"
    )
    record_check(
        f"the Mw = {PE_SHORT_MW_G_MOL:.0f} chain never leaves the old arithmetic", engaged == 0
    )
    record_check(
        f"and never had to: worst |ln W| = {worst_ln_capital_w:.1f} < {LN_W_WINDOW:.0f}",
        worst_ln_capital_w < LN_W_WINDOW,
    )

    states = [
        (
            "Peng-Robinson methane/ethane, 240 K, 3 MPa",
            ct.stability_tp(
                ct.Mixture.from_database(["Methane", "Ethane"], [0.5, 0.5], normalize=True),
                temperature_K=240.0,
                pressure_Pa=3.0e6,
                eos=ct.PengRobinsonEOS(),
            ),
        ),
        (
            "PC-SAFT water/n-hexane, 298.15 K, 1 atm",
            ct.stability_tp(
                ct.Mixture.from_database(["Water", "n-Hexane"], [0.5, 0.5], normalize=True),
                temperature_K=298.15,
                pressure_Pa=101325.0,
                eos=PCSAFTEOS(),
            ),
        ),
        (
            "NRTL n-butanol/water, 298.15 K (liquid-liquid)",
            ct.stability_tp(
                ct.Mixture.from_database(["n-Butanol", "Water"], [0.5, 0.5], normalize=True),
                temperature_K=298.15,
                pressure_Pa=101325.0,
                activity_model=ct.NRTL(
                    parameters=ct.NRTLParameters.from_pairs(
                        [("n-Butanol", "Water", 0.90047, 3.51307, 0.48, 0.48)]
                    )
                ),
            ),
        ),
    ]
    for label, result in states:
        record_check(
            f"{label}: no trial in log space",
            not any(trial.log_space for trial in result.trials)
            and "log_space_trial_count" not in result.diagnostics,
        )

    if not full:
        return

    path = REPO / "tests" / "test_stability_eos_surfaces.py"
    spec = importlib.util.spec_from_file_location("_stability_grid", path)
    assert spec is not None and spec.loader is not None
    grid = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(grid)

    eos = ct.PengRobinsonEOS()
    grid_states = 0
    grid_trials = 0
    grid_engaged = 0
    for names, z in grid.GRID_MIXTURES:
        for temperature_K in grid.GRID_T_K:
            for pressure_Pa in grid.GRID_P_PA:
                result = ct.stability_tp(
                    ct.Mixture.from_database(list(names), list(z), normalize=True),
                    temperature_K=temperature_K,
                    pressure_Pa=pressure_Pa,
                    eos=eos,
                )
                grid_states += 1
                grid_trials += len(result.trials)
                grid_engaged += sum(1 for trial in result.trials if trial.log_space)
    print(
        f"\n  (--full) 144-state Peng-Robinson stability grid:"
        f" {grid_engaged} of {grid_trials} trials over {grid_states} states in log space"
    )
    record_check(
        "the 144-state Peng-Robinson grid never leaves the old arithmetic", grid_engaged == 0
    )


def the_pressure_scan() -> None:
    section("6. (--full) The melt stationary point across the vapour-liquid region")
    print(
        "  Before ADR-0025 every one of these states raised or was seeded from a shallow\n"
        "  vapour-side point: 0.4 / 0.5 / 0.75 / 1.0 MPa raised `ConvergenceError` in the\n"
        "  log-space split stage, and 1.5 / 2.0 MPa converged from that shallow point to\n"
        "  the same answer this finds from the melt (they agree to 4e-13 relative)."
    )
    mixture = binary_mixture()
    eos = binary_eos()
    print(
        f"\n  {'P/MPa':>7} {'feed':>8} {'tpd_min':>14} {'ln W_PE':>11}"
        f" {'phases':>15} {'x_C5 (melt)':>14} {'beta_vapor':>13}"
    )
    found = 0
    for pressure_Pa in SCAN_PA:
        stability = ct.stability_tp(
            mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
        )
        result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
        assert stability.trial_ln_W is not None
        melt = "liquid" if "liquid" in result.phases else sorted(result.phases)[0]
        found += int(stability.tpd_min < -100.0)
        print(
            f"  {pressure_Pa / 1e6:7.2f} {str(stability.feed_branch):>8}"
            f" {stability.tpd_min:14.4f}"
            f" {stability.trial_ln_W[0]:11.2f} {'+'.join(sorted(result.phases)):>15}"
            f" {result.phases[melt].composition.fractions[1]:14.10f}"
            f" {result.vapor_fraction if result.vapor_fraction is not None else float('nan'):13.10f}"
        )
    record_check(
        f"the melt is the minimizer at all {len(SCAN_PA)} pressures", found == len(SCAN_PA)
    )

    print(
        "\n  Past the edge of that region the feed is no longer a vapour, there is no melt\n"
        "  stationary point distinct from it, and the answer is the one this slice did not\n"
        "  touch (printed, not asserted on):"
    )
    for pressure_Pa in SCAN_BEYOND_PA:
        stability = ct.stability_tp(
            mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
        )
        assert stability.trial_ln_W is not None
        print(
            f"  {pressure_Pa / 1e6:7.2f} {str(stability.feed_branch):>8}"
            f" {stability.tpd_min:14.4f} {stability.trial_ln_W[0]:11.2f}"
        )


# ---------------------------------------------------------------------------
# Route 6 (--full): the six refusals ADR-0026 retires (validation Case P-16)
# ---------------------------------------------------------------------------


def _ln_phi(pressure_Pa: float, x: Sequence[float], *, dense: bool) -> np.ndarray:
    """``ln phi`` on one density root, from the (T, rho, x) interface alone."""
    bound = PCSAFTEOS(
        components=("Polyethylene", "n-Pentane"), parameters=binary_parameters(), kij=KIJ
    )
    roots = bound.density_roots(
        temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=list(x)
    )
    density = roots[-1] if dense else roots[0]
    return np.asarray(
        bound.ln_fugacity_coefficients(
            temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=list(x)
        )
    )


def _independent_tie_line(pressure_Pa: float, start: tuple[float, float]) -> tuple[float, float]:
    """A liquid-liquid tie line from a two-equation Newton written here.

    The unknowns are ``ln x_polymer`` in each of the two liquids, and the
    equations are the two equal-fugacity conditions written in logarithms:

        ln x_P^I  + ln phi_P^I  =  ln x_P^II + ln phi_P^II
        ln(1-x_P^I) + ln phi_S^I = ln(1-x_P^II) + ln phi_S^II

    The feed never enters, so this knows nothing about Rachford-Rice, the phase
    count, the stability test or the split; it is the tie line of the *model*.
    The polymer-lean liquid holds ``exp(-94)`` polymer, which is why the
    unknowns are logarithms.
    """
    values = np.asarray(start, dtype=float)

    def equations(v: np.ndarray) -> np.ndarray:
        first, second = float(v[0]), float(v[1])
        x_i = [math.exp(first), -math.expm1(first)]
        x_ii = [math.exp(second), -math.expm1(second)]
        ln_phi_i = _ln_phi(pressure_Pa, x_i, dense=True)
        ln_phi_ii = _ln_phi(pressure_Pa, x_ii, dense=True)
        return np.array(
            [
                first + ln_phi_i[0] - second - ln_phi_ii[0],
                math.log1p(-math.exp(first))
                + ln_phi_i[1]
                - math.log1p(-math.exp(second))
                - ln_phi_ii[1],
            ]
        )

    for _ in range(80):
        residual = equations(values)
        size = float(np.max(np.abs(residual)))
        if size < NEWTON_TOL:
            break
        jacobian = np.zeros((2, 2))
        for column in range(2):
            step = 1e-6
            plus = values.copy()
            minus = values.copy()
            plus[column] += step
            minus[column] -= step
            jacobian[:, column] = (equations(plus) - equations(minus)) / (2.0 * step)
        direction = np.linalg.solve(jacobian, -residual)
        scale = 1.0
        while scale > 1e-12:
            if float(np.max(np.abs(equations(values + scale * direction)))) < size:
                break
            scale *= 0.5
        values = values + scale * direction
    return float(values[0]), float(values[1])


def the_six_refusals() -> None:
    section("7. (--full) The six states ADR-0026 retires (validation Case P-16)")
    print(
        "  At HEAD 584c508 a 0.3-3.6 MPa sweep at 0.1 MPa over both molar masses - 68\n"
        "  states - had six raising ConvergenceError, all on the Mw = 53 000 chain:\n"
        "    0.3, 2.8, 2.9 MPa   the log-space stage spent its budget next to the\n"
        "                        trivial solution (residual 3.9e-05 / 1.1e-08 / 6.0e+00)\n"
        "    3.0, 3.1, 3.2 MPa   'neither the stability-seeded nor the Wilson K-values\n"
        "                        bracket a Rachford-Rice root'\n"
    )
    mixture = binary_mixture()
    eos = binary_eos()

    # -- the iteration itself, before and after -----------------------------
    print("  The 0.3 MPa log-space iteration, residual by residual:\n")
    z = np.asarray(mixture.fractions, dtype=float)
    stability = ct.stability_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=P16_VLE_PA, eos=eos
    )
    assert stability.trial_composition is not None
    ln_capital_w = _detect._stationary_point_ln_capital_w(stability)
    _k, incipient = _detect._stability_k_seed(
        mixture,
        TEMPERATURE_K,
        P16_VLE_PA,
        z=z,
        w=np.asarray(stability.trial_composition, dtype=float),
        tpd_min=float(stability.tpd_min),
        ln_capital_w=ln_capital_w,
    )
    u0 = log_space_seed(
        z=z,
        w=np.asarray(stability.trial_composition, dtype=float),
        tpd_min=float(stability.tpd_min),
        incipient_vapor=incipient == "vapor",
        ln_capital_w=ln_capital_w,
    )
    roots_x, roots_y = _detect._phi_phi_roots(
        eos,
        mixture,
        TEMPERATURE_K,
        P16_VLE_PA,
        feed_branch=stability.feed_branch,
        incipient_branch=stability.phase_branch,
        incipient_phase=incipient,
        seed_label="stability-log",
    )

    def stage(iterations: int, safeguard: bool) -> float:
        return log_space_split(
            z=z,
            u0=u0,
            terms_i=roots_x.ln_fugacity_terms,
            terms_ii=roots_y.ln_fugacity_terms,
            settings=dataclasses.replace(
                ct.FlashSettings(), second_order_max_iter=max(iterations, 1)
            ),
            curvature_safeguard=safeguard,
        ).residual

    print(f"  {'iteration':>10} {'ADR-0024 rule':>18} {'with the safeguard':>22}")
    for iterations in (1, 2, 4, 6, 8, 10, 11, 12, 13, 20, 50, 100):
        print(
            f"  {iterations:10d} {stage(iterations, False):18.6e} {stage(iterations, True):22.6e}"
        )
    record_check(
        "the ADR-0024 rule is still on a residual of ~4e-05 after 100 iterations",
        3e-05 < stage(100, False) < 5e-05,
    )
    record_check(
        "the safeguarded rule is below 1e-12 within 13",
        stage(13, True) < 1e-12,
    )

    # -- the vapour-liquid state, against the 1-D solve ---------------------
    print("\n  0.3 MPa, vapour-liquid, against the route-3 one-dimensional solve:")
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=P16_VLE_PA, eos=eos)
    target = float(_ln_phi(P16_VLE_PA, [0.0, 1.0], dense=False)[1])

    def residual(ln_x: float) -> float:
        solvent = math.exp(ln_x)
        return ln_x + float(_ln_phi(P16_VLE_PA, [1.0 - solvent, solvent], dense=True)[1]) - target

    ln_x = math.log(float(result.phases["liquid"].composition.fractions[1])) - 1e-3
    for _ in range(80):
        value = residual(ln_x)
        step = 1e-7
        slope = (residual(ln_x + step) - residual(ln_x - step)) / (2.0 * step)
        correction = -value / slope
        ln_x += correction
        if abs(correction) < NEWTON_TOL:
            break
    ours = float(result.phases["liquid"].composition.fractions[1])
    theirs = math.exp(ln_x)
    print(
        f"    melt x_C5: flash {ours:.14f}   1-D {theirs:.14f}   |difference| {abs(ours - theirs):.2e}"
    )
    print(
        f"    beta_vapour {result.vapor_fraction:.12f}   ln y_PE"
        f" {float(result.diagnostics['log_space_ln_x_min']):.4f}"
        f"   dG/RT {float(result.diagnostics['delta_g_split_rt']):.4e}"
    )
    record_check(
        f"0.3 MPa: the melt matches the 1-D solve < {COMPOSITION_TOL:g}",
        abs(ours - theirs) < COMPOSITION_TOL,
    )
    record_check(
        "0.3 MPa: the curvature safeguard is what produced it",
        result.diagnostics["log_space_curvature_safeguard"] is True,
    )
    for label, ok in (
        ("a vapour and a liquid", sorted(result.phases) == ["liquid", "vapor"]),
        (
            f"mass balance < {MASS_BALANCE_TOL:g}",
            float(result.diagnostics["mass_balance_residual"]) < MASS_BALANCE_TOL,
        ),
        (
            f"equal fugacity < {LOG_FUGACITY_TOL:g}",
            float(result.diagnostics["log_space_residual"]) < LOG_FUGACITY_TOL,
        ),
        ("dG_split/RT < 0", float(result.diagnostics["delta_g_split_rt"]) < 0.0),
        ("post-split stable", result.diagnostics["post_split_status"] == "stable"),
    ):
        record_check(f"0.3 MPa: {label}", ok)

    # -- the five liquid-liquid states, against a 2-equation Newton ---------
    print("\n  2.8-3.2 MPa, liquid-liquid, against a two-equation Newton written here:")
    print(
        f"\n  {'P/MPa':>7} {'x_PE (rich, flash)':>21} {'x_PE (rich, Newton)':>21}"
        f" {'rel.':>9} {'ln x_PE (lean)':>15} {'beta_rich':>12}"
    )
    for pressure_Pa in P16_LLE_PA:
        result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
        rich = max(result.phases, key=lambda n: result.phases[n].composition.fractions[0])
        lean = min(result.phases, key=lambda n: result.phases[n].composition.fractions[0])
        ours_rich = float(result.phases[rich].composition.fractions[0])
        ours_lean = float(result.phases[lean].composition.fractions[0])
        # Started one percent away in ln x from the flash's answer, so the
        # Newton has real work to do and cannot be said to have been handed it.
        first, second = _independent_tie_line(
            pressure_Pa, (math.log(ours_rich) * 1.01, math.log(ours_lean) * 1.01)
        )
        theirs_rich = math.exp(first)
        relative = abs(theirs_rich - ours_rich) / ours_rich
        print(
            f"  {pressure_Pa / 1e6:7.1f} {ours_rich:21.15e} {theirs_rich:21.15e}"
            f" {relative:9.1e} {math.log(ours_lean):15.6f}"
            f" {result.phase_fractions[rich]:12.9f}"
        )
        record_check(
            f"{pressure_Pa / 1e6:.1f} MPa: the tie line matches the 2-equation Newton"
            f" < {TIE_LINE_REL_TOL:g} relative",
            relative < TIE_LINE_REL_TOL,
        )
        record_check(
            f"{pressure_Pa / 1e6:.1f} MPa: two liquids, verified",
            sorted(result.phases) == ["liquid1", "liquid2"]
            and float(result.diagnostics["mass_balance_residual"]) < MASS_BALANCE_TOL
            and float(result.diagnostics["fugacity_residual"]) < LOG_FUGACITY_TOL
            and float(result.diagnostics["delta_g_split_rt"]) < 0.0
            and result.diagnostics["post_split_status"] == "stable",
        )


def the_sweep() -> None:
    section("8. (--full) The whole 0.3-3.6 MPa sweep, both molar masses, 68 states")
    print(
        "  The claim ADR-0026 makes is not about six states but about the sweep they\n"
        "  came out of: no ConvergenceError anywhere in it, and a verdict sequence that\n"
        "  changes character once and only once for each chain."
    )
    print(
        f"\n  {'Mw':>7} {'states':>7} {'refused':>8} {'VLE':>5} {'LLE':>5}"
        f" {'single':>7} {'VLE -> LLE between':>22} {'worst resid':>12}"
    )
    failures_here = 0
    for mw_g_mol in (PE_SHORT_MW_G_MOL, PE_MW_G_MOL):
        mixture = binary_mixture(mw_g_mol)
        eos = binary_eos(KIJ, mw_g_mol)
        verdicts: list[tuple[float, str]] = []
        refused = 0
        worst = 0.0
        for pressure_Pa in SWEEP_PA:
            try:
                result = ct.flash_tp(
                    mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos
                )
            except ct.ConvergenceError:
                refused += 1
                verdicts.append((pressure_Pa, "REFUSED"))
                continue
            verdicts.append((pressure_Pa, str(result.diagnostics["phase_regime"])))
            if len(result.phases) > 1:
                worst = max(worst, float(result.diagnostics["fugacity_residual"]))
        counts = {name: sum(1 for _p, v in verdicts if v == name) for name in ("VLE", "LLE")}
        singles = sum(1 for _p, v in verdicts if v == "single-phase")
        changes = [
            (verdicts[index - 1][0], verdicts[index][0])
            for index in range(1, len(verdicts))
            if verdicts[index][1] != verdicts[index - 1][1]
        ]
        boundary = f"{changes[0][0] / 1e6:.1f}-{changes[0][1] / 1e6:.1f} MPa" if changes else "-"
        print(
            f"  {mw_g_mol:7.0f} {len(verdicts):7d} {refused:8d} {counts['VLE']:5d}"
            f" {counts['LLE']:5d} {singles:7d} {boundary:>22} {worst:12.2e}"
        )
        failures_here += refused
        record_check(f"Mw = {mw_g_mol:.0f}: no ConvergenceError in 34 states", refused == 0)
        record_check(
            f"Mw = {mw_g_mol:.0f}: the verdict changes character exactly once",
            len(changes) == 1,
        )
        record_check(
            f"Mw = {mw_g_mol:.0f}: every two-phase answer has residual < {LOG_FUGACITY_TOL:g}",
            worst < LOG_FUGACITY_TOL,
        )


# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help=(
            "also run the 144-state stability grid, the 0.4-2 MPa scan, the six states "
            "ADR-0026 retires and the whole 0.3-3.6 MPa sweep"
        ),
    )
    # `parse_known_args`, not `parse_args`: `tests/test_examples.py` runs this
    # script through `runpy`, so `sys.argv` is pytest's.
    args, _unknown = parser.parse_known_args()

    print("=" * 78)
    print("Tangent-plane stability in log mole numbers (ADR-0025, validation Case P-15)")
    print("=" * 78)

    the_stationary_points()
    results = the_splits()
    the_one_dimensional_solve(results)

    if si is None:
        section("4. FeOs chemical potentials")
        print("  feos is not installed; skipping route 4.")
        print("  Install it with:  pip install -e '.[validation]'")
    else:
        the_chemical_potentials()

    the_dormancy(args.full)
    if args.full:
        the_pressure_scan()
        the_six_refusals()
        the_sweep()

    section("Summary")
    if failures:
        print(f"  {len(failures)} FAILED:")
        for failure in failures:
            print(f"    - {failure}")
        raise SystemExit(1)
    print("  All checks passed.")
    if not args.full:
        print(
            "\n  (pass --full for the 144-state stability grid, the 0.4-2 MPa scan, the six"
            "\n   states ADR-0026 retires and the whole 0.3-3.6 MPa sweep)"
        )


if __name__ == "__main__":
    main()
