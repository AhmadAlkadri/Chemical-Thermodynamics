"""The fixed benchmark workload (ADR-0023).

Internal to :mod:`chemthermo.bench` (ADR-0001). The list is deliberately
*fixed*: a benchmark record is a comparison instrument, and a workload that
drifts between runs measures nothing. Adding a case is a normal change; editing
or removing one invalidates every committed baseline and has to be recorded in
``benchmarks/README.md``.

Each case names a reference path through the library:

======================================  =================================================
id                                      what it exercises
======================================  =================================================
``pr-flash-ternary``                    the phi-phi flash of a Peng-Robinson ternary
``pr-stability-ternary``                the tangent-plane stability test alone, same state
``pr-flash-grid-24``                    24 states of the ADR-0017 Peng-Robinson grid
``nrtl-lle-tessier-p1``                 the gamma-gamma liquid-liquid split
``modified-raoult-vle``                 the activity + ideal-gas vapour-liquid split
``vlle-364k``                           the ADR-0011 three-phase search
``pcsaft-vle-methane-hexane``           a PC-SAFT vapour-liquid flash
``pcsaft-lle-water-hexane``             the ADR-0019 associating liquid-liquid flash
``pcsaft-polymer-lle``                  the ADR-0022 polymer / solvent split
======================================  =================================================

Every model is built with **analytic** derivatives - no case anywhere in
chemthermo differentiates by finite difference - which is why
``derivative_mode`` is a constant in the record rather than a per-case choice.

``pcsaft-polymer-lle`` reads its parameters from
``tests/fixtures/pcsaft/martini2009_polymers.json``, because chemthermo
packages no polymer parameters (ADR-0022). Outside a source checkout that file
is absent and the case records ``status = "skipped"`` rather than failing the
run.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import chemthermo as ct

from ..eos import PCSAFTEOS
from ..models import PengRobinsonEOS
from ..parameters import PCSAFTParameters, PCSAFTRecord
from ._record import StateOutcome

#: Every model in this workload supplies exact analytic derivatives.
DERIVATIVE_MODE = "analytic"

_POLYMER_FIXTURE = (
    Path(__file__).resolve().parents[3]
    / "tests"
    / "fixtures"
    / "pcsaft"
    / "martini2009_polymers.json"
)
_PENTANE_MW_G_MOL = 72.146
_POLYMER_MW_G_MOL = 16400.0
_POLYMER_KIJ = -0.006


class CaseSkipped(Exception):
    """Raised by a case's ``prepare`` when its inputs are not available here."""


@dataclass(frozen=True)
class BenchCase:
    """One workload: how to build it once, and how to run it many times.

    Attributes:
        id: Stable identifier; the key every comparison is made on.
        description: One line of prose for the record and the README.
        model: The thermodynamic model the case runs on.
        route: The public entry point being measured.
        components: Component names, in the mixture's order.
        prepare: Builds mixtures, models and settings **once**, outside the
            timed region. Databank lookups and parameter resolution are not
            what this harness is measuring.
        invoke: Runs the workload against a prepared payload and returns one
            :class:`~chemthermo.bench._record.StateOutcome` per state. Called
            once per warm-up and once per timed repeat, so it must be free of
            state that would make the second call cheaper than the first.
        settings: The convergence criteria in force, recorded verbatim.
    """

    id: str
    description: str
    model: str
    route: str
    components: tuple[str, ...]
    prepare: Callable[[], Any]
    invoke: Callable[[Any], tuple[StateOutcome, ...]]
    settings: Mapping[str, Any]


# ---------------------------------------------------------------------------
# Shared readers
# ---------------------------------------------------------------------------


def _flash_settings(settings: ct.FlashSettings | None = None) -> dict[str, Any]:
    """The flash convergence criteria, as the record stores them."""
    active = settings or ct.FlashSettings()
    stability = active.stability_settings or ct.StabilitySettings()
    return {
        "flash_max_iter": active.max_iter,
        "flash_tol": active.tol,
        "flash_second_order": active.second_order,
        "flash_second_order_tol": active.second_order_tol,
        "flash_ssi_iterations": active.ssi_iterations,
        "flash_max_phases": active.max_phases,
        "post_split_stability": active.post_split_stability,
        "phase_detection": active.phase_detection,
        "stability_max_iter": stability.max_iter,
        "stability_tol": stability.tol,
        "stability_tpd_tol": stability.tpd_tol,
        "stability_trivial_tol": stability.trivial_tol,
    }


def _stability_settings(settings: ct.StabilitySettings | None = None) -> dict[str, Any]:
    active = settings or ct.StabilitySettings()
    return {
        "stability_max_iter": active.max_iter,
        "stability_tol": active.tol,
        "stability_tpd_tol": active.tpd_tol,
        "stability_trivial_tol": active.trivial_tol,
        "stability_second_order": active.second_order,
        "stability_ssi_iterations": active.ssi_iterations,
    }


_ITERATION_KEYS = (
    "iterations",
    "ssi_iterations",
    "second_order_iterations",
    "stability_trials",
    "phase_count",
    "search_iterations",
)
_INITIALIZATION_KEYS = ("k_seed", "phase_detection", "incipient_phase", "feed_branch")


def _flash_outcome(
    result: ct.FlashResult, *, temperature_K: float, pressure_Pa: float, z: Sequence[float]
) -> StateOutcome:
    names = tuple(result.phase_names())
    diagnostics = result.diagnostics
    return StateOutcome(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        composition=tuple(float(value) for value in z),
        phases=names,
        phase_compositions=tuple(
            tuple(float(value) for value in result.phases[name].composition.fractions)
            for name in names
        ),
        phase_fractions=tuple(float(result.phase_fractions[name]) for name in names),
        iterations={key: int(diagnostics[key]) for key in _ITERATION_KEYS if key in diagnostics},
        initialization={
            key: str(diagnostics[key]) for key in _INITIALIZATION_KEYS if key in diagnostics
        },
    )


def _refused(
    exc: Exception, *, temperature_K: float, pressure_Pa: float, z: Sequence[float]
) -> StateOutcome:
    return StateOutcome(
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        composition=tuple(float(value) for value in z),
        status=type(exc).__name__,
        error=str(exc),
    )


# ---------------------------------------------------------------------------
# Peng-Robinson
# ---------------------------------------------------------------------------

PR_TERNARY_NAMES = ("Methane", "Ethane", "Propane")
PR_TERNARY_FEED = (0.5, 0.3, 0.2)
PR_TERNARY_T_K = 240.0
PR_TERNARY_P_PA = 3.0e6

#: The ADR-0017 phase-detection grid, restricted to one temperature so the case
#: stays a few hundred milliseconds: six mixtures x four pressures = 24 states.
PR_GRID_MIXTURES: tuple[tuple[tuple[str, ...], tuple[float, ...]], ...] = (
    (("Methane", "Ethane"), (0.5, 0.5)),
    (("Methane", "Propane"), (0.7, 0.3)),
    (("Ethane", "n-Heptane"), (0.7, 0.3)),
    (("Methane", "n-Pentane"), (0.6, 0.4)),
    (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2)),
    (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3)),
)
PR_GRID_T_K = 240.0
PR_GRID_P_PA = (2.0e5, 1.0e6, 3.0e6, 8.0e6)


def _pr_ternary_prepare() -> Any:
    return (
        ct.Mixture.from_database(list(PR_TERNARY_NAMES), list(PR_TERNARY_FEED), normalize=True),
        PengRobinsonEOS(),
    )


def _pr_flash_invoke(payload: Any) -> tuple[StateOutcome, ...]:
    mixture, eos = payload
    result = ct.flash_tp(
        mixture, temperature_K=PR_TERNARY_T_K, pressure_Pa=PR_TERNARY_P_PA, eos=eos
    )
    return (
        _flash_outcome(
            result,
            temperature_K=PR_TERNARY_T_K,
            pressure_Pa=PR_TERNARY_P_PA,
            z=PR_TERNARY_FEED,
        ),
    )


def _pr_stability_invoke(payload: Any) -> tuple[StateOutcome, ...]:
    mixture, eos = payload
    result = ct.stability_tp(
        mixture, temperature_K=PR_TERNARY_T_K, pressure_Pa=PR_TERNARY_P_PA, eos=eos
    )
    # A stability verdict has no phase set, so the "answer" this case hashes is
    # the verdict itself, the composition the minimizing trial stopped at and
    # the tangent-plane distance there - carried in the `phases`,
    # `phase_compositions` and `phase_fractions` slots so that one hash rule
    # covers every case.
    return (
        StateOutcome(
            temperature_K=PR_TERNARY_T_K,
            pressure_Pa=PR_TERNARY_P_PA,
            composition=tuple(float(value) for value in PR_TERNARY_FEED),
            phases=(result.status,),
            phase_compositions=(tuple(float(value) for value in (result.trial_composition or ())),),
            phase_fractions=(float(result.tpd_min),),
            iterations={
                "stability_trials": len(result.trials),
                "trial_iterations": sum(int(trial.iterations) for trial in result.trials),
            },
            initialization={
                "phase_detection": "tangent-plane",
                "feed_branch": str(result.feed_branch),
            },
        ),
    )


def _pr_grid_prepare() -> Any:
    eos = PengRobinsonEOS()
    return [
        (names, z, ct.Mixture.from_database(list(names), list(z), normalize=True), eos)
        for names, z in PR_GRID_MIXTURES
    ]


def _pr_grid_invoke(payload: Any) -> tuple[StateOutcome, ...]:
    outcomes: list[StateOutcome] = []
    for _names, z, mixture, eos in payload:
        for pressure_Pa in PR_GRID_P_PA:
            try:
                result = ct.flash_tp(
                    mixture, temperature_K=PR_GRID_T_K, pressure_Pa=pressure_Pa, eos=eos
                )
            except ct.ChemThermoError as exc:
                outcomes.append(
                    _refused(exc, temperature_K=PR_GRID_T_K, pressure_Pa=pressure_Pa, z=z)
                )
                continue
            outcomes.append(
                _flash_outcome(result, temperature_K=PR_GRID_T_K, pressure_Pa=pressure_Pa, z=z)
            )
    return tuple(outcomes)


# ---------------------------------------------------------------------------
# Activity-model routes
# ---------------------------------------------------------------------------

#: Tessier, Brennecke & Stadtherr (2000) Table 1, Problem 1.
TESSIER_NAMES = ("1-Propanol", "n-Butanol", "Water")
TESSIER_TAU = (
    (0.0, -0.61259, -0.07149),
    (0.7164, 0.0, 0.90047),
    (2.7425, 3.51307, 0.0),
)
TESSIER_ALPHA = (
    (0.0, 0.3, 0.3),
    (0.3, 0.0, 0.48),
    (0.3, 0.48, 0.0),
)
#: The Problem 1 feed of the paper's Table 2.
TESSIER_FEED = (0.148, 0.052, 0.8)
TESSIER_T_K = 298.15
ATMOSPHERE_PA = 101325.0


def _tessier_model() -> ct.NRTL:
    pairs = [
        (
            TESSIER_NAMES[i],
            TESSIER_NAMES[j],
            TESSIER_TAU[i][j],
            TESSIER_TAU[j][i],
            TESSIER_ALPHA[i][j],
            TESSIER_ALPHA[j][i],
        )
        for i in range(len(TESSIER_NAMES))
        for j in range(i + 1, len(TESSIER_NAMES))
    ]
    return ct.NRTL(parameters=ct.NRTLParameters.from_pairs(pairs))


def _nrtl_lle_prepare() -> Any:
    return (
        ct.Mixture.from_database(list(TESSIER_NAMES), list(TESSIER_FEED), normalize=True),
        _tessier_model(),
    )


def _nrtl_lle_invoke(payload: Any) -> tuple[StateOutcome, ...]:
    mixture, model = payload
    result = ct.flash_tp(
        mixture,
        temperature_K=TESSIER_T_K,
        pressure_Pa=ATMOSPHERE_PA,
        activity_model=model,
    )
    return (
        _flash_outcome(
            result, temperature_K=TESSIER_T_K, pressure_Pa=ATMOSPHERE_PA, z=TESSIER_FEED
        ),
    )


#: The modified-Raoult vapour-liquid state of
#: ``examples/basic/flash_tp_modified_raoult_demo.py``.
RAOULT_NAMES = ("1-Propanol", "Water")
RAOULT_FEED = (0.5, 0.5)
RAOULT_T_K = 361.0


def _raoult_prepare() -> Any:
    parameters = ct.NRTLParameters.from_pairs(
        [(RAOULT_NAMES[0], RAOULT_NAMES[1], -0.07149, 2.7425, 0.3, 0.3)]
    )
    return (
        ct.Mixture.from_database(list(RAOULT_NAMES), list(RAOULT_FEED), normalize=True),
        ct.NRTL(parameters=parameters),
    )


def _raoult_invoke(payload: Any) -> tuple[StateOutcome, ...]:
    mixture, model = payload
    result = ct.flash_tp(
        mixture,
        temperature_K=RAOULT_T_K,
        pressure_Pa=ATMOSPHERE_PA,
        activity_model=model,
        flash_mode="modified-raoult",
    )
    return (
        _flash_outcome(result, temperature_K=RAOULT_T_K, pressure_Pa=ATMOSPHERE_PA, z=RAOULT_FEED),
    )


#: The three-phase feed of ``examples/basic/flash_tp_vlle_demo.py``: the
#: centroid of the 364 K tie triangle.
VLLE_T_K = 364.0
VLLE_FEED = (0.13418838, 0.08427618, 0.78153544)


def _vlle_prepare() -> Any:
    return (
        ct.Mixture.from_database(list(TESSIER_NAMES), list(VLLE_FEED), normalize=True),
        _tessier_model(),
    )


def _vlle_invoke(payload: Any) -> tuple[StateOutcome, ...]:
    mixture, model = payload
    result = ct.flash_tp(
        mixture,
        temperature_K=VLLE_T_K,
        pressure_Pa=ATMOSPHERE_PA,
        activity_model=model,
        flash_mode="modified-raoult",
    )
    return (_flash_outcome(result, temperature_K=VLLE_T_K, pressure_Pa=ATMOSPHERE_PA, z=VLLE_FEED),)


# ---------------------------------------------------------------------------
# PC-SAFT
# ---------------------------------------------------------------------------

PCSAFT_VLE_NAMES = ("Methane", "n-Hexane")
PCSAFT_VLE_FEED = (0.5, 0.5)
PCSAFT_VLE_T_K = 300.0
PCSAFT_VLE_P_PA = 3.0e6

PCSAFT_LLE_NAMES = ("Water", "n-Hexane")
PCSAFT_LLE_FEED = (0.5, 0.5)
PCSAFT_LLE_T_K = 298.15


def _pcsaft_vle_prepare() -> Any:
    return (
        ct.Mixture.from_database(list(PCSAFT_VLE_NAMES), list(PCSAFT_VLE_FEED), normalize=True),
        PCSAFTEOS(),
    )


def _pcsaft_vle_invoke(payload: Any) -> tuple[StateOutcome, ...]:
    mixture, eos = payload
    result = ct.flash_tp(
        mixture, temperature_K=PCSAFT_VLE_T_K, pressure_Pa=PCSAFT_VLE_P_PA, eos=eos
    )
    return (
        _flash_outcome(
            result,
            temperature_K=PCSAFT_VLE_T_K,
            pressure_Pa=PCSAFT_VLE_P_PA,
            z=PCSAFT_VLE_FEED,
        ),
    )


def _pcsaft_lle_prepare() -> Any:
    return (
        ct.Mixture.from_database(list(PCSAFT_LLE_NAMES), list(PCSAFT_LLE_FEED), normalize=True),
        PCSAFTEOS(),
    )


def _pcsaft_lle_invoke(payload: Any) -> tuple[StateOutcome, ...]:
    mixture, eos = payload
    result = ct.flash_tp(mixture, temperature_K=PCSAFT_LLE_T_K, pressure_Pa=ATMOSPHERE_PA, eos=eos)
    return (
        _flash_outcome(
            result,
            temperature_K=PCSAFT_LLE_T_K,
            pressure_Pa=ATMOSPHERE_PA,
            z=PCSAFT_LLE_FEED,
        ),
    )


POLYMER_T_K = 453.0
POLYMER_P_PA = 8.0e6
#: 5 wt% polymer, the feed of ``examples/basic/pcsaft_polymer_demo.py``.
POLYMER_WEIGHT_FRACTION = 0.05


def _polymer_prepare() -> Any:
    if not _POLYMER_FIXTURE.is_file():
        raise CaseSkipped(
            "polymer PC-SAFT parameters are not packaged (ADR-0022); this case needs "
            f"{_POLYMER_FIXTURE.name} from a source checkout"
        )
    with _POLYMER_FIXTURE.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    row = next(entry for entry in payload["polymers"] if entry["name"] == "Polyethylene")
    parameters = PCSAFTParameters.from_records(
        [
            PCSAFTRecord(
                name="Polyethylene",
                segments_per_g=row["segments_per_g"],
                MW_g_mol=_POLYMER_MW_G_MOL,
                sigma_A=row["sigma_A"],
                epsilon_k_K=row["epsilon_k_K"],
                source="Martini et al. 2009 Table 1 citing Gross & Sadowski 2002",
            ),
            PCSAFTRecord(
                name="n-Pentane",
                m=2.6896,
                sigma_A=3.7729,
                epsilon_k_K=231.20,
                MW_g_mol=_PENTANE_MW_G_MOL,
                source="Gross & Sadowski 2001 Table 1 (the packaged record)",
            ),
        ]
    )
    moles_polymer = POLYMER_WEIGHT_FRACTION / _POLYMER_MW_G_MOL
    moles_solvent = (1.0 - POLYMER_WEIGHT_FRACTION) / _PENTANE_MW_G_MOL
    total = moles_polymer + moles_solvent
    z = (moles_polymer / total, moles_solvent / total)
    mixture = ct.Mixture.from_components(
        [
            ct.Component.custom(
                "Polyethylene",
                mw_kg_per_mol=_POLYMER_MW_G_MOL / 1000.0,
                formula="(C2H4)n",
                volatile=False,
                source="see tests/fixtures/pcsaft/martini2009_polymers.json",
            ),
            ct.Component.from_database("n-Pentane"),
        ],
        list(z),
        normalize=True,
    )
    return mixture, PCSAFTEOS(parameters=parameters, kij=_POLYMER_KIJ), z


def _polymer_invoke(payload: Any) -> tuple[StateOutcome, ...]:
    mixture, eos, z = payload
    result = ct.flash_tp(mixture, temperature_K=POLYMER_T_K, pressure_Pa=POLYMER_P_PA, eos=eos)
    return (_flash_outcome(result, temperature_K=POLYMER_T_K, pressure_Pa=POLYMER_P_PA, z=z),)


# ---------------------------------------------------------------------------
# The workload
# ---------------------------------------------------------------------------

CASES: tuple[BenchCase, ...] = (
    BenchCase(
        id="pr-flash-ternary",
        description="Peng-Robinson phi-phi flash, methane/ethane/propane at 240 K and 3 MPa",
        model="Peng-Robinson",
        route="flash_tp (phi-phi)",
        components=PR_TERNARY_NAMES,
        prepare=_pr_ternary_prepare,
        invoke=_pr_flash_invoke,
        settings=_flash_settings(),
    ),
    BenchCase(
        id="pr-stability-ternary",
        description="Tangent-plane stability test alone, at the pr-flash-ternary state",
        model="Peng-Robinson",
        route="stability_tp",
        components=PR_TERNARY_NAMES,
        prepare=_pr_ternary_prepare,
        invoke=_pr_stability_invoke,
        settings=_stability_settings(),
    ),
    BenchCase(
        id="pr-flash-grid-24",
        description="24 states of the ADR-0017 Peng-Robinson grid (six mixtures, 240 K, four pressures)",
        model="Peng-Robinson",
        route="flash_tp (phi-phi)",
        components=("mixed",),
        prepare=_pr_grid_prepare,
        invoke=_pr_grid_invoke,
        settings=_flash_settings(),
    ),
    BenchCase(
        id="nrtl-lle-tessier-p1",
        description="NRTL gamma-gamma liquid-liquid split, Tessier (2000) Problem 1 feed",
        model="NRTL",
        route="flash_tp (gamma-gamma)",
        components=TESSIER_NAMES,
        prepare=_nrtl_lle_prepare,
        invoke=_nrtl_lle_invoke,
        settings=_flash_settings(),
    ),
    BenchCase(
        id="modified-raoult-vle",
        description="Modified-Raoult vapour-liquid split, 1-propanol/water at 361 K and 1 atm",
        model="NRTL + ideal gas (Antoine reference)",
        route="flash_tp (modified-raoult)",
        components=RAOULT_NAMES,
        prepare=_raoult_prepare,
        invoke=_raoult_invoke,
        settings=_flash_settings(),
    ),
    BenchCase(
        id="vlle-364k",
        description="Three-phase modified-Raoult search, 364 K tie-triangle centroid feed",
        model="NRTL + ideal gas (Antoine reference)",
        route="flash_tp (modified-raoult, phase addition)",
        components=TESSIER_NAMES,
        prepare=_vlle_prepare,
        invoke=_vlle_invoke,
        settings=_flash_settings(),
    ),
    BenchCase(
        id="pcsaft-vle-methane-hexane",
        description="PC-SAFT phi-phi vapour-liquid flash, methane/n-hexane at 300 K and 3 MPa",
        model="PC-SAFT",
        route="flash_tp (phi-phi)",
        components=PCSAFT_VLE_NAMES,
        prepare=_pcsaft_vle_prepare,
        invoke=_pcsaft_vle_invoke,
        settings=_flash_settings(),
    ),
    BenchCase(
        id="pcsaft-lle-water-hexane",
        description="PC-SAFT associating liquid-liquid flash, water/n-hexane at 298.15 K and 1 atm",
        model="PC-SAFT (2B association)",
        route="flash_tp (phi-phi)",
        components=PCSAFT_LLE_NAMES,
        prepare=_pcsaft_lle_prepare,
        invoke=_pcsaft_lle_invoke,
        settings=_flash_settings(),
    ),
    BenchCase(
        id="pcsaft-polymer-lle",
        description="PC-SAFT polymer/solvent split, polyethylene(16400)/n-pentane at 453 K and 8 MPa",
        model="PC-SAFT (polymer, k_ij = -0.006)",
        route="flash_tp (phi-phi)",
        components=("Polyethylene", "n-Pentane"),
        prepare=_polymer_prepare,
        invoke=_polymer_invoke,
        settings=_flash_settings(),
    ),
)

CASES_BY_ID: Mapping[str, BenchCase] = {case.id: case for case in CASES}
