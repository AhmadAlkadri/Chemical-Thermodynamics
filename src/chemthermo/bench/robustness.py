"""The robustness map: what every model family refuses, and why (ADR-0027).

**Internal module.** Nothing here is re-exported from ``chemthermo`` and
nothing here is part of the public API (ADR-0001). Like the rest of
:mod:`chemthermo.bench` it is a maintainer's instrument, and it is run the same
way::

    python -m chemthermo.bench robustness --out benchmarks/robustness_<sha>.json

What it is for
--------------
ADR-0023 gave this repository an instrument for *speed*. This one is the
instrument for *coverage*: it sweeps every model family over fixed state and
composition grids, classifies each state into a phase verdict or a refusal
class, checks the invariants of every converged answer, and writes one record.
The point is to replace "the polymer path sometimes refuses" with a count, a
class, an exact state and the exact message - so the next solver slice is
picked from evidence rather than from the last thing someone happened to hit.

What it is **not**
------------------
- It is not a correctness check. Nothing here is compared against a published
  number or an independent implementation; the ledger cases are where that
  lives (F-1..F-5, L-1..L-4, P-0..P-16, R-1..R-4, V-1..V-5). What is checked
  here is *internal consistency* - mass balance, the equilibrium residual the
  solver itself reports, a negative Gibbs-energy change, phase fractions in
  range, the post-split verdict - which catches a wrong answer only when it is
  wrong in one of those ways.
- It is not a bit-identity fixture. Phase *compositions* are deliberately not
  recorded: that is what ``benchmarks/*.json`` (ADR-0023) and
  ``refactor_bit_identity_v3.json`` are for, and duplicating them here would
  make a coverage map into a second baseline that has to be regenerated
  whenever a last bit moves.
- Its wall times are not portable between machines, for the reasons
  ``benchmarks/README.md`` gives. They are recorded so the *shape* of the cost
  (which family, which state) is visible, not so two machines can be compared.
- A verdict here is never better than the stability test that produced it
  (brain.md section 10). A state counted "single" is a state the deterministic
  trial set found no negative tangent-plane distance at.

The grids
---------
Ten families, every one of them fixed in this module and documented at its
definition below. The first six are the ADR-0027 grid; the last four
(ADR-0027 amendment, "robustness-map-coverage") fill the gaps that ADR-0027 /
ADR-0028 named as unranked: EOS three-phase windows, the legacy ``gamma-phi``
mode, near-critical PR states and the CO2/n-decane / methane/n-pentane windows
that produced past defects, and associating ternaries.

===============================  =======  =============================================
family                            states   what it sweeps
===============================  =======  =============================================
``pr-phi-phi``                     1270   8 Peng-Robinson mixtures over the Case F-1 wide
                                          scan (11 T x 13 P), plus a second feed for each
                                          ternary on a coarser grid
``pcsaft``                          204   the 188-state Case F-4 grid, plus two nonzero
                                          ``kij`` binaries at 8 states each
``pcsaft-associating``               260   water / n-hexane, water / ethanol, water /
                                          1-propanol / n-hexane, methanol / n-hexane over
                                          4 T x 5 P x 3-4 feeds
``modified-raoult``                  108   Tessier P1 ternary and P2 quaternary at 1 atm
                                          over 5 T x 8 feeds, plus water / 1-butanol
``gamma-gamma``                       16   the Tessier P1 / P2 feeds, activity-only
``polymer``                          252   PE 16400 and PE 53000 in n-pentane at 453 K,
                                          0.3-12 MPa, 3 weight fractions, plus a ternary
``eos-three-phase``                  117   PC-SAFT water/n-hexane around T3 (Case P-9),
                                          PC-SAFT and Peng-Robinson water/ethanol/n-hexane
                                          tie-triangles (Case P-10)
``gamma-phi-legacy``                  30   the deprecated ``flash_mode="gamma-phi"`` path,
                                          NRTL-packaged Methane/Ethane, over the CLI's
                                          state and a small T/P grid
``pr-near-critical``                 104   Peng-Robinson near the two-phase boundary of
                                          two mixtures, plus Case F-2 and the CO2/n-decane
                                          / methane/n-pentane defect windows
``pcsaft-associating-ternary``       144   water/ethanol/n-hexane and water/1-propanol/
                                          n-hexane PC-SAFT over 3 T x 3 P x 8 feeds
===============================  =======  =============================================

``--family NAME`` runs one family, so the sweep partitions and resumes.
``--quick`` runs a fixed cost-bounded subset across all ten, cheap enough for
the default test suite; see :data:`QUICK_NOTE`.

Classification
--------------
Every state ends in exactly one bucket:

- a **verdict** - ``single-liquid``, ``single-vapor``, ``VLE``, ``LLE``,
  ``VLLE``, ``LLL`` or ``other`` - when ``flash_tp`` returned and every
  invariant held;
- ``converged-invariant-violated`` when it returned and one did not (the
  details are recorded per state);
- one of the **refusal classes** of :data:`REFUSAL_CLASSES` when it raised,
  read off the exception type and message by :func:`classify_refusal`.

The message and the state are always recorded, so a class is a starting point
for a diagnosis and never a substitute for one.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from collections import Counter
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

import chemthermo as ct

from ..eos import PCSAFTEOS
from ..models import PengRobinsonEOS
from ..parameters import PCSAFTParameters, PCSAFTRecord
from ._record import environment, git_state

__all__ = [
    "FAMILIES",
    "REFUSAL_CLASSES",
    "VERDICTS",
    "classify_refusal",
    "main",
    "run_sweep",
    "summary_markdown",
    "systems",
]

#: Schema tag written into every record, so a later field change is detectable
#: rather than silently mis-read. Separate from ``chemthermo-bench/1``: this is
#: a different artefact with different fields.
SCHEMA = "chemthermo-robustness/1"

#: Mass-balance tolerance. Recomputed here from the returned phase fractions
#: and compositions rather than read out of ``diagnostics``, so the check does
#: not rest on the solver's own bookkeeping.
MASS_BALANCE_TOL = 1e-10

#: Equilibrium-residual tolerance. This one *is* the solver's own number
#: (``fugacity_residual`` / ``equilibrium_residual``): recomputing equal
#: fugacities here would need per-model code this module deliberately does not
#: carry. It is the same 1e-06 the Case F-3 / F-4 invariants assert.
EQUILIBRIUM_TOL = 1e-6

#: Composition-sum tolerance for a returned phase.
COMPOSITION_TOL = 1e-10

#: How many slowest states the record lists explicitly.
WORST_N = 20

FAMILIES: tuple[str, ...] = (
    "pr-phi-phi",
    "pcsaft",
    "pcsaft-associating",
    "modified-raoult",
    "gamma-gamma",
    "polymer",
    "eos-three-phase",
    "gamma-phi-legacy",
    "pr-near-critical",
    "pcsaft-associating-ternary",
)

VERDICTS: tuple[str, ...] = (
    "single-liquid",
    "single-vapor",
    "VLE",
    "LLE",
    "VLLE",
    "LLL",
    "other",
)

REFUSAL_CLASSES: tuple[str, ...] = (
    "stability-inconclusive",
    "rr-no-bracket",
    "density-root-failure",
    "post-split-third-phase",
    "multiphase-solver-failure",
    "split-non-convergence",
    "model-error",
    "other-refusal",
)

#: The bucket a state lands in when ``flash_tp`` returned an answer that fails
#: one of this module's invariant checks. Counted separately from the refusal
#: classes because it is a strictly worse outcome than a refusal: a wrong
#: answer nobody is told about.
INVARIANT_VIOLATED = "converged-invariant-violated"

QUICK_NOTE = (
    "The quick subset is a fixed offset/stride sample of each system's state list "
    "plus explicitly pinned states. It is COST-BOUNDED, not representative: the "
    "offsets and strides of the costly families (the three original PC-SAFT ones "
    "plus the eos-three-phase and pcsaft-associating-ternary families slice "
    "robustness-map-coverage added) were chosen from the full run's per-state wall "
    "times so that the default test suite stays inside its budget, and the pinned "
    "indices are refusal states, one per refusal stage the full sweep found that "
    "costs under ~1.5 s (a PC-SAFT three-phase state can cost 13-40 s, so those "
    "families' quick bucket counts are not representative of their full-sweep "
    "verdict mix). Read it as a smoke test of the instrument plus those pinned "
    "defects - 224 states in about 14 s, against 2505 states in about 37 minutes "
    "for the full map."
)


# ---------------------------------------------------------------------------
# Record shapes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StateSpec:
    """One swept state: the feed and the (T, P) it is flashed at."""

    composition: tuple[float, ...]
    temperature_K: float
    pressure_Pa: float


@dataclass(frozen=True)
class RobustnessSystem:
    """One mixture + model, and the states it is swept over.

    Attributes:
        family: Which of :data:`FAMILIES` this belongs to.
        name: Stable identifier, unique across all families.
        description: One line of prose for the record and the summary table.
        components: Component names in the mixture's order.
        model: The thermodynamic model, as the record names it.
        route: The public entry point being swept.
        prepare: Builds whatever the flash needs, **once**. May raise
            :class:`SystemSkipped` when its inputs are not available here.
        flash: ``(context, spec) -> FlashResult``.
        states: Every state of this system, in a fixed order.
        quick_offset: First ``states`` index the ``--quick`` sample starts at.
        quick_stride: Stride into ``states`` for the ``--quick`` subset.
        quick_limit: Cap on the strided sample (``None`` for no cap).
        quick_pinned: Extra ``states`` indices always in the quick subset.
    """

    family: str
    name: str
    description: str
    components: tuple[str, ...]
    model: str
    route: str
    prepare: Callable[[], Any]
    flash: Callable[[Any, StateSpec], ct.FlashResult]
    states: tuple[StateSpec, ...]
    quick_offset: int = 0
    quick_stride: int = 1
    quick_limit: int | None = None
    quick_pinned: tuple[int, ...] = ()

    def selected(self, *, quick: bool) -> tuple[tuple[int, StateSpec], ...]:
        """Return ``(index, state)`` pairs for a full or quick sweep."""
        if not quick:
            return tuple(enumerate(self.states))
        sampled = list(range(self.quick_offset, len(self.states), max(1, self.quick_stride)))
        if self.quick_limit is not None:
            sampled = sampled[: self.quick_limit]
        chosen = sorted(set(sampled) | {i for i in self.quick_pinned if i < len(self.states)})
        return tuple((index, self.states[index]) for index in chosen)


class SystemSkipped(Exception):
    """Raised by a system's ``prepare`` when its inputs are not available here."""


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

#: ``(needle, class, stage)`` in priority order; the first needle found in the
#: refusal message decides. Order matters where a message could match twice: a
#: density-root failure inside a split is a density-root failure, and a
#: post-split refusal that names ``max_phases`` is a max-phases refusal rather
#: than a generic stability one.
_REFUSAL_RULES: tuple[tuple[str, str, str], ...] = (
    ("Tangent-plane stability analysis was inconclusive", "stability-inconclusive", "feed"),
    ("A post-split stability test was inconclusive", "stability-inconclusive", "post-split"),
    (
        "Tangent-plane stability reported an unstable feed without a minimizing",
        "stability-inconclusive",
        "no-minimizer",
    ),
    ("failed to bracket a vapor fraction", "rr-no-bracket", "phi-phi"),
    ("bracket a Rachford-Rice root", "rr-no-bracket", "seeded"),
    ("PC-SAFT found no density root", "density-root-failure", "pcsaft-no-root"),
    ("only mechanically unstable density roots", "density-root-failure", "pcsaft-unstable-root"),
    ("No real compressibility roots", "density-root-failure", "cubic-no-root"),
    ("No usable fugacity-coefficient root", "density-root-failure", "split-branch"),
    ("Z <= B for Peng-Robinson", "density-root-failure", "cubic-degenerate"),
    ("packing fraction eta", "density-root-failure", "pcsaft-eta"),
    ("mechanically unstable", "density-root-failure", "unstable-root"),
    ("A third phase is required", "post-split-third-phase", "max-phases"),
    ("A further phase is required", "post-split-third-phase", "max-phases"),
    ("Multiphase Rachford-Rice", "multiphase-solver-failure", "rachford-rice"),
    ("multiphase split", "multiphase-solver-failure", "split"),
    ("multiphase flash", "multiphase-solver-failure", "split"),
    ("phase addition/removal search did not settle", "multiphase-solver-failure", "search"),
    ("multiphase Rachford-Rice feasible region", "multiphase-solver-failure", "feasible-region"),
    (
        "A two-phase set converged to a non-positive phase fraction",
        "multiphase-solver-failure",
        "collapsed",
    ),
    ("phi-phi split in log mole numbers", "split-non-convergence", "log-space"),
    ("log-space split cannot start from this seed", "split-non-convergence", "log-space-seed"),
    ("did not converge the phi-phi split", "split-non-convergence", "phi-phi"),
    ("did not converge the liquid-liquid split", "split-non-convergence", "gamma-gamma"),
    ("did not converge the modified-Raoult split", "split-non-convergence", "modified-raoult"),
    ("did not converge within the iteration limit", "split-non-convergence", "legacy"),
    ("converged to a vapor fraction outside", "split-non-convergence", "beta-outside-window"),
    ("is not a stable phase set", "post-split-third-phase", "unstable-phase-set"),
)


def classify_refusal(exc: BaseException) -> tuple[str, str]:
    """Return ``(class, stage)`` for one refusal.

    The class comes from the message where the message is specific enough, and
    from the exception type otherwise. ``stage`` narrows it - which solver
    stage, or which seed - and is ``""`` when the rule that matched carries no
    stage.

    A class is a *bucket*, not a diagnosis. The exact message and state are
    recorded next to it in the record for exactly that reason.
    """
    message = str(exc)
    for needle, refusal_class, stage in _REFUSAL_RULES:
        if needle in message:
            return refusal_class, stage
    if isinstance(exc, (ct.ModelError, ct.PropertyNotFoundError, ct.PCSAFTParameterError)):
        return "model-error", type(exc).__name__
    if isinstance(exc, (ct.InputRangeError, ct.CompositionError)):
        return "model-error", type(exc).__name__
    return "other-refusal", type(exc).__name__


def _verdict(result: ct.FlashResult) -> str:
    """The phase verdict of a converged result."""
    names = result.phase_names()
    if len(names) == 1:
        return "single-vapor" if names[0].startswith("vapor") else "single-liquid"
    vapors = sum(1 for name in names if name.startswith("vapor"))
    liquids = len(names) - vapors
    if len(names) == 2:
        return "VLE" if vapors else "LLE"
    if len(names) == 3:
        if vapors == 1 and liquids == 2:
            return "VLLE"
        if vapors == 0:
            return "LLL"
    return "other"


def _invariant_violations(
    result: ct.FlashResult, spec: StateSpec
) -> tuple[list[str], dict[str, float]]:
    """Check every converged answer; return the violations and the residuals.

    Mass balance and the composition sums are **recomputed** from the returned
    phases. The equilibrium residual and the Gibbs-energy change are read from
    the solver's own diagnostics, which is stated rather than hidden: a bug
    that made the solver mis-report its own residual would not be caught here.
    """
    diagnostics = result.diagnostics
    names = result.phase_names()
    violations: list[str] = []
    residuals: dict[str, float] = {}

    # Compositions sum to one, phase by phase.
    worst_sum = 0.0
    for name in names:
        fractions = np.asarray(result.phases[name].composition.fractions, dtype=float)
        worst_sum = max(worst_sum, abs(float(fractions.sum()) - 1.0))
    residuals["composition_sum"] = worst_sum
    if worst_sum > COMPOSITION_TOL:
        violations.append(f"composition sum off by {worst_sum:.3e} > {COMPOSITION_TOL:.0e}")

    # Mass balance, recomputed.
    feed = np.asarray(spec.composition, dtype=float)
    feed = feed / feed.sum()
    recombined = np.zeros_like(feed)
    for name in names:
        fraction = float(result.phase_fractions.get(name, 0.0))
        recombined += fraction * np.asarray(result.phases[name].composition.fractions, dtype=float)
    mass_balance = float(np.max(np.abs(recombined - feed)))
    residuals["mass_balance"] = mass_balance
    if mass_balance > MASS_BALANCE_TOL:
        violations.append(f"mass balance {mass_balance:.3e} > {MASS_BALANCE_TOL:.0e}")

    # Phase fractions in (0, 1) and summing to one.
    if len(names) > 1:
        fractions = [float(result.phase_fractions.get(name, float("nan"))) for name in names]
        total = sum(fractions)
        residuals["phase_fraction_sum"] = abs(total - 1.0)
        if abs(total - 1.0) > COMPOSITION_TOL:
            violations.append(f"phase fractions sum to {total!r}")
        for name, fraction in zip(names, fractions):
            if not (0.0 < fraction < 1.0):
                violations.append(f"phase fraction {name}={fraction!r} is not in (0, 1)")

    # The equilibrium residual the solver reports.
    for key in ("fugacity_residual", "equilibrium_residual"):
        if key in diagnostics:
            value = float(diagnostics[key])
            residuals["equilibrium"] = value
            if value > EQUILIBRIUM_TOL:
                violations.append(f"{key} {value:.3e} > {EQUILIBRIUM_TOL:.0e}")
            break

    # The split must lower the Gibbs energy.
    if len(names) > 1 and "delta_g_split_rt" in diagnostics:
        delta_g = float(diagnostics["delta_g_split_rt"])
        residuals["delta_g_split_rt"] = delta_g
        if not (delta_g < 0.0):
            violations.append(f"delta_g_split_rt = {delta_g!r} is not negative")

    # The post-split stability test must have passed where it ran. The legacy
    # gamma-phi path never sets `post_split_checked` (ADR-0007, ADR-0010: it
    # has no stability test to check against), so this is silently dormant
    # there by construction - not by a family-specific carve-out - and the
    # `gamma-phi-legacy` family's summary says so explicitly rather than
    # relying on that absence to be read correctly.
    if diagnostics.get("post_split_checked") and diagnostics.get("post_split_status") != "stable":
        violations.append(f"post_split_status = {diagnostics.get('post_split_status')!r}")

    # A three-phase answer must be the lower-Gibbs choice against the
    # two-phase pair the ADR-0011/ADR-0020 phase addition/removal search
    # started from. Present only on a result that entered the search and
    # returned three phases; absent everywhere else (a two-phase answer, or a
    # result the search never touched), so this is dormant on every family
    # but `eos-three-phase` and the rare state elsewhere that lands on three
    # phases by chance.
    if len(names) == 3 and "delta_g_vs_two_phase_rt" in diagnostics:
        delta_g_vs_two_phase = float(diagnostics["delta_g_vs_two_phase_rt"])
        residuals["delta_g_vs_two_phase_rt"] = delta_g_vs_two_phase
        if not (delta_g_vs_two_phase < 0.0):
            violations.append(f"delta_g_vs_two_phase_rt = {delta_g_vs_two_phase!r} is not negative")

    return violations, residuals


# ---------------------------------------------------------------------------
# Family 1: Peng-Robinson phi-phi
# ---------------------------------------------------------------------------

#: The wide Peng-Robinson scan of validation Case F-1: 11 temperatures and 13
#: geometric pressures. 8 mixtures x 11 x 13 = 1144 states, which is the scan
#: size Case F-1 and ADR-0008 record.
#:
#: **Reconstruction caveat.** Case F-1 names the grid ("8 databank mixtures,
#: T 150-450 K in 11 steps, P 1e5-3.2e7 Pa in 13 geometric steps") and names 6
#: of the mixtures (the in-repo grid) plus the Methane/Propane/n-Decane feed of
#: the adjudicated disagreements; no file in this repository lists all 8. The
#: eighth here is Nitrogen/Methane, chosen because it is the one inorganic /
#: hydrocarbon pair the PC-SAFT validation grid also uses. So this is *a*
#: 1144-state Case F-1-shaped scan, not provably the identical one.
PR_T_K: tuple[float, ...] = tuple(float(150.0 + 30.0 * i) for i in range(11))
PR_P_PA: tuple[float, ...] = tuple(float(value) for value in np.geomspace(1.0e5, 3.2e7, 13))

PR_MIXTURES: tuple[tuple[tuple[str, ...], tuple[float, ...]], ...] = (
    (("Methane", "Ethane"), (0.5, 0.5)),
    (("Methane", "Propane"), (0.7, 0.3)),
    (("Ethane", "n-Heptane"), (0.7, 0.3)),
    (("Methane", "n-Pentane"), (0.6, 0.4)),
    (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2)),
    (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3)),
    (("Methane", "Propane", "n-Decane"), (0.7, 0.2, 0.1)),
    (("Nitrogen", "Methane"), (0.3, 0.7)),
)

#: A second feed for each of the three ternaries, on the coarser grid of every
#: other temperature and every other pressure (6 x 7 = 42 states each). The
#: composition axis is the one the Case F-1 scan does not have, and a ternary
#: is where a second feed can change the answer most.
PR_TERNARY_SECOND_FEEDS: tuple[tuple[tuple[str, ...], tuple[float, ...]], ...] = (
    (("Methane", "Ethane", "Propane"), (0.2, 0.3, 0.5)),
    (("Propane", "n-Butane", "n-Pentane"), (0.1, 0.2, 0.7)),
    (("Methane", "Propane", "n-Decane"), (0.3, 0.3, 0.4)),
)
PR_COARSE_T_K: tuple[float, ...] = PR_T_K[::2]
PR_COARSE_P_PA: tuple[float, ...] = PR_P_PA[::2]


def _pr_prepare(names: Sequence[str], z: Sequence[float]) -> Callable[[], Any]:
    def prepare() -> Any:
        mixture = ct.Mixture.from_database(list(names), list(z), normalize=True)
        return mixture, PengRobinsonEOS()

    return prepare


def _pr_flash(context: Any, spec: StateSpec) -> ct.FlashResult:
    mixture, eos = context
    return ct.flash_tp(
        mixture, temperature_K=spec.temperature_K, pressure_Pa=spec.pressure_Pa, eos=eos
    )


def _pr_composition_prepare(names: Sequence[str]) -> Callable[[], Any]:
    """Like :func:`_pr_prepare`, but the composition varies within one system.

    ``_pr_prepare`` fixes ``z`` at prepare time and ``_pr_flash`` ignores
    ``spec.composition`` entirely, which is fine where every state of a
    system shares one feed (every ``pr-phi-phi`` and ``pr-near-critical``
    system does). The ``eos-three-phase`` PR ternary sweeps 12 feeds within
    one system, so it needs the mixture rebuilt per state instead - the same
    shape :func:`_pcsaft_prepare` / :func:`_pcsaft_flash` already use.
    """

    def prepare() -> Any:
        return tuple(names), PengRobinsonEOS()

    return prepare


def _pr_composition_flash(context: Any, spec: StateSpec) -> ct.FlashResult:
    names, eos = context
    mixture = ct.Mixture.from_database(list(names), list(spec.composition), normalize=True)
    return ct.flash_tp(
        mixture, temperature_K=spec.temperature_K, pressure_Pa=spec.pressure_Pa, eos=eos
    )


def _pr_systems() -> list[RobustnessSystem]:
    built: list[RobustnessSystem] = []
    for names, z in PR_MIXTURES:
        states = tuple(
            StateSpec(tuple(z), temperature_K, pressure_Pa)
            for temperature_K in PR_T_K
            for pressure_Pa in PR_P_PA
        )
        built.append(
            RobustnessSystem(
                family="pr-phi-phi",
                name="pr-" + "-".join(name.lower().replace(" ", "") for name in names),
                description=f"{'/'.join(names)} z={z}, 11 T x 13 P (Case F-1 wide scan)",
                components=tuple(names),
                model="Peng-Robinson (kij = 0)",
                route="flash_tp (phi-phi)",
                prepare=_pr_prepare(names, z),
                flash=_pr_flash,
                states=states,
            )
        )
    for names, z in PR_TERNARY_SECOND_FEEDS:
        states = tuple(
            StateSpec(tuple(z), temperature_K, pressure_Pa)
            for temperature_K in PR_COARSE_T_K
            for pressure_Pa in PR_COARSE_P_PA
        )
        built.append(
            RobustnessSystem(
                family="pr-phi-phi",
                name="pr-" + "-".join(name.lower().replace(" ", "") for name in names) + "-feed2",
                description=f"{'/'.join(names)} z={z}, second feed, 6 T x 7 P",
                components=tuple(names),
                model="Peng-Robinson (kij = 0)",
                route="flash_tp (phi-phi)",
                prepare=_pr_prepare(names, z),
                flash=_pr_flash,
                states=states,
            )
        )
    return built


# ---------------------------------------------------------------------------
# Family 2: PC-SAFT, non-associating
# ---------------------------------------------------------------------------

#: The Case F-4 grid verbatim (``tests/validation/test_flash_split_robustness_pcsaft.py``):
#: 3 x 4 x 4 + 4 x 5 x 7 = 188 states.
PCSAFT_F4_GRID: tuple[
    tuple[tuple[str, str], tuple[float, ...], tuple[float, ...], tuple[float, ...]], ...
] = (
    (
        ("Carbon dioxide", "n-Decane"),
        (0.6, 0.8, 0.9),
        (230.0, 240.0, 250.0, 260.0),
        (1.0e6, 1.5e6, 2.0e6, 2.5e6),
    ),
    (
        ("Methane", "n-Hexane"),
        (0.5, 0.8, 0.9, 0.95),
        (170.0, 180.0, 190.0, 195.0, 200.0),
        (0.5e6, 1.0e6, 1.5e6, 2.0e6, 2.5e6, 3.0e6, 3.5e6),
    ),
)

#: Two binaries run with a **nonzero** ``kij``, 8 states each. chemthermo
#: packages no ``kij`` dataset (ADR-0014), so these two values are
#: ILLUSTRATIVE - they are of the right order for these pairs and are not
#: fitted, not cited and not asserted against anything. They are here because
#: every other PC-SAFT state in this repository runs at ``kij = 0``, and a
#: robustness map that never leaves ``kij = 0`` cannot see a ``kij``-driven
#: refusal.
PCSAFT_KIJ_SYSTEMS: tuple[
    tuple[tuple[str, str], float, tuple[tuple[float, ...], ...], tuple[tuple[float, float], ...]],
    ...,
] = (
    (
        ("Carbon dioxide", "n-Decane"),
        0.12,
        ((0.8, 0.2), (0.5, 0.5)),
        ((250.0, 2.0e6), (280.0, 4.0e6), (310.0, 6.0e6), (340.0, 8.0e6)),
    ),
    (
        ("Methane", "n-Decane"),
        0.045,
        ((0.8, 0.2), (0.5, 0.5)),
        ((250.0, 2.0e6), (280.0, 5.0e6), (310.0, 1.0e7), (340.0, 2.0e7)),
    ),
)


def _pcsaft_prepare(names: Sequence[str], kij: float = 0.0) -> Callable[[], Any]:
    def prepare() -> Any:
        return tuple(names), PCSAFTEOS(kij=kij)

    return prepare


def _pcsaft_flash(context: Any, spec: StateSpec) -> ct.FlashResult:
    names, eos = context
    mixture = ct.Mixture.from_database(list(names), list(spec.composition), normalize=True)
    return ct.flash_tp(
        mixture, temperature_K=spec.temperature_K, pressure_Pa=spec.pressure_Pa, eos=eos
    )


def _pcsaft_systems() -> list[RobustnessSystem]:
    built: list[RobustnessSystem] = []
    for names, feeds, temperatures, pressures in PCSAFT_F4_GRID:
        states = tuple(
            StateSpec((z1, 1.0 - z1), temperature_K, pressure_Pa)
            for z1 in feeds
            for temperature_K in temperatures
            for pressure_Pa in pressures
        )
        built.append(
            RobustnessSystem(
                family="pcsaft",
                name="pcsaft-f4-" + "-".join(name.lower().replace(" ", "") for name in names),
                description=f"{'/'.join(names)}, Case F-4 grid ({len(states)} states)",
                components=tuple(names),
                model="PC-SAFT (kij = 0)",
                route="flash_tp (phi-phi)",
                prepare=_pcsaft_prepare(names),
                flash=_pcsaft_flash,
                states=states,
            )
        )
    for names, kij, feeds, state_pairs in PCSAFT_KIJ_SYSTEMS:
        states = tuple(
            StateSpec(tuple(z), temperature_K, pressure_Pa)
            for z in feeds
            for temperature_K, pressure_Pa in state_pairs
        )
        built.append(
            RobustnessSystem(
                family="pcsaft",
                name="pcsaft-kij-" + "-".join(name.lower().replace(" ", "") for name in names),
                description=f"{'/'.join(names)} at kij = {kij} (illustrative, not fitted)",
                components=tuple(names),
                model=f"PC-SAFT (kij = {kij}, illustrative)",
                route="flash_tp (phi-phi)",
                prepare=_pcsaft_prepare(names, kij),
                flash=_pcsaft_flash,
                states=states,
            )
        )
    return built


# ---------------------------------------------------------------------------
# Family 3: PC-SAFT with association
# ---------------------------------------------------------------------------

#: 4 temperatures x 5 pressures, the same for every associating system.
ASSOCIATING_T_K: tuple[float, ...] = (290.0, 320.0, 350.0, 380.0)
ASSOCIATING_P_PA: tuple[float, ...] = (0.1e6, 0.5e6, 1.0e6, 2.0e6, 5.0e6)

#: The four associating systems and their feeds. Water and the C1-C3 1-alkanols
#: are the 2B associating records ADR-0018 packages; the pairings are the ones
#: this repository already has validated states in (water / n-hexane, Cases
#: P-7..P-9) plus three that it does not, which is where a robustness map earns
#: its keep.
ASSOCIATING_SYSTEMS: tuple[tuple[tuple[str, ...], tuple[tuple[float, ...], ...]], ...] = (
    (("Water", "n-Hexane"), ((0.2, 0.8), (0.5, 0.5), (0.8, 0.2), (0.95, 0.05))),
    (("Water", "Ethanol"), ((0.2, 0.8), (0.5, 0.5), (0.8, 0.2))),
    (
        ("Water", "1-Propanol", "n-Hexane"),
        ((0.4, 0.2, 0.4), (0.2, 0.1, 0.7), (0.7, 0.1, 0.2)),
    ),
    (("Methanol", "n-Hexane"), ((0.2, 0.8), (0.5, 0.5), (0.8, 0.2))),
)


def _associating_systems() -> list[RobustnessSystem]:
    built: list[RobustnessSystem] = []
    for names, feeds in ASSOCIATING_SYSTEMS:
        states = tuple(
            StateSpec(tuple(z), temperature_K, pressure_Pa)
            for z in feeds
            for temperature_K in ASSOCIATING_T_K
            for pressure_Pa in ASSOCIATING_P_PA
        )
        built.append(
            RobustnessSystem(
                family="pcsaft-associating",
                name="assoc-" + "-".join(name.lower().replace(" ", "") for name in names),
                description=f"{'/'.join(names)}, {len(feeds)} feeds x 4 T x 5 P",
                components=tuple(names),
                model="PC-SAFT (2B association, kij = 0)",
                route="flash_tp (phi-phi)",
                prepare=_pcsaft_prepare(names),
                flash=_pcsaft_flash,
                states=states,
            )
        )
    return built


# ---------------------------------------------------------------------------
# Families 4 and 5: the NRTL activity routes
# ---------------------------------------------------------------------------

ATMOSPHERE_PA = 101325.0

#: Tessier, Brennecke & Stadtherr (2000) Problem 1
#: (``tests/fixtures/nrtl/tessier2000_problem1.json``).
TESSIER_P1_NAMES: tuple[str, ...] = ("1-Propanol", "n-Butanol", "Water")
TESSIER_P1_TAU: tuple[tuple[float, ...], ...] = (
    (0.0, -0.61259, -0.07149),
    (0.7164, 0.0, 0.90047),
    (2.7425, 3.51307, 0.0),
)
TESSIER_P1_ALPHA: tuple[tuple[float, ...], ...] = (
    (0.0, 0.3, 0.3),
    (0.3, 0.0, 0.48),
    (0.3, 0.48, 0.0),
)
#: The four Table 2 feeds, then four midpoints between consecutive ones. The
#: midpoints are this module's own, not the paper's, and are labelled so in the
#: record: they exist because a robustness map wants composition coverage and
#: four published feeds are not a grid.
TESSIER_P1_FEEDS: tuple[tuple[float, ...], ...] = (
    (0.148, 0.052, 0.8),
    (0.12, 0.08, 0.8),
    (0.13, 0.07, 0.8),
    (0.12, 0.05, 0.83),
    (0.134, 0.066, 0.8),
    (0.125, 0.075, 0.8),
    (0.125, 0.06, 0.815),
    (0.134, 0.051, 0.815),
)

#: Problem 2 (``tests/fixtures/nrtl/tessier2000_problem2.json``).
TESSIER_P2_NAMES: tuple[str, ...] = ("1-Propanol", "n-Butanol", "Benzene", "Water")
TESSIER_P2_TAU: tuple[tuple[float, ...], ...] = (
    (0.0, 2.16486, 0.23689, 0.1306),
    (-1.2007, 0.0, -0.0973, 0.19154),
    (2.01911, 1.73912, 0.0, 4.01932),
    (2.31985, 4.31706, 4.09334, 0.0),
)
TESSIER_P2_ALPHA: tuple[tuple[float, ...], ...] = (
    (0.0, 0.494, 0.286, 0.282),
    (0.494, 0.0, 0.297, 0.344),
    (0.286, 0.297, 0.0, 0.281),
    (0.282, 0.344, 0.281, 0.0),
)
#: The five Table 5 feeds, then three midpoints (this module's own).
TESSIER_P2_FEEDS: tuple[tuple[float, ...], ...] = (
    (0.148, 0.052, 0.6, 0.2),
    (0.25, 0.25, 0.25, 0.25),
    (0.148, 0.052, 0.7, 0.1),
    (0.25, 0.15, 0.4, 0.2),
    (0.25, 0.15, 0.35, 0.25),
    (0.199, 0.151, 0.425, 0.225),
    (0.199, 0.151, 0.475, 0.175),
    (0.25, 0.15, 0.375, 0.225),
)

#: The modified-Raoult temperature axis: two below the three-phase region, then
#: the 350/364/366 K band validation Cases R-3 and V-1..V-5 live in.
RAOULT_T_K: tuple[float, ...] = (300.0, 330.0, 350.0, 364.0, 366.0)

#: The water / 1-butanol binary of Case R-3 / Case L-3, swept across 320-380 K.
BUTANOL_WATER_NAMES: tuple[str, ...] = ("n-Butanol", "Water")
BUTANOL_WATER_TAU: tuple[float, float] = (0.90047, 3.51307)
BUTANOL_WATER_ALPHA = 0.48
BUTANOL_WATER_T_K: tuple[float, ...] = tuple(float(320.0 + 10.0 * i) for i in range(7))
BUTANOL_WATER_FEEDS: tuple[tuple[float, ...], ...] = (
    (0.05, 0.95),
    (0.2, 0.8),
    (0.5, 0.5),
    (0.8, 0.2),
)


def _nrtl_model(
    names: Sequence[str], tau: Sequence[Sequence[float]], alpha: Sequence[Sequence[float]]
) -> ct.NRTL:
    pairs = [
        (names[i], names[j], tau[i][j], tau[j][i], alpha[i][j], alpha[j][i])
        for i in range(len(names))
        for j in range(i + 1, len(names))
    ]
    return ct.NRTL(parameters=ct.NRTLParameters.from_pairs(pairs))


def _activity_prepare(
    names: Sequence[str], tau: Sequence[Sequence[float]], alpha: Sequence[Sequence[float]]
) -> Callable[[], Any]:
    def prepare() -> Any:
        return tuple(names), _nrtl_model(names, tau, alpha)

    return prepare


def _modified_raoult_flash(context: Any, spec: StateSpec) -> ct.FlashResult:
    names, model = context
    mixture = ct.Mixture.from_database(list(names), list(spec.composition), normalize=True)
    return ct.flash_tp(
        mixture,
        temperature_K=spec.temperature_K,
        pressure_Pa=spec.pressure_Pa,
        activity_model=model,
        flash_mode="modified-raoult",
    )


def _gamma_gamma_flash(context: Any, spec: StateSpec) -> ct.FlashResult:
    names, model = context
    mixture = ct.Mixture.from_database(list(names), list(spec.composition), normalize=True)
    return ct.flash_tp(
        mixture,
        temperature_K=spec.temperature_K,
        pressure_Pa=spec.pressure_Pa,
        activity_model=model,
    )


def _modified_raoult_systems() -> list[RobustnessSystem]:
    built: list[RobustnessSystem] = []
    for label, names, tau, alpha, feeds in (
        ("p1", TESSIER_P1_NAMES, TESSIER_P1_TAU, TESSIER_P1_ALPHA, TESSIER_P1_FEEDS),
        ("p2", TESSIER_P2_NAMES, TESSIER_P2_TAU, TESSIER_P2_ALPHA, TESSIER_P2_FEEDS),
    ):
        states = tuple(
            StateSpec(tuple(z), temperature_K, ATMOSPHERE_PA)
            for z in feeds
            for temperature_K in RAOULT_T_K
        )
        built.append(
            RobustnessSystem(
                family="modified-raoult",
                name=f"raoult-tessier-{label}",
                description=(
                    f"Tessier (2000) Problem {label.upper()} {'/'.join(names)} at 1 atm, "
                    f"{len(feeds)} feeds x {len(RAOULT_T_K)} T"
                ),
                components=tuple(names),
                model="NRTL + ideal gas (Antoine reference)",
                route="flash_tp (modified-raoult)",
                prepare=_activity_prepare(names, tau, alpha),
                flash=_modified_raoult_flash,
                states=states,
            )
        )

    butanol_tau = (
        (0.0, BUTANOL_WATER_TAU[0]),
        (BUTANOL_WATER_TAU[1], 0.0),
    )
    butanol_alpha = ((0.0, BUTANOL_WATER_ALPHA), (BUTANOL_WATER_ALPHA, 0.0))
    states = tuple(
        StateSpec(tuple(z), temperature_K, ATMOSPHERE_PA)
        for z in BUTANOL_WATER_FEEDS
        for temperature_K in BUTANOL_WATER_T_K
    )
    built.append(
        RobustnessSystem(
            family="modified-raoult",
            name="raoult-butanol-water",
            description="n-Butanol/Water at 1 atm, 4 feeds x 320-380 K in 10 K steps",
            components=tuple(BUTANOL_WATER_NAMES),
            model="NRTL + ideal gas (Antoine reference)",
            route="flash_tp (modified-raoult)",
            prepare=_activity_prepare(BUTANOL_WATER_NAMES, butanol_tau, butanol_alpha),
            flash=_modified_raoult_flash,
            states=states,
        )
    )
    return built


def _gamma_gamma_systems() -> list[RobustnessSystem]:
    built: list[RobustnessSystem] = []
    for label, names, tau, alpha, feeds in (
        ("p1", TESSIER_P1_NAMES, TESSIER_P1_TAU, TESSIER_P1_ALPHA, TESSIER_P1_FEEDS),
        ("p2", TESSIER_P2_NAMES, TESSIER_P2_TAU, TESSIER_P2_ALPHA, TESSIER_P2_FEEDS),
    ):
        # The Tessier tau are dimensionless and temperature-independent, so
        # sweeping T here would be a repetition, not a state. One temperature,
        # every feed: 8 + 8 = 16 states.
        states = tuple(StateSpec(tuple(z), 298.15, ATMOSPHERE_PA) for z in feeds)
        built.append(
            RobustnessSystem(
                family="gamma-gamma",
                name=f"gamma-tessier-{label}",
                description=(
                    f"Tessier (2000) Problem {label.upper()} {'/'.join(names)}, "
                    f"{len(feeds)} feeds, activity-only"
                ),
                components=tuple(names),
                model="NRTL (gamma-gamma)",
                route="flash_tp (gamma-gamma)",
                prepare=_activity_prepare(names, tau, alpha),
                flash=_gamma_gamma_flash,
                states=states,
            )
        )
    return built


# ---------------------------------------------------------------------------
# Family 6: polymer / solvent
# ---------------------------------------------------------------------------

_POLYMER_FIXTURE = (
    Path(__file__).resolve().parents[3]
    / "tests"
    / "fixtures"
    / "pcsaft"
    / "martini2009_polymers.json"
)
_PENTANE_MW_G_MOL = 72.146
_HEXANE_MW_G_MOL = 86.177

#: ``k_ij`` for polyethylene / n-pentane. Martini et al. (2009) Table 3 lists
#: -0.006 for the Mw = 16400 sample only; the Mw = 53000 sample is not in that
#: table and this repository has used -0.006 for it throughout (ADR-0024,
#: ADR-0026, ``tests/test_pcsaft_polymer.py``). Kept the same here so the sweep
#: is comparable with those cases.
POLYMER_KIJ = -0.006
POLYMER_T_K = 453.0
POLYMER_MW_G_MOL: tuple[float, ...] = (16400.0, 53000.0)
#: 0.3 to 12.0 MPa in 0.3 MPa steps: 40 pressures.
POLYMER_P_PA: tuple[float, ...] = tuple(float(0.3e6 * (i + 1)) for i in range(40))
#: Weight fractions of polymer: dilute, the ADR-0022 demo feed, concentrated.
POLYMER_WEIGHT_FRACTIONS: tuple[float, ...] = (0.01, 0.05, 0.15)
#: The ternary: the same polymer in a 50/50 (by mass) pentane / hexane solvent.
POLYMER_TERNARY_P_PA: tuple[float, ...] = (1.0e6, 3.0e6, 5.0e6, 7.0e6, 9.0e6, 1.1e7)


def _polymer_parameters(mw_g_mol: float, *, with_hexane: bool) -> PCSAFTParameters:
    with _POLYMER_FIXTURE.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    row = next(entry for entry in payload["polymers"] if entry["name"] == "Polyethylene")
    records = [
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
            MW_g_mol=_PENTANE_MW_G_MOL,
            source="Gross & Sadowski 2001 Table 1 (the packaged record)",
        ),
    ]
    if with_hexane:
        records.append(
            PCSAFTRecord(
                name="n-Hexane",
                m=3.0576,
                sigma_A=3.7983,
                epsilon_k_K=236.77,
                MW_g_mol=_HEXANE_MW_G_MOL,
                source="Gross & Sadowski 2001 Table 1 (the packaged record)",
            )
        )
    return PCSAFTParameters.from_records(records)


def _weight_to_mole(weights: Sequence[float], molar_masses: Sequence[float]) -> tuple[float, ...]:
    moles = [w / m for w, m in zip(weights, molar_masses)]
    total = sum(moles)
    return tuple(value / total for value in moles)


def _polymer_prepare(mw_g_mol: float, *, with_hexane: bool) -> Callable[[], Any]:
    def prepare() -> Any:
        if not _POLYMER_FIXTURE.is_file():
            raise SystemSkipped(
                "polymer PC-SAFT parameters are not packaged (ADR-0022); this family "
                f"needs {_POLYMER_FIXTURE.name} from a source checkout"
            )
        parameters = _polymer_parameters(mw_g_mol, with_hexane=with_hexane)
        components = [
            ct.Component.custom(
                "Polyethylene",
                mw_kg_per_mol=mw_g_mol / 1000.0,
                formula="(C2H4)n",
                volatile=False,
                source="see tests/fixtures/pcsaft/martini2009_polymers.json",
            ),
            ct.Component.from_database("n-Pentane"),
        ]
        if with_hexane:
            components.append(ct.Component.from_database("n-Hexane"))
        eos = PCSAFTEOS(parameters=parameters, kij={("Polyethylene", "n-Pentane"): POLYMER_KIJ})
        return components, eos

    return prepare


def _polymer_flash(context: Any, spec: StateSpec) -> ct.FlashResult:
    components, eos = context
    mixture = ct.Mixture.from_components(components, list(spec.composition), normalize=True)
    return ct.flash_tp(
        mixture, temperature_K=spec.temperature_K, pressure_Pa=spec.pressure_Pa, eos=eos
    )


def _polymer_systems() -> list[RobustnessSystem]:
    built: list[RobustnessSystem] = []
    for mw_g_mol in POLYMER_MW_G_MOL:
        states = tuple(
            StateSpec(
                _weight_to_mole((weight, 1.0 - weight), (mw_g_mol, _PENTANE_MW_G_MOL)),
                POLYMER_T_K,
                pressure_Pa,
            )
            for weight in POLYMER_WEIGHT_FRACTIONS
            for pressure_Pa in POLYMER_P_PA
        )
        built.append(
            RobustnessSystem(
                family="polymer",
                name=f"polymer-pe{int(mw_g_mol)}-pentane",
                description=(
                    f"Polyethylene(Mw={int(mw_g_mol)})/n-Pentane at 453 K, "
                    "0.3-12 MPa in 0.3 MPa steps, 1/5/15 wt% polymer"
                ),
                components=("Polyethylene", "n-Pentane"),
                model=f"PC-SAFT (polymer, kij = {POLYMER_KIJ})",
                route="flash_tp (phi-phi)",
                prepare=_polymer_prepare(mw_g_mol, with_hexane=False),
                flash=_polymer_flash,
                states=states,
            )
        )
    ternary_states = tuple(
        StateSpec(
            _weight_to_mole((0.05, 0.475, 0.475), (mw_g_mol, _PENTANE_MW_G_MOL, _HEXANE_MW_G_MOL)),
            POLYMER_T_K,
            pressure_Pa,
        )
        for mw_g_mol in POLYMER_MW_G_MOL
        for pressure_Pa in POLYMER_TERNARY_P_PA
    )
    # Both molar masses share one system record; the molar mass is part of the
    # model, so the states are split into two systems instead.
    for index, mw_g_mol in enumerate(POLYMER_MW_G_MOL):
        built.append(
            RobustnessSystem(
                family="polymer",
                name=f"polymer-pe{int(mw_g_mol)}-pentane-hexane",
                description=(
                    f"Polyethylene(Mw={int(mw_g_mol)})/n-Pentane/n-Hexane at 453 K, "
                    "5 wt% polymer in a 50/50 by mass solvent, 6 pressures"
                ),
                components=("Polyethylene", "n-Pentane", "n-Hexane"),
                model=f"PC-SAFT (polymer, kij = {POLYMER_KIJ} on the pentane pair only)",
                route="flash_tp (phi-phi)",
                prepare=_polymer_prepare(mw_g_mol, with_hexane=True),
                flash=_polymer_flash,
                states=ternary_states[
                    index * len(POLYMER_TERNARY_P_PA) : (index + 1) * len(POLYMER_TERNARY_P_PA)
                ],
            )
        )
    return built


# ---------------------------------------------------------------------------
# Family 7: EOS three-phase windows (slice robustness-map-coverage)
# ---------------------------------------------------------------------------
#
# ADR-0027's "What remains" and brain.md roadmap item 1 named this the biggest
# gap: "no three-phase EOS window is in the grid (the only VLLE states are 4
# modified-Raoult ones)". ADR-0011 / ADR-0020 gave `flash_tp` a phase
# addition/removal search that serves both the activity path and an equation
# of state; this family sweeps the three windows that search is validated on
# (Cases P-9, P-10) instead of the single feed each ledger case happens to
# pin.

#: The water / n-hexane three-phase temperature, from the independent
#: 4-equation Newton of validation Case P-9 (`tests/test_flash_vlle_eos.py`,
#: residual 1.74e-12, cross-checked against FeOs to 6.8e-11). Cited here
#: rather than re-derived, so building the grid needs no solve.
EOS3P_T3_K = 334.807826336

#: Offsets around T3: four below (where the search runs V -> LV -> LLV -> LL,
#: Case P-9 (i)), the point itself (Case P-9 (iii), no three-phase
#: `FlashResult` claimed - Gibbs' phase rule), four above (Case P-9 (ii),
#: straight VLE, the search is not entered).
EOS3P_T3_OFFSETS_K: tuple[float, ...] = (-1.0, -0.5, -0.1, -0.01, 0.0, 0.01, 0.1, 0.5, 1.0)
EOS3P_WATER_HEXANE_Z: tuple[float, ...] = (0.05, 0.3, 0.5, 0.7, 0.95)

#: PC-SAFT water / ethanol / n-hexane tie-triangle temperatures: 333 K is
#: Case P-10 (i)'s validated VLLE point, 331/335/337 K bracket it (Case P-10
#: (i) itself scans 328-337 K and finds the region closes below ~331 K).
EOS3P_TERNARY_T_K: tuple[float, ...] = (331.0, 333.0, 335.0, 337.0)

#: A coarse simplex grid, every mole fraction strictly positive as
#: `test_every_swept_composition_is_a_normalizable_feed` requires: the three
#: corners, the three edge midpoints, three interior points, the centroid and
#: two of the feeds Case P-10 itself already uses (kept for continuity, not
#: because they are special).
EOS3P_TERNARY_FEEDS: tuple[tuple[float, float, float], ...] = (
    (0.8, 0.1, 0.1),
    (0.1, 0.8, 0.1),
    (0.1, 0.1, 0.8),
    (0.6, 0.2, 0.2),
    (0.2, 0.6, 0.2),
    (0.2, 0.2, 0.6),
    (0.4, 0.4, 0.2),
    (0.4, 0.2, 0.4),
    (0.2, 0.4, 0.4),
    (0.34, 0.33, 0.33),
    (0.5, 0.3, 0.2),
    (0.3, 0.3, 0.4),
)

#: The Peng-Robinson three-*liquid* temperature of Case P-10 (ii) (280 K) and
#: a second temperature on the same feed grid (300 K), `kij = 0` throughout.
EOS3P_PR_T_K: tuple[float, ...] = (280.0, 300.0)


def _eos_three_phase_systems() -> list[RobustnessSystem]:
    built: list[RobustnessSystem] = []

    t3_states = tuple(
        StateSpec((z1, 1.0 - z1), EOS3P_T3_K + offset, ATMOSPHERE_PA)
        for offset in EOS3P_T3_OFFSETS_K
        for z1 in EOS3P_WATER_HEXANE_Z
    )
    built.append(
        RobustnessSystem(
            family="eos-three-phase",
            name="eos3p-pcsaft-water-n-hexane-t3-scan",
            description=(
                "Water/n-Hexane PC-SAFT (2B water, kij = 0) at 1 atm, T3 +/- 1 K "
                f"({len(EOS3P_T3_OFFSETS_K)} offsets) x {len(EOS3P_WATER_HEXANE_Z)} "
                "z_water feeds (Case P-9)"
            ),
            components=("Water", "n-Hexane"),
            model="PC-SAFT (2B association, kij = 0)",
            route="flash_tp (phi-phi, phase addition/removal, ADR-0020)",
            prepare=_pcsaft_prepare(("Water", "n-Hexane")),
            flash=_pcsaft_flash,
            states=t3_states,
        )
    )

    ternary_states = tuple(
        StateSpec(feed, temperature_K, ATMOSPHERE_PA)
        for temperature_K in EOS3P_TERNARY_T_K
        for feed in EOS3P_TERNARY_FEEDS
    )
    built.append(
        RobustnessSystem(
            family="eos-three-phase",
            name="eos3p-pcsaft-water-ethanol-n-hexane",
            description=(
                "Water/Ethanol/n-Hexane PC-SAFT (2B water, 2B ethanol, kij = 0) at "
                f"1 atm, {len(EOS3P_TERNARY_T_K)} T x {len(EOS3P_TERNARY_FEEDS)} feeds "
                "(Case P-10 tie-triangle)"
            ),
            components=("Water", "Ethanol", "n-Hexane"),
            model="PC-SAFT (2B association, kij = 0)",
            route="flash_tp (phi-phi, phase addition/removal, ADR-0020)",
            prepare=_pcsaft_prepare(("Water", "Ethanol", "n-Hexane")),
            flash=_pcsaft_flash,
            states=ternary_states,
        )
    )

    pr_states = tuple(
        StateSpec(feed, temperature_K, ATMOSPHERE_PA)
        for temperature_K in EOS3P_PR_T_K
        for feed in EOS3P_TERNARY_FEEDS
    )
    built.append(
        RobustnessSystem(
            family="eos-three-phase",
            name="eos3p-pr-water-ethanol-n-hexane",
            description=(
                "Water/Ethanol/n-Hexane Peng-Robinson (kij = 0) at 1 atm, "
                f"{len(EOS3P_PR_T_K)} T x {len(EOS3P_TERNARY_FEEDS)} feeds "
                "(Case P-10 three-liquid state)"
            ),
            components=("Water", "Ethanol", "n-Hexane"),
            model="Peng-Robinson (kij = 0)",
            route="flash_tp (phi-phi, phase addition/removal, ADR-0020)",
            prepare=_pr_composition_prepare(("Water", "Ethanol", "n-Hexane")),
            flash=_pr_composition_flash,
            states=pr_states,
        )
    )
    return built


# ---------------------------------------------------------------------------
# Family 8: the legacy gamma-phi path
# ---------------------------------------------------------------------------
#
# ADR-0027's roadmap said it plainly: "`gamma-phi` is not swept at all".
# `flash_mode="gamma-phi"` has been DEPRECATED since ADR-0010 - no stability
# test (ADR-0007), Wilson-heuristic phase detection only (ADR-0008 decision
# 4) - but it is still public and still reachable from the CLI
# (`--flash-mode gamma-phi`), whose only packaged activity parameters are the
# synthetic Methane/Ethane NRTL pair (`tests/fixtures` has no gamma-phi
# fixture; the pair lives in `chemthermo/parameters/data/activity/nrtl.json`).
# This family sweeps it on that pair with Peng-Robinson, over a grid that
# contains the CLI contract test's exact state (240 K, 3 MPa,
# `tests/test_cli_tp_flash.py`).

GAMMA_PHI_LEGACY_NAMES: tuple[str, ...] = ("Methane", "Ethane")
GAMMA_PHI_LEGACY_Z: tuple[float, float] = (0.5, 0.5)
GAMMA_PHI_LEGACY_T_K: tuple[float, ...] = (200.0, 220.0, 240.0, 260.0, 280.0, 300.0)
GAMMA_PHI_LEGACY_P_PA: tuple[float, ...] = (1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6)


def _gamma_phi_legacy_prepare() -> Callable[[], Any]:
    def prepare() -> Any:
        mixture = ct.Mixture.from_database(
            list(GAMMA_PHI_LEGACY_NAMES), list(GAMMA_PHI_LEGACY_Z), normalize=True
        )
        return mixture, PengRobinsonEOS(), ct.NRTL()

    return prepare


def _gamma_phi_legacy_flash(context: Any, spec: StateSpec) -> ct.FlashResult:
    mixture, eos, model = context
    return ct.flash_tp(
        mixture,
        temperature_K=spec.temperature_K,
        pressure_Pa=spec.pressure_Pa,
        eos=eos,
        activity_model=model,
        flash_mode="gamma-phi",
    )


def _gamma_phi_legacy_systems() -> list[RobustnessSystem]:
    states = tuple(
        StateSpec(GAMMA_PHI_LEGACY_Z, temperature_K, pressure_Pa)
        for temperature_K in GAMMA_PHI_LEGACY_T_K
        for pressure_Pa in GAMMA_PHI_LEGACY_P_PA
    )
    return [
        RobustnessSystem(
            family="gamma-phi-legacy",
            name="gammaphi-methane-ethane",
            description=(
                "Methane/Ethane z=(0.5, 0.5), NRTL (packaged synthetic pair) + "
                f"Peng-Robinson, DEPRECATED gamma-phi mode, {len(GAMMA_PHI_LEGACY_T_K)} T "
                f"x {len(GAMMA_PHI_LEGACY_P_PA)} P including the CLI contract state "
                "(240 K, 3 MPa)"
            ),
            components=GAMMA_PHI_LEGACY_NAMES,
            model="Peng-Robinson (kij = 0) vapor / NRTL (synthetic pair) liquid",
            route="flash_tp (gamma-phi, DEPRECATED, ADR-0010)",
            prepare=_gamma_phi_legacy_prepare(),
            flash=_gamma_phi_legacy_flash,
            states=states,
        )
    ]


# ---------------------------------------------------------------------------
# Family 9: Peng-Robinson near-critical and past-defect windows
# ---------------------------------------------------------------------------
#
# Two mixtures swept fine around a two-phase boundary located here by
# bisecting `stability_tp`'s verdict (stable / unstable) - a plain binary
# search on pressure at a fixed anchor temperature, run once and its result
# cited as a constant (like `EOS3P_T3_K` above), so building the grid needs no
# solve. Both anchors are states this repository already names: the
# Methane/Ethane/Propane one is ADR-0016's own words, "the weakly unstable
# near-critical methane/ethane/propane at 290 K / 8 MPa"; the Methane/n-Pentane
# one is Case F-2's disagreement state. Plus Case F-2's three states verbatim,
# and two more windows - CO2/n-decane around the ADR-0016 negative-flash
# defect state (PC-SAFT there, Peng-Robinson here: new coverage of the same
# binary and state) and a second Methane/n-Pentane window - both "windows that
# produced past defects" per the slice declaration.
NEARCRIT_MEP_NAMES: tuple[str, ...] = ("Methane", "Ethane", "Propane")
NEARCRIT_MEP_Z: tuple[float, ...] = (0.5, 0.3, 0.2)
NEARCRIT_MEP_T0_K = 290.0
#: Bisected on stability_tp's stable/unstable verdict at T0, 60 iterations
#: halving a (8.0e6, 8.2e6) Pa bracket: stable above, unstable below, to a
#: final bracket width under 2e-9 Pa.
NEARCRIT_MEP_P0_PA = 8_076_427.991645763

NEARCRIT_MP_NAMES: tuple[str, ...] = ("Methane", "n-Pentane")
NEARCRIT_MP_Z: tuple[float, ...] = (0.6, 0.4)
NEARCRIT_MP_T0_K = 175.0
#: Bisected the same way over a (1.778e6, 2.0e6) Pa bracket.
NEARCRIT_MP_P0_PA = 1_884_640.4242924503

#: Case F-2's disagreement state and its two further states (all
#: Methane(0.6)/n-Pentane(0.4)): the legacy Wilson heuristic calls each one
#: single-phase, the tangent plane finds a two-phase split, and both agree
#: with `thermo`'s `FlashVL` and a Gibbs-energy comparison.
F2_STATES_K_PA: tuple[tuple[float, float], ...] = (
    (175.0, 1.778e6),
    (150.0, 6.8399e5),
    (210.0, 4.6784e6),
)

#: CO2/n-decane around the ADR-0016 negative-flash defect state (successive
#: substitution oscillates without a curvature-safeguarded stage; PC-SAFT
#: there, `kij = 0`).
CO2_DECANE_WINDOW_NAMES: tuple[str, ...] = ("Carbon dioxide", "n-Decane")
CO2_DECANE_WINDOW_Z: tuple[float, float] = (0.9, 0.1)
CO2_DECANE_WINDOW_T_K: tuple[float, ...] = (230.0, 240.0, 250.0)
CO2_DECANE_WINDOW_P_PA: tuple[float, ...] = (0.5e6, 1.0e6, 1.5e6, 2.0e6)

#: A second Methane/n-Pentane window, away from Case F-2's own three states,
#: bracketing the Rachford-Rice region found while bisecting
#: `NEARCRIT_MP_P0_PA` above.
METHANE_PENTANE_WINDOW_T_K: tuple[float, ...] = (185.0, 195.0, 205.0)
METHANE_PENTANE_WINDOW_P_PA: tuple[float, ...] = (1.5e6, 2.5e6, 3.5e6)


def _nearcrit_grid_states(
    composition: tuple[float, ...], t0_k: float, p0_pa: float
) -> tuple[StateSpec, ...]:
    """8 T x 5 P states: T0 +/- 5 K, P0 +/- 10 %, around a bisected boundary."""
    temperatures = tuple(t0_k - 5.0 + 10.0 * i / 7.0 for i in range(8))
    pressures = tuple(p0_pa * (0.9 + 0.2 * i / 4.0) for i in range(5))
    return tuple(
        StateSpec(composition, temperature_K, pressure_Pa)
        for temperature_K in temperatures
        for pressure_Pa in pressures
    )


def _nearcrit_systems() -> list[RobustnessSystem]:
    return [
        RobustnessSystem(
            family="pr-near-critical",
            name="nearcrit-methane-ethane-propane",
            description=(
                "Methane/Ethane/Propane z=(0.5, 0.3, 0.2), 8 T x 5 P around the "
                f"boundary at {NEARCRIT_MEP_T0_K:g} K located by bisecting "
                "stability_tp (ADR-0016's near-critical state)"
            ),
            components=NEARCRIT_MEP_NAMES,
            model="Peng-Robinson (kij = 0)",
            route="flash_tp (phi-phi)",
            prepare=_pr_prepare(NEARCRIT_MEP_NAMES, NEARCRIT_MEP_Z),
            flash=_pr_flash,
            states=_nearcrit_grid_states(NEARCRIT_MEP_Z, NEARCRIT_MEP_T0_K, NEARCRIT_MEP_P0_PA),
        ),
        RobustnessSystem(
            family="pr-near-critical",
            name="nearcrit-methane-n-pentane",
            description=(
                "Methane/n-Pentane z=(0.6, 0.4), 8 T x 5 P around the boundary at "
                f"{NEARCRIT_MP_T0_K:g} K located by bisecting stability_tp"
            ),
            components=NEARCRIT_MP_NAMES,
            model="Peng-Robinson (kij = 0)",
            route="flash_tp (phi-phi)",
            prepare=_pr_prepare(NEARCRIT_MP_NAMES, NEARCRIT_MP_Z),
            flash=_pr_flash,
            states=_nearcrit_grid_states(NEARCRIT_MP_Z, NEARCRIT_MP_T0_K, NEARCRIT_MP_P0_PA),
        ),
        RobustnessSystem(
            family="pr-near-critical",
            name="nearcrit-f2-methane-n-pentane",
            description="Methane/n-Pentane z=(0.6, 0.4), Case F-2's 3 states",
            components=NEARCRIT_MP_NAMES,
            model="Peng-Robinson (kij = 0)",
            route="flash_tp (phi-phi)",
            prepare=_pr_prepare(NEARCRIT_MP_NAMES, NEARCRIT_MP_Z),
            flash=_pr_flash,
            states=tuple(
                StateSpec(NEARCRIT_MP_Z, temperature_K, pressure_Pa)
                for temperature_K, pressure_Pa in F2_STATES_K_PA
            ),
        ),
        RobustnessSystem(
            family="pr-near-critical",
            name="nearcrit-window-co2-n-decane",
            description=(
                "Carbon dioxide/n-Decane z=(0.9, 0.1), 3 T x 4 P around the "
                "ADR-0016 negative-flash defect state (240 K, 1 MPa; PC-SAFT "
                "there, Peng-Robinson here)"
            ),
            components=CO2_DECANE_WINDOW_NAMES,
            model="Peng-Robinson (kij = 0)",
            route="flash_tp (phi-phi)",
            prepare=_pr_prepare(CO2_DECANE_WINDOW_NAMES, CO2_DECANE_WINDOW_Z),
            flash=_pr_flash,
            states=tuple(
                StateSpec(CO2_DECANE_WINDOW_Z, temperature_K, pressure_Pa)
                for temperature_K in CO2_DECANE_WINDOW_T_K
                for pressure_Pa in CO2_DECANE_WINDOW_P_PA
            ),
        ),
        RobustnessSystem(
            family="pr-near-critical",
            name="nearcrit-window-methane-n-pentane",
            description=(
                "Methane/n-Pentane z=(0.6, 0.4), 3 T x 3 P, a second window away "
                "from Case F-2's own three states"
            ),
            components=NEARCRIT_MP_NAMES,
            model="Peng-Robinson (kij = 0)",
            route="flash_tp (phi-phi)",
            prepare=_pr_prepare(NEARCRIT_MP_NAMES, NEARCRIT_MP_Z),
            flash=_pr_flash,
            states=tuple(
                StateSpec(NEARCRIT_MP_Z, temperature_K, pressure_Pa)
                for temperature_K in METHANE_PENTANE_WINDOW_T_K
                for pressure_Pa in METHANE_PENTANE_WINDOW_P_PA
            ),
        ),
    ]


# ---------------------------------------------------------------------------
# Family 10: PC-SAFT associating ternaries near their cloud points
# ---------------------------------------------------------------------------
#
# `pcsaft-associating` (family 3) sweeps four binary/ternary associating
# systems but only one of them - water/1-propanol/n-hexane - is a ternary.
# This family adds a second associating ternary and a finer feed grid on both,
# aimed at their cloud points (the LLE/VLE boundary a ternary tie-line sweep
# crosses) rather than at the four-feed coverage sample family 3 already has.

ASSOC3_T_K: tuple[float, ...] = (320.0, 340.0, 360.0)
ASSOC3_P_PA: tuple[float, ...] = (0.1e6, 0.5e6, 1.0e6)
#: A coarse simplex grid, every mole fraction strictly positive: the three
#: near-corner points, the three near-edge midpoints, one near-centroid and
#: the ADR-0022-style 40/30/30 feed.
ASSOC3_FEEDS: tuple[tuple[float, float, float], ...] = (
    (0.7, 0.15, 0.15),
    (0.15, 0.7, 0.15),
    (0.15, 0.15, 0.7),
    (0.5, 0.25, 0.25),
    (0.25, 0.5, 0.25),
    (0.25, 0.25, 0.5),
    (0.4, 0.3, 0.3),
    (0.34, 0.33, 0.33),
)
ASSOC3_SYSTEMS: tuple[tuple[str, ...], ...] = (
    ("Water", "Ethanol", "n-Hexane"),
    ("Water", "1-Propanol", "n-Hexane"),
)


def _associating_ternary_systems() -> list[RobustnessSystem]:
    built: list[RobustnessSystem] = []
    for names in ASSOC3_SYSTEMS:
        states = tuple(
            StateSpec(feed, temperature_K, pressure_Pa)
            for temperature_K in ASSOC3_T_K
            for pressure_Pa in ASSOC3_P_PA
            for feed in ASSOC3_FEEDS
        )
        built.append(
            RobustnessSystem(
                family="pcsaft-associating-ternary",
                name="assoc3-" + "-".join(name.lower().replace(" ", "") for name in names),
                description=(
                    f"{'/'.join(names)}, {len(ASSOC3_T_K)} T x {len(ASSOC3_P_PA)} P x "
                    f"{len(ASSOC3_FEEDS)} feeds, near their cloud points"
                ),
                components=names,
                model="PC-SAFT (2B association, kij = 0)",
                route="flash_tp (phi-phi)",
                prepare=_pcsaft_prepare(names),
                flash=_pcsaft_flash,
                states=states,
            )
        )
    return built


# ---------------------------------------------------------------------------
# The sweep
# ---------------------------------------------------------------------------

#: How the ``--quick`` subset is drawn from each system:
#: ``(offset, stride, limit, pinned indices)``.
#:
#: **These numbers are cost-bounded, not representative.** The default test
#: suite has a standing wall-time budget (brain.md section 10), so the offsets
#: and strides of the three PC-SAFT families were picked from the *full run's*
#: per-state wall times to land on cheap states; a uniform ``stride`` sample of
#: the same systems costs about four times as much and says nothing more about
#: the instrument. The Peng-Robinson and activity families are a plain stride
#: from index 0, because every state of theirs costs milliseconds.
#:
#: The pinned indices are refusal states, one per refusal *stage* the full
#: sweep found that costs under ~1.5 s: they are what makes
#: ``tests/test_robustness_map.py``'s pinned-defect assertions non-vacuous.
#: ``stability-inconclusive`` has no pinned entry - its two states cost ~4 s
#: each - and is covered by the ``slow`` full-sweep test instead.
QUICK_SAMPLING: Mapping[str, tuple[int, int, int | None, tuple[int, ...]]] = {
    "pcsaft-f4-carbondioxide-n-decane": (5, 5, 3, ()),
    "pcsaft-f4-methane-n-hexane": (3, 8, 4, ()),
    "pcsaft-kij-carbondioxide-n-decane": (5, 1, 2, ()),
    "pcsaft-kij-methane-n-decane": (4, 3, 2, ()),
    "assoc-water-n-hexane": (15, 40, 2, ()),
    "assoc-water-ethanol": (13, 6, 2, ()),
    "assoc-water-1-propanol-n-hexane": (35, 1, 1, ()),
    "assoc-methanol-n-hexane": (19, 40, 2, ()),
    "raoult-tessier-p1": (0, 5, None, ()),
    "raoult-tessier-p2": (0, 5, None, ()),
    "raoult-butanol-water": (0, 7, None, ()),
    "gamma-tessier-p1": (0, 2, None, ()),
    "gamma-tessier-p2": (0, 2, None, ()),
    # index 0 = 1 wt%, 0.3 MPa: split-non-convergence / beta-outside-window.
    "polymer-pe16400-pentane": (39, 40, 2, (0,)),
    # index 0 = 1 wt%, 0.3 MPa: .../log-space.
    # index 106 = 15 wt%, 8.1 MPa: .../phi-phi, the 30-state band.
    "polymer-pe53000-pentane": (39, 73, 2, (0, 106)),
    "polymer-pe16400-pentane-hexane": (4, 1, 1, ()),
    "polymer-pe53000-pentane-hexane": (4, 1, 1, ()),
    # index 40 = z_water=0.05, T3+1 K: multiphase-solver-failure / collapsed
    # (a converged two-phase set with a non-positive phase fraction).
    "eos3p-pcsaft-water-n-hexane-t3-scan": (0, 1, 0, (40,)),
    # index 13 = feed (0.1, 0.8, 0.1), 333 K: cheap single-phase state (~0.4 s);
    # the VLLE vertices themselves cost 13-40 s each (ADR-0020) and are left to
    # the full sweep and the dedicated slow-marked verdict test below.
    "eos3p-pcsaft-water-ethanol-n-hexane": (0, 1, 0, (13,)),
    # Every state of this system shares one pressure and one of two
    # temperatures, so a stride sample risks two feeds landing on the same
    # (T, P) - ambiguous for the (system, T, P) refusal pin below. Pinned
    # only: index 4 = feed (0.2, 0.6, 0.2), 280 K (multiphase-solver-failure /
    # split); index 18 = feed (0.4, 0.4, 0.2), 300 K (same class).
    "eos3p-pr-water-ethanol-n-hexane": (0, 1, 0, (4, 18)),
    "nearcrit-methane-ethane-propane": (0, 10, 4, ()),
    "nearcrit-methane-n-pentane": (0, 10, 4, ()),
    "nearcrit-window-co2-n-decane": (0, 4, 3, ()),
    "nearcrit-window-methane-n-pentane": (0, 3, 3, ()),
    # index 9 = feed (0.15, 0.7, 0.15), 320 K, 0.5 MPa: single-phase, ~0.6 s;
    # the LLE states cost several seconds each and are left to the full sweep.
    "assoc3-water-ethanol-n-hexane": (0, 1, 0, (9,)),
    "assoc3-water-1-propanol-n-hexane": (0, 1, 0, (9,)),
}

#: Every Peng-Robinson system takes this stride from index 0.
PR_QUICK_STRIDE = 11


def _with_quick_sampling(system: RobustnessSystem) -> RobustnessSystem:
    default = (0, PR_QUICK_STRIDE, None, ()) if system.family == "pr-phi-phi" else (0, 1, None, ())
    offset, stride, limit, pinned = QUICK_SAMPLING.get(system.name, default)
    return replace(
        system,
        quick_offset=offset,
        quick_stride=stride,
        quick_limit=limit,
        quick_pinned=pinned,
    )


def systems(family: str | None = None) -> tuple[RobustnessSystem, ...]:
    """Every system, or the systems of one family, in a fixed order."""
    built: list[RobustnessSystem] = []
    built.extend(_pr_systems())
    built.extend(_pcsaft_systems())
    built.extend(_associating_systems())
    built.extend(_modified_raoult_systems())
    built.extend(_gamma_gamma_systems())
    built.extend(_polymer_systems())
    built.extend(_eos_three_phase_systems())
    built.extend(_gamma_phi_legacy_systems())
    built.extend(_nearcrit_systems())
    built.extend(_associating_ternary_systems())
    built = [_with_quick_sampling(system) for system in built]
    if family is None:
        return tuple(built)
    if family not in FAMILIES:
        raise KeyError(f"unknown family {family!r}; known families are {', '.join(FAMILIES)}")
    return tuple(system for system in built if system.family == family)


def _state_record(
    system: RobustnessSystem, index: int, spec: StateSpec, *, quick: bool
) -> dict[str, Any]:
    return {
        "family": system.family,
        "system": system.name,
        "state_index": index,
        "components": list(system.components),
        "composition": [float(value) for value in spec.composition],
        "temperature_K": float(spec.temperature_K),
        "pressure_Pa": float(spec.pressure_Pa),
        "quick": quick,
    }


def run_system(system: RobustnessSystem, *, quick: bool = False) -> list[dict[str, Any]]:
    """Run one system's states and return one classified record each."""
    selected = system.selected(quick=quick)
    try:
        context = system.prepare()
    except SystemSkipped as exc:
        return [
            {
                **_state_record(system, index, spec, quick=quick),
                "outcome": "skipped",
                "bucket": "skipped",
                "reason": str(exc),
                "wall_time_s": 0.0,
            }
            for index, spec in selected
        ]

    records: list[dict[str, Any]] = []
    for index, spec in selected:
        entry = _state_record(system, index, spec, quick=quick)
        started = time.perf_counter()
        try:
            result = system.flash(context, spec)
        except Exception as exc:  # noqa: BLE001 - the classifier is the point
            entry["wall_time_s"] = time.perf_counter() - started
            refusal_class, stage = classify_refusal(exc)
            entry.update(
                {
                    "outcome": "refused",
                    "bucket": refusal_class,
                    "refusal_class": refusal_class,
                    "refusal_stage": stage,
                    "exception": type(exc).__name__,
                    "message": str(exc),
                }
            )
            records.append(entry)
            continue
        entry["wall_time_s"] = time.perf_counter() - started
        violations, residuals = _invariant_violations(result, spec)
        verdict = _verdict(result)
        entry.update(
            {
                "outcome": "converged",
                "bucket": INVARIANT_VIOLATED if violations else verdict,
                "verdict": verdict,
                "phases": result.phase_names(),
                "phase_fractions": [
                    float(result.phase_fractions.get(name, float("nan")))
                    for name in result.phase_names()
                ],
                "residuals": {key: float(value) for key, value in residuals.items()},
                "invariant_violations": violations,
                "converged_stage": str(result.diagnostics.get("converged_stage", "")),
                "phase_set_history": str(result.diagnostics.get("phase_set_history", "")),
            }
        )
        records.append(entry)
    return records


def _summarize(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Roll the per-state records up into counts, worst residuals and times."""
    by_family: dict[str, dict[str, Any]] = {}
    for record in records:
        family = str(record["family"])
        bucket = by_family.setdefault(
            family,
            {
                "states": 0,
                "converged": 0,
                "refused": 0,
                "skipped": 0,
                "verdicts": Counter(),
                "refusal_classes": Counter(),
                "invariant_violations": 0,
                "worst_mass_balance": 0.0,
                "worst_equilibrium": 0.0,
                "worst_delta_g_split_rt": -math.inf,
                "wall_time_s": 0.0,
            },
        )
        bucket["states"] += 1
        bucket["wall_time_s"] += float(record.get("wall_time_s", 0.0))
        outcome = record.get("outcome")
        if outcome == "skipped":
            bucket["skipped"] += 1
            continue
        if outcome == "refused":
            bucket["refused"] += 1
            bucket["refusal_classes"][str(record["refusal_class"])] += 1
            continue
        bucket["converged"] += 1
        bucket["verdicts"][str(record["verdict"])] += 1
        if record.get("invariant_violations"):
            bucket["invariant_violations"] += 1
        residuals = record.get("residuals") or {}
        bucket["worst_mass_balance"] = max(
            bucket["worst_mass_balance"], float(residuals.get("mass_balance", 0.0))
        )
        bucket["worst_equilibrium"] = max(
            bucket["worst_equilibrium"], float(residuals.get("equilibrium", 0.0))
        )
        if "delta_g_split_rt" in residuals:
            bucket["worst_delta_g_split_rt"] = max(
                bucket["worst_delta_g_split_rt"], float(residuals["delta_g_split_rt"])
            )

    families: dict[str, Any] = {}
    for family, bucket in by_family.items():
        worst_delta_g = bucket["worst_delta_g_split_rt"]
        families[family] = {
            **{key: value for key, value in bucket.items() if not isinstance(value, Counter)},
            "verdicts": dict(sorted(bucket["verdicts"].items())),
            "refusal_classes": dict(sorted(bucket["refusal_classes"].items())),
            "worst_delta_g_split_rt": None if worst_delta_g == -math.inf else worst_delta_g,
        }

    totals = {
        "states": len(records),
        "converged": sum(1 for r in records if r.get("outcome") == "converged"),
        "refused": sum(1 for r in records if r.get("outcome") == "refused"),
        "skipped": sum(1 for r in records if r.get("outcome") == "skipped"),
        "invariant_violations": sum(1 for r in records if r.get("invariant_violations")),
        "verdicts": dict(
            sorted(Counter(str(r["verdict"]) for r in records if r.get("verdict")).items())
        ),
        "refusal_classes": dict(
            sorted(
                Counter(str(r["refusal_class"]) for r in records if r.get("refusal_class")).items()
            )
        ),
        "buckets": dict(sorted(Counter(str(r["bucket"]) for r in records).items())),
        "wall_time_s": sum(float(r.get("wall_time_s", 0.0)) for r in records),
    }
    return {"totals": totals, "families": families}


def _refusal_examples(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """One representative state per (class, stage, system), with the message."""
    seen: dict[tuple[str, str, str], dict[str, Any]] = {}
    counts: Counter[tuple[str, str, str]] = Counter()
    for record in records:
        if record.get("outcome") != "refused":
            continue
        key = (
            str(record["refusal_class"]),
            str(record.get("refusal_stage", "")),
            str(record["system"]),
        )
        counts[key] += 1
        seen.setdefault(
            key,
            {
                "refusal_class": key[0],
                "refusal_stage": key[1],
                "system": key[2],
                "family": record["family"],
                "components": record["components"],
                "composition": record["composition"],
                "temperature_K": record["temperature_K"],
                "pressure_Pa": record["pressure_Pa"],
                "exception": record["exception"],
                "message": record["message"],
            },
        )
    examples = []
    for key, entry in seen.items():
        examples.append({**entry, "count": counts[key]})
    examples.sort(key=lambda entry: (-int(entry["count"]), str(entry["refusal_class"])))
    return examples


def run_sweep(
    *, quick: bool = False, family: str | None = None, progress: bool = False
) -> dict[str, Any]:
    """Run the sweep and return a complete record.

    Args:
        quick: Run the ~220-state subset instead of the whole map.
        family: Restrict to one of :data:`FAMILIES`; ``None`` runs all ten.
        progress: Print one line per system as it finishes.

    Returns:
        A JSON-ready record: schema tag, timestamp, git state, environment,
        the grid description, the per-state records, the rolled-up summary,
        one example per refusal class and the slowest states.
    """
    selected = systems(family)
    started = time.perf_counter()
    records: list[dict[str, Any]] = []
    for system in selected:
        system_started = time.perf_counter()
        system_records = run_system(system, quick=quick)
        records.extend(system_records)
        if progress:
            elapsed = time.perf_counter() - system_started
            refused = sum(1 for entry in system_records if entry.get("outcome") == "refused")
            print(
                f"  {system.name:<44} {len(system_records):5d} states "
                f"{refused:4d} refused  {elapsed:8.2f} s",
                flush=True,
            )
    total = time.perf_counter() - started

    worst = sorted(records, key=lambda entry: -float(entry.get("wall_time_s", 0.0)))[:WORST_N]
    summary = _summarize(records)
    return {
        "schema": SCHEMA,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git": git_state(),
        "environment": environment(),
        "quick": quick,
        "quick_note": QUICK_NOTE,
        "family": family,
        "tolerances": {
            "mass_balance": MASS_BALANCE_TOL,
            "equilibrium": EQUILIBRIUM_TOL,
            "composition_sum": COMPOSITION_TOL,
        },
        "systems": [
            {
                "name": system.name,
                "family": system.family,
                "description": system.description,
                "components": list(system.components),
                "model": system.model,
                "route": system.route,
                "states": len(system.states),
                "quick_states": len(system.selected(quick=True)),
            }
            for system in selected
        ],
        "totals": summary["totals"],
        "families": summary["families"],
        "refusal_examples": _refusal_examples(records),
        "slowest_states": [
            {
                "system": entry["system"],
                "temperature_K": entry["temperature_K"],
                "pressure_Pa": entry["pressure_Pa"],
                "composition": entry["composition"],
                "bucket": entry["bucket"],
                "wall_time_s": entry["wall_time_s"],
            }
            for entry in worst
        ],
        "wall_time_s": total,
        "states": records,
    }


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def _format_counter(counts: Mapping[str, int]) -> str:
    if not counts:
        return "-"
    return ", ".join(f"{key} {value}" for key, value in counts.items())


def summary_lines(record: Mapping[str, Any]) -> list[str]:
    """The plain-text table the CLI prints."""
    lines = [
        f"{'family':<22} {'states':>7} {'conv':>6} {'ref':>5} {'inv':>4} "
        f"{'worst |mb|':>11} {'worst eq':>10} {'time / s':>9}"
    ]
    lines.append("-" * 82)
    families: Mapping[str, Any] = record["families"]
    for family in FAMILIES:
        entry = families.get(family)
        if entry is None:
            continue
        lines.append(
            f"{family:<22} {entry['states']:7d} {entry['converged']:6d} "
            f"{entry['refused']:5d} {entry['invariant_violations']:4d} "
            f"{entry['worst_mass_balance']:11.2e} {entry['worst_equilibrium']:10.2e} "
            f"{entry['wall_time_s']:9.2f}"
        )
    totals = record["totals"]
    lines.append("-" * 82)
    lines.append(
        f"{'TOTAL':<22} {totals['states']:7d} {totals['converged']:6d} "
        f"{totals['refused']:5d} {totals['invariant_violations']:4d} "
        f"{'':>11} {'':>10} {record['wall_time_s']:9.2f}"
    )
    lines.append("")
    lines.append("verdicts:        " + _format_counter(totals["verdicts"]))
    lines.append("refusal classes: " + _format_counter(totals["refusal_classes"]))
    return lines


def _state_label(entry: Mapping[str, Any]) -> str:
    """``components``, feed and state, written so a trace mole fraction survives."""
    composition = ", ".join(f"{float(value):.4g}" for value in entry["composition"])
    return (
        f"`{'/'.join(entry['components'])}` z=({composition}) "
        f"T={float(entry['temperature_K']):g} K P={float(entry['pressure_Pa']):.4g} Pa"
    )


def summary_markdown(record: Mapping[str, Any]) -> str:
    """The committed ``benchmarks/robustness_<sha>.md`` summary."""
    git = record.get("git") or {}
    sha = str(git.get("sha") or "unknown")
    env = record.get("environment") or {}
    lines = [
        f"# Robustness map at `{sha[:7]}`",
        "",
        f"- Generated: {record['generated_utc']}",
        f"- Scope: {'quick subset' if record['quick'] else 'full sweep'}"
        + (f", family `{record['family']}`" if record.get("family") else ", all families"),
        f"- Machine: {env.get('cpu')} / Python {env.get('python')} / numpy {env.get('numpy')}",
        f"- Total wall time: {record['wall_time_s']:.1f} s",
        f"- Tree dirty at measurement: {git.get('dirty')}",
        "",
        "What this measures, and what it does not: ADR-0027.",
        "",
        "## Family x outcome",
        "",
        "| family | states | verdicts | refusal classes | invariant violations "
        "| worst mass balance | worst equilibrium residual | worst dG_split/RT | time / s |",
        "| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    families: Mapping[str, Any] = record["families"]
    for family in FAMILIES:
        entry = families.get(family)
        if entry is None:
            continue
        worst_delta_g = entry["worst_delta_g_split_rt"]
        lines.append(
            f"| `{family}` | {entry['states']} | {_format_counter(entry['verdicts'])} "
            f"| {_format_counter(entry['refusal_classes'])} | {entry['invariant_violations']} "
            f"| {entry['worst_mass_balance']:.2e} | {entry['worst_equilibrium']:.2e} "
            f"| {'-' if worst_delta_g is None else f'{worst_delta_g:.2e}'} "
            f"| {entry['wall_time_s']:.1f} |"
        )
    totals = record["totals"]
    lines += [
        "",
        f"**Totals:** {totals['states']} states, {totals['converged']} converged, "
        f"{totals['refused']} refused, {totals['invariant_violations']} converged-but-violating, "
        f"{totals['skipped']} skipped.",
        "",
        "## Refusal classes, ranked by count",
        "",
        "| refusal class | count | families | example state | message (first line) |",
        "| --- | ---: | --- | --- | --- |",
    ]
    by_class: dict[str, dict[str, Any]] = {}
    for entry in record["refusal_examples"]:
        bucket = by_class.setdefault(
            str(entry["refusal_class"]),
            {"count": 0, "families": set(), "example": entry},
        )
        bucket["count"] += int(entry["count"])
        bucket["families"].add(str(entry["family"]))
    for refusal_class, bucket in sorted(by_class.items(), key=lambda kv: -int(kv[1]["count"])):
        example = bucket["example"]
        message = str(example["message"]).split(". ")[0].replace("|", "\\|")
        lines.append(
            f"| `{refusal_class}` | {bucket['count']} | "
            f"{', '.join(sorted(bucket['families']))} | {_state_label(example)} | {message} |"
        )
    lines += [
        "",
        "### ... and by stage, which is what names the defect",
        "",
        "| refusal class | stage | count | system | example state |",
        "| --- | --- | ---: | --- | --- |",
    ]
    for entry in record["refusal_examples"]:
        lines.append(
            f"| `{entry['refusal_class']}` | `{entry['refusal_stage']}` | {entry['count']} "
            f"| `{entry['system']}` | {_state_label(entry)} |"
        )
    lines += [
        "",
        "## Slowest states",
        "",
        "| system | T / K | P / Pa | outcome | time / s |",
        "| --- | ---: | ---: | --- | ---: |",
    ]
    for entry in record["slowest_states"][:10]:
        lines.append(
            f"| `{entry['system']}` | {entry['temperature_K']:g} | {entry['pressure_Pa']:.3g} "
            f"| {entry['bucket']} | {entry['wall_time_s']:.2f} |"
        )
    lines.append("")
    return "\n".join(lines)


def _write_json(path: Path, record: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=False)
        handle.write("\n")


def _iter_state_counts(selected: Iterable[RobustnessSystem]) -> tuple[int, int]:
    full = 0
    quick = 0
    for system in selected:
        full += len(system.states)
        quick += len(system.selected(quick=True))
    return full, quick


def main(argv: Sequence[str] | None = None) -> int:
    """``python -m chemthermo.bench robustness --help``."""
    parser = argparse.ArgumentParser(
        prog="python -m chemthermo.bench robustness",
        description=(
            "Sweep every model family over fixed state and composition grids and write a "
            "classified robustness record (ADR-0027)."
        ),
    )
    parser.add_argument("--out", type=Path, default=None, help="write the JSON record here")
    parser.add_argument(
        "--summary-out",
        type=Path,
        default=None,
        help="write the Markdown summary table here",
    )
    parser.add_argument(
        "--family",
        default=None,
        choices=FAMILIES,
        help="run one family only, so the sweep partitions and resumes",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="run the ~150-state subset across all families instead of the whole map",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list the families and systems with their state counts, and exit",
    )
    parser.add_argument(
        "--quiet", action="store_true", help="do not print a line per system while running"
    )
    args = parser.parse_args(argv)

    if args.list:
        for family in FAMILIES:
            selected = systems(family)
            full, quick = _iter_state_counts(selected)
            print(f"{family:<22} {full:5d} states ({quick} in --quick)")
            for system in selected:
                print(
                    f"    {system.name:<44} {len(system.states):5d}  {system.description}",
                )
        total_full, total_quick = _iter_state_counts(systems())
        print(f"\n{'ALL':<22} {total_full:5d} states ({total_quick} in --quick)")
        return 0

    record = run_sweep(quick=args.quick, family=args.family, progress=not args.quiet)
    print("\n".join(summary_lines(record)))
    if args.out is not None:
        _write_json(args.out, record)
        print(f"\nwrote {args.out}")
    if args.summary_out is not None:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        args.summary_out.write_text(summary_markdown(record), encoding="utf-8")
        print(f"wrote {args.summary_out}")
    return 0
