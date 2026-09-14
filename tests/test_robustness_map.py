"""The robustness map runs, classifies, and pins today's refusals (ADR-0027).

Three different things are checked here and it is worth keeping them apart.

**The instrument.** The classifier is a pure function of an exception, so it is
tested against the *message strings the solvers actually raise*, with no flash
run at all. That part is fast, deterministic, and the part that would break
first if a message were reworded.

**The map.** The ``--quick`` subset - 224 states, a cost-bounded sample of
every family (see ``QUICK_SAMPLING`` in the module under test) - is run once per
module and its per-family bucket counts are compared against the committed
expectation below. Those counts were measured at the commit named in
``COMMITTED_RECORD``; a change to any of them is a change in what chemthermo
can answer, which is exactly what this file exists to notice.

**The refusals.** Every refusal the quick subset contains is pinned by name,
state and class, with a test that **expects the failure**. If a later slice
fixes one of them, this file fails - deliberately. The fix is to move the state
out of :data:`PINNED_REFUSALS` and into the record, not to loosen the test.
The six original families are still empty since ADR-0028 - the full sweep of
*that* grid refuses nothing (ledger Cases R-MAP-1, P-17) - but slice
``robustness-map-coverage`` (ADR-0027 amendment) adds four families that were
named as coverage gaps in ADR-0027's own roadmap, and two of them refuse:
see ledger Case R-MAP-2.

**The three-phase verdicts.** A handful of states the ledger already pins by
an independent solve (Case P-9: the water/n-hexane three-phase temperature
T3; Case P-10: the water/ethanol/n-hexane tie-triangle at 333 K) are re-tested
here against the same constants the ``eos-three-phase`` family sweeps, so a
regression in the phase-addition/removal search shows up two ways: a moved
bucket count above, and a wrong verdict below.

The full sweep is marked ``slow``: the quick subset runs the same code over
every family by default, so this is the "full grid whose representative
subset runs by default" case of `.agents/dev-contract.md`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

import pytest

import chemthermo as ct
from chemthermo.bench import main as bench_main
from chemthermo.bench.robustness import (
    EOS3P_T3_K,
    EOS3P_TERNARY_FEEDS,
    FAMILIES,
    INVARIANT_VIOLATED,
    REFUSAL_CLASSES,
    classify_refusal,
    run_sweep,
    summary_markdown,
    systems,
)
from chemthermo.eos import PCSAFTEOS

# ---------------------------------------------------------------------------
# The committed expectation for the quick subset (measured at 9adf390).
# ---------------------------------------------------------------------------

#: The committed full sweep this module checks the code against. Regenerating
#: it at a new commit means a new file (``benchmarks/README.md``), so the name
#: lives in one place.
COMMITTED_RECORD = Path(__file__).resolve().parents[1] / "benchmarks" / "robustness_74820b8.json"

#: Total states in the ``--quick`` subset.
EXPECTED_QUICK_STATES = 224

#: Per family: (states, converged, refused, invariant violations).
EXPECTED_QUICK_FAMILIES: dict[str, tuple[int, int, int, int]] = {
    "pr-phi-phi": (116, 116, 0, 0),
    "pcsaft": (11, 11, 0, 0),
    "pcsaft-associating": (7, 7, 0, 0),
    "modified-raoult": (20, 20, 0, 0),
    "gamma-gamma": (8, 8, 0, 0),
    "polymer": (9, 9, 0, 0),
    "eos-three-phase": (4, 1, 3, 0),
    "gamma-phi-legacy": (30, 25, 5, 0),
    "pr-near-critical": (17, 17, 0, 0),
    "pcsaft-associating-ternary": (2, 2, 0, 0),
}

#: Every refusal the quick subset contains, pinned. Each entry is
#: ``(system, T / K, P / Pa, refusal class, refusal stage)``.
#:
#: The six original families are still empty since ADR-0028 - the three
#: polyethylene / n-pentane states pinned here at ``87f0820`` now converge to
#: verified two-phase results (ledger Cases R-MAP-1 and P-17), and nothing in
#: those six families refuses in the quick subset either. Slice
#: ``robustness-map-coverage`` (ADR-0027 amendment) adds four new families and
#: finds real refusals in two of them - see ledger Case R-MAP-2 for the
#: diagnosis of each.
PINNED_REFUSALS: tuple[tuple[str, float, float, str, str], ...] = (
    # eos-three-phase: z_water = 0.05, T3 + 1 K. A converged two-phase set
    # with a non-positive phase fraction, one step past the phase-addition
    # search's usual add-then-remove repair (Case R-MAP-2).
    (
        "eos3p-pcsaft-water-n-hexane-t3-scan",
        335.807826336,
        101325.0,
        "multiphase-solver-failure",
        "collapsed",
    ),
    # eos-three-phase: the Peng-Robinson water/ethanol/n-hexane ternary,
    # feed (0.2, 0.6, 0.2), at both swept temperatures. The multiphase split
    # itself does not converge (Case R-MAP-2), not the search around it.
    (
        "eos3p-pr-water-ethanol-n-hexane",
        280.0,
        101325.0,
        "multiphase-solver-failure",
        "split",
    ),
    (
        "eos3p-pr-water-ethanol-n-hexane",
        300.0,
        101325.0,
        "multiphase-solver-failure",
        "split",
    ),
    # gamma-phi-legacy: the deprecated path's own known failure mode
    # (ADR-0008/ADR-0016 fixed this on the tangent-plane path; the legacy
    # Wilson-heuristic path is deliberately unchanged and still has it).
    ("gammaphi-methane-ethane", 200.0, 3.0e6, "rr-no-bracket", "phi-phi"),
    ("gammaphi-methane-ethane", 220.0, 1.0e6, "rr-no-bracket", "phi-phi"),
    ("gammaphi-methane-ethane", 240.0, 2.0e6, "rr-no-bracket", "phi-phi"),
    ("gammaphi-methane-ethane", 260.0, 4.0e6, "rr-no-bracket", "phi-phi"),
    ("gammaphi-methane-ethane", 260.0, 5.0e6, "rr-no-bracket", "phi-phi"),
)

#: States that converge but violate an invariant. Empty is the claim - the full
#: 2110-state map found none - and a new entry here has to be argued for, not
#: added to make a test pass.
PINNED_INVARIANT_VIOLATIONS: tuple[tuple[str, float, float], ...] = ()


@pytest.fixture(scope="module")
def quick_record() -> Mapping[str, Any]:
    """The quick sweep, run once for the whole module."""
    return run_sweep(quick=True)


# ---------------------------------------------------------------------------
# The instrument: the classifier, with no flash run.
# ---------------------------------------------------------------------------

#: Message text taken verbatim from the ``raise`` sites in
#: ``chemthermo/flash/`` and ``chemthermo/eos/``, so this test fails when a
#: message is reworded rather than silently reclassifying its states.
_MESSAGES: tuple[tuple[str, type[Exception], str, str], ...] = (
    (
        "Tangent-plane stability analysis was inconclusive (no trial converged), so "
        "flash_tp cannot decide whether the feed is one phase or two.",
        ct.ConvergenceError,
        "stability-inconclusive",
        "feed",
    ),
    (
        "A post-split stability test was inconclusive for phase(s) liquid, so flash_tp "
        "cannot decide whether the two-phase set is the answer.",
        ct.ConvergenceError,
        "stability-inconclusive",
        "post-split",
    ),
    (
        "Rachford-Rice failed to bracket a vapor fraction.",
        ct.ConvergenceError,
        "rr-no-bracket",
        "phi-phi",
    ),
    (
        "Feed is unstable (tpd_min=-1.000000e-03) but the stability-seeded K-values do "
        "not bracket a Rachford-Rice root, so no split can be started.",
        ct.ConvergenceError,
        "rr-no-bracket",
        "seeded",
    ),
    (
        "PC-SAFT found no density root: the isotherm never crosses P = 100000.0 Pa.",
        ct.ModelError,
        "density-root-failure",
        "pcsaft-no-root",
    ),
    (
        "No real compressibility roots found for Peng-Robinson EOS.",
        ct.ModelError,
        "density-root-failure",
        "cubic-no-root",
    ),
    (
        "The converged two-phase solution is not a stable phase set: the post-split "
        "stability test reports 'unstable' for phase(s) liquid (most negative "
        "post-split tpd = -1.000000e-03). A third phase is required, and "
        "FlashSettings.max_phases = 2 forbids it.",
        ct.ConvergenceError,
        "post-split-third-phase",
        "max-phases",
    ),
    (
        "Multiphase Rachford-Rice did not converge; max |f_j| = 1.000e-03.",
        ct.ConvergenceError,
        "multiphase-solver-failure",
        "rachford-rice",
    ),
    (
        "The phase addition/removal search did not settle on a stable phase set within "
        "8 rounds; the sets visited were L -> LV -> LLV.",
        ct.ConvergenceError,
        "multiphase-solver-failure",
        "search",
    ),
    (
        "The multiphase split did not converge; equal-fugacity residual=1.836e-08 after "
        "50 successive-substitution and 1 second-order iterations.",
        ct.ConvergenceError,
        "multiphase-solver-failure",
        "split",
    ),
    (
        "A two-phase set converged to a non-positive phase fraction, which would leave "
        "no split at all. This is a solver failure, not a phase count: the tangent-plane "
        "test had already proved the feed unstable.",
        ct.ConvergenceError,
        "multiphase-solver-failure",
        "collapsed",
    ),
    (
        "flash_tp did not converge the phi-phi split in log mole numbers; "
        "equal-fugacity residual=1.000e-04 after 100 log-space Newton iterations "
        "from the stability seed.",
        ct.ConvergenceError,
        "split-non-convergence",
        "log-space",
    ),
    (
        "flash_tp did not converge the phi-phi split; equal-fugacity residual=5.768e-04 "
        "after 100 successive-substitution (max_delta_k=4.093e+106) and 101 "
        "second-order iterations.",
        ct.ConvergenceError,
        "split-non-convergence",
        "phi-phi",
    ),
    (
        "flash_tp did not converge the liquid-liquid split; equal-activity "
        "residual=1.000e-04 after 50 successive-substitution and 100 second-order "
        "iterations.",
        ct.ConvergenceError,
        "split-non-convergence",
        "gamma-gamma",
    ),
    (
        "flash_tp did not converge the modified-Raoult split; equilibrium "
        "residual=1.000e-04 after 50 successive-substitution and 100 second-order "
        "iterations.",
        ct.ConvergenceError,
        "split-non-convergence",
        "modified-raoult",
    ),
    (
        "The phi-phi split converged to a vapor fraction outside (0, 1) "
        "(beta=1.500000e+00), i.e. to a single phase.",
        ct.ConvergenceError,
        "split-non-convergence",
        "beta-outside-window",
    ),
    (
        "Antoine correlation for 'Water' is valid over [288.00, 400.00] K; got 450.0000 K.",
        ct.InputRangeError,
        "model-error",
        "InputRangeError",
    ),
    (
        "Component 'Polyethylene' has no Antoine vapor-pressure record.",
        ct.PropertyNotFoundError,
        "model-error",
        "PropertyNotFoundError",
    ),
)


@pytest.mark.parametrize(
    ("message", "exception_type", "expected_class", "expected_stage"), _MESSAGES
)
def test_the_classifier_maps_each_raised_message_to_its_class(
    message: str, exception_type: type[Exception], expected_class: str, expected_stage: str
) -> None:
    """Every message the solvers raise lands in the class ADR-0027 names for it."""
    assert classify_refusal(exception_type(message)) == (expected_class, expected_stage)


def test_an_unrecognised_refusal_falls_into_other_rather_than_a_wrong_class() -> None:
    """The catch-all is a bucket, not a silent miscount."""
    refusal_class, stage = classify_refusal(RuntimeError("something nobody wrote a rule for"))
    assert refusal_class == "other-refusal"
    assert stage == "RuntimeError"


def test_every_class_the_rules_can_emit_is_declared() -> None:
    """No rule may name a class that :data:`REFUSAL_CLASSES` does not list."""
    from chemthermo.bench.robustness import _REFUSAL_RULES

    emitted = {refusal_class for _needle, refusal_class, _stage in _REFUSAL_RULES}
    assert emitted <= set(REFUSAL_CLASSES)


# ---------------------------------------------------------------------------
# The grid itself.
# ---------------------------------------------------------------------------


def test_every_system_name_is_unique_and_belongs_to_a_declared_family() -> None:
    names = [system.name for system in systems()]
    assert len(names) == len(set(names))
    assert {system.family for system in systems()} == set(FAMILIES)


def test_every_family_is_reachable_on_its_own_so_the_sweep_partitions() -> None:
    """``--family`` partitions the map: the parts sum to the whole, exactly."""
    total = sum(len(systems(family)) for family in FAMILIES)
    assert total == len(systems())
    with pytest.raises(KeyError):
        systems("not-a-family")


def test_every_swept_composition_is_a_normalizable_feed() -> None:
    for system in systems():
        for spec in system.states:
            assert len(spec.composition) == len(system.components), system.name
            assert all(value > 0.0 for value in spec.composition), system.name
            assert spec.temperature_K > 0.0 and spec.pressure_Pa > 0.0, system.name


# ---------------------------------------------------------------------------
# The map: the quick subset against its committed counts.
# ---------------------------------------------------------------------------


def test_the_quick_subset_is_the_size_it_was_committed_at(
    quick_record: Mapping[str, Any],
) -> None:
    assert quick_record["totals"]["states"] == EXPECTED_QUICK_STATES
    assert quick_record["quick"] is True
    assert set(quick_record["families"]) == set(FAMILIES)


@pytest.mark.parametrize("family", FAMILIES)
def test_each_family_matches_its_committed_bucket_counts(
    quick_record: Mapping[str, Any], family: str
) -> None:
    """The counts are the measurement; a move in any of them is a finding."""
    entry = quick_record["families"][family]
    expected = EXPECTED_QUICK_FAMILIES[family]
    assert (
        entry["states"],
        entry["converged"],
        entry["refused"],
        entry["invariant_violations"],
    ) == expected


def test_the_totals_are_the_sum_of_the_families(quick_record: Mapping[str, Any]) -> None:
    totals = quick_record["totals"]
    assert totals["states"] == sum(entry["states"] for entry in quick_record["families"].values())
    assert totals["converged"] + totals["refused"] + totals["skipped"] == totals["states"]


# ---------------------------------------------------------------------------
# The refusals: pinned, expecting today's failure.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("system_name", "temperature_K", "pressure_Pa", "refusal_class", "refusal_stage"),
    PINNED_REFUSALS,
)
def test_each_pinned_refusal_still_refuses_with_the_same_class(
    quick_record: Mapping[str, Any],
    system_name: str,
    temperature_K: float,
    pressure_Pa: float,
    refusal_class: str,
    refusal_stage: str,
) -> None:
    """A pinned defect that gets fixed **fails this test**, which is the point.

    When a slice repairs one of these, delete its entry here and record the
    repair in the ledger; do not relax the assertion.
    """
    matches = [
        state
        for state in quick_record["states"]
        if state["system"] == system_name
        and state["temperature_K"] == pytest.approx(temperature_K)
        and state["pressure_Pa"] == pytest.approx(pressure_Pa)
    ]
    assert len(matches) == 1, (system_name, temperature_K, pressure_Pa)
    state = matches[0]
    assert state["outcome"] == "refused"
    assert (state["refusal_class"], state["refusal_stage"]) == (refusal_class, refusal_stage)


def test_the_quick_subset_contains_exactly_the_pinned_refusals(
    quick_record: Mapping[str, Any],
) -> None:
    """No unpinned refusal, and no pinned refusal that no longer happens."""
    found = {
        (
            state["system"],
            round(float(state["temperature_K"]), 6),
            round(float(state["pressure_Pa"]), 6),
        )
        for state in quick_record["states"]
        if state["outcome"] == "refused"
    }
    pinned = {
        (name, round(float(temperature_K), 6), round(float(pressure_Pa), 6))
        for name, temperature_K, pressure_Pa, _cls, _stage in PINNED_REFUSALS
    }
    assert found == pinned


def test_no_unpinned_invariant_violation_appears(quick_record: Mapping[str, Any]) -> None:
    """A converged answer that breaks mass balance, dG or the post-split test.

    This is the outcome that outranks a refusal, so the expectation is exact
    rather than an upper bound.
    """
    found = {
        (
            state["system"],
            round(float(state["temperature_K"]), 6),
            round(float(state["pressure_Pa"]), 6),
        )
        for state in quick_record["states"]
        if state.get("invariant_violations")
    }
    pinned = {
        (name, round(float(temperature_K), 6), round(float(pressure_Pa), 6))
        for name, temperature_K, pressure_Pa in PINNED_INVARIANT_VIOLATIONS
    }
    assert found == pinned
    assert all(
        state["bucket"] != INVARIANT_VIOLATED
        for state in quick_record["states"]
        if not state.get("invariant_violations")
    )


# ---------------------------------------------------------------------------
# Three-phase verdict expectations (ADR-0027/0028 coverage gap, Cases P-9/P-10)
# ---------------------------------------------------------------------------


#: `tests/test_flash_vlle_eos.py` already flashes exactly these two states by
#: default (`below_t3` / `above_t3`, `OFFSET_K = 0.05`) and is where the
#: independent 4-equation Newton that derived `EOS3P_T3_K` lives. This is a
#: repetition of that same capability through the robustness module's own
#: constant, so it is marked slow rather than paid for twice in the default
#: suite (`.agents/dev-contract.md`).
@pytest.mark.slow
def test_the_water_n_hexane_scan_brackets_t3_with_the_expected_verdicts() -> None:
    """T3 - 0.05 K -> LLE, T3 + 0.05 K -> VLE (Case P-9), from EOS3P_T3_K."""
    eos = PCSAFTEOS()
    mixture = ct.Mixture.from_database(["Water", "n-Hexane"], [0.5, 0.5], normalize=True)

    below = ct.flash_tp(mixture, temperature_K=EOS3P_T3_K - 0.05, pressure_Pa=101325.0, eos=eos)
    assert sorted(below.phases) == ["liquid1", "liquid2"]
    assert below.diagnostics["phase_regime"] == "LLE"
    assert float(below.diagnostics["delta_g_vs_two_phase_rt"]) < 0.0

    above = ct.flash_tp(mixture, temperature_K=EOS3P_T3_K + 0.05, pressure_Pa=101325.0, eos=eos)
    assert sorted(above.phases) == ["liquid", "vapor"]
    assert above.diagnostics["phase_regime"] == "VLE"


#: A PC-SAFT ternary VLLE flash costs 13-40 s (ADR-0020's own measurement is
#: "~35 s"); two feeds only, chosen from the twelve
#: `eos3p-pcsaft-water-ethanol-n-hexane` sweeps at their cheaper end. Marked
#: slow for the same reason `test_flash_vlle_eos.py`'s own ternary VLLE check
#: is: the ``eos-three-phase`` family's quick pin and the (cheap, default-run)
#: Peng-Robinson three-liquid family already exercise the phase
#: addition/removal search machinery by default; this is the PC-SAFT-specific
#: numeric confirmation on top of it.
@pytest.mark.slow
@pytest.mark.parametrize("feed", [EOS3P_TERNARY_FEEDS[3], EOS3P_TERNARY_FEEDS[9]])
def test_the_ternary_feeds_are_a_vlle_tie_triangle_at_333k(
    feed: tuple[float, float, float],
) -> None:
    """Every feed inside the 333 K tie-triangle returns VLLE (Case P-10)."""
    eos = PCSAFTEOS()
    mixture = ct.Mixture.from_database(["Water", "Ethanol", "n-Hexane"], list(feed), normalize=True)
    result = ct.flash_tp(mixture, temperature_K=333.0, pressure_Pa=101325.0, eos=eos)
    assert sorted(result.phases) == ["liquid1", "liquid2", "vapor"]
    assert result.diagnostics["phase_regime"] == "VLLE"
    assert float(result.diagnostics["delta_g_vs_two_phase_rt"]) < 0.0


# ---------------------------------------------------------------------------
# The record and the CLI.
# ---------------------------------------------------------------------------


def test_the_record_is_json_and_carries_its_provenance(
    quick_record: Mapping[str, Any],
) -> None:
    assert quick_record["schema"] == "chemthermo-robustness/1"
    assert quick_record["environment"]["python"]
    assert "sha" in quick_record["git"]
    assert quick_record["tolerances"]["mass_balance"] == 1e-10
    json.dumps(quick_record)


def test_the_markdown_summary_names_every_family(quick_record: Mapping[str, Any]) -> None:
    markdown = summary_markdown(quick_record)
    for family in FAMILIES:
        assert f"`{family}`" in markdown
    assert "ADR-0027" in markdown


def test_the_cli_runs_one_family_and_writes_both_artefacts(tmp_path: Path) -> None:
    """The subcommand, end to end, on the cheapest family."""
    out = tmp_path / "robustness.json"
    summary = tmp_path / "robustness.md"

    status = bench_main(
        [
            "robustness",
            "--quick",
            "--quiet",
            "--family",
            "gamma-gamma",
            "--out",
            str(out),
            "--summary-out",
            str(summary),
        ]
    )

    assert status == 0
    record = json.loads(out.read_text(encoding="utf-8"))
    assert record["family"] == "gamma-gamma"
    assert record["totals"]["states"] > 0
    assert summary.read_text(encoding="utf-8").startswith("# Robustness map at")


def test_the_timing_cli_still_means_what_it_meant(tmp_path: Path) -> None:
    """Adding a subcommand may not move the ADR-0023 flags (ADR-0027 decision 1)."""
    out = tmp_path / "record.json"
    assert bench_main(["--repeats", "1", "--case", "pr-stability-ternary", "--out", str(out)]) == 0
    record = json.loads(out.read_text(encoding="utf-8"))
    assert record["schema"] == "chemthermo-bench/1"


def test_the_committed_record_matches_the_grid_this_module_defines() -> None:
    """The committed full record describes the same 2110-state map as the code."""
    path = COMMITTED_RECORD
    if not path.is_file():  # pragma: no cover - a wheel has no benchmarks/ directory
        pytest.skip("the committed record lives in a source checkout only")
    with path.open("r", encoding="utf-8") as handle:
        record = json.load(handle)
    assert record["schema"] == "chemthermo-robustness/1"
    assert record["quick"] is False
    by_name = {system.name: len(system.states) for system in systems()}
    assert {entry["name"]: entry["states"] for entry in record["systems"]} == by_name


# ---------------------------------------------------------------------------
# The full sweep.
# ---------------------------------------------------------------------------


@pytest.mark.slow  # minutes; the quick subset above runs the same code over every family
def test_the_full_sweep_reproduces_the_committed_totals() -> None:
    """The whole 2110-state map, against the committed record's counts."""
    with COMMITTED_RECORD.open("r", encoding="utf-8") as handle:
        committed = json.load(handle)

    record = run_sweep()

    assert record["totals"]["states"] == committed["totals"]["states"]
    assert record["totals"]["verdicts"] == committed["totals"]["verdicts"]
    assert record["totals"]["refusal_classes"] == committed["totals"]["refusal_classes"]
    assert record["totals"]["invariant_violations"] == committed["totals"]["invariant_violations"]
