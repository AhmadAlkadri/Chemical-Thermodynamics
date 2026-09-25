"""Fixed density-root surfaces in the equation-of-state trials (ADR-0021, Case P-11).

ADR-0012 pinned each *modified-Raoult* trial to one phase candidate and
deliberately left the equation-of-state family re-selecting the lowest-Gibbs
root at every iterate. ADR-0021 finishes the job, because the failure ADR-0012
described turned out to have an exact analogue here: from the hexane-rich
liquid of a water / n-hexane feed above the three-phase temperature, every
trial slid onto the liquid root and stopped at the trivial solution or at the
partner liquid, while a vapour stationary point with ``tpd = -6.5e-03`` sat
unvisited, and ``flash_tp`` returned two liquids that are metastable by
2.5e-03 RT (validation Case P-9 (iv)).

What is under test here:

1. the trial set - which root each start is pinned to, and the two cases that
   are deliberately *not* pinned;
2. the repaired search itself, at the state Case P-9 (iv) pinned;
3. the fallback: a density root that does not exist is not an error, because
   both ``phase`` labels name the only root the model has there, and the trial
   records that it had no choice of surface;
4. what did **not** change: the reported distance is still the lowest-Gibbs one
   and the reported branch is still the lowest-Gibbs candidate, so
   ``chemthermo.flash._detect``'s pinning semantics are untouched.

The whole-result version of (4) is ``tests/test_flash_refactor_bit_identity.py``
(155 states against ``refactor_bit_identity_v3.json``), and the flash-level
consequence is ``tests/test_flash_vlle_eos.py``.
"""

from __future__ import annotations

import numpy as np
import pytest
from _capture_identity import on_capture_platform

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.stability._evaluator import _EOSTangentPlane

PRESSURE_PA = 101325.0
BINARY = ("Water", "n-Hexane")

#: Validation Case P-9 (iv): 335 K, and the hexane-rich vertex of the
#: two-liquid pair that `flash_tp` used to return at `z_water = 0.7` there.
LEDGER_T_K = 335.0
LEDGER_HEXANE_RICH = (0.02273, 0.97727)

#: A Peng-Robinson state whose liquid-surface trials stop where the cubic has a
#: single real root, so the fallback of ADR-0021 is exercised.
SINGLE_ROOT_STATE = (("Methane", "Ethane"), (0.5, 0.5), 240.0, 3.0e6)

#: The 144-state Peng-Robinson grid of `tests/test_flash_phase_detection.py`.
GRID_MIXTURES: tuple[tuple[tuple[str, ...], tuple[float, ...]], ...] = (
    (("Methane", "Ethane"), (0.5, 0.5)),
    (("Methane", "Propane"), (0.7, 0.3)),
    (("Ethane", "n-Heptane"), (0.7, 0.3)),
    (("Methane", "n-Pentane"), (0.6, 0.4)),
    (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2)),
    (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3)),
)
GRID_T_K = (170.0, 200.0, 240.0, 280.0, 320.0, 360.0)
GRID_P_PA = (2.0e5, 1.0e6, 3.0e6, 8.0e6)

#: Two surfaces whose best `tpd` differ by no more than this reached the same
#: stationary point (ADR-0032): 5 orders above the largest tie measured on the
#: grid (5e-16), 8 below the smallest real gap (1.7e-02), and 100x below
#: `StabilitySettings.tpd_tol`, so it cannot hide a difference the verdict
#: logic would see.
TIE_MARGIN = 1e-10
#: Per-surface count of the grid states whose minimizing surface is decisive,
#: and the ties (ledger Case P-11, "cross-platform").
DECISIVE_SURFACES = {"vapor": 19, "liquid": 38, "tie": 20}


def _mixture(names, z) -> ct.Mixture:
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


# ---------------------------------------------------------------------------
# 1. The trial set
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "eos, names, z, temperature_K, pressure_Pa",
    [
        (ct.PengRobinsonEOS(), ("Methane", "Ethane"), (0.5, 0.5), 240.0, 3.0e6),
        (
            ct.PengRobinsonEOS(),
            ("Methane", "Ethane", "Propane"),
            (0.5, 0.3, 0.2),
            240.0,
            1.0e6,
        ),
        (PCSAFTEOS(), BINARY, (0.7, 0.3), LEDGER_T_K, PRESSURE_PA),
    ],
    ids=["pr-binary", "pr-ternary", "pcsaft-binary"],
)
def test_every_multicomponent_eos_trial_names_the_root_its_start_estimates(
    eos, names, z, temperature_K: float, pressure_Pa: float
) -> None:
    """`wilson-vapor` -> vapour root, `wilson-liquid` and `pure-<name>` -> liquid root.

    One vapour-root trial and ``n + 1`` liquid-root ones, for ``n`` components:
    the count is what it always was, only the surfaces are new.
    """
    result = ct.stability_tp(
        _mixture(names, z), temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=eos
    )
    surfaces = {trial.label: trial.surface for trial in result.trials}
    assert surfaces["wilson-vapor"] == "vapor"
    assert surfaces["wilson-liquid"] == "liquid"
    for name in names:
        assert surfaces[f"pure-{name}"] == "liquid"
    assert len(result.trials) == len(names) + 2
    assert result.diagnostics["trial_surfaces"] == f"vapor:1,liquid:{len(names) + 1}"


def test_the_pure_component_starts_are_not_repeated_on_the_vapor_root() -> None:
    """Measured, not assumed: they find nothing the four shipped trials do not.

    Running the ``n`` pure-component-dominant estimates on the vapour root as
    well was tried over the 144-state Peng-Robinson grid, the 188-state PC-SAFT
    grid of validation Case F-4 and the two water / n-hexane states of Case
    P-9 (iv) - **334 states, zero verdict changes and not one state where
    ``tpd_min`` fell by more than 2.6e-13**, which is the same stationary point
    reached by a different trial, not a new one. So they are not run, and this
    test pins the trial count that decision produces (validation Case P-11;
    ADR-0012 decision 3 recorded the same measurement for the ideal-gas
    surface, for a different reason).
    """
    result = ct.stability_tp(
        _mixture(("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2)),
        temperature_K=240.0,
        pressure_Pa=1.0e6,
        eos=ct.PengRobinsonEOS(),
    )
    vapor_trials = [trial for trial in result.trials if trial.surface == "vapor"]
    assert [trial.label for trial in vapor_trials] == ["wilson-vapor"]


def test_a_single_active_component_feed_keeps_min_gibbs_selection() -> None:
    """The one case ADR-0021 leaves unpinned, for ADR-0012's reason.

    With one active component there is no composition degree of freedom: both
    Wilson estimates *are* the feed, and pinning them would test a root the
    feed is not on. The degenerate trials therefore name no surface, and the
    per-surface diagnostics stay absent.
    """
    mixture = ct.Mixture.from_database(["Methane", "Ethane"], [1.0, 0.0], normalize=False)
    result = ct.stability_tp(
        mixture, temperature_K=170.0, pressure_Pa=2.0e5, eos=ct.PengRobinsonEOS()
    )
    assert [(trial.label, trial.surface) for trial in result.trials] == [
        ("wilson-vapor", None),
        ("wilson-liquid", None),
    ]
    assert "trial_surfaces" not in result.diagnostics
    assert "surface_fallback_trial_count" not in result.diagnostics


# ---------------------------------------------------------------------------
# 2. The repaired search, at the state Case P-9 (iv) pinned
# ---------------------------------------------------------------------------


def test_the_hexane_rich_liquid_above_t3_now_reaches_the_vapor_stationary_point() -> None:
    """Case P-9 (iv), at the stability layer: the miss, and where it went.

    Before ADR-0021 every one of these four trials came back trivial or on the
    partner liquid, `tpd_min = -3.0e-09`, verdict `"stable"`. The numbers below
    are the repaired ones; the stationary point they find is the vapour the
    ledger recorded as unvisited, and it is found by the trial that is now
    pinned to the vapour root.
    """
    result = ct.stability_tp(
        _mixture(BINARY, LEDGER_HEXANE_RICH),
        temperature_K=LEDGER_T_K,
        pressure_Pa=PRESSURE_PA,
        eos=PCSAFTEOS(),
    )
    assert result.status == "unstable"
    assert result.tpd_min == pytest.approx(-6.5237494612e-03, abs=1e-12)
    assert result.feed_branch == "liquid"
    assert result.phase_branch == "vapor"
    assert result.trial_composition is not None
    assert result.trial_composition[0] == pytest.approx(0.2135593273947, abs=1e-12)
    assert result.diagnostics["minimizing_trial"] == "wilson-vapor"
    assert result.diagnostics["minimizing_trial_surface"] == "vapor"

    # The other three trials are exactly the ones that used to be the whole
    # answer: two trivial, one on a shallow liquid-surface stationary point.
    others = {trial.label: trial for trial in result.trials if trial.label != "wilson-vapor"}
    assert others["wilson-liquid"].trivial is True
    assert others["pure-n-Hexane"].trivial is True
    assert others["pure-Water"].tpd > 0.0


@pytest.mark.slow
def test_the_vapor_liquid_pair_it_leads_to_is_the_lower_gibbs_answer() -> None:
    """The flash-level consequence, with the arithmetic that decides it.

    `slow` only because the flash at this state runs the phase-addition search:
    the default run covers the same claim through
    `examples/validation/19_eos_stability_surfaces.py`, which `tests/test_examples.py`
    executes, and through the verdict scans of `tests/test_flash_vlle_eos.py`.
    """
    eos = PCSAFTEOS()
    mixture = _mixture(BINARY, (0.5, 0.5))
    result = ct.flash_tp(
        _mixture(BINARY, (0.7, 0.3)),
        temperature_K=LEDGER_T_K,
        pressure_Pa=PRESSURE_PA,
        eos=eos,
    )
    assert sorted(result.phases) == ["liquid", "vapor"]

    def reduced_g(composition, branch: str) -> float:
        x = np.asarray(composition, dtype=float)
        x = x / float(np.sum(x))
        phi = np.asarray(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=LEDGER_T_K,
                pressure_Pa=PRESSURE_PA,
                composition=x.tolist(),
                phase=branch,
            ),
            dtype=float,
        )
        return float(np.sum(x * (np.log(x) + np.log(phi))))

    g_returned = sum(
        result.phase_fractions[name]
        * reduced_g(
            result.phases[name].composition.fractions,
            "vapor" if name == "vapor" else "liquid",
        )
        for name in result.phases
    )
    # The two-liquid pair that used to be returned, from a 2-equation Newton on
    # the liquid branch (validation Case P-9 (iv)'s own numbers).
    water_rich = (0.99993557, 1.0 - 0.99993557)
    hexane_rich = (0.02273239, 1.0 - 0.02273239)
    beta = (0.7 - hexane_rich[0]) / (water_rich[0] - hexane_rich[0])
    g_liquids = beta * reduced_g(water_rich, "liquid") + (1.0 - beta) * reduced_g(
        hexane_rich, "liquid"
    )
    assert g_returned == pytest.approx(-1.1643059308, abs=1e-9)
    assert g_liquids == pytest.approx(-1.1618107137, abs=1e-8)
    assert g_liquids - g_returned == pytest.approx(2.4952e-03, abs=1e-6)


# ---------------------------------------------------------------------------
# 3. The fallback
# ---------------------------------------------------------------------------


def test_the_fallback_is_recorded_where_the_model_has_a_single_root() -> None:
    """A density root that does not exist is not an error, and is not silent.

    ``PengRobinsonEOS`` and ``PCSAFTEOS`` both answer ``phase="vapor"`` and
    ``phase="liquid"`` with the *same* numbers where they have one admissible
    root, so a trial pinned to either label simply walks the only surface there
    - no exception, no wrong phase. That is the documented difference from the
    modified-Raoult pair of ADR-0012, whose two candidates exist at every
    composition, and it is reported rather than hidden.
    """
    names, z, temperature_K, pressure_Pa = SINGLE_ROOT_STATE
    eos = ct.PengRobinsonEOS()
    mixture = _mixture(names, z)
    result = ct.stability_tp(mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=eos)

    fell_back = [trial for trial in result.trials if trial.surface_fallback]
    assert [trial.label for trial in fell_back] == [
        "wilson-liquid",
        "pure-Methane",
        "pure-Ethane",
    ]
    assert result.diagnostics["surface_fallback_trial_count"] == 3
    assert result.diagnostics["surface_fallback_evaluation_count"] == 3

    def roots(composition) -> tuple[float, float]:
        values = [
            eos.compressibility_factor(
                mixture=mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                composition=list(composition),
                phase=phase,
            )
            for phase in ("vapor", "liquid")
        ]
        return values[0], values[1]

    for trial in result.trials:
        assert trial.composition is not None
        vapor_z, liquid_z = roots(trial.composition)
        # The flag says exactly this and nothing else: one root, or two.
        assert trial.surface_fallback == (vapor_z == liquid_z)
        assert trial.surface_fallback == (trial.surface_fallback_count > 0)


def test_a_trial_that_did_not_fall_back_reports_zero() -> None:
    names, z, temperature_K, pressure_Pa = SINGLE_ROOT_STATE
    result = ct.stability_tp(
        _mixture(names, z),
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=ct.PengRobinsonEOS(),
    )
    vapor_trial = next(trial for trial in result.trials if trial.label == "wilson-vapor")
    assert vapor_trial.surface_fallback is False
    assert vapor_trial.surface_fallback_count == 0


def test_an_unknown_root_surface_is_a_model_error() -> None:
    evaluator = _EOSTangentPlane(
        ct.PengRobinsonEOS(),
        mixture=_mixture(("Methane", "Ethane"), (0.5, 0.5)),
        temperature=240.0,
        pressure=3.0e6,
    )
    with pytest.raises(ct.ModelError, match="Unknown phase-candidate surface"):
        evaluator.ln_terms_on_surface(np.array([0.5, 0.5]), "solid")
    with pytest.raises(ct.ModelError, match="Unknown phase-candidate surface"):
        evaluator.ln_report_terms(np.array([0.5, 0.5]), "solid")


# ---------------------------------------------------------------------------
# 4. What did not change
# ---------------------------------------------------------------------------


def test_the_reported_distance_and_branch_are_the_lowest_gibbs_ones() -> None:
    """ADR-0012 decision 4, now for density roots: `tpd` is off the lower envelope.

    A trial pinned to the vapour root must not report a vapour distance where
    the liquid root lies lower: the tangent-plane distance is the distance to
    the *minimum* over the roots. `surface` says what the trial walked,
    `phase_branch` says which root was lowest where it stopped, and
    `chemthermo.flash._detect` reads the second - which is why the incipient
    phase's pinned branch is unchanged by this slice.
    """
    eos = PCSAFTEOS()
    mixture = _mixture(BINARY, LEDGER_HEXANE_RICH)
    result = ct.stability_tp(mixture, temperature_K=LEDGER_T_K, pressure_Pa=PRESSURE_PA, eos=eos)

    def ln_phi(composition, branch: str) -> np.ndarray:
        return np.log(
            np.asarray(
                eos.fugacity_coefficients(
                    mixture=mixture,
                    temperature_K=LEDGER_T_K,
                    pressure_Pa=PRESSURE_PA,
                    composition=list(composition),
                    phase=branch,
                ),
                dtype=float,
            )
        )

    z = np.asarray(result.feed_composition, dtype=float)
    feed_terms = min(
        (ln_phi(z, branch) for branch in ("vapor", "liquid")),
        key=lambda terms: float(np.sum(z * terms)),
    )
    d = np.log(z) + feed_terms

    for trial in result.trials:
        assert trial.composition is not None
        w = np.asarray(trial.composition, dtype=float)
        terms = min(
            (ln_phi(w, branch) for branch in ("vapor", "liquid")),
            key=lambda values: float(np.sum(w * values)),
        )
        assert trial.tpd == pytest.approx(float(np.sum(w * (np.log(w) + terms - d))), abs=1e-12)


def test_the_peng_robinson_grid_still_reaches_a_verdict_everywhere() -> None:
    """144 states, no `"inconclusive"`, and both root surfaces find minimizers.

    The verdicts and every downstream number of this grid are pinned
    bit-for-bit by `tests/test_flash_refactor_bit_identity.py`; what this adds
    is the statement about the *trials*, which that fixture does not carry:
    neither surface is decorative. Measured over the grid, the vapour root
    supplies the minimizing trial on 32 states and the liquid root on 45
    (validation Case P-11).

    On 20 of those 77 states both surfaces reach the *same* stationary point
    (the two best trials' `tpd` agree to 5e-16 and their compositions to
    6e-13), so which surface "wins" is decided by the last bit of the
    arithmetic and differs between machines (45/32 on macOS arm64, 46/31 and
    47/30 on two Linux x86_64 hosts). The platform-independent statement,
    asserted everywhere, counts a state for a surface only when that surface's
    best `tpd` beats the other's by more than `TIE_MARGIN`, and counts the
    rest as ties; the real gaps on this grid are 1.7e-02 and 1.0, the ties at
    most 5e-16. That gives 19 decisive vapour, 38 decisive liquid, 20 ties:
    still neither surface decorative. The exact 45/32 split stays pinned on
    the capture platform (ADR-0032).
    """
    eos = ct.PengRobinsonEOS()
    surfaces: dict[str, int] = {}
    decisive: dict[str, int] = {}
    statuses: dict[str, int] = {}
    for names, z in GRID_MIXTURES:
        for temperature_K in GRID_T_K:
            for pressure_Pa in GRID_P_PA:
                result = ct.stability_tp(
                    _mixture(names, z),
                    temperature_K=temperature_K,
                    pressure_Pa=pressure_Pa,
                    eos=eos,
                )
                statuses[result.status] = statuses.get(result.status, 0) + 1
                surface = result.diagnostics.get("minimizing_trial_surface")
                if surface is not None:
                    surfaces[str(surface)] = surfaces.get(str(surface), 0) + 1
                best_by_surface: dict[str, float] = {}
                for trial in result.trials:
                    if (
                        trial.converged
                        and not trial.trivial
                        and np.isfinite(trial.tpd)
                        and trial.surface is not None
                    ):
                        previous = best_by_surface.get(trial.surface)
                        if previous is None or trial.tpd < previous:
                            best_by_surface[trial.surface] = trial.tpd
                if not best_by_surface:
                    continue
                ranked = sorted(best_by_surface.items(), key=lambda item: item[1])
                if len(ranked) > 1 and ranked[1][1] - ranked[0][1] <= TIE_MARGIN:
                    decisive["tie"] = decisive.get("tie", 0) + 1
                else:
                    # A decisive winner must be the surface the result reports.
                    assert surface == ranked[0][0], (names, temperature_K, pressure_Pa)
                    decisive[ranked[0][0]] = decisive.get(ranked[0][0], 0) + 1
    assert sum(statuses.values()) == 144
    assert "inconclusive" not in statuses
    assert statuses == {"unstable": 47, "stable": 97}
    assert sum(surfaces.values()) == 77
    assert decisive == DECISIVE_SURFACES
    if on_capture_platform():
        assert surfaces == {"vapor": 32, "liquid": 45}
