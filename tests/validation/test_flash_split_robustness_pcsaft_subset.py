"""Validation Case F-4, default-CI subset (slice `flash-phase-labels-by-compressibility`).

`test_flash_split_robustness_pcsaft.py::test_the_whole_grid_answers_and_every_answer_is_verified`
is the full 188-state Case F-4 grid; it is marked `@pytest.mark.slow` and
deselected from the default `pytest -q` run (`pyproject.toml` `addopts`)
because running it there, on top of the same grid running again as the
`examples/validation/15_flash_split_robustness.py` smoke test, took the full
suite from ~116 s to ~389 s. This module is what keeps default `pytest -q`
exercising the PC-SAFT phi-phi path at all: a fixed 16-state subset of the
same grid, including all four states that need the ADR-0016 second-order
stage, checked against the same per-state invariants the full grid checks
(mass balance, equal fugacity, a negative Gibbs-energy change, a stable
post-split phase set). It does not repeat the full grid's exact phase-count
assertions (`two_phase == 123`, `single_phase == 65`), since those are
properties of the whole grid, not of an arbitrary subset. 16 states, not the
24 first considered, to leave headroom in the ~180 s default-run target;
validation Case F-5 (`.agents/brain/validation-cases.md`) records the
before/after timing and why the two other PC-SAFT-adjacent slow examples
(`12_vlle_verdict_map.py`, `14_pcsaft_flash_vs_teqp.py`, ~33 s together) were
left alone - they predate and are unrelated to the regression this slice
fixes.

Run the full grid with `pytest -q -m slow`.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

import chemthermo as ct

#: The four states that need the ADR-0016 second-order stage (verbatim from
#: `test_flash_split_robustness_pcsaft.py`).
PREVIOUSLY_FAILING: tuple[tuple[tuple[str, str], float, float, float], ...] = (
    (("Carbon dioxide", "n-Decane"), 0.8, 250.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 240.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 250.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 260.0, 1.5e6),
)

#: 16 states spread over both binaries, the full feed/temperature/pressure
#: range of the Case F-4 grid, and both single- and two-phase verdicts -
#: always including the four states above. Every tuple is a valid member of
#: the grid in `test_flash_split_robustness_pcsaft.py::GRID`; membership is
#: asserted below rather than trusted.
SUBSET: tuple[tuple[tuple[str, str], float, float, float], ...] = PREVIOUSLY_FAILING + (
    (("Carbon dioxide", "n-Decane"), 0.6, 230.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.6, 260.0, 2.5e6),
    (("Carbon dioxide", "n-Decane"), 0.8, 230.0, 1.5e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 230.0, 2.0e6),
    (("Methane", "n-Hexane"), 0.5, 170.0, 0.5e6),
    (("Methane", "n-Hexane"), 0.5, 200.0, 3.5e6),
    (("Methane", "n-Hexane"), 0.8, 195.0, 2.5e6),
    (("Methane", "n-Hexane"), 0.8, 200.0, 3.0e6),
    (("Methane", "n-Hexane"), 0.9, 190.0, 3.5e6),
    (("Methane", "n-Hexane"), 0.95, 170.0, 1.0e6),
    (("Methane", "n-Hexane"), 0.95, 195.0, 2.5e6),
    (("Methane", "n-Hexane"), 0.95, 200.0, 3.5e6),
)


def _full_grid_states() -> set[tuple[tuple[str, str], float, float, float]]:
    # Loaded by path, not as ``tests.validation....``: whether the repository
    # root is on ``sys.path`` depends on how pytest was invoked
    # (``python -m pytest`` puts it there, the ``pytest`` console script does
    # not), and this assertion should not depend on that.
    module_path = Path(__file__).with_name("test_flash_split_robustness_pcsaft.py")
    spec = importlib.util.spec_from_file_location("_pcsaft_full_grid", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    grid = module.GRID

    return {
        (components, z1, temperature_K, pressure_Pa)
        for components, feeds, temperatures, pressures in grid
        for z1 in feeds
        for temperature_K in temperatures
        for pressure_Pa in pressures
    }


def test_subset_is_well_formed() -> None:
    assert len(SUBSET) == 16, len(SUBSET)
    assert len(set(SUBSET)) == 16, "SUBSET must not contain duplicates"
    assert set(PREVIOUSLY_FAILING) <= set(SUBSET)
    assert set(SUBSET) <= _full_grid_states(), "every SUBSET state must belong to the full grid"


def _flash(components: tuple[str, str], z1: float, temperature_K: float, pressure_Pa: float):
    mixture = ct.Mixture.from_database(list(components), [z1, 1.0 - z1], normalize=True)
    return ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=ct.PCSAFTEOS(),
    )


def test_the_subset_answers_and_every_answer_is_verified() -> None:
    """The 16-state subset: 0 `ConvergenceError`s, every answer verified."""
    failures: list[str] = []
    two_phase = 0
    single_phase = 0
    rescued = 0
    worst_mass_balance = 0.0
    worst_fugacity = 0.0
    worst_delta_g = -np.inf

    for components, z1, temperature_K, pressure_Pa in SUBSET:
        label = f"{'/'.join(components)} z1={z1} T={temperature_K} P={pressure_Pa}"
        try:
            result = _flash(components, z1, temperature_K, pressure_Pa)
        except ct.ConvergenceError as error:
            failures.append(f"{label}: {error}")
            continue

        diagnostics = result.diagnostics
        if diagnostics["phase_count"] == 1:
            single_phase += 1
            assert diagnostics["stability_status"] == "stable", label
            assert diagnostics["termination_reason"] == "feed_stable_tangent_plane", label
            continue

        two_phase += 1
        beta = result.vapor_fraction
        assert beta is not None and 0.0 < beta < 1.0, label
        mass_balance = float(diagnostics["mass_balance_residual"])
        fugacity = float(diagnostics["fugacity_residual"])
        delta_g = float(diagnostics["delta_g_split_rt"])
        assert mass_balance < 1e-12, label
        assert fugacity < 1e-6, label
        assert delta_g < 0.0, label
        assert diagnostics["post_split_status"] == "stable", label
        worst_mass_balance = max(worst_mass_balance, mass_balance)
        worst_fugacity = max(worst_fugacity, fugacity)
        worst_delta_g = max(worst_delta_g, delta_g)
        if diagnostics.get("converged_stage") == "second-order":
            rescued += 1

    assert not failures, "\n".join(failures)
    assert two_phase + single_phase == len(SUBSET), (two_phase, single_phase)
    assert rescued == len(PREVIOUSLY_FAILING), rescued
    assert worst_mass_balance < 1e-12
    assert worst_fugacity < 1e-6
    assert worst_delta_g < 0.0
