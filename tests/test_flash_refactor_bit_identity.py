"""Bit-identity regression for the ``flash-module-split`` refactor.

This is not a behavior test in the usual sense: it exists purely to prove that
splitting ``src/chemthermo/flash/tp.py`` into single-purpose internal modules
moved code without changing a single floating-point operation.

Every field a caller of ``flash_tp`` can observe - phase names, every phase's
composition, phase fractions, ``vapor_fraction`` and the full ``diagnostics``
mapping - is captured for a fixed set of states and pinned, bit-for-bit
(floats compared with ``==``, ints/bools/strings exact), against
``tests/fixtures/flash/refactor_bit_identity_v1.json``. That fixture was
generated from this repository at HEAD ``e927623`` (the commit immediately
before the module split) and is committed in its own commit, before the
refactor commit, so the git history itself proves the fixture predates the
code move: see the ``test(flash): capture pre-refactor bit-identity fixture``
commit versus ``refactor(flash): split tp.py into single-purpose internal
modules``.

States captured (see :func:`_build_states`):

- The phi-phi grid of ``tests/test_flash_phase_detection.py``
  (``GRID_MIXTURES x GRID_T_K x GRID_P_PA``, default tangent-plane settings).
- The gamma-gamma n-Butanol/Water binary feeds and the single-component feed
  of ``tests/test_flash_lle.py``, plus the Tessier (2000) Problem 1
  near-plait feed exercised there via ``tests/conftest.py``'s
  ``tessier2000_*`` fixtures.
- The two gamma-phi cases of ``tests/test_flash_gamma_phi.py``.
- The two ``phase_detection="wilson-heuristic"`` cases pinned in
  ``tests/test_flash_phase_detection.py::test_legacy_path_reproduces_the_pre_slice_numbers_exactly``.

If this fixture ever legitimately needs to change (a deliberate, reviewed
numerical change - never to make a refactor "pass"), regenerate it with a
throwaway script that imports :func:`_build_states` and :func:`_encode` from
this module, builds the Tessier ``(names, model)`` pair the way
``tests/conftest.py`` does, runs every state, and writes the JSON from the
commit *before* the change.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Sequence

import chemthermo as ct

FIXTURE_PATH = (
    Path(__file__).resolve().parent / "fixtures" / "flash" / "refactor_bit_identity_v1.json"
)

EOS = ct.PengRobinsonEOS()
LEGACY = ct.FlashSettings(phase_detection="wilson-heuristic")

# Mirrors tests/test_flash_phase_detection.py's GRID_MIXTURES / GRID_T_K /
# GRID_P_PA. Duplicated rather than imported so this module stays
# self-contained; keep the two grids in sync if the phase-detection grid
# ever changes.
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

# n-Butanol / Water NRTL parameters, tests/test_flash_lle.py.
LLE_NAMES = ("n-Butanol", "Water")
LLE_TAU_12 = 0.90047
LLE_TAU_21 = 3.51307
LLE_ALPHA = 0.48
LLE_TEMPERATURE_K = 298.15
LLE_PRESSURE_PA = 101325.0
# Union of the feeds named across tests/test_flash_lle.py (stable feed,
# two-phase feed, permutation-invariance feed, and the binodal-loop feeds).
LLE_BINARY_FEEDS = (0.05, 0.10, 0.20, 0.30, 0.45)


def _mixture(names: Sequence[str], z: Sequence[float]) -> ct.Mixture:
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


def _lle_model() -> ct.NRTL:
    return ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [(LLE_NAMES[0], LLE_NAMES[1], LLE_TAU_12, LLE_TAU_21, LLE_ALPHA, LLE_ALPHA)]
        )
    )


def _encode(result: ct.FlashResult) -> dict[str, Any]:
    """Every user-visible field of a `FlashResult`, JSON-ready."""
    return {
        "phase_names": result.phase_names(),
        "phases": {
            name: list(phase.composition.fractions) for name, phase in result.phases.items()
        },
        "phase_fractions": dict(result.phase_fractions),
        "vapor_fraction": result.vapor_fraction,
        "diagnostics": dict(result.diagnostics),
    }


def _build_states(
    tessier_names: list[str], tessier_model: ct.NRTL
) -> list[tuple[str, Callable[[], ct.FlashResult]]]:
    """Every ``(label, thunk)`` state this regression pins, in a fixed order."""
    states: list[tuple[str, Callable[[], ct.FlashResult]]] = []

    # -- phi-phi grid (tangent-plane, default settings) ---------------------
    for names, z in GRID_MIXTURES:
        for temperature_K in GRID_T_K:
            for pressure_Pa in GRID_P_PA:
                label = (
                    f"phi-phi-grid|{'-'.join(names)}|z={z!r}|T={temperature_K!r}|P={pressure_Pa!r}"
                )

                def thunk(
                    names: tuple[str, ...] = names,
                    z: tuple[float, ...] = z,
                    temperature_K: float = temperature_K,
                    pressure_Pa: float = pressure_Pa,
                ) -> ct.FlashResult:
                    return ct.flash_tp(
                        _mixture(names, z),
                        temperature_K=temperature_K,
                        pressure_Pa=pressure_Pa,
                        eos=EOS,
                    )

                states.append((label, thunk))

    # -- gamma-gamma: n-Butanol/Water binary feeds ---------------------------
    for x1 in LLE_BINARY_FEEDS:
        label = f"gamma-gamma-binary|z1={x1!r}"

        def binary_thunk(x1: float = x1) -> ct.FlashResult:
            return ct.flash_tp(
                _mixture(LLE_NAMES, (x1, 1.0 - x1)),
                temperature_K=LLE_TEMPERATURE_K,
                pressure_Pa=LLE_PRESSURE_PA,
                activity_model=_lle_model(),
            )

        states.append((label, binary_thunk))

    # -- gamma-gamma: single-component feed ----------------------------------
    def water_thunk() -> ct.FlashResult:
        return ct.flash_tp(
            _mixture(("Water",), (1.0,)),
            temperature_K=LLE_TEMPERATURE_K,
            pressure_Pa=LLE_PRESSURE_PA,
            activity_model=_lle_model(),
        )

    states.append(("gamma-gamma-single-component-water", water_thunk))

    # -- gamma-gamma: Tessier (2000) Problem 1 near-plait feed ---------------
    def tessier_thunk() -> ct.FlashResult:
        return ct.flash_tp(
            _mixture(tessier_names, (0.148, 0.052, 0.80)),
            temperature_K=LLE_TEMPERATURE_K,
            pressure_Pa=LLE_PRESSURE_PA,
            activity_model=tessier_model,
        )

    states.append(("gamma-gamma-tessier2000-near-plait", tessier_thunk))

    # -- gamma-phi (tests/test_flash_gamma_phi.py) ---------------------------
    def gamma_phi_two_phase() -> ct.FlashResult:
        components = tuple(ct.Component.from_database(name) for name in ("Methane", "Ethane"))
        mixture = ct.Mixture(
            components=components,
            composition=ct.Composition(fractions=(0.5, 0.5), basis="mole", normalize=False),
        )
        return ct.flash_tp(
            mixture,
            temperature_K=240.0,
            pressure_Pa=3.0e6,
            eos=EOS,
            activity_model=ct.NRTL(),
            flash_mode="gamma-phi",
        )

    states.append(("gamma-phi-two-phase", gamma_phi_two_phase))

    def gamma_phi_single_phase() -> ct.FlashResult:
        component = ct.Component.from_database("Methane")
        mixture = ct.Mixture(
            components=(component,),
            composition=ct.Composition(fractions=(1.0,), basis="mole", normalize=False),
        )
        return ct.flash_tp(
            mixture,
            temperature_K=350.0,
            pressure_Pa=101325.0,
            eos=EOS,
            activity_model=ct.NRTL(),
            flash_mode="gamma-phi",
        )

    states.append(("gamma-phi-single-phase", gamma_phi_single_phase))

    # -- legacy wilson-heuristic: both pinned cases --------------------------
    def legacy_binary() -> ct.FlashResult:
        return ct.flash_tp(
            _mixture(("Methane", "Ethane"), (0.5, 0.5)),
            temperature_K=240.0,
            pressure_Pa=3.0e6,
            eos=EOS,
            settings=LEGACY,
        )

    states.append(("legacy-wilson-heuristic-binary", legacy_binary))

    def legacy_ternary() -> ct.FlashResult:
        return ct.flash_tp(
            _mixture(("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2)),
            temperature_K=240.0,
            pressure_Pa=3.0e6,
            eos=EOS,
            settings=LEGACY,
        )

    states.append(("legacy-wilson-heuristic-ternary", legacy_ternary))

    return states


def test_flash_tp_is_bit_identical_to_the_pre_refactor_capture(
    tessier2000_names: list[str], tessier2000_model: ct.NRTL
) -> None:
    with FIXTURE_PATH.open("r", encoding="utf-8") as handle:
        fixture: dict[str, dict[str, Any]] = json.load(handle)

    states = _build_states(tessier2000_names, tessier2000_model)
    computed: dict[str, dict[str, Any]] = {}
    skipped: list[str] = []
    for label, thunk in states:
        try:
            result = thunk()
        except ct.ConvergenceError:
            skipped.append(label)
            continue
        computed[label] = _encode(result)

    # The label sets must match exactly, not just intersect: a state that
    # silently stopped converging (or started converging when it used not
    # to) must fail this test rather than be quietly skipped.
    assert set(computed) == set(fixture), (
        f"State set changed since capture. Missing from a fresh run: "
        f"{sorted(set(fixture) - set(computed))}; new/unexpected: "
        f"{sorted(set(computed) - set(fixture))}; skipped this run: {skipped}"
    )
    # 144 phi-phi grid states (all 6 x 6 x 4 combinations converge; none is
    # skipped) + 5 gamma-gamma binary feeds + 1 single-component feed + 1
    # Tessier near-plait feed + 2 gamma-phi cases + 2 legacy cases = 155.
    assert len(fixture) == 155, len(fixture)
    assert not skipped, skipped

    for label, expected in fixture.items():
        assert computed[label] == expected, label
