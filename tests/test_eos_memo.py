"""The call-local equation-of-state solve memo answers exactly what re-solving answers.

ADR-0030 gives every ``flash_tp`` / ``stability_tp`` call a bounded memo of the
density/compressibility root solves it makes, so that a state solved once in a
call is not solved again in the same call. That is a **performance** change and
its whole gate is bit-identity, so this module compares with ``==`` and never
with a tolerance:

1. every user-visible field of every result - phase names, compositions,
   fractions, ``vapor_fraction`` and the whole ``diagnostics`` mapping - is
   identical with the memo and with the memo taken away, over the Case F-4
   default subset, the water / n-hexane liquid-liquid states and the
   Peng-Robinson grid slice;
2. the memo is **bounded**: it never holds more than ``max_entries`` and
   evicting is invisible in the answer, because an evicted entry can only be
   recomputed into the same doubles - a flash run with a one-entry memo is
   bit-identical to one run with the default bound;
3. the memo is **call-local**: nothing survives a call, and two models that
   differ only in a parameter never read each other's entries even when they
   are deliberately run inside one shared scope;
4. the memo actually hits. A cache that never fires would pass (1)-(3) and
   deliver nothing, so the hit counts are asserted too.

"Without the memo" is produced by making the models' ``active_memo`` lookup
return ``None``, which is exactly the code path they take outside a
``flash_tp`` / ``stability_tp`` call and exactly the pre-ADR-0030 one.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo import _eos_memo
from chemthermo.eos import PCSAFTEOS
from chemthermo.models import PengRobinsonEOS

#: Verbatim from ``tests/test_eos_branch_reuse.py``: the default-run subset of
#: the Case F-4 PC-SAFT grid, including the four states that need the ADR-0016
#: second-order stage.
F4_SUBSET: tuple[tuple[tuple[str, str], float, float, float], ...] = (
    (("Carbon dioxide", "n-Decane"), 0.8, 250.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 240.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 250.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 260.0, 1.5e6),
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

#: The ADR-0019 associating water / n-hexane liquid-liquid states.
WATER_HEXANE_STATES: tuple[tuple[float, float, tuple[float, float]], ...] = (
    (298.15, 101325.0, (0.5, 0.5)),
    (298.15, 101325.0, (0.2, 0.8)),
    (298.15, 1.0e6, (0.5, 0.5)),
    (335.0, 101325.0, (0.7, 0.3)),
)

#: A slice of the ADR-0017 Peng-Robinson grid, at the pressures and
#: temperatures where the cubic has one real root and where it has three.
PR_STATES: tuple[tuple[tuple[str, ...], tuple[float, ...], float, float], ...] = (
    (("Methane", "Ethane"), (0.5, 0.5), 170.0, 2.0e5),
    (("Methane", "Ethane"), (0.5, 0.5), 240.0, 3.0e6),
    (("Methane", "Ethane"), (0.5, 0.5), 360.0, 8.0e6),
    (("Ethane", "n-Heptane"), (0.7, 0.3), 280.0, 1.0e6),
    (("Methane", "n-Pentane"), (0.6, 0.4), 200.0, 8.0e6),
    (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2), 240.0, 3.0e6),
    (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3), 320.0, 2.0e5),
)

POLYMER_FIXTURE = (
    Path(__file__).resolve().parent / "fixtures" / "pcsaft" / "martini2009_polymers.json"
)
POLYMER_MW_G_MOL = 16400.0
PENTANE_MW_G_MOL = 72.146


# ---------------------------------------------------------------------------
# Running with and without the memo
# ---------------------------------------------------------------------------


def _no_memo() -> Any:
    """A context manager that disables the memo inside a single ``with``."""
    return pytest.MonkeyPatch().context()


def _run_without_memo(thunk: Callable[[], Any]) -> Any:
    from chemthermo.eos import pcsaft as pcsaft_module
    from chemthermo.models import peng_robinson as pr_module

    with _no_memo() as patch:
        patch.setattr(pcsaft_module, "active_memo", lambda: None)
        patch.setattr(pr_module, "active_memo", lambda: None)
        return thunk()


def _encode_flash(result: ct.FlashResult) -> dict[str, Any]:
    """Every user-visible field of a ``FlashResult``, JSON-ready."""
    return {
        "phase_names": list(result.phase_names()),
        "phases": {
            name: list(phase.composition.fractions) for name, phase in result.phases.items()
        },
        "phase_fractions": dict(result.phase_fractions),
        "vapor_fraction": result.vapor_fraction,
        "diagnostics": dict(result.diagnostics),
    }


def _encode_stability(result: ct.StabilityResult) -> dict[str, Any]:
    return {
        "status": result.status,
        "stable": result.stable,
        "tpd_min": result.tpd_min,
        "trial_composition": (
            None if result.trial_composition is None else list(result.trial_composition)
        ),
        "trial_ln_W": None if result.trial_ln_W is None else list(result.trial_ln_W),
        "diagnostics": dict(result.diagnostics),
    }


def _identical(left: Any, right: Any, path: str = "") -> None:
    """Assert two encoded results agree field for field, floats with ``==``."""
    assert type(left) is type(right), f"{path}: type {type(left)} vs {type(right)}"
    if isinstance(left, dict):
        assert left.keys() == right.keys(), f"{path}: keys {sorted(left)} vs {sorted(right)}"
        for key in left:
            _identical(left[key], right[key], f"{path}.{key}")
    elif isinstance(left, list):
        assert len(left) == len(right), f"{path}: length {len(left)} vs {len(right)}"
        for index, (a, b) in enumerate(zip(left, right)):
            _identical(a, b, f"{path}[{index}]")
    elif isinstance(left, float):
        if math.isnan(left) and math.isnan(right):
            return
        assert left == right, f"{path}: {left!r} != {right!r}"
    else:
        assert left == right, f"{path}: {left!r} != {right!r}"


def _mixture(names: Sequence[str], z: Sequence[float]) -> ct.Mixture:
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


def _polymer_eos() -> tuple[ct.Mixture, PCSAFTEOS]:
    """The ADR-0022 polyethylene(16400) / n-pentane split of the bench workload.

    Built exactly the way ``chemthermo.bench._cases._polymer_prepare`` builds
    it, from the same cited fixture: the state where ``exp(ln phi)`` underflows
    and the ADR-0022 logarithmic route runs.
    """
    payload = json.loads(POLYMER_FIXTURE.read_text(encoding="utf-8"))
    row = next(entry for entry in payload["polymers"] if entry["name"] == "Polyethylene")
    parameters = ct.PCSAFTParameters.from_records(
        [
            ct.PCSAFTRecord(
                name="Polyethylene",
                segments_per_g=row["segments_per_g"],
                MW_g_mol=POLYMER_MW_G_MOL,
                sigma_A=row["sigma_A"],
                epsilon_k_K=row["epsilon_k_K"],
                source="Martini et al. 2009 Table 1",
            ),
            ct.PCSAFTRecord(
                name="n-Pentane",
                m=2.6896,
                sigma_A=3.7729,
                epsilon_k_K=231.20,
                MW_g_mol=PENTANE_MW_G_MOL,
                source="Gross & Sadowski 2001 Table 1",
            ),
        ]
    )
    moles_polymer = 0.05 / POLYMER_MW_G_MOL
    moles_solvent = 0.95 / PENTANE_MW_G_MOL
    total = moles_polymer + moles_solvent
    mixture = ct.Mixture.from_components(
        [
            ct.Component.custom(
                "Polyethylene",
                mw_kg_per_mol=POLYMER_MW_G_MOL / 1000.0,
                formula="(C2H4)n",
                volatile=False,
                source="see tests/fixtures/pcsaft/martini2009_polymers.json",
            ),
            ct.Component.from_database("n-Pentane"),
        ],
        [moles_polymer / total, moles_solvent / total],
        normalize=True,
    )
    return mixture, PCSAFTEOS(parameters=parameters, kij=-0.006)


# ---------------------------------------------------------------------------
# 1) The answers do not move
# ---------------------------------------------------------------------------


def test_the_f4_subset_is_bit_identical_with_and_without_the_memo() -> None:
    """Every PC-SAFT Case F-4 default-subset flash, both routes, field for field."""
    eos = PCSAFTEOS()
    checked = 0
    for names, z1, temperature, pressure in F4_SUBSET:
        mixture = _mixture(names, (z1, 1.0 - z1))
        label = f"{names} z1={z1} T={temperature} P={pressure}"

        def thunk(
            mixture: ct.Mixture = mixture, T: float = temperature, P: float = pressure
        ) -> Any:
            return _encode_flash(ct.flash_tp(mixture, temperature_K=T, pressure_Pa=P, eos=eos))

        with_memo = thunk()
        without = _run_without_memo(thunk)
        _identical(with_memo, without, label)
        checked += 1
    assert checked == len(F4_SUBSET)


def test_the_associating_liquid_liquid_states_are_bit_identical() -> None:
    """Water / n-hexane, the ADR-0019 states, flash and stability alike."""
    eos = PCSAFTEOS(components=("Water", "n-Hexane"))
    for temperature, pressure, z in WATER_HEXANE_STATES:
        mixture = _mixture(("Water", "n-Hexane"), z)
        label = f"water/hexane T={temperature} P={pressure} z={z}"

        def flash(
            mixture: ct.Mixture = mixture, T: float = temperature, P: float = pressure
        ) -> Any:
            return _encode_flash(ct.flash_tp(mixture, temperature_K=T, pressure_Pa=P, eos=eos))

        def stability(
            mixture: ct.Mixture = mixture, T: float = temperature, P: float = pressure
        ) -> Any:
            return _encode_stability(
                ct.stability_tp(mixture, temperature_K=T, pressure_Pa=P, eos=eos)
            )

        _identical(flash(), _run_without_memo(flash), f"{label}/flash")
        _identical(stability(), _run_without_memo(stability), f"{label}/stability")


def test_the_peng_robinson_grid_slice_is_bit_identical() -> None:
    """The cubic's memo is keyed on ``(A, B)`` alone; one and three real roots."""
    eos = PengRobinsonEOS()
    for names, z, temperature, pressure in PR_STATES:
        mixture = _mixture(names, z)
        label = f"{names} T={temperature} P={pressure}"

        def flash(
            mixture: ct.Mixture = mixture, T: float = temperature, P: float = pressure
        ) -> Any:
            return _encode_flash(ct.flash_tp(mixture, temperature_K=T, pressure_Pa=P, eos=eos))

        def stability(
            mixture: ct.Mixture = mixture, T: float = temperature, P: float = pressure
        ) -> Any:
            return _encode_stability(
                ct.stability_tp(mixture, temperature_K=T, pressure_Pa=P, eos=eos)
            )

        _identical(flash(), _run_without_memo(flash), f"{label}/flash")
        _identical(stability(), _run_without_memo(stability), f"{label}/stability")


@pytest.mark.skipif(not POLYMER_FIXTURE.exists(), reason="polymer fixture not in this checkout")
def test_the_polymer_split_is_bit_identical() -> None:
    """The ADR-0022 state where ``exp(ln phi)`` underflows and the log route runs."""
    mixture, eos = _polymer_eos()

    def flash() -> Any:
        return _encode_flash(ct.flash_tp(mixture, temperature_K=453.15, pressure_Pa=8.0e6, eos=eos))

    _identical(flash(), _run_without_memo(flash), "polymer")


# ---------------------------------------------------------------------------
# 2) The memo is bounded, and the bound cannot move a number
# ---------------------------------------------------------------------------


def test_the_memo_evicts_in_insertion_order_and_never_exceeds_its_bound() -> None:
    memo = _eos_memo.EosSolveMemo(max_entries=3)
    for index in range(10):
        memo.store(index, index * 10)
        assert len(memo) <= 3
    assert len(memo) == 3
    # The three most recent survive; the seven oldest were evicted.
    assert memo.lookup(9) == 90
    assert memo.lookup(8) == 80
    assert memo.lookup(7) == 70
    assert memo.lookup(6) is _eos_memo.MISS
    assert memo.lookup(0) is _eos_memo.MISS
    # Re-storing an existing key replaces it rather than growing the memo.
    memo.store(9, 999)
    assert len(memo) == 3
    assert memo.lookup(9) == 999


def test_a_memo_bound_of_one_gives_the_same_answer_as_the_default_bound() -> None:
    """Eviction can only cost a re-solve, never a different double."""
    eos = PCSAFTEOS(components=("Water", "n-Hexane"))
    mixture = _mixture(("Water", "n-Hexane"), (0.5, 0.5))

    default = _encode_flash(
        ct.flash_tp(mixture, temperature_K=298.15, pressure_Pa=101325.0, eos=eos)
    )
    tiny = _eos_memo.EosSolveMemo(max_entries=1)
    with _eos_memo.activated(tiny):
        bounded = _encode_flash(
            ct.flash_tp(mixture, temperature_K=298.15, pressure_Pa=101325.0, eos=eos)
        )
    assert len(tiny) == 1
    _identical(default, bounded, "water/hexane one-entry memo")


def test_max_entries_must_be_at_least_one() -> None:
    with pytest.raises(ValueError):
        _eos_memo.EosSolveMemo(max_entries=0)


# ---------------------------------------------------------------------------
# 3) The memo is call-local
# ---------------------------------------------------------------------------


def test_no_memo_is_active_outside_a_call() -> None:
    assert _eos_memo.active_memo() is None
    eos = PengRobinsonEOS()
    mixture = _mixture(("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2))
    ct.flash_tp(mixture, temperature_K=240.0, pressure_Pa=3.0e6, eos=eos)
    assert _eos_memo.active_memo() is None, "the scope must not outlive the call"


def test_a_nested_stability_call_shares_the_flashs_memo() -> None:
    """One flash is one memo: the nested stability tests do not each open one."""
    seen: list[_eos_memo.EosSolveMemo] = []
    eos = PengRobinsonEOS()
    mixture = _mixture(("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2))
    with _eos_memo.activated() as outer:
        ct.flash_tp(mixture, temperature_K=240.0, pressure_Pa=3.0e6, eos=eos)
        seen.append(outer)
        inner_observed = _eos_memo.active_memo()
    assert inner_observed is outer
    assert outer.hits > 0


def test_two_models_never_read_each_others_entries() -> None:
    """Different ``kij`` at one ``(T, P, x)``: same arguments, different answers.

    Both flashes run inside **one** shared memo - the worst case, wider than
    anything the library itself opens - so if the key did not separate the two
    model instances the second call would be served the first one's roots.
    """
    names = ("Methane", "n-Hexane")
    mixture = _mixture(names, (0.5, 0.5))
    plain = PCSAFTEOS(components=names)
    shifted = PCSAFTEOS(components=names, kij=0.03)

    def run(eos: PCSAFTEOS) -> Any:
        return _encode_flash(ct.flash_tp(mixture, temperature_K=200.0, pressure_Pa=3.5e6, eos=eos))

    alone_plain = run(plain)
    alone_shifted = run(shifted)

    with _eos_memo.activated() as shared:
        together_plain = run(plain)
        together_shifted = run(shifted)
    assert shared.hits > 0, "the shared scope must actually be memoizing"

    _identical(alone_plain, together_plain, "plain kij")
    _identical(alone_shifted, together_shifted, "shifted kij")
    assert alone_plain["phases"] != alone_shifted["phases"], (
        "the two models must disagree, or this test proves nothing"
    )


def test_the_peng_robinson_key_separates_two_states_that_share_a_composition() -> None:
    """``(A, B)`` carries temperature, pressure and the mixing rule; nothing else has to."""
    eos = PengRobinsonEOS()
    mixture = _mixture(("Methane", "Ethane"), (0.5, 0.5))
    with _eos_memo.activated():
        cold = eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=170.0,
            pressure_Pa=2.0e5,
            composition=[0.5, 0.5],
            phase="liquid",
        )
        hot = eos.fugacity_coefficients(
            mixture=mixture,
            temperature_K=360.0,
            pressure_Pa=2.0e5,
            composition=[0.5, 0.5],
            phase="liquid",
        )
    assert not np.array_equal(np.asarray(cold), np.asarray(hot))


# ---------------------------------------------------------------------------
# 4) The memo hits, and the key is exact
# ---------------------------------------------------------------------------


def test_the_memo_hits_on_the_reference_paths() -> None:
    """A memo that never fires would pass every test above and deliver nothing."""
    cases: list[tuple[str, Callable[[], Any]]] = []

    pcsaft = PCSAFTEOS(components=("Water", "n-Hexane"))
    water_hexane = _mixture(("Water", "n-Hexane"), (0.5, 0.5))
    cases.append(
        (
            "pcsaft-lle",
            lambda: ct.flash_tp(
                water_hexane, temperature_K=298.15, pressure_Pa=101325.0, eos=pcsaft
            ),
        )
    )

    pr = PengRobinsonEOS()
    ternary = _mixture(("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2))
    cases.append(
        ("pr-flash", lambda: ct.flash_tp(ternary, temperature_K=240.0, pressure_Pa=3.0e6, eos=pr))
    )
    cases.append(
        (
            "pr-stability",
            lambda: ct.stability_tp(ternary, temperature_K=240.0, pressure_Pa=3.0e6, eos=pr),
        )
    )

    for label, thunk in cases:
        with _eos_memo.activated() as memo:
            thunk()
        assert memo.hits > 0, f"{label}: the memo never hit"
        assert memo.misses > 0, f"{label}: nothing was solved at all"


def test_the_composition_key_separates_the_two_zeros() -> None:
    """``-0.0 == 0.0`` and they hash alike; the key must not merge them."""
    positive = _eos_memo.composition_key([0.0, 0.5, 0.5])
    negative = _eos_memo.composition_key([-0.0, 0.5, 0.5])
    assert positive != negative
    assert hash(positive) != hash(negative) or positive != negative
    # Ordinary values pass through unchanged, so the key is the arguments.
    assert _eos_memo.composition_key([0.25, 0.75]) == (0.25, 0.75)


def test_a_not_a_number_composition_simply_misses() -> None:
    """``nan != nan``, so a ``nan`` argument can never read another state's entry."""
    memo = _eos_memo.EosSolveMemo()
    key = _eos_memo.composition_key([float("nan"), 1.0])
    memo.store(key, "solved")
    assert memo.lookup(_eos_memo.composition_key([float("nan"), 1.0])) is _eos_memo.MISS
