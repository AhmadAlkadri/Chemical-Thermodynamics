"""The three-phase neighbourhood of water / n-hexane against FeOs (Case P-9).

ADR-0020 wires the ADR-0011 phase addition/removal search to the phi-phi path,
so the band of temperatures just below the water / n-hexane three-phase
temperature `T3` - which raised `ConvergenceError` before - now returns the two
conjugate liquids, reached by *adding* a second liquid to a vapour-liquid pair
and then *removing* the vapour. This file checks those answers against an
implementation that shares no code with them.

**FeOs** (feos-org/feos, MIT OR Apache-2.0) is that implementation: the same
Gross & Sadowski model in Rust, every derivative by automatic differentiation.
Two independent things are compared:

1. FeOs's own **chemical potentials evaluated at chemthermo's converged phases
   and densities** - the reference model's equilibrium condition, at
   chemthermo's answer, which does not depend on FeOs's flash converging;
2. FeOs's own **two-phase flash** on the two-liquid state below `T3`.

As in Cases P-6, P-7 and P-8, one input is deliberately not shared: the 42
universal constants of the 2001 dispersion term (chemthermo packages the ten
figures as printed, FeOs hard-codes fourteen). Every chemical-potential
comparison is therefore run twice, as shipped and with FeOs's constants
substituted, and both numbers are recorded in validation Case P-9. Skipped
when `feos` is not installed (`pip install -e ".[validation]"`).
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.eos import pcsaft as pcsaft_module
from chemthermo.eos.pcsaft import R_J_PER_MOL_K

feos = pytest.importorskip("feos")
si = pytest.importorskip("si_units")

from feos import (  # noqa: E402
    Contributions,
    EquationOfState,
    Parameters,
    PureRecord,
    State,
)

BINARY = ("Water", "n-Hexane")
TERNARY = ("Water", "Ethanol", "n-Hexane")
ATMOSPHERE_PA = 101325.0
OFFSET_K = 0.05
TERNARY_T_K = 333.0

#: Gross & Sadowski (2002) Table 1 for water and ethanol, (2001) Table 1 for
#: n-hexane, written out here so the reference model is built from this file
#: and not from the package databank.
#: ``name -> (MW, m, sigma/A, eps/k, kappa^AB, eps^AB/k)``.
PARAMETERS: dict[str, tuple[float, float, float, float, float | None, float]] = {
    "Water": (18.015, 1.0656, 3.0007, 366.51, 0.034868, 2500.7),
    "Ethanol": (46.069, 2.3827, 3.1771, 198.24, 0.032384, 2653.4),
    "n-Hexane": (86.177, 3.0576, 3.7983, 236.77, None, 0.0),
}

_MOL_PER_M3 = si.MOL / si.METER**3

#: Asserted tolerances (validation Case P-9).
POTENTIAL_TOL = 1e-8
COMPOSITION_TOL = 1e-8
DENSITY_TOL = 1e-6
FRACTION_TOL = 1e-8


def _feos_universal_constants() -> tuple[np.ndarray, np.ndarray]:
    """FeOs's fourteen-figure ``a`` and ``b`` tables, from the Case P-6 module.

    Loaded by path, not with ``from tests.validation... import``: CI runs the
    ``pytest`` console script, which does not put the working directory on
    ``sys.path`` (see ``.agents/dev-contract.md``).
    """
    path = Path(__file__).with_name("test_pcsaft_association_vs_feos.py")
    spec = importlib.util.spec_from_file_location("_pcsaft_feos_constants", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.FEOS_A_UNIVERSAL, module.FEOS_B_UNIVERSAL


@pytest.fixture
def matched_constants(monkeypatch: pytest.MonkeyPatch) -> None:
    a_universal, b_universal = _feos_universal_constants()
    monkeypatch.setattr(pcsaft_module, "A_UNIVERSAL", a_universal)
    monkeypatch.setattr(pcsaft_module, "B_UNIVERSAL", b_universal)


def _pure_record(name: str) -> PureRecord:
    mw, m, sigma, epsilon, kappa, epsilon_ab = PARAMETERS[name]
    payload: dict[str, object] = {
        "identifier": {"name": name},
        "molarweight": mw,
        "m": m,
        "sigma": sigma,
        "epsilon_k": epsilon,
    }
    if kappa is not None:
        payload["association_sites"] = [
            {"kappa_ab": kappa, "epsilon_k_ab": epsilon_ab, "na": 1.0, "nb": 1.0}
        ]
    return PureRecord.from_json_str(json.dumps(payload))


def _feos_eos(names: Sequence[str]) -> EquationOfState:
    """The reference model, always with ``k_ij = 0`` (FeOs's own default)."""
    return EquationOfState.pcsaft(Parameters.from_records([_pure_record(n) for n in names]))


def _feos_reduced_potentials(
    names: Sequence[str], temperature_K: float, density: float, x: Sequence[float] | np.ndarray
) -> np.ndarray:
    """``mu_i / RT`` from **FeOs** at chemthermo's ``(T, rho, x)``, up to a constant.

    ``mu_i / RT = ln(x_i phi_i) + ln P + const``, and ``ln phi_i`` follows from
    FeOs's residual chemical potential as ``mu_i^res / RT - ln Z``. All phases
    of one flash are at the same ``T`` and ``P``, so the omitted constant
    cancels in the differences taken below.
    """
    values = np.asarray(x, dtype=float)
    state = State(
        _feos_eos(names),
        temperature=temperature_K * si.KELVIN,
        density=density * _MOL_PER_M3,
        composition=values,
    )
    factor = R_J_PER_MOL_K * temperature_K
    z_factor = (state.pressure() / si.PASCAL) / (density * factor)
    mu_residual = np.atleast_1d(
        np.asarray(state.chemical_potential(Contributions.Residual) / si.JOULE * si.MOL)
    )
    return mu_residual / factor - math.log(z_factor) + np.log(values)


def _three_phase_temperature() -> float:
    """The binary three-phase temperature from a 4-equation Newton (this file's own).

    Deliberately re-derived here rather than imported from
    ``tests/test_flash_vlle_eos.py``: see the ``sys.path`` note above, and the
    two solves are independent evidence of the same number.
    """
    eos = PCSAFTEOS()
    mixture = ct.Mixture.from_database(list(BINARY), [0.5, 0.5], normalize=True)

    def ln_f(x1: float, temperature: float, branch: str) -> np.ndarray:
        x = np.array([x1, 1.0 - x1])
        phi = np.asarray(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature,
                pressure_Pa=ATMOSPHERE_PA,
                composition=x.tolist(),
                phase=branch,
            ),
            dtype=float,
        )
        return np.log(x) + np.log(phi)

    def residual(u: np.ndarray) -> np.ndarray:
        first = ln_f(u[0], u[3], "liquid")
        second = ln_f(u[1], u[3], "liquid")
        vapor = ln_f(u[2], u[3], "vapor")
        return np.concatenate([first - second, first - vapor])

    u = np.array([0.9999, 0.02, 0.20, 334.5])
    steps = (1e-8, 1e-8, 1e-8, 3.345e-5)
    for _iteration in range(40):
        f = residual(u)
        worst = float(np.max(np.abs(f)))
        if worst < 1e-11:
            break
        jacobian = np.zeros((4, 4))
        for column in range(4):
            shifted = u.copy()
            shifted[column] += steps[column]
            jacobian[:, column] = (residual(shifted) - f) / steps[column]
        direction = np.linalg.solve(jacobian, -f)
        scale = 1.0
        while scale > 1e-10:
            candidate = u + scale * direction
            inside = bool(np.all(candidate[:3] > 0.0) and np.all(candidate[:3] < 1.0))
            if inside and float(np.max(np.abs(residual(candidate)))) < worst:
                break
            scale *= 0.5
        u = u + scale * direction
    assert float(np.max(np.abs(residual(u)))) < 1e-10
    return float(u[3])


@pytest.fixture(scope="module")
def three_phase_temperature() -> float:
    return _three_phase_temperature()


def _our_phases(
    names: Sequence[str], temperature_K: float, z: Sequence[float]
) -> tuple[ct.FlashResult, dict[str, tuple[np.ndarray, float, float]]]:
    """``flash_tp`` plus, per phase, ``(composition, density, phase fraction)``.

    The density reported for a phase is the root its measured identity names -
    the lowest for a vapour, the highest for a liquid - which is the root the
    split put that phase on (ADR-0019).
    """
    mixture = ct.Mixture.from_database(list(names), list(z), normalize=True)
    eos = PCSAFTEOS()
    result = ct.flash_tp(mixture, temperature_K=temperature_K, pressure_Pa=ATMOSPHERE_PA, eos=eos)
    phases: dict[str, tuple[np.ndarray, float, float]] = {}
    for name, phase in result.phases.items():
        fractions = np.asarray(phase.composition.fractions, dtype=float)
        roots = eos.density_roots(
            mixture=mixture,
            temperature_K=temperature_K,
            pressure_Pa=ATMOSPHERE_PA,
            composition=fractions.tolist(),
        )
        density = roots[0] if name == "vapor" else roots[-1]
        phases[name] = (fractions, density, result.phase_fractions[name])
    return result, phases


def _worst_potential_difference(
    names: Sequence[str],
    temperature_K: float,
    phases: dict[str, tuple[np.ndarray, float, float]],
) -> float:
    potentials = [
        _feos_reduced_potentials(names, temperature_K, density, composition)
        for composition, density, _fraction in phases.values()
    ]
    return max(
        float(np.max(np.abs(potentials[i] - potentials[j])))
        for i in range(len(potentials))
        for j in range(i + 1, len(potentials))
    )


# ---------------------------------------------------------------------------
# Case P-9: the binary window
# ---------------------------------------------------------------------------


def test_feos_own_flash_returns_the_metastable_pair_below_t3(three_phase_temperature) -> None:
    """Observed behaviour of the reference, recorded rather than worked around.

    Just below `T3` the model has **two** two-phase stationary states: the
    vapour-liquid pair and the pair of conjugate liquids. FeOs's own
    `State.tp_flash` converges on the vapour-liquid one (feos 0.10.1, measured
    here). chemthermo converges on the same pair first - the assertions below
    compare them digit for digit, with `post_split_stability=False` to see it -
    and then its post-split stability test refuses it and the ADR-0020 search
    replaces it with the two liquids, which have the **lower** reduced Gibbs
    energy.

    So this is not a disagreement about the model: both answers are stationary
    states of it, and the Gibbs comparison at the end says which is the
    equilibrium. Every quantity in that comparison is computed here from the
    public `fugacity_coefficients`, not read out of `flash_tp`.
    """
    temperature = three_phase_temperature - OFFSET_K

    state = State(
        _feos_eos(BINARY),
        temperature=temperature * si.KELVIN,
        pressure=ATMOSPHERE_PA * si.PASCAL,
        composition=np.array([0.5, 0.5]),
    )
    equilibrium = state.tp_flash()
    reference_liquid = np.atleast_1d(np.asarray(equilibrium.liquid.molefracs, dtype=float))
    reference_vapor = np.atleast_1d(np.asarray(equilibrium.vapor.molefracs, dtype=float))
    assert float(equilibrium.liquid.density / _MOL_PER_M3) > 5000.0
    assert float(equilibrium.vapor.density / _MOL_PER_M3) < 200.0

    unpoliced = ct.flash_tp(
        ct.Mixture.from_database(list(BINARY), [0.5, 0.5], normalize=True),
        temperature_K=temperature,
        pressure_Pa=ATMOSPHERE_PA,
        eos=PCSAFTEOS(),
        settings=ct.FlashSettings(post_split_stability=False),
    )
    assert sorted(unpoliced.phases) == ["liquid", "vapor"]
    assert unpoliced.diagnostics["post_split_stable"] is False
    assert float(unpoliced.phases["liquid"].composition.fractions[0]) == pytest.approx(
        reference_liquid[0], abs=COMPOSITION_TOL
    )
    assert float(unpoliced.phases["vapor"].composition.fractions[0]) == pytest.approx(
        reference_vapor[0], abs=COMPOSITION_TOL
    )
    assert unpoliced.phase_fractions["vapor"] == pytest.approx(
        float(equilibrium.vapor_phase_fraction), abs=FRACTION_TOL
    )

    # And the answer the search returns instead is the lower-Gibbs one.
    result, _ours = _our_phases(BINARY, temperature, (0.5, 0.5))
    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert result.diagnostics["phase_set_history"] == "V -> LV -> LLV -> LL"

    mixture = ct.Mixture.from_database(list(BINARY), [0.5, 0.5], normalize=True)
    eos = PCSAFTEOS()

    def reduced_g(composition: Sequence[float], branch: str) -> float:
        x = np.asarray(composition, dtype=float)
        phi = np.asarray(
            eos.fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature,
                pressure_Pa=ATMOSPHERE_PA,
                composition=x.tolist(),
                phase=branch,
            ),
            dtype=float,
        )
        return float(np.sum(x * (np.log(x) + np.log(phi))))

    vapor_fraction = float(unpoliced.phase_fractions["vapor"])
    g_vapor_liquid = (1.0 - vapor_fraction) * reduced_g(
        unpoliced.phases["liquid"].composition.fractions, "liquid"
    ) + vapor_fraction * reduced_g(unpoliced.phases["vapor"].composition.fractions, "vapor")
    g_liquids = sum(
        result.phase_fractions[name]
        * reduced_g(result.phases[name].composition.fractions, "liquid")
        for name in ("liquid1", "liquid2")
    )
    assert g_liquids < g_vapor_liquid
    assert g_vapor_liquid - g_liquids > 1e-4


@pytest.mark.parametrize(
    "offset_K",
    # The below-T3 state is the one this slice made reachable, so it runs by
    # default; the above-T3 side is the pre-existing vapour-liquid answer and
    # is `slow` (ADR-0020 runtime trim).
    [-OFFSET_K, pytest.param(OFFSET_K, marks=pytest.mark.slow)],
    ids=["below-T3", "above-T3"],
)
@pytest.mark.usefixtures("matched_constants")
def test_feos_chemical_potentials_are_equal_across_the_window(
    three_phase_temperature, offset_K: float
) -> None:
    """The reference model's equilibrium condition at chemthermo's answers.

    Below `T3` those are the two liquids the search reaches; above it, the
    vapour-liquid pair. Run with matched constants, as Cases P-7 and P-8 did
    and for the same reason - FeOs is evaluating chemthermo's own densities, so
    the universal-constants difference would otherwise floor the residual near
    2e-06. Case P-9 records both numbers.
    """
    temperature = three_phase_temperature + offset_K
    phases = _our_phases(BINARY, temperature, (0.5, 0.5))[1]
    expected = ["liquid1", "liquid2"] if offset_K < 0 else ["liquid", "vapor"]
    assert sorted(phases) == expected
    assert _worst_potential_difference(BINARY, temperature, phases) < POTENTIAL_TOL


@pytest.mark.slow
def test_the_agreement_below_t3_is_not_vacuous(three_phase_temperature) -> None:
    """Perturbing water's association energy by 1 % must move the tie line."""
    temperature = three_phase_temperature - OFFSET_K
    _result, ours = _our_phases(BINARY, temperature, (0.5, 0.5))
    perturbed = ct.PCSAFTParameters.from_records(
        [
            {
                "name": "Water",
                "m": 1.0656,
                "sigma_A": 3.0007,
                "epsilon_k_K": 366.51,
                "association": {"kappa_ab": 0.034868, "epsilon_ab_k_K": 2500.7 * 1.01},
            },
            {
                "name": "n-Hexane",
                "m": 3.0576,
                "sigma_A": 3.7983,
                "epsilon_k_K": 236.77,
            },
        ]
    )
    shifted = ct.flash_tp(
        ct.Mixture.from_database(list(BINARY), [0.5, 0.5], normalize=True),
        temperature_K=temperature,
        pressure_Pa=ATMOSPHERE_PA,
        eos=PCSAFTEOS(components=BINARY, parameters=perturbed),
    )
    moved = abs(
        float(shifted.phases["liquid2"].composition.fractions[0]) - float(ours["liquid2"][0][0])
    )
    assert moved > 1e-4, moved


# ---------------------------------------------------------------------------
# Case P-10: the ternary tie-triangle
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.usefixtures("matched_constants")
def test_feos_chemical_potentials_are_equal_across_the_ternary_triangle() -> None:
    """Case P-10: three phases, three chemical potentials, one reference model."""
    result, ours = _our_phases(TERNARY, TERNARY_T_K, (0.4, 0.3, 0.3))
    assert sorted(result.phases) == ["liquid1", "liquid2", "vapor"]
    assert result.diagnostics["phase_regime"] == "VLLE"
    assert _worst_potential_difference(TERNARY, TERNARY_T_K, ours) < POTENTIAL_TOL


# ---------------------------------------------------------------------------
# Case P-18: the states ADR-0029 repaired
# ---------------------------------------------------------------------------
#
# Both are `multiphase-solver-failure` states of the robustness map at
# `74820b8` that ADR-0029 turned into verified vapour-liquid answers, and both
# are checked here the same way everything above is: FeOs's own chemical
# potentials at chemthermo's converged phases and densities, with matched
# universal constants. They are `slow`-marked as *repetitions* - the default
# run already makes this comparison on the two liquids below `T3`, which is the
# same measurement on the same binary - and the solver-side verification of
# these states (an independent Newton, the Gibbs ordering, the mass balance)
# lives in `tests/test_flash_eos_multiphase_robustness.py` and runs by default.


#: The water-lean feed of Case P-18 (i): 4 states above `T3` that used to raise
#: "a two-phase set converged to a non-positive phase fraction".
WATER_LEAN_FEED = (0.05, 0.95)


@pytest.mark.slow
@pytest.mark.parametrize("offset_K", [0.01, 1.0], ids=["T3+0.01", "T3+1.0"])
@pytest.mark.usefixtures("matched_constants")
def test_feos_chemical_potentials_agree_on_the_water_lean_band_above_t3(
    three_phase_temperature, offset_K: float
) -> None:
    """Case P-18 (i): the hexane-rich liquid against the vapour, per FeOs."""
    temperature = three_phase_temperature + offset_K
    result, ours = _our_phases(BINARY, temperature, WATER_LEAN_FEED)
    assert sorted(result.phases) == ["liquid", "vapor"]
    assert result.diagnostics["phase_regime"] == "VLE"
    assert float(ours["liquid"][0][0]) < 0.05
    assert _worst_potential_difference(BINARY, temperature, ours) < POTENTIAL_TOL


@pytest.mark.slow
@pytest.mark.usefixtures("matched_constants")
def test_feos_chemical_potentials_agree_at_the_hexane_rich_ternary_corner() -> None:
    """Case P-18 (i), ternary: `(0.1, 0.1, 0.8)` at 333 K, two phases not three."""
    result, ours = _our_phases(TERNARY, TERNARY_T_K, (0.1, 0.1, 0.8))
    assert sorted(result.phases) == ["liquid", "vapor"]
    assert _worst_potential_difference(TERNARY, TERNARY_T_K, ours) < POTENTIAL_TOL
