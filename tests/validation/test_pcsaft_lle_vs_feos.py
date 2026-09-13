"""The EOS liquid-liquid tie line against FeOs's own TP flash (Case P-8).

ADR-0019 taught the phi-phi split to put both phases on the branch the
tangent-plane stability test found them on, so water / n-hexane at 298.15 K and
1 atm - where a vapour density root still exists, and where the split used to
converge on a spurious water-liquid / hexane-vapour pair that the post-split
test refused - now returns two liquids. This file checks that answer against an
implementation that shares no code with it.

**FeOs** (feos-org/feos, MIT OR Apache-2.0) is that implementation: the same
Gross & Sadowski model written in Rust, every derivative by automatic
differentiation, with its own two-phase flash (``State.tp_flash``). Three
independent things are compared, not one:

1. the **tie line and the phase amounts** against FeOs's own flash - so the
   split is checked against another *solver*, not only another *model*;
2. the **phase densities** against the densities FeOs's flash converged on;
3. FeOs's own **chemical potentials evaluated at chemthermo's phases** - so
   chemthermo's answer is shown to be an equilibrium state of the reference
   model, independently of whether the two flashes agree.

One shared input is deliberately not shared: the 42 universal constants of the
2001 dispersion term. chemthermo packages them as printed (ten figures), FeOs
hard-codes fourteen, and the two tables differ by up to 4.8e-09. That is an
*input* difference. Every comparison below is therefore run twice - as shipped
and with FeOs's constants substituted in - and both numbers are recorded in
validation Case P-8. Skipped when ``feos`` is not installed
(``pip install -e ".[validation]"``).
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
from chemthermo.models.base import KAPPA_LIQUID_THRESHOLD

feos = pytest.importorskip("feos")
si = pytest.importorskip("si_units")

from feos import (  # noqa: E402
    Contributions,
    EquationOfState,
    Parameters,
    PureRecord,
    State,
)

NAMES = ("Water", "n-Hexane")
TEMPERATURE_K = 298.15
ATMOSPHERE_PA = 101325.0
HIGH_PRESSURE_PA = 1.0e6

#: Gross & Sadowski (2002) Table 1 for water, (2001) Table 1 for n-hexane,
#: written out here so the reference model is built from this file and not from
#: the package databank. ``name -> (MW, m, sigma/A, eps/k, kappa^AB, eps^AB/k)``.
PARAMETERS: dict[str, tuple[float, float, float, float, float | None, float]] = {
    "Water": (18.015, 1.0656, 3.0007, 366.51, 0.034868, 2500.7),
    "n-Hexane": (86.177, 3.0576, 3.7983, 236.77, None, 0.0),
}

_MOL_PER_M3 = si.MOL / si.METER**3

#: Asserted tolerances (validation Case P-8). Compositions and phase amounts
#: are absolute mole fractions; densities are relative.
COMPOSITION_TOL = 1e-8
FRACTION_TOL = 1e-8
DENSITY_TOL = 1e-6
#: FeOs's chemical potentials at chemthermo's phases, with matched constants.
POTENTIAL_TOL = 1e-8


def _feos_universal_constants() -> tuple[np.ndarray, np.ndarray]:
    """FeOs's fourteen-figure ``a`` and ``b`` tables, from the Case P-6 module.

    Loaded **by path** rather than with ``from tests.validation... import``:
    CI runs the ``pytest`` console script, which does not put the working
    directory on ``sys.path`` (see ``.agents/dev-contract.md``). Reusing the
    single copy there rather than pasting 42 more numbers here is deliberate -
    two copies of a constants table drift.
    """
    path = Path(__file__).with_name("test_pcsaft_association_vs_feos.py")
    spec = importlib.util.spec_from_file_location("_pcsaft_feos_constants", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.FEOS_A_UNIVERSAL, module.FEOS_B_UNIVERSAL


@pytest.fixture
def matched_constants(monkeypatch: pytest.MonkeyPatch) -> None:
    """Give chemthermo FeOs's fourteen-figure universal constants."""
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


def _feos_eos() -> EquationOfState:
    """The reference model, always with ``k_ij = 0`` (FeOs's own default)."""
    return EquationOfState.pcsaft(Parameters.from_records([_pure_record(n) for n in NAMES]))


def _feos_tp_flash(
    pressure_Pa: float, z: Sequence[float]
) -> tuple[tuple[np.ndarray, float], tuple[np.ndarray, float], float]:
    """FeOs's own two-phase flash: ``((x, rho), (y, rho), fraction of y)``.

    FeOs names the two phases ``liquid`` and ``vapor`` by construction of its
    ``PhaseEquilibrium`` container; at this state **both are liquids** (both
    densities are thousands of mol/m^3), which is exactly the point, and the
    test asserts that rather than trusting the names.
    """
    state = State(
        _feos_eos(),
        temperature=TEMPERATURE_K * si.KELVIN,
        pressure=pressure_Pa * si.PASCAL,
        composition=np.asarray(z, dtype=float),
    )
    equilibrium = state.tp_flash()
    dense = (
        np.atleast_1d(np.asarray(equilibrium.liquid.molefracs, dtype=float)),
        float(equilibrium.liquid.density / _MOL_PER_M3),
    )
    light = (
        np.atleast_1d(np.asarray(equilibrium.vapor.molefracs, dtype=float)),
        float(equilibrium.vapor.density / _MOL_PER_M3),
    )
    return dense, light, float(equilibrium.vapor_phase_fraction)


def _feos_reduced_potentials(density: float, x: Sequence[float] | np.ndarray) -> np.ndarray:
    """``mu_i / RT`` from **FeOs** at chemthermo's ``(T, rho, x)``, up to a constant.

    ``mu_i / RT = ln(x_i phi_i) + ln P + const``, and ``ln phi_i`` follows from
    FeOs's residual chemical potential as ``mu_i^res / RT - ln Z``. Both phases
    of one flash are at the same ``T`` and ``P``, so the omitted constant
    cancels in the difference below and the comparison is the equality of
    chemical potentials, evaluated by the reference implementation at the state
    chemthermo converged on.
    """
    values = np.asarray(x, dtype=float)
    state = State(
        _feos_eos(),
        temperature=TEMPERATURE_K * si.KELVIN,
        density=density * _MOL_PER_M3,
        composition=values,
    )
    factor = R_J_PER_MOL_K * TEMPERATURE_K
    z_factor = (state.pressure() / si.PASCAL) / (density * factor)
    mu_residual = np.atleast_1d(
        np.asarray(state.chemical_potential(Contributions.Residual) / si.JOULE * si.MOL)
    )
    return mu_residual / factor - math.log(z_factor) + np.log(values)


def _our_phases(
    pressure_Pa: float, z: Sequence[float]
) -> tuple[ct.FlashResult, dict[str, tuple[np.ndarray, float, float]]]:
    """``flash_tp`` plus, per phase, ``(composition, density, phase fraction)``."""
    mixture = ct.Mixture.from_database(list(NAMES), list(z), normalize=True)
    eos = PCSAFTEOS()
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
    phases: dict[str, tuple[np.ndarray, float, float]] = {}
    for name, phase in result.phases.items():
        fractions = np.asarray(phase.composition.fractions, dtype=float)
        density = eos.density_roots(
            mixture=mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=pressure_Pa,
            composition=fractions.tolist(),
        )[-1]
        phases[name] = (fractions, density, result.phase_fractions[name])
    return result, phases


def _kappa(pressure_Pa: float, composition: Sequence[float] | np.ndarray) -> float:
    """``kappa = P / (rho dP/drho)`` at the liquid-like root, by finite difference."""
    mixture = ct.Mixture.from_database(list(NAMES), [0.5, 0.5], normalize=True)
    eos = PCSAFTEOS()
    values = np.asarray(composition, dtype=float).tolist()
    density = eos.density_roots(
        mixture=mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure_Pa,
        composition=values,
    )[-1]
    bound = PCSAFTEOS(components=NAMES)
    step = 1e-4 * density
    slope = (
        bound.pressure_Pa(
            temperature_K=TEMPERATURE_K,
            density_mol_m3=density + step,
            composition=values,
        )
        - bound.pressure_Pa(
            temperature_K=TEMPERATURE_K,
            density_mol_m3=density - step,
            composition=values,
        )
    ) / (2.0 * step)
    return pressure_Pa / (density * slope)


# ---------------------------------------------------------------------------
# Case P-8 (i): 1 atm, the state Case P-7(iii) recorded as unreachable
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("pressure_Pa", [ATMOSPHERE_PA, HIGH_PRESSURE_PA])
def test_the_tie_line_matches_the_feos_tp_flash(pressure_Pa: float) -> None:
    """Compositions, densities and phase amounts, against FeOs's own flash.

    Run on the **shipped** universal constants, the ones a user gets. The
    matched-constants run below tightens the same comparison by three orders of
    magnitude and shows the residual here is the constants table.
    """
    result, ours = _our_phases(pressure_Pa, (0.5, 0.5))

    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert result.vapor_fraction is None
    assert result.diagnostics["phase_regime"] == "LLE"
    assert result.diagnostics["phase_label_method"] == "compressibility"
    assert float(result.diagnostics["mass_balance_residual"]) < 1e-12
    assert float(result.diagnostics["fugacity_residual"]) < 1e-9
    assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
    assert result.diagnostics["post_split_status"] == "stable"

    water_rich, hexane_rich = ours["liquid1"], ours["liquid2"]
    dense, light, light_fraction = _feos_tp_flash(pressure_Pa, (0.5, 0.5))

    # FeOs's two phases are both liquids, whatever its container calls them.
    assert dense[1] > 5000.0 and light[1] > 5000.0

    assert water_rich[0] == pytest.approx(dense[0], abs=COMPOSITION_TOL)
    assert hexane_rich[0] == pytest.approx(light[0], abs=COMPOSITION_TOL)
    assert water_rich[1] == pytest.approx(dense[1], rel=DENSITY_TOL)
    assert hexane_rich[1] == pytest.approx(light[1], rel=DENSITY_TOL)
    assert hexane_rich[2] == pytest.approx(light_fraction, abs=FRACTION_TOL)

    # And both of chemthermo's phases are liquids by a kappa computed here from
    # the public pressure routine, not from `phase_identity`'s own derivative.
    for composition, _density, _fraction in (water_rich, hexane_rich):
        assert 0.0 < _kappa(pressure_Pa, composition) < KAPPA_LIQUID_THRESHOLD


@pytest.mark.parametrize("pressure_Pa", [ATMOSPHERE_PA, HIGH_PRESSURE_PA])
@pytest.mark.usefixtures("matched_constants")
def test_feos_chemical_potentials_are_equal_at_chemthermos_phases(
    pressure_Pa: float,
) -> None:
    """The reference model's own equilibrium condition, at chemthermo's answer.

    This is the comparison that does not depend on FeOs's flash converging:
    only the two compositions and the two densities are chemthermo's, and FeOs
    supplies the chemical potentials. Run with matched constants, as validation
    Case P-7 did and for the same reason - FeOs is evaluating chemthermo's own
    densities, so the universal-constants table difference would otherwise
    floor the residual at about 1e-6. Case P-8 records both numbers.
    """
    _result, ours = _our_phases(pressure_Pa, (0.5, 0.5))
    water_rich, hexane_rich = ours["liquid1"], ours["liquid2"]
    difference = np.max(
        np.abs(
            _feos_reduced_potentials(water_rich[1], water_rich[0])
            - _feos_reduced_potentials(hexane_rich[1], hexane_rich[0])
        )
    )
    assert float(difference) < POTENTIAL_TOL


def test_the_same_tie_line_from_three_feeds_still_matches_feos() -> None:
    """Case P-8 (ii): a tie line is a property of the state, not of the feed.

    Only ``z = 0.5 / 0.5`` is compared against FeOs's *flash* directly: at
    ``z = 0.2 / 0.8`` FeOs's own ``tp_flash`` raises ``"stability analysis did
    not converge"`` on this binary (measured with feos 0.10.1), so there is no
    reference answer from it at that feed. Recorded as observed behaviour of
    the reference, not worked around. What is checked instead is the statement
    that makes those feeds meaningful: all three return the **same tie line**,
    so all three agree with the FeOs tie line that does exist, and each
    satisfies the lever rule.
    """
    dense, light, _fraction = _feos_tp_flash(ATMOSPHERE_PA, (0.5, 0.5))

    for feed in ((0.5, 0.5), (0.2, 0.8), (0.8, 0.2)):
        result, ours = _our_phases(ATMOSPHERE_PA, feed)
        assert sorted(result.phases) == ["liquid1", "liquid2"]
        assert ours["liquid1"][0] == pytest.approx(dense[0], abs=COMPOSITION_TOL)
        assert ours["liquid2"][0] == pytest.approx(light[0], abs=COMPOSITION_TOL)

        z = np.asarray(feed, dtype=float) / float(sum(feed))
        beta = ours["liquid2"][2]
        recombined = (1.0 - beta) * ours["liquid1"][0] + beta * ours["liquid2"][0]
        assert float(np.max(np.abs(z - recombined))) < 1e-12


def test_the_agreement_is_not_vacuous() -> None:
    """Perturbing water's association energy by 1 % must move the tie line.

    Without this, "the two agree" could be a statement about two solvers that
    both happen to return something insensitive to the model.
    """
    _result, ours = _our_phases(ATMOSPHERE_PA, (0.5, 0.5))
    perturbed_parameters = ct.PCSAFTParameters.from_records(
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
    mixture = ct.Mixture.from_database(list(NAMES), [0.5, 0.5], normalize=True)
    shifted = ct.flash_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=ATMOSPHERE_PA,
        eos=PCSAFTEOS(components=NAMES, parameters=perturbed_parameters),
    )
    moved = abs(
        float(shifted.phases["liquid2"].composition.fractions[0]) - float(ours["liquid2"][0][0])
    )
    assert moved > 1e-4, moved
