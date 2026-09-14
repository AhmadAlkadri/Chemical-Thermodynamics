"""Polymer/solvent PC-SAFT against FeOs (validation Cases P-12 and P-13).

**FeOs** (feos-org/feos, MIT OR Apache-2.0) is the independent route, as it has
been since ADR-0018: the same Gross & Sadowski model written in Rust with every
derivative by automatic differentiation. It carries no polymer parameter file,
but it takes a segment number directly, so the same
``m = (m/M) * Mw = 431.32`` this package derives from the segments-per-mass
record can be handed to it as a plain ``m``.

Three comparisons, and they are deliberately different in kind:

1. **Properties** at states chemthermo's own density-root solver found -
   ``A^res/RT``, ``Z`` and ``ln phi_i`` for the pure melt at two molar masses
   and for three polymer concentrations in n-pentane (Case P-12).
2. **FeOs's own TP flash**, where it converges: a tie line from another
   *solver*, not only another model (Case P-13).
3. **FeOs's chemical potentials evaluated at chemthermo's converged phases**
   (the Case P-8 route), which does not depend on FeOs's flash converging at
   all.

Two honest limitations are recorded here rather than worked around.

*FeOs's flash mostly fails on this system.* At 453 K with a 5 wt% feed it
raises ``RuntimeError: `rachford_rice` encountered illegal values during the
iteration`` at 5 and 8 MPa, returns a degenerate pair (both phases equal, a
vapour fraction of exactly 0.5) at 3 MPa, and says "No phase split according to
stability analysis" at and above 11 MPa. It converges at 10 MPa, and that is
the one state comparison 2 uses. This is a fact about the reference
implementation at this extreme size ratio; it is not a disagreement.

*``k_ij`` cannot be given to FeOs's PC-SAFT here.* ``BinaryRecord`` and
``Parameters.from_records`` both accept one, but ``EquationOfState.pcsaft``
then raises ``RuntimeError: missing field `k_ij``` in feos 0.10.1 for every
serialization tried (a bare float, ``{"k_ij": x}``, ``{"k_ij": [x]}``, with and
without ``l_ij``). Every comparison below is therefore run at ``k_ij = 0`` on
**both** sides - a code cross-check of the same model, at a state of the same
kind. chemthermo's headline answer for this system uses the fitted
``k_ij = -0.006``; that one is checked against an independently written
two-equation Newton solve in ``tests/test_pcsaft_polymer.py``.

The 42 universal constants of the 2001 dispersion term are the one input that
is genuinely not shared (chemthermo packages them as printed, ten figures; FeOs
hard-codes fourteen), so every comparison is run twice, as shipped and with
FeOs's table substituted in, and both numbers go into the ledger.

Skipped when ``feos`` is not installed (``pip install -e ".[validation]"``).
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
from chemthermo.parameters import PCSAFTParameters, PCSAFTRecord

feos = pytest.importorskip("feos")
si = pytest.importorskip("si_units")

from feos import (  # noqa: E402
    Contributions,
    EquationOfState,
    Parameters,
    PureRecord,
    State,
)

TEMPERATURE_K = 453.0
PENTANE_MW_G_MOL = 72.146
PE_MW_G_MOL = 16400.0

#: The polymer row of ``tests/fixtures/pcsaft/martini2009_polymers.json``,
#: written out here so the reference model is built from this file. **As
#: tabulated by Martini et al. (2009) citing Gross & Sadowski (2002); not
#: verified against the primary table** - see the fixture's provenance block.
PE_SEGMENTS_PER_G = 0.0263
PE_SIGMA_A = 4.0217
PE_EPSILON_K_K = 247.5
#: n-pentane, Gross & Sadowski (2001) Table 1 - the packaged record's values.
PENTANE = (2.6896, 3.7729, 231.20)

_MOL_PER_M3 = si.MOL / si.METER**3

#: Asserted tolerances (validation Cases P-12, P-13), with matched constants.
PROPERTY_TOL = 1e-10
POTENTIAL_TOL = 1e-8
#: The FeOs-flash comparison runs at 10 MPa, which is within 8% of this
#: system's cloud point at ``k_ij = 0`` (10.77 MPa), so both solvers are
#: ill-conditioned there; see
#: ``test_the_feos_flash_tie_line_matches_where_feos_converges``.
FLASH_COMPOSITION_TOL = 1e-8


def _feos_universal_constants() -> tuple[np.ndarray, np.ndarray]:
    """FeOs's fourteen-figure ``a`` and ``b`` tables, from the Case P-6 module.

    Loaded **by path**: CI runs the ``pytest`` console script, which does not
    put the working directory on ``sys.path`` (see ``.agents/dev-contract.md``).
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


# ---------------------------------------------------------------------------
# The two models, built from the same numbers
# ---------------------------------------------------------------------------


def _our_parameters(mw_g_mol: float) -> PCSAFTParameters:
    return PCSAFTParameters.from_records(
        [
            PCSAFTRecord(
                name="Polyethylene",
                segments_per_g=PE_SEGMENTS_PER_G,
                MW_g_mol=mw_g_mol,
                sigma_A=PE_SIGMA_A,
                epsilon_k_K=PE_EPSILON_K_K,
            ),
            PCSAFTRecord(
                name="n-Pentane",
                m=PENTANE[0],
                sigma_A=PENTANE[1],
                epsilon_k_K=PENTANE[2],
                MW_g_mol=PENTANE_MW_G_MOL,
            ),
        ]
    )


def _feos_eos(names: Sequence[str], mw_g_mol: float) -> EquationOfState:
    """The reference model, always with ``k_ij = 0`` - see the module docstring."""
    payloads = {
        "Polyethylene": {
            "identifier": {"name": "Polyethylene"},
            "molarweight": mw_g_mol,
            # FeOs is given the *derived* segment number, which is what makes
            # this a check of the segments-per-mass convention as well.
            "m": PE_SEGMENTS_PER_G * mw_g_mol,
            "sigma": PE_SIGMA_A,
            "epsilon_k": PE_EPSILON_K_K,
        },
        "n-Pentane": {
            "identifier": {"name": "n-Pentane"},
            "molarweight": PENTANE_MW_G_MOL,
            "m": PENTANE[0],
            "sigma": PENTANE[1],
            "epsilon_k": PENTANE[2],
        },
    }
    records = [PureRecord.from_json_str(json.dumps(payloads[name])) for name in names]
    return EquationOfState.pcsaft(Parameters.from_records(records))


def _polymer_component(mw_g_mol: float) -> ct.Component:
    return ct.Component.custom(
        "Polyethylene",
        mw_kg_per_mol=mw_g_mol / 1000.0,
        formula="(C2H4)n",
        volatile=False,
    )


def _weight_to_mole_fractions(weight_fraction: float, mw_g_mol: float) -> list[float]:
    polymer = weight_fraction / mw_g_mol
    solvent = (1.0 - weight_fraction) / PENTANE_MW_G_MOL
    total = polymer + solvent
    return [polymer / total, solvent / total]


def _mixture(weight_fraction: float, mw_g_mol: float = PE_MW_G_MOL) -> ct.Mixture:
    return ct.Mixture.from_components(
        [_polymer_component(mw_g_mol), ct.Component.from_database("n-Pentane")],
        _weight_to_mole_fractions(weight_fraction, mw_g_mol),
        normalize=True,
    )


def _feos_report(
    names: Sequence[str], mw_g_mol: float, density: float, x: Sequence[float] | np.ndarray
) -> tuple[float, float, np.ndarray]:
    """``(A^res/RT, Z, ln phi)`` from FeOs at ``(T, rho, x)``."""
    state = State(
        _feos_eos(names, mw_g_mol),
        temperature=TEMPERATURE_K * si.KELVIN,
        density=density * _MOL_PER_M3,
        composition=np.asarray(x, dtype=float),
    )
    factor = R_J_PER_MOL_K * TEMPERATURE_K
    a_res = sum(
        (value / si.JOULE * si.MOL) / factor
        for _label, value in state.residual_molar_helmholtz_energy_contributions()
    )
    z_factor = (state.pressure() / si.PASCAL) / (density * factor)
    mu_res = np.atleast_1d(
        np.asarray(state.chemical_potential(Contributions.Residual) / si.JOULE * si.MOL)
    )
    return float(a_res), float(z_factor), mu_res / factor - math.log(z_factor)


def _feos_reduced_potentials(
    mw_g_mol: float, density: float, x: Sequence[float] | np.ndarray
) -> np.ndarray:
    """``mu_i / RT`` from FeOs at chemthermo's ``(T, rho, x)``, up to a constant.

    The omitted constant is the same for both phases of one flash, so it
    cancels in the difference and the comparison is the equality of chemical
    potentials evaluated by the reference implementation at the state
    chemthermo converged on.
    """
    values = np.asarray(x, dtype=float)
    _a_res, z_factor, ln_phi = _feos_report(
        ("Polyethylene", "n-Pentane"), mw_g_mol, density, values
    )
    del z_factor
    return ln_phi + np.log(values)


# ---------------------------------------------------------------------------
# Case P-12: properties
# ---------------------------------------------------------------------------

#: ``(label, names, mw, weight fraction or None for the pure melt, P/Pa)``.
#: Densities are not written down: each state's density is the **model's own**
#: liquid root at that pressure, which is what makes this a check of the root
#: solver as well as of the residual Helmholtz energy.
PROPERTY_STATES: tuple[tuple[str, tuple[str, ...], float, float | None, float], ...] = (
    ("melt-16400-1MPa", ("Polyethylene",), 16400.0, None, 1.0e6),
    ("melt-16400-10MPa", ("Polyethylene",), 16400.0, None, 1.0e7),
    ("melt-16400-30MPa", ("Polyethylene",), 16400.0, None, 3.0e7),
    ("melt-53000-1MPa", ("Polyethylene",), 53000.0, None, 1.0e6),
    ("melt-53000-10MPa", ("Polyethylene",), 53000.0, None, 1.0e7),
    ("melt-53000-30MPa", ("Polyethylene",), 53000.0, None, 3.0e7),
    ("mix-05wt-10MPa", ("Polyethylene", "n-Pentane"), 16400.0, 0.05, 1.0e7),
    ("mix-10wt-10MPa", ("Polyethylene", "n-Pentane"), 16400.0, 0.10, 1.0e7),
    ("mix-15wt-10MPa", ("Polyethylene", "n-Pentane"), 16400.0, 0.15, 1.0e7),
    ("mix-05wt-15MPa", ("Polyethylene", "n-Pentane"), 16400.0, 0.05, 1.5e7),
    ("mix-15wt-15MPa", ("Polyethylene", "n-Pentane"), 16400.0, 0.15, 1.5e7),
)


def _state_composition(
    names: Sequence[str], mw_g_mol: float, weight_fraction: float | None
) -> list[float]:
    if weight_fraction is None:
        return [1.0]
    del names
    return _weight_to_mole_fractions(weight_fraction, mw_g_mol)


@pytest.mark.usefixtures("matched_constants")
def test_the_polymer_properties_match_feos_with_matched_constants() -> None:
    """Case P-12: ``A^res/RT``, ``Z`` and ``ln phi`` over eleven states, to 1e-10."""
    worst = [0.0, 0.0, 0.0]
    for label, names, mw_g_mol, weight_fraction, pressure_Pa in PROPERTY_STATES:
        eos = PCSAFTEOS(components=names, parameters=_our_parameters(mw_g_mol))
        x = _state_composition(names, mw_g_mol, weight_fraction)
        density = eos.density_roots(
            temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=x
        )[-1]

        ours = (
            eos.residual_helmholtz(
                temperature_K=TEMPERATURE_K, volume_m3=1.0 / density, composition=x
            ),
            eos.compressibility_factor(
                temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=x
            ),
            np.asarray(
                eos.ln_fugacity_coefficients(
                    temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=x
                )
            ),
        )
        theirs = _feos_report(names, mw_g_mol, density, x)

        differences = (
            abs(ours[0] - theirs[0]),
            abs(ours[1] - theirs[1]),
            float(np.max(np.abs(ours[2] - theirs[2]))),
        )
        for index, difference in enumerate(differences):
            assert difference < PROPERTY_TOL, (label, index, difference)
            worst[index] = max(worst[index], difference)

    # Measured, recorded in validation Case P-12; as shipped the same three are
    # 2.9e-07, 1.1e-06 and 1.4e-06, which is the universal-constants table.
    assert worst[0] < 1e-11
    assert worst[1] < 1e-11
    assert worst[2] < 1e-11


def test_the_properties_still_agree_with_the_shipped_constants() -> None:
    """The user-visible number: agreement is floored by the constants table, not the model."""
    names = ("Polyethylene", "n-Pentane")
    eos = PCSAFTEOS(components=names, parameters=_our_parameters(PE_MW_G_MOL))
    x = _weight_to_mole_fractions(0.05, PE_MW_G_MOL)
    density = eos.density_roots(temperature_K=TEMPERATURE_K, pressure_Pa=1.0e7, composition=x)[-1]
    ours = np.asarray(
        eos.ln_fugacity_coefficients(
            temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=x
        )
    )
    _a_res, _z, theirs = _feos_report(names, PE_MW_G_MOL, density, x)
    assert float(np.max(np.abs(ours - theirs))) < 1e-5


def test_the_melt_density_is_the_same_state_in_both_implementations() -> None:
    """chemthermo's liquid root reproduces the pressure FeOs computes there."""
    for mw_g_mol in (16400.0, 53000.0):
        eos = PCSAFTEOS(components=("Polyethylene",), parameters=_our_parameters(mw_g_mol))
        for pressure_Pa in (1.0e6, 1.0e7, 3.0e7):
            density = eos.density_roots(
                temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=[1.0]
            )[-1]
            state = State(
                _feos_eos(("Polyethylene",), mw_g_mol),
                temperature=TEMPERATURE_K * si.KELVIN,
                density=density * _MOL_PER_M3,
                composition=np.asarray([1.0]),
            )
            assert float(state.pressure() / si.PASCAL) == pytest.approx(pressure_Pa, rel=1e-5)


# ---------------------------------------------------------------------------
# Case P-13: equilibrium
# ---------------------------------------------------------------------------


def _our_phases(
    pressure_Pa: float, kij: float
) -> tuple[ct.FlashResult, dict[str, tuple[np.ndarray, float, float]]]:
    """``flash_tp`` plus, per phase, ``(composition, density, phase fraction)``."""
    mixture = _mixture(0.05)
    eos = PCSAFTEOS(parameters=_our_parameters(PE_MW_G_MOL), kij=kij)
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


@pytest.mark.usefixtures("matched_constants")
def test_feos_chemical_potentials_are_equal_at_chemthermos_phases() -> None:
    """Case P-13: the check that does not need FeOs's flash to converge.

    Only the two compositions and the two densities are chemthermo's; FeOs
    supplies the chemical potentials. Run at 8 MPa, where FeOs's own flash
    raises, and at ``k_ij = 0`` because FeOs cannot be given one here.
    """
    _result, ours = _our_phases(8.0e6, 0.0)
    (first, second) = list(ours)
    difference = np.max(
        np.abs(
            _feos_reduced_potentials(PE_MW_G_MOL, ours[first][1], ours[first][0])
            - _feos_reduced_potentials(PE_MW_G_MOL, ours[second][1], ours[second][0])
        )
    )
    assert float(difference) < POTENTIAL_TOL


def test_the_feos_flash_tie_line_matches_where_feos_converges() -> None:
    """Case P-13: FeOs's own ``tp_flash``, at the one pressure where it converges.

    10 MPa is within 8% of this system's ``k_ij = 0`` cloud point (10.77 MPa),
    so both solvers are working near a plait point: chemthermo's converged pair
    has an equal-fugacity residual of 1.06e-07 in chemthermo's own model and
    FeOs's has 9.07e-06 in the same model, i.e. chemthermo's is the tighter
    stationary point by about 86x. The tolerances below are set by that
    conditioning, not by a model difference - the property comparison above
    agrees to 1e-11 at the same kind of state.
    """
    result, ours = _our_phases(1.0e7, 0.0)
    assert sorted(result.phases) == ["liquid1", "liquid2"]
    assert result.vapor_fraction is None
    assert result.diagnostics["phase_regime"] == "LLE"

    state = State(
        _feos_eos(("Polyethylene", "n-Pentane"), PE_MW_G_MOL),
        temperature=TEMPERATURE_K * si.KELVIN,
        pressure=1.0e7 * si.PASCAL,
        composition=np.asarray(_weight_to_mole_fractions(0.05, PE_MW_G_MOL)),
    )
    equilibrium = state.tp_flash()
    feos_phases = [
        (
            np.atleast_1d(np.asarray(equilibrium.liquid.molefracs, dtype=float)),
            float(equilibrium.liquid.density / _MOL_PER_M3),
            1.0 - float(equilibrium.vapor_phase_fraction),
        ),
        (
            np.atleast_1d(np.asarray(equilibrium.vapor.molefracs, dtype=float)),
            float(equilibrium.vapor.density / _MOL_PER_M3),
            float(equilibrium.vapor_phase_fraction),
        ),
    ]
    # Both of FeOs's phases are liquids here, whatever its container calls them.
    assert all(density > 5000.0 for _x, density, _f in feos_phases)

    by_polymer = sorted(ours, key=lambda name: ours[name][0][0])
    feos_sorted = sorted(feos_phases, key=lambda entry: entry[0][0])
    for name, (their_x, their_density, their_fraction) in zip(by_polymer, feos_sorted):
        our_x, our_density, our_fraction = ours[name]
        assert float(np.max(np.abs(our_x - their_x))) < FLASH_COMPOSITION_TOL
        assert our_density == pytest.approx(their_density, rel=1e-6)
        assert our_fraction == pytest.approx(their_fraction, abs=1e-5)

    # The polymer partition itself, which is the number the asymmetry threatens:
    # measured 5.8e-13 in the polymer-rich phase and 1.2e-09 in the solvent-rich
    # one (matched constants), both well above 1e-12 so both are real digits.
    polymer_rich, solvent_rich = by_polymer[-1], by_polymer[0]
    assert ours[polymer_rich][0][0] == pytest.approx(feos_sorted[-1][0][0], rel=1e-8)
    assert ours[solvent_rich][0][0] == pytest.approx(feos_sorted[0][0][0], rel=1e-4)


@pytest.mark.slow  # the same comparison at the other end of the two-phase window
def test_feos_refuses_this_flash_where_chemthermo_returns_a_verified_split() -> None:
    """Recorded, not asserted against: FeOs's flash is not a reference here.

    At 5 and 8 MPa FeOs's ``tp_flash`` raises; at 3 MPa it returns a degenerate
    pair (the two phases equal to nine figures, vapour fraction exactly 0.5).
    chemthermo returns a verified two-liquid split at all three. This test
    exists so that a future FeOs release changing any of it is noticed.
    """
    composition = np.asarray(_weight_to_mole_fractions(0.05, PE_MW_G_MOL))
    raised = 0
    degenerate = 0
    for pressure_Pa in (3.0e6, 5.0e6, 8.0e6):
        state = State(
            _feos_eos(("Polyethylene", "n-Pentane"), PE_MW_G_MOL),
            temperature=TEMPERATURE_K * si.KELVIN,
            pressure=pressure_Pa * si.PASCAL,
            composition=composition,
        )
        try:
            equilibrium = state.tp_flash()
        except RuntimeError:
            raised += 1
        else:
            liquid = np.atleast_1d(np.asarray(equilibrium.liquid.molefracs, dtype=float))
            vapor = np.atleast_1d(np.asarray(equilibrium.vapor.molefracs, dtype=float))
            if float(np.max(np.abs(liquid - vapor))) < 1e-9:
                degenerate += 1

        result, _ours = _our_phases(pressure_Pa, 0.0)
        assert sorted(result.phases) == ["liquid1", "liquid2"]
        assert float(result.diagnostics["delta_g_split_rt"]) < 0.0
        assert result.diagnostics["post_split_status"] == "stable"

    assert raised == 2
    assert degenerate == 1
