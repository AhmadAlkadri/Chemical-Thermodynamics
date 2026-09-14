"""Polymer/solvent **vapour-liquid** equilibrium against FeOs (validation Case P-14).

The companion of ``test_pcsaft_polymer_vs_feos.py``, one regime lower in
pressure: below n-pentane's saturation pressure at 453 K the polyethylene /
n-pentane equilibrium is vapour-liquid, the vapour's polymer mole fraction is
``exp(-450)``, and ``flash_tp`` reaches it only because ADR-0024 carries the
split in log mole numbers.

Three comparisons:

1. **FeOs's chemical potentials at chemthermo's converged phases and
   densities** - the reference implementation's own equilibrium condition,
   evaluated at chemthermo's answer. It does not depend on FeOs's flash
   converging, which on this system it does not (recorded in the ADR-0022
   module above: it raises or returns a degenerate pair at every pressure tried
   but one, all of them *higher* than these).
2. **A one-dimensional equal-fugacity solve written here**, with the vapour
   taken as exactly pure solvent, which is the limit the log-space answer is
   supposed to approach. One unknown, no part of the flash used.
3. **The polymer's own equal-fugacity condition in logarithms**, the equation
   that cannot be written in linear mole numbers at all.

``k_ij`` cannot be given to FeOs's PC-SAFT here (feos 0.10.1;
``EquationOfState.pcsaft`` raises "missing field ``k_ij``" for every
serialization tried), so comparison 1 runs at ``k_ij = 0`` on **both** sides
and the flash it checks is re-run at ``k_ij = 0``. Comparisons 2 and 3 use the
fitted ``k_ij = -0.006``, the value the headline numbers use.

The 42 universal constants of the 2001 dispersion term are the one input that
is genuinely not shared (chemthermo packages them as printed, ten figures; FeOs
hard-codes fourteen), so comparison 1 is run twice, as shipped and with FeOs's
table substituted in, and both numbers go into the ledger.

**The polymer parameters are not verified data** - see the provenance block of
``tests/fixtures/pcsaft/martini2009_polymers.json``.

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
PE_MW_G_MOL = 16400.0
PENTANE_MW_G_MOL = 72.146
KIJ = -0.006

#: The polymer row of ``tests/fixtures/pcsaft/martini2009_polymers.json``.
PE_SEGMENTS_PER_G = 0.0263
PE_SIGMA_A = 4.0217
PE_EPSILON_K_K = 247.5
#: n-pentane, Gross & Sadowski (2001) Table 1 - the packaged record's values.
PENTANE = (2.6896, 3.7729, 231.20)

#: The three vapour-liquid states, all below the solvent's saturation pressure.
PRESSURES_PA = (5.0e5, 1.0e6, 2.0e6)
#: The same three for the per-state comparisons, with the two that are *another
#: pressure on the same map* marked ``slow`` (see ``.agents/dev-contract.md``).
#: 1 MPa runs by default and is the state the ledger pins.
PRESSURE_PARAMS = (
    pytest.param(5.0e5, marks=pytest.mark.slow),
    pytest.param(1.0e6),
    pytest.param(2.0e6, marks=pytest.mark.slow),
)

_MOL_PER_M3 = si.MOL / si.METER**3

POTENTIAL_TOL = 1e-8
#: The melt's solvent mole fraction, flash against the one-dimensional solve.
COMPOSITION_TOL = 1e-10
LOG_FUGACITY_TOL = 1e-8


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


def _our_parameters() -> PCSAFTParameters:
    return PCSAFTParameters.from_records(
        [
            PCSAFTRecord(
                name="Polyethylene",
                segments_per_g=PE_SEGMENTS_PER_G,
                MW_g_mol=PE_MW_G_MOL,
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


def _feos_eos() -> EquationOfState:
    """The reference model, always with ``k_ij = 0`` - see the module docstring."""
    payloads = [
        {
            "identifier": {"name": "Polyethylene"},
            "molarweight": PE_MW_G_MOL,
            "m": PE_SEGMENTS_PER_G * PE_MW_G_MOL,
            "sigma": PE_SIGMA_A,
            "epsilon_k": PE_EPSILON_K_K,
        },
        {
            "identifier": {"name": "n-Pentane"},
            "molarweight": PENTANE_MW_G_MOL,
            "m": PENTANE[0],
            "sigma": PENTANE[1],
            "epsilon_k": PENTANE[2],
        },
    ]
    records = [PureRecord.from_json_str(json.dumps(payload)) for payload in payloads]
    return EquationOfState.pcsaft(Parameters.from_records(records))


def _mixture(weight_fraction: float = 0.05) -> ct.Mixture:
    polymer = weight_fraction / PE_MW_G_MOL
    solvent = (1.0 - weight_fraction) / PENTANE_MW_G_MOL
    total = polymer + solvent
    return ct.Mixture.from_components(
        [
            ct.Component.custom(
                "Polyethylene",
                mw_kg_per_mol=PE_MW_G_MOL / 1000.0,
                formula="(C2H4)n",
                volatile=False,
            ),
            ct.Component.from_database("n-Pentane"),
        ],
        [polymer / total, solvent / total],
        normalize=True,
    )


def _our_phases(pressure_Pa: float, kij: float) -> tuple[ct.FlashResult, dict[str, tuple]]:
    """``flash_tp`` plus, per phase, ``(composition, density)`` on that phase's root."""
    mixture = _mixture()
    eos = PCSAFTEOS(parameters=_our_parameters(), kij=kij)
    result = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, eos=eos)
    assert sorted(result.phases) == ["liquid", "vapor"]
    phases: dict[str, tuple] = {}
    for name, phase in result.phases.items():
        fractions = np.asarray(phase.composition.fractions, dtype=float)
        roots = eos.density_roots(
            mixture=mixture,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=pressure_Pa,
            composition=fractions.tolist(),
        )
        phases[name] = (fractions, roots[-1] if name == "liquid" else roots[0])
    return result, phases


def _feos_reduced_potentials(density: float, x: Sequence[float]) -> np.ndarray:
    """``mu_i / RT`` from FeOs at chemthermo's ``(T, rho, x)``, up to one constant.

    The omitted constant is the same for both phases of one flash, so it
    cancels in the difference and the comparison is the equality of chemical
    potentials evaluated by the reference implementation at the state
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
    mu_res = np.atleast_1d(
        np.asarray(state.chemical_potential(Contributions.Residual) / si.JOULE * si.MOL)
    )
    with np.errstate(divide="ignore"):
        return mu_res / factor - math.log(z_factor) + np.log(values)


# ---------------------------------------------------------------------------
# 1. FeOs at chemthermo's phases
# ---------------------------------------------------------------------------


@pytest.mark.usefixtures("matched_constants")
def test_feos_chemical_potentials_are_equal_at_the_vapour_liquid_phases() -> None:
    """Case P-14: the check that does not need FeOs's flash to converge.

    Only the two compositions and the two densities are chemthermo's. A polymer
    mole fraction of ``1e-196`` is an ordinary double and FeOs takes it as
    given; what is not an ordinary double is the *path* chemthermo took to it.
    """
    worst = 0.0
    for pressure_Pa in PRESSURES_PA:
        _result, ours = _our_phases(pressure_Pa, 0.0)
        difference = np.max(
            np.abs(
                _feos_reduced_potentials(*reversed(ours["liquid"]))
                - _feos_reduced_potentials(*reversed(ours["vapor"]))
            )
        )
        worst = max(worst, float(difference))
    # Measured 2.6e-12 with matched constants; 2.5e-07 as shipped, which is the
    # universal-constants table and not the model.
    assert worst < POTENTIAL_TOL


def test_the_feos_comparison_is_floored_by_the_constants_table() -> None:
    """The user-visible number, with the shipped ten-figure constants."""
    _result, ours = _our_phases(1.0e6, 0.0)
    difference = float(
        np.max(
            np.abs(
                _feos_reduced_potentials(*reversed(ours["liquid"]))
                - _feos_reduced_potentials(*reversed(ours["vapor"]))
            )
        )
    )
    assert difference < 1e-5


def test_the_flash_densities_reproduce_the_pressure_in_feos() -> None:
    """Both converged phases really are states of the reference model at ``P``.

    At 1 MPa only: the chemical-potential comparison above visits all three
    states, and it is evaluated at these same densities.
    """
    for pressure_Pa in (1.0e6,):
        _result, ours = _our_phases(pressure_Pa, 0.0)
        for _name, (x, density) in ours.items():
            state = State(
                _feos_eos(),
                temperature=TEMPERATURE_K * si.KELVIN,
                density=density * _MOL_PER_M3,
                composition=np.asarray(x, dtype=float),
            )
            assert float(state.pressure() / si.PASCAL) == pytest.approx(pressure_Pa, rel=1e-6)


# ---------------------------------------------------------------------------
# 2 and 3. An independent one-dimensional solve, at the fitted k_ij
# ---------------------------------------------------------------------------


def _ln_phi(bound: PCSAFTEOS, pressure_Pa: float, x: Sequence[float], *, vapor: bool) -> np.ndarray:
    roots = bound.density_roots(
        temperature_K=TEMPERATURE_K, pressure_Pa=pressure_Pa, composition=list(x)
    )
    density = roots[0] if vapor else roots[-1]
    return np.asarray(
        bound.ln_fugacity_coefficients(
            temperature_K=TEMPERATURE_K, density_mol_m3=density, composition=list(x)
        )
    )


@pytest.mark.parametrize("pressure_Pa", PRESSURE_PARAMS)
def test_the_melt_matches_a_one_dimensional_equal_fugacity_solve(pressure_Pa: float) -> None:
    """Case P-14: the same tie line from a solver written here, not reused.

    The flash's vapour is pure solvent to ``1e-196``, so the equilibrium melt
    must satisfy the single equation

        ln phi_s^V(pure solvent vapour) = ln x_s + ln phi_s^L(x)

    exactly. One unknown, carried as ``ln x_s``, solved by a finite-difference
    Newton started a thousandth away from the flash's answer so that it has to
    do real work. Nothing of the flash's stability test, Rachford-Rice,
    log-space stage or phase-count logic is used; only ``PCSAFTEOS``'s own
    ``(T, rho, x)`` interface and ``density_roots``.
    """
    result, _ours = _our_phases(pressure_Pa, KIJ)
    bound = PCSAFTEOS(
        components=("Polyethylene", "n-Pentane"), parameters=_our_parameters(), kij=KIJ
    )
    target = float(_ln_phi(bound, pressure_Pa, [0.0, 1.0], vapor=True)[1])

    def residual(ln_x: float) -> float:
        solvent = math.exp(ln_x)
        return (
            ln_x
            + float(_ln_phi(bound, pressure_Pa, [1.0 - solvent, solvent], vapor=False)[1])
            - target
        )

    ln_x = math.log(float(result.phases["liquid"].composition.fractions[1])) - 1e-3
    for _ in range(60):
        value = residual(ln_x)
        step = 1e-7
        slope = (residual(ln_x + step) - residual(ln_x - step)) / (2.0 * step)
        correction = -value / slope
        ln_x += correction
        if abs(correction) < 1e-14:
            break

    assert abs(residual(ln_x)) < LOG_FUGACITY_TOL
    ours = float(result.phases["liquid"].composition.fractions[1])
    # Measured: 2.5e-14 at 0.5 MPa, 5.6e-15 at 1 MPa, 1.2e-15 at 2 MPa.
    assert abs(ours - math.exp(ln_x)) < COMPOSITION_TOL

    # The flash's own vapour-phase solvent fugacity is the pure-solvent one.
    y = list(result.phases["vapor"].composition.fractions)
    flash_vapor = math.log(y[1]) + float(_ln_phi(bound, pressure_Pa, y, vapor=True)[1])
    assert abs(flash_vapor - target) < LOG_FUGACITY_TOL


@pytest.mark.parametrize("pressure_Pa", PRESSURE_PARAMS)
def test_the_polymers_own_equal_fugacity_condition_holds_in_logarithms(
    pressure_Pa: float,
) -> None:
    """``ln x_PE + ln phi_PE^L = ln y_PE + ln phi_PE^V`` with ``ln y_PE ~ -450``.

    This is the equation the linear split cannot write down: the right-hand
    side is ``-492`` to ``-530`` and its first term is not a representable mole
    fraction's logarithm by accident - it is one because the stage kept it as a
    logarithm throughout.
    """
    result, _ours = _our_phases(pressure_Pa, KIJ)
    bound = PCSAFTEOS(
        components=("Polyethylene", "n-Pentane"), parameters=_our_parameters(), kij=KIJ
    )
    x = list(result.phases["liquid"].composition.fractions)
    y = list(result.phases["vapor"].composition.fractions)
    ln_y_polymer = math.log(y[0]) if y[0] > 0.0 else float(result.diagnostics["log_space_ln_x_min"])
    melt = math.log(x[0]) + float(_ln_phi(bound, pressure_Pa, x, vapor=False)[0])
    vapor = ln_y_polymer + float(_ln_phi(bound, pressure_Pa, y, vapor=True)[0])
    assert melt < -400.0  # the number itself, so a sign or scale change is visible
    assert abs(melt - vapor) < LOG_FUGACITY_TOL
