"""PC-SAFT association against FeOs, validation Cases P-6 and P-7 (ADR-0018).

FeOs (feos-org/feos, MIT OR Apache-2.0) is an independent Rust implementation
of the same Gross & Sadowski model. Like teqp it obtains **every** derivative
by automatic differentiation (``num-dual``) of one hand-written Helmholtz
energy, where chemthermo writes them analytically, and unlike teqp it
implements the **association** term - teqp's ``PCSAFT`` kind does not, which
is why this file exists alongside ``test_pcsaft_vs_teqp.py`` instead of
replacing it.

Skipped when ``feos`` is not installed (``pip install -e ".[validation]"``).

One shared input is deliberately not shared: the 42 universal constants
--------------------------------------------------------------------------
chemthermo packages Gross & Sadowski (2001) Table 1 as **printed**, to ten
significant figures, which is also what teqp uses - the two agree on
``A^res/RT`` to 4e-15. FeOs hard-codes the same constants to **fourteen**
figures (``crates/feos/src/pcsaft/eos/dispersion.rs``), i.e. an extended-
precision variant of the published table; the two tables differ by up to
4.8e-09 termwise. That difference is an input, not an implementation
difference, and it floors any comparison of the *dispersion* term at about
1e-9 relative. So every dispersion-dependent quantity is compared twice:

* as shipped, where the residual difference is the constants and is recorded
  rather than asserted tightly, and
* with FeOs's own constants substituted into chemthermo, where the whole
  model - association included - comes back to 1.3e-11 or better.

The **association** term does not depend on those constants and matches to
round-off either way, which is the comparison this slice is actually about.
"""

from __future__ import annotations

import json
import math
from typing import Iterable, Mapping, Sequence

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
    PhaseEquilibrium,
    PureRecord,
    State,
)

# ---------------------------------------------------------------------------
# Parameters, written out here so the reference model is built from this file
# ---------------------------------------------------------------------------

#: ``name -> (MW, m, sigma/A, eps/k in K, kappa^AB or None, eps^AB/k in K)``.
#: The associating five are Gross & Sadowski (2002) Table 1 (2B scheme,
#: na = nb = 1); n-hexane is the 2001 non-associating record.
PARAMETERS: dict[str, tuple[float, float, float, float, float | None, float]] = {
    "Water": (18.015, 1.0656, 3.0007, 366.51, 0.034868, 2500.7),
    "Ethanol": (46.069, 2.3827, 3.1771, 198.24, 0.032384, 2653.4),
    "Methanol": (32.042, 1.5255, 3.2300, 188.90, 0.035176, 2899.5),
    "n-Hexane": (86.177, 3.0576, 3.7983, 236.77, None, 0.0),
}

#: FeOs's own universal constants (its ``dispersion.rs``), fourteen figures.
FEOS_A_UNIVERSAL = np.array(
    [
        [
            0.91056314451539,
            0.63612814494991,
            2.68613478913903,
            -26.5473624914884,
            97.7592087835073,
            -159.591540865600,
            91.2977740839123,
        ],
        [
            -0.30840169182720,
            0.18605311591713,
            -2.50300472586548,
            21.4197936296668,
            -65.2558853303492,
            83.3186804808856,
            -33.7469229297323,
        ],
        [
            -0.09061483509767,
            0.45278428063920,
            0.59627007280101,
            -1.72418291311787,
            -4.13021125311661,
            13.7766318697211,
            -8.67284703679646,
        ],
    ],
    dtype=float,
)
FEOS_B_UNIVERSAL = np.array(
    [
        [
            0.72409469413165,
            2.23827918609380,
            -4.00258494846342,
            -21.00357681484648,
            26.8556413626615,
            206.5513384066188,
            -355.60235612207947,
        ],
        [
            -0.57554980753450,
            0.69950955214436,
            3.89256733895307,
            -17.21547164777212,
            192.6722644652495,
            -161.8264616487648,
            -165.2076934555607,
        ],
        [
            0.09768831158356,
            -0.25575749816100,
            -9.15585615297321,
            20.64207597439724,
            -38.80443005206285,
            93.6267740770146,
            -29.66690558514725,
        ],
    ],
    dtype=float,
)

#: Asserted tolerance for a comparison in which the universal constants match,
#: and for the terms that do not depend on them at all. Measured worst with
#: matched constants: 1.3e-11 (in ln phi) over all eighteen states below.
TIGHT = 1e-10
#: Asserted tolerance for an energy-like dispersion-dependent quantity as
#: shipped, where the two tables of universal constants differ. Measured worst:
#: 7.1e-10 in ``a_disp`` and ``A^res/RT``, 3.9e-09 in ``Z``.
DISPERSION_LIMITED = 1e-8
#: The same, for ``ln phi``. A dilute component's ``ln phi`` amplifies the
#: constants difference; measured worst 1.7e-06 (water/n-hexane 0.9/0.1 at
#: 298 K, where the hexane mole fraction is 0.1 in a water-like liquid).
LN_PHI_LIMITED = 1e-5


@pytest.fixture
def matched_constants(monkeypatch: pytest.MonkeyPatch) -> None:
    """Give chemthermo FeOs's universal constants for the duration of a test."""
    monkeypatch.setattr(pcsaft_module, "A_UNIVERSAL", FEOS_A_UNIVERSAL)
    monkeypatch.setattr(pcsaft_module, "B_UNIVERSAL", FEOS_B_UNIVERSAL)


# ---------------------------------------------------------------------------
# The FeOs side
# ---------------------------------------------------------------------------


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
    """Build the FeOs PC-SAFT model, always with ``k_ij = 0``.

    ``Parameters.from_records`` takes an optional list of ``BinaryRecord``s and
    defaults it to an empty list, which FeOs reads as no binary interaction -
    i.e. exactly ``k_ij = 0``, the same default ``PCSAFTEOS`` has.
    """
    records = [_pure_record(name) for name in names]
    if len(records) == 1:
        return EquationOfState.pcsaft(Parameters.new_pure(records[0]))
    return EquationOfState.pcsaft(Parameters.from_records(records))


def _feos_state(names: Sequence[str], temperature: float, density: float, x: Sequence[float]):
    eos = _feos_eos(names)
    if len(names) == 1:
        return State(eos, temperature=temperature * si.KELVIN, density=density * _MOL_PER_M3)
    return State(
        eos,
        temperature=temperature * si.KELVIN,
        density=density * _MOL_PER_M3,
        composition=np.asarray(x, dtype=float),
    )


_MOL_PER_M3 = si.MOL / si.METER**3


def _feos_report(
    names: Sequence[str], temperature: float, density: float, x: Sequence[float]
) -> tuple[dict[str, float], float, np.ndarray]:
    """Return ``(contributions in A/RT, Z, ln phi)`` from FeOs at ``(T, rho, x)``."""
    state = _feos_state(names, temperature, density, x)
    factor = R_J_PER_MOL_K * temperature
    contributions = {
        label: (value / si.JOULE * si.MOL) / factor
        for label, value in state.residual_molar_helmholtz_energy_contributions()
    }
    pressure = state.pressure() / si.PASCAL
    z_factor = pressure / (density * factor)
    mu_res = np.atleast_1d(state.chemical_potential(Contributions.Residual) / si.JOULE * si.MOL)
    return contributions, float(z_factor), mu_res / factor - math.log(z_factor)


# ---------------------------------------------------------------------------
# Case P-6: term by term
# ---------------------------------------------------------------------------

#: ``(label, components, x, T/K, rho/(mol/m^3))``. Four states per pure
#: associating fluid - gas-like, two liquid-like, high temperature - plus the
#: mixtures of Case P-6's second half.
PURE_STATES: list[tuple[str, tuple[str, ...], list[float], float, float]] = [
    ("water 300 K / 55000 (liquid-like)", ("Water",), [1.0], 300.0, 55000.0),
    ("water 350 K / 50000 (liquid-like)", ("Water",), [1.0], 350.0, 50000.0),
    ("water 373.15 K / 100 (gas-like)", ("Water",), [1.0], 373.15, 100.0),
    ("water 550 K / 41241 (high T, liquid-like)", ("Water",), [1.0], 550.0, 41241.1907),
    ("water 550 K / 1250 (high T, gas-like)", ("Water",), [1.0], 550.0, 1250.012),
    ("ethanol 300 K / 17000 (liquid-like)", ("Ethanol",), [1.0], 300.0, 17000.0),
    ("ethanol 350 K / 15928 (liquid-like)", ("Ethanol",), [1.0], 350.0, 15928.6188),
    ("ethanol 450 K / 691 (gas-like)", ("Ethanol",), [1.0], 450.0, 691.1054),
    ("ethanol 500 K / 9929 (high T)", ("Ethanol",), [1.0], 500.0, 9929.2725),
]

MIXTURE_STATES: list[tuple[str, tuple[str, ...], list[float], float, float]] = [
    ("water/ethanol 0.2/0.8 320 K / 20000", ("Water", "Ethanol"), [0.2, 0.8], 320.0, 20000.0),
    ("water/ethanol 0.5/0.5 320 K / 25401", ("Water", "Ethanol"), [0.5, 0.5], 320.0, 25401.3411),
    ("water/ethanol 0.8/0.2 320 K / 36508", ("Water", "Ethanol"), [0.8, 0.2], 320.0, 36507.9501),
    (
        "water/ethanol 0.5/0.5 400 K / 94.6 (gas-like)",
        ("Water", "Ethanol"),
        [0.5, 0.5],
        400.0,
        94.6239,
    ),
    ("water/ethanol 0.2/0.8 351 K / 18575", ("Water", "Ethanol"), [0.2, 0.8], 351.0, 18574.9573),
    ("water/n-hexane 0.3/0.7 298 K / 10229", ("Water", "n-Hexane"), [0.3, 0.7], 298.15, 10229.0998),
    ("water/n-hexane 0.9/0.1 298 K / 33473", ("Water", "n-Hexane"), [0.9, 0.1], 298.15, 33472.7843),
    (
        "water/n-hexane 0.5/0.5 400 K / 93.9 (gas-like)",
        ("Water", "n-Hexane"),
        [0.5, 0.5],
        400.0,
        93.9327,
    ),
    (
        "methanol/water/n-hexane 320 K / 16685",
        ("Methanol", "Water", "n-Hexane"),
        [0.3, 0.4, 0.3],
        320.0,
        16684.835,
    ),
]

ALL_STATES = PURE_STATES + MIXTURE_STATES


_StateSpec = tuple[str, tuple[str, ...], list[float], float, float]


def _ids(states: Iterable[_StateSpec]) -> list[str]:
    return [state[0] for state in states]


def _number(diagnostics: Mapping[str, object], key: str) -> float:
    """Return one diagnostics entry as a float.

    ``FlashResult.diagnostics`` is a free-form mapping, so a bare
    ``diagnostics[key] < 1e-8`` is untyped; this narrows it in one place.
    """
    value = diagnostics[key]
    assert isinstance(value, (int, float)) and not isinstance(value, bool), (key, value)
    return float(value)


@pytest.mark.parametrize("state", ALL_STATES, ids=_ids(ALL_STATES))
@pytest.mark.usefixtures("matched_constants")
def test_every_term_matches_feos(
    state: tuple[str, tuple[str, ...], list[float], float, float],
) -> None:
    """Hard sphere, hard chain, dispersion, association, total, Z and ln phi.

    Run with FeOs's universal constants, so a disagreement here is a
    disagreement about the *model*, not about a shared table's last digits.
    """
    _, names, x, temperature, density = state
    eos = PCSAFTEOS(components=names)
    contributions, z_reference, ln_phi_reference = _feos_report(names, temperature, density, x)

    internal = eos._state(temperature, density, x)
    # FeOs splits the hard-chain reference into 'Hard Sphere' + 'Hard Chain';
    # chemthermo reports their sum as a_hc (Eq. A.4).
    assert internal.a_hc == pytest.approx(
        contributions["Hard Sphere"] + contributions["Hard Chain"], rel=TIGHT, abs=TIGHT
    )
    assert internal.a_disp == pytest.approx(contributions["Dispersion"], rel=TIGHT, abs=TIGHT)
    assert internal.a_assoc == pytest.approx(contributions["Association"], rel=TIGHT, abs=TIGHT)

    a_res = eos.residual_helmholtz(
        temperature_K=temperature, volume_m3=1.0 / density, composition=x
    )
    assert a_res == pytest.approx(sum(contributions.values()), rel=TIGHT, abs=TIGHT)

    z_factor = eos.compressibility_factor(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    assert z_factor == pytest.approx(z_reference, rel=TIGHT, abs=TIGHT)

    ln_phi = eos.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    np.testing.assert_allclose(ln_phi, ln_phi_reference, rtol=TIGHT, atol=TIGHT)


@pytest.mark.parametrize("state", ALL_STATES, ids=_ids(ALL_STATES))
def test_the_association_term_matches_feos_with_the_published_constants(
    state: tuple[str, tuple[str, ...], list[float], float, float],
) -> None:
    """As shipped. The association term is independent of the a/b table.

    So it must still match to ``TIGHT`` here, while the dispersion-dependent
    quantities are only asserted at :data:`CONSTANTS_LIMITED` - the floor the
    two tables of universal constants impose.
    """
    _, names, x, temperature, density = state
    eos = PCSAFTEOS(components=names)
    contributions, z_reference, ln_phi_reference = _feos_report(names, temperature, density, x)

    internal = eos._state(temperature, density, x)
    assert internal.a_assoc == pytest.approx(contributions["Association"], rel=TIGHT, abs=TIGHT)
    assert internal.a_hc == pytest.approx(
        contributions["Hard Sphere"] + contributions["Hard Chain"], rel=TIGHT, abs=TIGHT
    )
    assert internal.a_disp == pytest.approx(
        contributions["Dispersion"], rel=DISPERSION_LIMITED, abs=DISPERSION_LIMITED
    )
    a_res = eos.residual_helmholtz(
        temperature_K=temperature, volume_m3=1.0 / density, composition=x
    )
    assert a_res == pytest.approx(
        sum(contributions.values()), rel=DISPERSION_LIMITED, abs=DISPERSION_LIMITED
    )
    z_factor = eos.compressibility_factor(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    assert z_factor == pytest.approx(z_reference, rel=DISPERSION_LIMITED, abs=DISPERSION_LIMITED)
    ln_phi = eos.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    np.testing.assert_allclose(ln_phi, ln_phi_reference, rtol=LN_PHI_LIMITED, atol=LN_PHI_LIMITED)


def test_the_two_tables_of_universal_constants_really_do_differ() -> None:
    """Guards the explanation above: if they agreed, the two tests are the same."""
    difference = max(
        float(np.max(np.abs(FEOS_A_UNIVERSAL - pcsaft_module.A_UNIVERSAL))),
        float(np.max(np.abs(FEOS_B_UNIVERSAL - pcsaft_module.B_UNIVERSAL))),
    )
    assert 1e-10 < difference < 1e-7
    # ... and that they agree to the ten figures the paper prints.
    np.testing.assert_allclose(FEOS_A_UNIVERSAL, pcsaft_module.A_UNIVERSAL, rtol=0.0, atol=5e-9)
    np.testing.assert_allclose(FEOS_B_UNIVERSAL, pcsaft_module.B_UNIVERSAL, rtol=0.0, atol=5e-9)


def test_the_association_strength_uses_sigma_cubed_and_not_d_cubed() -> None:
    """Settles Eq. (3)'s prefactor numerically; see ``_pcsaft_association``.

    Both spellings are in circulation. This test recomputes the pure-water 2B
    association term from scratch, in this file, both ways, and asserts that
    ``sigma^3`` is the one that reproduces FeOs while ``d^3`` does not - by
    eleven orders of magnitude, so the verdict cannot be ambiguous.
    """
    temperature, density = 300.0, 55000.0
    _, m, sigma, epsilon_k, kappa, epsilon_ab = PARAMETERS["Water"]
    assert kappa is not None
    rho = density * 6.02214076e23 * 1e-30
    d = sigma * (1.0 - 0.12 * math.exp(-3.0 * epsilon_k / temperature))
    zeta_2 = math.pi / 6.0 * rho * m * d**2
    eta = math.pi / 6.0 * rho * m * d**3
    u = 1.0 - eta
    contact = 0.5 * d
    g = 1.0 / u + 3.0 * contact * zeta_2 / u**2 + 2.0 * contact**2 * zeta_2**2 / u**3

    def association(prefactor: float) -> float:
        strength = rho * prefactor * g * kappa * math.expm1(epsilon_ab / temperature)
        site = (-1.0 + math.sqrt(1.0 + 4.0 * strength)) / (2.0 * strength)
        return 2.0 * math.log(site) - site + 1.0

    contributions, _, _ = _feos_report(("Water",), temperature, density, [1.0])
    reference = contributions["Association"]
    assert association(sigma**3) == pytest.approx(reference, rel=1e-12, abs=1e-12)
    assert abs(association(d**3) - reference) > 1e-3

    # And the shipped model agrees with the sigma^3 form.
    internal = PCSAFTEOS(components=("Water",))._state(temperature, density, [1.0])
    assert internal.a_assoc == pytest.approx(association(sigma**3), rel=1e-13, abs=1e-13)


def test_the_agreement_is_not_vacuous() -> None:
    """Negative control: change one association parameter and it must break."""
    temperature, density = 300.0, 55000.0
    contributions, _, _ = _feos_report(("Water",), temperature, density, [1.0])
    perturbed = ct.PCSAFTParameters.from_records(
        [
            {
                "name": "Water",
                "m": 1.0656,
                "sigma_A": 3.0007,
                "epsilon_k_K": 366.51,
                "association": {"kappa_ab": 0.034868 * 1.01, "epsilon_ab_k_K": 2500.7},
            }
        ]
    )
    internal = PCSAFTEOS(components=("Water",), parameters=perturbed)._state(
        temperature, density, [1.0]
    )
    assert abs(internal.a_assoc - contributions["Association"]) > 1e-3


def test_non_associating_hexane_is_unchanged() -> None:
    """Case P-1's pinned n-hexane numbers, reproduced bit for bit."""
    eos = PCSAFTEOS(components=("n-Hexane",))
    assert not eos.associates()
    assert (
        eos.residual_helmholtz(temperature_K=300.0, volume_m3=1.0 / 7700.0, composition=[1.0])
        == -5.783742760059239
    )
    assert (
        eos.compressibility_factor(temperature_K=300.0, density_mol_m3=7700.0, composition=[1.0])
        == 0.661534529144653
    )
    assert (
        eos.ln_fugacity_coefficients(temperature_K=300.0, density_mol_m3=7700.0, composition=[1.0])[
            0
        ]
        == -5.709015132378622
    )


# ---------------------------------------------------------------------------
# Case P-7: the equilibrium exam, on the shipped constants
# ---------------------------------------------------------------------------


def _bisect_saturation(
    eos: PCSAFTEOS, temperature: float, low: float, high: float, iterations: int = 90
) -> tuple[float, float, float]:
    """Return ``(Psat, rho_liquid, rho_vapour)`` for a pure fluid.

    A **test-only** solver, as in Case P-2: bisection on
    ``ln phi(liquid root) - ln phi(vapour root) = 0``, with chemthermo's own
    density roots supplying the two branches. FeOs reaches the same state by a
    Newton iteration of its own, so only the model is shared.
    """

    def gap(pressure: float) -> float:
        roots = eos.density_roots(
            temperature_K=temperature, pressure_Pa=pressure, composition=[1.0]
        )
        assert len(roots) >= 2, f"only one density root at {pressure!r} Pa"
        ln_phi = [
            eos.ln_fugacity_coefficients(
                temperature_K=temperature, density_mol_m3=root, composition=[1.0]
            )[0]
            for root in (roots[0], roots[-1])
        ]
        return ln_phi[1] - ln_phi[0]

    low_value = gap(low)
    assert low_value * gap(high) < 0.0, "the saturation pressure is not bracketed"
    for _ in range(iterations):
        middle = 0.5 * (low + high)
        value = gap(middle)
        if low_value * value <= 0.0:
            high = middle
        else:
            low, low_value = middle, value
    pressure = 0.5 * (low + high)
    roots = eos.density_roots(temperature_K=temperature, pressure_Pa=pressure, composition=[1.0])
    return pressure, roots[-1], roots[0]


def test_pure_water_saturation_matches_feos() -> None:
    """Case P-7 (i): 373.15 K, against FeOs's own ``PhaseEquilibrium.pure``."""
    temperature = 373.15
    eos = PCSAFTEOS(components=("Water",))
    pressure, rho_liquid, rho_vapour = _bisect_saturation(eos, temperature, 5.0e4, 2.0e5)

    equilibrium = PhaseEquilibrium.pure(_feos_eos(("Water",)), temperature * si.KELVIN)
    pressure_reference = equilibrium.liquid.pressure() / si.PASCAL
    rho_liquid_reference = equilibrium.liquid.density / _MOL_PER_M3
    rho_vapour_reference = equilibrium.vapor.density / _MOL_PER_M3

    assert pressure == pytest.approx(pressure_reference, rel=1e-6)
    assert rho_liquid == pytest.approx(rho_liquid_reference, rel=1e-6)
    assert rho_vapour == pytest.approx(rho_vapour_reference, rel=1e-6)

    # The saturation condition restated on chemthermo's own numbers.
    ln_phi_liquid = eos.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=rho_liquid, composition=[1.0]
    )[0]
    ln_phi_vapour = eos.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=rho_vapour, composition=[1.0]
    )[0]
    assert abs(ln_phi_liquid - ln_phi_vapour) < 1e-9

    # Model versus experiment, a remark and not an assertion: water boils at
    # 101.325 kPa at 373.15 K by definition of the normal boiling point, and
    # PC-SAFT with the 2002 parameters puts its saturation pressure ~0.43 %
    # below that. Nothing here asserts the experimental number.
    assert 0.99e5 < pressure < 1.02e5


def _feos_ln_fugacity(
    names: Sequence[str], temperature: float, density: float, x: Sequence[float]
) -> np.ndarray:
    _, _, ln_phi = _feos_report(names, temperature, density, x)
    return np.atleast_1d(ln_phi)


def _equal_fugacity_residual(
    names: Sequence[str],
    temperature: float,
    pressure: float,
    phases: Sequence[Sequence[float]],
    densities: Sequence[float],
) -> float:
    """Return ``max_i |ln(x_i^I phi_i^I) - ln(x_i^II phi_i^II)|`` from FeOs.

    The compositions and the densities are chemthermo's; only the fugacity
    coefficients are FeOs's. A small residual therefore says that the state
    chemthermo converged on is an equilibrium state *of the reference model*.
    """
    logs = []
    for x, density in zip(phases, densities):
        ln_phi = _feos_ln_fugacity(names, temperature, density, x)
        logs.append(np.log(np.asarray(x, dtype=float)) + ln_phi)
    return float(np.max(np.abs(logs[0] - logs[1])))


@pytest.mark.usefixtures("matched_constants")
def test_water_ethanol_vapor_liquid_flash_matches_feos() -> None:
    """Case P-7 (ii): a two-phase VLE feed at 351 K, verified against FeOs.

    Run with FeOs's universal constants for the same reason as
    :func:`test_every_term_matches_feos`: the equal-fugacity residual below is
    evaluated by FeOs at *chemthermo's* compositions and densities, so a
    difference in the shared a/b table shows up in it directly. Measured: the
    tie line is the same to eight decimal places either way, and the residual
    is 6.1e-10 with FeOs's constants against 2.5e-06 with the published ones.
    """
    names = ("Water", "Ethanol")
    temperature, pressure = 351.0, 80.0e3
    mixture = ct.Mixture.from_database(list(names), [0.7, 0.3])
    eos = PCSAFTEOS()

    stability = ct.stability_tp(mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos)
    assert stability.status == "unstable"

    result = ct.flash_tp(mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos)
    assert set(result.phases) == {"liquid", "vapor"}
    assert result.vapor_fraction is not None
    assert 0.0 < result.vapor_fraction < 1.0
    assert result.diagnostics["phase_regime"] == "VLE"
    assert _number(result.diagnostics, "delta_g_split_rt") < 0.0
    assert _number(result.diagnostics, "mass_balance_residual") < 1e-10
    assert _number(result.diagnostics, "fugacity_residual") < 1e-8
    assert result.diagnostics["post_split_status"] == "stable"

    liquid = list(result.phases["liquid"].composition.fractions)
    vapor = list(result.phases["vapor"].composition.fractions)
    # The liquid is the water-rich phase and the vapor the ethanol-enriched one.
    assert liquid[0] > 0.9 and vapor[1] > vapor[0] * 0.5

    rho_liquid = eos.density_roots(
        temperature_K=temperature, pressure_Pa=pressure, composition=liquid, mixture=mixture
    )[-1]
    rho_vapor = eos.density_roots(
        temperature_K=temperature, pressure_Pa=pressure, composition=vapor, mixture=mixture
    )[0]
    residual = _equal_fugacity_residual(
        names, temperature, pressure, (liquid, vapor), (rho_liquid, rho_vapor)
    )
    assert residual < 1e-8


def test_water_hexane_is_unstable_and_the_phi_phi_split_cannot_express_two_liquids() -> None:
    """Case P-7 (iii), part 1: 298.15 K, 1 atm, z = 0.5 / 0.5.

    The tangent-plane test is right (the feed is unstable) and the phi-phi
    split then fails, **honestly**: it pairs a vapour-root phase with a
    liquid-root one, which is the only pairing ``_split`` can express, and the
    post-split stability test refuses the result. At 298 K and 1 atm the true
    answer is two liquids - water's and n-hexane's vapour pressures sum to
    about 23 kPa, far below 1 atm, so there is no vapour phase - and FeOs's own
    two-phase flash finds exactly that. Recorded as a limitation, not worked
    around; see ADR-0018 and the ``flash-phase-addition-eos`` slice.
    """
    names = ("Water", "n-Hexane")
    temperature, pressure = 298.15, 101325.0
    mixture = ct.Mixture.from_database(list(names), [0.5, 0.5])
    eos = PCSAFTEOS()

    stability = ct.stability_tp(mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos)
    assert stability.status == "unstable"
    assert stability.tpd_min is not None and stability.tpd_min < -0.5

    with pytest.raises(ct.ConvergenceError, match="third phase"):
        ct.flash_tp(mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos)

    # What the split actually converges on, with the guard switched off: a
    # water-rich liquid against a hexane-rich *vapour* whose Gibbs energy is
    # above the feed's, which is why the guard fires.
    unguarded = ct.flash_tp(
        mixture,
        temperature_K=temperature,
        pressure_Pa=pressure,
        eos=eos,
        settings=ct.FlashSettings(post_split_stability=False),
    )
    assert _number(unguarded.diagnostics, "delta_g_split_rt") > 0.0
    assert unguarded.diagnostics["post_split_stable"] is False

    # FeOs's two-phase flash at the same state returns the liquid-liquid pair.
    reference = State(
        _feos_eos(names),
        temperature=temperature * si.KELVIN,
        pressure=pressure * si.PASCAL,
        composition=np.array([0.5, 0.5]),
    ).tp_flash()
    compositions = sorted(
        float(np.atleast_1d(phase.molefracs)[0]) for phase in (reference.liquid, reference.vapor)
    )
    assert compositions[0] < 0.02 and compositions[1] > 0.98
    for phase in (reference.liquid, reference.vapor):
        assert float(phase.density / _MOL_PER_M3) > 5000.0  # both are liquids


@pytest.mark.usefixtures("matched_constants")
def test_water_hexane_liquid_liquid_split_above_the_vapour_root() -> None:
    """Case P-7 (iii), part 2: the same tie line at 1 MPa, where it is reachable.

    Raising the pressure above the vapour branch's existence limit leaves the
    isotherm with a single density root at every composition, so ``"vapor"``
    and ``"liquid"`` name the *same* (liquid) root and the phi-phi machinery
    expresses a genuine liquid-liquid split. Both converged phases are liquids
    by the ADR-0017 compressibility criterion; ADR-0017's naming then falls
    back to the Wilson ranking and calls them ``"liquid"`` / ``"vapor"``, and
    ``vapor_fraction`` is the hexane-rich *liquid*'s fraction. That mislabel is
    recorded here, not hidden.

    Run with FeOs's universal constants, as the VLE test above and for the same
    reason. Measured: the tie line is the same to eight decimal places with the
    published constants, and the FeOs equal-fugacity residual is 1.2e-11 here
    against 2.6e-07 there.
    """
    names = ("Water", "n-Hexane")
    temperature, pressure = 298.15, 1.0e6
    mixture = ct.Mixture.from_database(list(names), [0.5, 0.5])
    eos = PCSAFTEOS()

    for composition in ([0.5, 0.5], [0.999, 0.001], [0.001, 0.999]):
        roots = eos.density_roots(
            temperature_K=temperature,
            pressure_Pa=pressure,
            composition=composition,
            mixture=mixture,
        )
        assert len(roots) == 1, "the vapour root must be gone for this to be a clean LL split"

    assert (
        ct.stability_tp(mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos).status
        == "unstable"
    )
    result = ct.flash_tp(mixture, temperature_K=temperature, pressure_Pa=pressure, eos=eos)
    assert _number(result.diagnostics, "delta_g_split_rt") < 0.0
    assert _number(result.diagnostics, "fugacity_residual") < 1e-10
    assert result.diagnostics["post_split_status"] == "stable"
    # The labels actually produced (ADR-0017 fallback), recorded as fact:
    assert set(result.phases) == {"liquid", "vapor"}
    assert result.diagnostics["phase_label_method"] == "wilson-ranking"

    phases = {name: list(phase.composition.fractions) for name, phase in result.phases.items()}
    densities = {}
    for name, composition in phases.items():
        assert (
            eos.phase_identity(
                mixture=mixture,
                temperature_K=temperature,
                pressure_Pa=pressure,
                composition=composition,
                phase="liquid",
            )
            == "liquid"
        ), f"phase {name!r} is not liquid-like, so this is not a liquid-liquid split"
        densities[name] = eos.density_roots(
            temperature_K=temperature,
            pressure_Pa=pressure,
            composition=composition,
            mixture=mixture,
        )[0]

    water_rich = phases["liquid"]
    hexane_rich = phases["vapor"]
    assert water_rich[0] > 0.99 and hexane_rich[1] > 0.99

    residual = _equal_fugacity_residual(
        names,
        temperature,
        pressure,
        (water_rich, hexane_rich),
        (densities["liquid"], densities["vapor"]),
    )
    assert residual < 1e-8

    # Mutual solubilities. PC-SAFT with k_ij = 0 is known to be poor for
    # water/hydrocarbon systems, so these are a **code** check against FeOs,
    # not a model validation: the commonly quoted experimental values are
    # ~5e-4 mole fraction water in hexane and ~2e-6 hexane in water at 298 K
    # (tabulated figures, not verified against a primary source here), and the
    # model is orders of magnitude away on both. Nothing asserts them.
    assert 0.0 < water_rich[1] < 1e-3
    assert 0.0 < hexane_rich[0] < 1e-1
