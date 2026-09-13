"""PC-SAFT (Gross & Sadowski 2001, non-associating) unit tests.

Replaces the former ``tests/test_pcsaft_sanity.py``, whose only expectations
were ``"missing_parameters"`` / ``"not_implemented"``.

Nothing here needs an optional dependency; the teqp cross-check lives in
``tests/validation/test_pcsaft_vs_teqp.py``. These tests are the *internal*
route: independently written pure-component formulas, finite differences, and
thermodynamic identities derived in the test docstrings.

See validation Cases P-0, P-1, P-2 in ``.agents/brain/validation-cases.md``.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.eos.pcsaft import (
    A_UNIVERSAL,
    AVOGADRO_PER_MOL,
    B_UNIVERSAL,
    BOLTZMANN_J_PER_K,
    PCSAFTEOS,
    R_J_PER_MOL_K,
    _c1_terms,
)


def _z(eos: PCSAFTEOS, temperature_K: float, density_mol_m3: float, x: Sequence[float]) -> float:
    """Typed shorthand; keeps the tests free of untyped ``**kwargs`` dicts."""
    return eos.compressibility_factor(
        temperature_K=temperature_K, density_mol_m3=density_mol_m3, composition=x
    )


def _ln_phi(
    eos: PCSAFTEOS, temperature_K: float, density_mol_m3: float, x: Sequence[float]
) -> list[float]:
    return eos.ln_fugacity_coefficients(
        temperature_K=temperature_K, density_mol_m3=density_mol_m3, composition=x
    )


# Gross & Sadowski (2001) Table 1 values used by the hand-written reference
# formulas below: (m, sigma / Angstrom, epsilon/k / K).
METHANE = (1.0000, 3.7039, 150.03)
HEXANE = (3.0576, 3.7983, 236.77)
DECANE = (4.6627, 3.8384, 243.87)
NITROGEN = (1.2053, 3.3130, 90.96)


# ---------------------------------------------------------------------------
# A. Universal constants and the typo-corrected C1
# ---------------------------------------------------------------------------


def test_universal_constants_have_the_published_shape_and_checksums() -> None:
    """Transcription guard for the 21 + 21 constants of G&S (2001) Table 1.

    The paper is paywalled and was not read directly. The values in
    ``chemthermo.eos.pcsaft`` were taken from two independent sources that
    agree digit for digit:

    1. teqp (NIST, MIT), ``src/data/PCSAFT.cpp``, namespace
       ``teqp::saft::PCSAFT::PCSAFTMatrices::GrossSadowski2001``;
    2. the table in the Wikipedia article "PC-SAFT", section "Dispersion Term",
       which cites the same paper.

    This test does not re-derive them - nothing can, they are regressed
    constants - it pins the shape, a few individually quoted entries, and row /
    column checksums so that a single mistyped digit cannot pass unnoticed.
    """
    assert A_UNIVERSAL.shape == (3, 7)
    assert B_UNIVERSAL.shape == (3, 7)

    # Individually quoted corners, exactly as both sources print them.
    assert A_UNIVERSAL[0, 0] == 0.9105631445
    assert A_UNIVERSAL[2, 6] == -8.6728470368
    assert B_UNIVERSAL[0, 0] == 0.7240946941
    assert B_UNIVERSAL[2, 6] == -29.666905585

    np.testing.assert_allclose(
        A_UNIVERSAL.sum(axis=1),
        [7.1509055854999986, 3.1103125482000067, 0.20783018529999886],
        rtol=0.0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        B_UNIVERSAL.sum(axis=1),
        [-144.23916423030002, -147.5608352044, 36.4835890744],
        rtol=0.0,
        atol=1e-10,
    )
    assert A_UNIVERSAL.sum() == pytest.approx(10.469048318999999, abs=1e-12)
    assert B_UNIVERSAL.sum() == pytest.approx(-255.31641036030004, abs=1e-10)


def test_c1_is_the_reciprocal_of_the_bracket_not_the_bracket() -> None:
    """Eq. (A.11) as printed in the paper drops the outer exponent -1.

    NIST TRC's PC-SAFT page states the erratum: the equation is printed as
    ``C1 = (...)^-1 = (...)`` and the right-hand side should also carry the
    ``-1``. The two readings are only equal where the bracket equals 1, i.e.
    at ``eta = 0``; everywhere else the typo form is the reciprocal of the
    right answer. This test pins the reciprocal reading.
    """
    for mbar in (1.0, 2.0, 3.0576, 4.6627):
        for eta in (1e-8, 0.01, 0.1, 0.3, 0.45):
            c1, _c1_deta, _c1_dmbar = _c1_terms(eta, mbar)
            one_minus = 1.0 - eta
            gap = one_minus * (2.0 - eta)
            bracket = (
                1.0
                + mbar * (8.0 * eta - 2.0 * eta**2) / one_minus**4
                + (1.0 - mbar)
                * (20.0 * eta - 27.0 * eta**2 + 12.0 * eta**3 - 2.0 * eta**4)
                / gap**2
            )
            assert c1 * bracket == pytest.approx(1.0, abs=1e-14)
            if mbar >= 1.0 and eta > 1e-6:
                # The bracket grows with eta for a chain fluid, so the correct
                # C1 is strictly below 1 while the typo form is above it.
                assert bracket > 1.0
                assert 0.0 < c1 < 1.0


def test_c1_tends_to_one_as_density_tends_to_zero() -> None:
    for mbar in (1.0, 2.5, 4.6627):
        assert _c1_terms(0.0, mbar)[0] == pytest.approx(1.0, abs=0.0)
        assert _c1_terms(1e-12, mbar)[0] == pytest.approx(1.0, abs=1e-10)


def test_c1_eta_derivative_matches_finite_differences() -> None:
    step = 1e-7
    for mbar in (1.0, 2.0, 3.0576):
        for eta in (0.01, 0.1, 0.3, 0.45):
            analytic = _c1_terms(eta, mbar)[1]
            numeric = (_c1_terms(eta + step, mbar)[0] - _c1_terms(eta - step, mbar)[0]) / (
                2.0 * step
            )
            assert analytic == pytest.approx(numeric, rel=1e-6)


def test_c1_mbar_derivative_matches_finite_differences() -> None:
    step = 1e-7
    for mbar in (1.5, 2.0, 3.0576):
        for eta in (0.01, 0.1, 0.3, 0.45):
            analytic = _c1_terms(eta, mbar)[2]
            numeric = (_c1_terms(eta, mbar + step)[0] - _c1_terms(eta, mbar - step)[0]) / (
                2.0 * step
            )
            assert analytic == pytest.approx(numeric, rel=1e-6)


def test_dispersion_density_derivative_matches_finite_differences() -> None:
    """``rho d a_disp / d rho`` is the eta-derivative, since eta is linear in rho.

    Exercises the dispersion term's derivative on its own, separately from the
    hard-chain term that the total ``Z`` mixes it with.
    """
    eos = PCSAFTEOS(components=("Methane", "n-Hexane"))
    x = [0.4, 0.6]
    for temperature, density in ((300.0, 200.0), (300.0, 5000.0), (450.0, 8000.0)):
        step = density * 1e-6
        state = eos._state(temperature, density, x)
        upper = eos._state(temperature, density + step, x)
        lower = eos._state(temperature, density - step, x)
        numeric = density * (upper.a_disp - lower.a_disp) / (2.0 * step)
        assert state.z_disp == pytest.approx(numeric, rel=1e-7)

        numeric_hc = density * (upper.a_hc - lower.a_hc) / (2.0 * step)
        assert state.z_hc == pytest.approx(numeric_hc, rel=1e-7)


# ---------------------------------------------------------------------------
# B. Ideal-gas limit
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("components", "x"),
    [
        (("n-Hexane",), [1.0]),
        (("Methane", "n-Hexane"), [0.5, 0.5]),
        (("Methane", "n-Hexane", "Nitrogen"), [0.2, 0.5, 0.3]),
    ],
)
def test_ideal_gas_limit_at_vanishing_density(components: tuple[str, ...], x: list[float]) -> None:
    eos = PCSAFTEOS(components=components)
    density = 1e-6  # mol/m^3
    a_res = eos.residual_helmholtz(temperature_K=300.0, volume_m3=1.0 / density, composition=x)
    z_factor = eos.compressibility_factor(
        temperature_K=300.0, density_mol_m3=density, composition=x
    )
    ln_phi = eos.ln_fugacity_coefficients(
        temperature_K=300.0, density_mol_m3=density, composition=x
    )
    assert abs(a_res) < 1e-8
    assert abs(z_factor - 1.0) < 1e-8
    assert max(abs(value) for value in ln_phi) < 1e-8


def test_pressure_approaches_the_ideal_gas_law_at_vanishing_density() -> None:
    eos = PCSAFTEOS(components=("Methane",))
    density = 1e-6
    pressure = eos.pressure_Pa(temperature_K=300.0, density_mol_m3=density, composition=[1.0])
    ideal = density * R_J_PER_MOL_K * 300.0
    assert pressure == pytest.approx(ideal, rel=1e-8)


# ---------------------------------------------------------------------------
# C. Pure-component consistency against independently written formulas
# ---------------------------------------------------------------------------


def _pure_a_res(
    temperature_K: float,
    density_mol_m3: float,
    m: float,
    sigma_A: float,
    epsilon_k_K: float,
) -> float:
    """Reduced residual Helmholtz energy of a *pure* PC-SAFT fluid.

    Written from the pure-component form of the equations, not from the
    mixture code under test: Carnahan-Starling for the hard-sphere term,
    ``g = (1 - eta/2) / (1 - eta)^3`` for the contact value, and the
    single-component dispersion sums. Only the 42 universal constants and the
    exact SI constants are shared with the implementation.
    """
    rho = density_mol_m3 * AVOGADRO_PER_MOL * 1e-30  # molecules / Angstrom^3
    d = sigma_A * (1.0 - 0.12 * math.exp(-3.0 * epsilon_k_K / temperature_K))
    eta = math.pi / 6.0 * rho * m * d**3

    a_hs = (4.0 * eta - 3.0 * eta**2) / (1.0 - eta) ** 2
    g_hs = (1.0 - eta / 2.0) / (1.0 - eta) ** 3
    a_hc = m * a_hs - (m - 1.0) * math.log(g_hs)

    i1 = 0.0
    i2 = 0.0
    for n in range(7):
        a_n = (
            A_UNIVERSAL[0, n]
            + (m - 1.0) / m * A_UNIVERSAL[1, n]
            + (m - 1.0) / m * (m - 2.0) / m * A_UNIVERSAL[2, n]
        )
        b_n = (
            B_UNIVERSAL[0, n]
            + (m - 1.0) / m * B_UNIVERSAL[1, n]
            + (m - 1.0) / m * (m - 2.0) / m * B_UNIVERSAL[2, n]
        )
        i1 += a_n * eta**n
        i2 += b_n * eta**n

    c1 = 1.0 / (
        1.0
        + m * (8.0 * eta - 2.0 * eta**2) / (1.0 - eta) ** 4
        + (1.0 - m)
        * (20.0 * eta - 27.0 * eta**2 + 12.0 * eta**3 - 2.0 * eta**4)
        / ((1.0 - eta) * (2.0 - eta)) ** 2
    )

    eps_over_kt = epsilon_k_K / temperature_K
    m2es3 = m**2 * eps_over_kt * sigma_A**3
    m2e2s3 = m**2 * eps_over_kt**2 * sigma_A**3
    a_disp = -2.0 * math.pi * rho * i1 * m2es3 - math.pi * rho * m * c1 * i2 * m2e2s3

    return a_hc + a_disp


@pytest.mark.parametrize(
    ("name", "params"),
    [("Methane", METHANE), ("n-Hexane", HEXANE), ("n-Decane", DECANE), ("Nitrogen", NITROGEN)],
)
@pytest.mark.parametrize(
    ("temperature", "density"), [(200.0, 50.0), (300.0, 2000.0), (450.0, 6000.0)]
)
def test_mixture_code_reduces_to_the_pure_component_formulas(
    name: str, params: tuple[float, float, float], temperature: float, density: float
) -> None:
    expected = _pure_a_res(temperature, density, *params)
    actual = PCSAFTEOS(components=(name,)).residual_helmholtz(
        temperature_K=temperature, volume_m3=1.0 / density, composition=[1.0]
    )
    assert actual == pytest.approx(expected, rel=1e-13, abs=1e-15)


def test_a_binary_of_one_component_with_itself_equals_the_pure_fluid() -> None:
    """Two labels, one set of parameters: the mixture must collapse to the pure.

    A genuine internal invariant of the mixing rules (it is not imposed
    anywhere in the code), and it is blind to nothing: every quadratic sum,
    every ``zeta`` moment and the ``mbar`` mixing all have to be right.
    """
    parameters = ct.PCSAFTParameters.from_records(
        [
            {"name": "A", "m": HEXANE[0], "sigma_A": HEXANE[1], "epsilon_k_K": HEXANE[2]},
            {"name": "B", "m": HEXANE[0], "sigma_A": HEXANE[1], "epsilon_k_K": HEXANE[2]},
        ]
    )
    binary = PCSAFTEOS(components=("A", "B"), parameters=parameters)
    pure = PCSAFTEOS(components=("n-Hexane",))
    for x1 in (0.0 + 1e-12, 0.25, 0.5, 0.9):
        assert _z(binary, 300.0, 7700.0, [x1, 1.0 - x1]) == pytest.approx(
            _z(pure, 300.0, 7700.0, [1.0]), rel=1e-14
        )
        ln_phi = _ln_phi(binary, 300.0, 7700.0, [x1, 1.0 - x1])
        pure_ln_phi = _ln_phi(pure, 300.0, 7700.0, [1.0])[0]
        assert ln_phi[0] == pytest.approx(pure_ln_phi, rel=1e-12)
        assert ln_phi[1] == pytest.approx(pure_ln_phi, rel=1e-12)


# ---------------------------------------------------------------------------
# D. Derivative discipline
# ---------------------------------------------------------------------------

_DERIVATIVE_STATES = [
    (("n-Hexane",), [1.0], 300.0, 100.0),
    (("n-Hexane",), [1.0], 300.0, 7700.0),
    (("Methane", "n-Hexane"), [0.5, 0.5], 300.0, 200.0),
    (("Methane", "n-Hexane"), [0.2, 0.8], 450.0, 7000.0),
    (("Methane", "n-Hexane", "Nitrogen"), [0.1, 0.7, 0.2], 250.0, 10000.0),
]


@pytest.mark.parametrize(("components", "x", "temperature", "density"), _DERIVATIVE_STATES)
def test_compressibility_matches_the_density_finite_difference(
    components: tuple[str, ...], x: list[float], temperature: float, density: float
) -> None:
    """``Z = 1 + rho (d a_res / d rho)_{T,x}`` against a central difference."""
    eos = PCSAFTEOS(components=components)
    step = density * 1e-6

    def a_res(rho: float) -> float:
        return eos.residual_helmholtz(temperature_K=temperature, volume_m3=1.0 / rho, composition=x)

    numeric = 1.0 + density * (a_res(density + step) - a_res(density - step)) / (2.0 * step)
    analytic = eos.compressibility_factor(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    assert analytic == pytest.approx(numeric, rel=1e-8)


@pytest.mark.parametrize(("components", "x", "temperature", "density"), _DERIVATIVE_STATES[2:])
def test_residual_chemical_potentials_match_mole_number_finite_differences(
    components: tuple[str, ...], x: list[float], temperature: float, density: float
) -> None:
    """``mu_i^res/RT = d(n a_res)/dn_i`` at fixed T and V, by definition.

    The volume is held fixed while one mole number moves, so both the total
    density and the composition change - which is exactly the derivative
    ``ln phi_i + ln Z`` is supposed to be.
    """
    eos = PCSAFTEOS(components=components)
    volume_m3 = 1.0  # fixed container; n_total = density * volume
    moles = np.array(x, dtype=float) * density * volume_m3

    def extensive(n: np.ndarray) -> float:
        total = float(n.sum())
        return total * eos.residual_helmholtz(
            temperature_K=temperature,
            volume_m3=volume_m3 / total,
            composition=(n / total).tolist(),
        )

    z_factor = eos.compressibility_factor(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    analytic = [
        value + math.log(z_factor)
        for value in eos.ln_fugacity_coefficients(
            temperature_K=temperature, density_mol_m3=density, composition=x
        )
    ]

    for index in range(len(x)):
        step = moles[index] * 1e-6
        upper = moles.copy()
        lower = moles.copy()
        upper[index] += step
        lower[index] -= step
        numeric = (extensive(upper) - extensive(lower)) / (2.0 * step)
        assert analytic[index] == pytest.approx(numeric, rel=1e-7, abs=1e-9)


@pytest.mark.parametrize(("components", "x", "temperature", "density"), _DERIVATIVE_STATES)
def test_euler_identity_between_ln_phi_a_res_and_z(
    components: tuple[str, ...], x: list[float], temperature: float, density: float
) -> None:
    """``sum_i x_i ln phi_i = a_res + Z - 1 - ln Z``.

    Derivation: ``n a_res(T, n/V, x)`` is homogeneous of degree one in
    ``(V, n)`` jointly, so Euler's theorem gives
    ``n a_res = sum_i n_i mu_i^res/RT + V (d(n a_res)/dV)`` and the volume
    derivative is ``-(n/V)(Z - 1)``. Dividing by ``n`` and substituting
    ``mu_i^res/RT = ln phi_i + ln Z`` gives the identity. It is an exact
    consequence of the definitions, so it must hold to round-off.
    """
    eos = PCSAFTEOS(components=components)
    a_res = eos.residual_helmholtz(
        temperature_K=temperature, volume_m3=1.0 / density, composition=x
    )
    z_factor = eos.compressibility_factor(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    ln_phi = eos.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    left = float(np.array(x) @ np.array(ln_phi))
    right = a_res + z_factor - 1.0 - math.log(z_factor)
    assert left == pytest.approx(right, rel=1e-13, abs=1e-13)


@pytest.mark.parametrize(("components", "x", "temperature", "density"), _DERIVATIVE_STATES[2:])
def test_gibbs_duhem_at_fixed_temperature_and_density(
    components: tuple[str, ...], x: list[float], temperature: float, density: float
) -> None:
    """The fixed-(T, rho) form of Gibbs-Duhem for residual properties.

    The familiar ``sum_i x_i d ln phi_i = 0`` holds at fixed **T and P**. These
    ``ln phi_i`` are evaluated at fixed **T and molar density**, where the
    pressure moves with the composition, so that form is the wrong test.

    Differentiating the Euler relation above at fixed T and V while holding the
    total mole number fixed (so ``rho`` is fixed and ``sum_i dx_i = 0``) leaves

        sum_i x_i d(mu_i^res/RT) = dZ

    and with ``mu_i^res/RT = ln phi_i + ln Z`` that becomes

        sum_i x_i d ln phi_i = (1 - 1/Z) dZ .

    Both sides are taken along a direction ``s`` in the composition simplex by
    central differences.
    """
    eos = PCSAFTEOS(components=components)
    n = len(x)
    base = np.array(x, dtype=float)
    step = 1e-6

    directions = []
    for i in range(n - 1):
        s = np.zeros(n)
        s[i] = 1.0
        s[i + 1] = -1.0
        directions.append(s)

    for s in directions:

        def ln_phi(t: float, s: np.ndarray = s) -> np.ndarray:
            return np.array(
                eos.ln_fugacity_coefficients(
                    temperature_K=temperature,
                    density_mol_m3=density,
                    composition=(base + t * s).tolist(),
                )
            )

        def z_factor(t: float, s: np.ndarray = s) -> float:
            return eos.compressibility_factor(
                temperature_K=temperature,
                density_mol_m3=density,
                composition=(base + t * s).tolist(),
            )

        d_ln_phi = (ln_phi(step) - ln_phi(-step)) / (2.0 * step)
        d_z = (z_factor(step) - z_factor(-step)) / (2.0 * step)
        z_here = z_factor(0.0)
        left = float(base @ d_ln_phi)
        right = (1.0 - 1.0 / z_here) * d_z
        assert left == pytest.approx(right, rel=1e-6, abs=1e-6)


# ---------------------------------------------------------------------------
# G. Registry, parameters and API contract
# ---------------------------------------------------------------------------


def test_registry_exposes_pcsaft_and_builds_a_working_model() -> None:
    assert "pcsaft" in ct.list_eos()
    eos = ct.get_eos("pcsaft", components=["Methane", "Ethane"])
    assert eos.num_components() == 2
    assert eos.name == "PC-SAFT"
    value = eos.residual_helmholtz(
        temperature_K=280.0, volume_m3=1.0 / 500.0, composition=[0.5, 0.5]
    )
    assert math.isfinite(value)
    assert value < 0.0


def test_packaged_parameters_cover_the_eleven_published_compounds() -> None:
    names = ct.PCSAFTParameters.load().names()
    expected = {
        "methane",
        "ethane",
        "propane",
        "n-butane",
        "n-pentane",
        "n-hexane",
        "n-heptane",
        "n-octane",
        "n-decane",
        "nitrogen",
        "carbon dioxide",
    }
    assert set(names) == expected


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("Methane", METHANE),
        ("n-Hexane", HEXANE),
        ("n-Decane", DECANE),
        ("Nitrogen", NITROGEN),
        ("Carbon dioxide", (2.0729, 2.7852, 169.21)),
    ],
)
def test_packaged_parameter_values(name: str, expected: tuple[float, float, float]) -> None:
    """Pins the packaged values against Gross & Sadowski (2001) Table 1.

    Verified against two independent secondary sources that both cite the
    paper: FeOs ``parameters/pcsaft/gross2001.json`` and Clapeyron.jl
    ``database/SAFT/PCSAFT/PCSAFT_like.csv`` (whose ``source`` column is the
    DOI 10.1021/ie0003887). See validation Case P-0.
    """
    m, sigma_A, epsilon_k_K = ct.PCSAFTParameters.load().for_components([name])
    assert (float(m[0]), float(sigma_A[0]), float(epsilon_k_K[0])) == expected


def test_packaged_parameter_names_resolve_in_the_component_databank() -> None:
    """Every packaged PC-SAFT name must also name a databank component."""
    for name in ct.PCSAFTParameters.load().names():
        record = ct.Component.from_database(name)
        assert record.name


def test_missing_parameters_raise_pcsaft_parameter_error() -> None:
    eos = PCSAFTEOS(components=("Water",))
    with pytest.raises(ct.PCSAFTParameterError, match="Missing PC-SAFT parameters"):
        eos.residual_helmholtz(temperature_K=350.0, volume_m3=1.0, composition=[1.0])


def test_user_supplied_parameters_override_the_packaged_set() -> None:
    overridden = ct.PCSAFTParameters.from_records(
        [
            {
                "name": "n-Hexane",
                "m": 3.5,
                "sigma_A": 3.9,
                "epsilon_k_K": 250.0,
                "MW_g_mol": 86.177,
                "source": "made up for this test",
            }
        ]
    )
    default_eos = PCSAFTEOS(components=("n-Hexane",))
    custom_eos = PCSAFTEOS(components=("n-Hexane",), parameters=overridden)
    assert _z(custom_eos, 320.0, 4000.0, [1.0]) != _z(default_eos, 320.0, 4000.0, [1.0])
    expected = _pure_a_res(320.0, 4000.0, 3.5, 3.9, 250.0)
    assert custom_eos.residual_helmholtz(
        temperature_K=320.0, volume_m3=1.0 / 4000.0, composition=[1.0]
    ) == pytest.approx(expected, rel=1e-13)
    assert overridden.record("n-hexane").MW_g_mol == 86.177


def test_parameter_records_reject_non_positive_values() -> None:
    for bad in ({"m": 0.0}, {"sigma_A": -1.0}, {"epsilon_k_K": 0.0}, {"MW_g_mol": -2.0}):
        record = {"name": "X", "m": 1.0, "sigma_A": 3.0, "epsilon_k_K": 200.0} | bad
        with pytest.raises(ct.PCSAFTParameterError):
            ct.PCSAFTParameters.from_records([record])


def test_parameter_records_reject_duplicates_and_empty_sets() -> None:
    record = {"name": "X", "m": 1.0, "sigma_A": 3.0, "epsilon_k_K": 200.0}
    with pytest.raises(ct.PCSAFTParameterError, match="Duplicate"):
        ct.PCSAFTParameters.from_records([record, dict(record)])
    with pytest.raises(ct.PCSAFTParameterError, match="at least one record"):
        ct.PCSAFTParameters.from_records([])


def test_kij_contract_matches_peng_robinson() -> None:
    """Scalar and per-pair mapping are the same contract as ADR-0006."""
    scalar = PCSAFTEOS(components=("Methane", "n-Decane"), kij=0.03)
    mapping = PCSAFTEOS(components=("Methane", "n-Decane"), kij={("Methane", "n-Decane"): 0.03})
    recased = PCSAFTEOS(components=("Methane", "n-Decane"), kij={("n-decane", "METHANE"): 0.03})
    x = [0.3, 0.7]
    assert mapping.kij == ((("methane", "n-decane"), 0.03),)
    assert _z(mapping, 350.0, 100.0, x) == _z(scalar, 350.0, 100.0, x)
    assert _z(recased, 350.0, 100.0, x) == _z(scalar, 350.0, 100.0, x)

    zero = PCSAFTEOS(components=("Methane", "n-Decane"))
    assert _z(zero, 350.0, 100.0, x) != _z(scalar, 350.0, 100.0, x)

    with pytest.raises(ct.ModelError):
        PCSAFTEOS(components=("Methane", "n-Decane"), kij={("Methane", "Methane"): 0.01})
    with pytest.raises(ct.ModelError):
        PCSAFTEOS(
            components=("Methane", "n-Decane"),
            kij={("Methane", "n-Decane"): 0.03, ("n-Decane", "Methane"): 0.05},
        )


def test_kij_never_touches_the_diagonal() -> None:
    """A pure fluid must be independent of kij, as for Peng-Robinson."""
    baseline = _z(PCSAFTEOS(components=("n-Hexane",)), 320.0, 4000.0, [1.0])
    for kij in (0.3, {("Methane", "n-Hexane"): 0.1}):
        candidate = PCSAFTEOS(components=("n-Hexane",), kij=kij)  # type: ignore[arg-type]
        assert _z(candidate, 320.0, 4000.0, [1.0]) == baseline
    matrix = PCSAFTEOS(components=("Methane", "n-Decane"), kij=0.3).kij_matrix()
    assert matrix[0, 0] == 0.0
    assert matrix[1, 1] == 0.0
    assert matrix[0, 1] == matrix[1, 0] == 0.3


def test_ln_fugacity_coefficients_are_permutation_invariant() -> None:
    forward = PCSAFTEOS(
        components=("Methane", "n-Hexane", "Nitrogen"),
        kij={("Methane", "n-Hexane"): 0.02},
    )
    reversed_eos = PCSAFTEOS(
        components=("Nitrogen", "n-Hexane", "Methane"),
        kij={("n-Hexane", "Methane"): 0.02},
    )
    x = [0.2, 0.5, 0.3]
    a = _ln_phi(forward, 280.0, 400.0, x)
    b = _ln_phi(reversed_eos, 280.0, 400.0, list(reversed(x)))
    np.testing.assert_allclose(a, list(reversed(b)), rtol=0.0, atol=1e-15)
    assert _z(forward, 280.0, 400.0, x) == pytest.approx(
        _z(reversed_eos, 280.0, 400.0, list(reversed(x))), rel=1e-15
    )


def test_results_are_deterministic() -> None:
    eos = PCSAFTEOS(components=("Methane", "n-Hexane"))
    first = _ln_phi(eos, 300.0, 200.0, [0.5, 0.5])
    for _ in range(3):
        assert _ln_phi(eos, 300.0, 200.0, [0.5, 0.5]) == first


def test_invalid_inputs_raise_the_documented_errors() -> None:
    eos = PCSAFTEOS(components=("Methane", "n-Hexane"))
    with pytest.raises(ct.CompositionError):
        eos.compressibility_factor(temperature_K=300.0, density_mol_m3=100.0, composition=[1.0])
    with pytest.raises(ct.CompositionError):
        eos.compressibility_factor(
            temperature_K=300.0, density_mol_m3=100.0, composition=[0.5, 0.6]
        )
    with pytest.raises(ct.InputRangeError):
        eos.compressibility_factor(temperature_K=300.0, density_mol_m3=0.0, composition=[0.5, 0.5])
    with pytest.raises(ct.InputRangeError):
        eos.residual_helmholtz(temperature_K=300.0, volume_m3=-1.0, composition=[0.5, 0.5])
    with pytest.raises(ct.InputRangeError):
        eos.compressibility_factor(temperature_K=-1.0, density_mol_m3=100.0, composition=[0.5, 0.5])
    # Since ADR-0015 an empty ``components`` is legal at construction - it means
    # "take the order from the Mixture" - so the failure moved to the point of
    # use, where the order is genuinely needed and none is available.
    with pytest.raises(ct.ModelError, match="at least one component"):
        PCSAFTEOS().compressibility_factor(
            temperature_K=300.0, density_mol_m3=100.0, composition=[1.0]
        )
    with pytest.raises(ct.ModelError, match="no Mixture was supplied"):
        PCSAFTEOS().density_roots(temperature_K=300.0, pressure_Pa=1.0e5, composition=[1.0])


def test_ln_fugacity_coefficients_refuse_a_non_positive_compressibility() -> None:
    """Inside the spinodal ``Z`` can be negative; ln phi does not exist there."""
    eos = PCSAFTEOS(components=("n-Hexane",))
    z_factor = eos.compressibility_factor(
        temperature_K=300.0, density_mol_m3=1000.0, composition=[1.0]
    )
    assert z_factor < 0.0
    with pytest.raises(ct.ModelError, match="not positive"):
        eos.ln_fugacity_coefficients(temperature_K=300.0, density_mol_m3=1000.0, composition=[1.0])


def test_packing_fraction_outside_the_unit_interval_raises() -> None:
    eos = PCSAFTEOS(components=("n-Hexane",))
    with pytest.raises(ct.ModelError, match="packing fraction"):
        eos.compressibility_factor(temperature_K=300.0, density_mol_m3=5.0e4, composition=[1.0])


def test_pressure_uses_the_exact_si_gas_constant() -> None:
    assert R_J_PER_MOL_K == AVOGADRO_PER_MOL * BOLTZMANN_J_PER_K
    assert R_J_PER_MOL_K == pytest.approx(ct.units.R_J_PER_MOL_K, rel=1e-10)
    eos = PCSAFTEOS(components=("n-Hexane",))
    density, temperature = 100.0, 300.0
    z_factor = eos.compressibility_factor(
        temperature_K=temperature, density_mol_m3=density, composition=[1.0]
    )
    assert eos.pressure_Pa(
        temperature_K=temperature, density_mol_m3=density, composition=[1.0]
    ) == pytest.approx(z_factor * density * R_J_PER_MOL_K * temperature, rel=1e-15)
