"""``PCSAFTEOS.residual_helmholtz_temperature_derivative`` (ADR-0034, Case P-19).

Default-suite checks that need no optional dependency: central finite
differences of ``residual_helmholtz`` in temperature, a negative control, the
refusal for associating mixtures, and the input contract. The external check
against teqp's ``get_Ar10`` lives in ``tests/validation/test_pcsaft_vs_teqp.py``.
"""

from __future__ import annotations

import pytest

from chemthermo.eos import PCSAFTEOS
from chemthermo.exceptions import InputRangeError, ModelError
from chemthermo.parameters import PCSAFTParameters, PCSAFTRecord

#: ``(components, kij, x, T/K, rho/(mol/m^3))``: gas-like and liquid-like, a
#: per-pair ``kij``, a ternary, and a state inside the spinodal (``Z < 0``).
STATES = [
    (("n-Hexane",), 0.0, [1.0], 300.0, 100.0),
    (("n-Hexane",), 0.0, [1.0], 300.0, 7700.0),
    (("Methane", "n-Hexane"), 0.0, [0.5, 0.5], 300.0, 11000.0),
    (("Methane", "n-Hexane"), 0.0, [0.2, 0.8], 300.0, 8000.0),
    (("Methane", "n-Decane"), {("Methane", "n-Decane"): 0.03}, [0.3, 0.7], 350.0, 6500.0),
    (("Methane", "n-Hexane", "Nitrogen"), 0.0, [0.1, 0.7, 0.2], 250.0, 10000.0),
]


def _central_difference(eos: PCSAFTEOS, x: list[float], t: float, rho: float) -> float:
    h = 1e-4 * t

    def a(temperature: float) -> float:
        return eos.residual_helmholtz(temperature_K=temperature, volume_m3=1.0 / rho, composition=x)

    # Fourth-order stencil: truncation ~h^4, so ~1e-13 relative here.
    return (-a(t + 2 * h) + 8 * a(t + h) - 8 * a(t - h) + a(t - 2 * h)) / (12 * h)


@pytest.mark.parametrize(("components", "kij", "x", "temperature", "density"), STATES)
def test_matches_a_central_difference(
    components: tuple[str, ...], kij: object, x: list[float], temperature: float, density: float
) -> None:
    eos = PCSAFTEOS(components=components, kij=kij)  # type: ignore[arg-type]
    derivative = eos.residual_helmholtz_temperature_derivative(
        temperature_K=temperature, volume_m3=1.0 / density, composition=x
    )
    assert derivative == pytest.approx(
        _central_difference(eos, x, temperature, density), rel=1e-8, abs=1e-14
    )


def test_is_not_vacuous() -> None:
    """A 1 % change in one epsilon/k moves the derivative far outside the FD tolerance."""
    base = PCSAFTEOS(components=("n-Hexane",))
    moved = PCSAFTEOS(
        components=("n-Hexane",),
        parameters=PCSAFTParameters.from_records(
            [PCSAFTRecord(name="n-Hexane", m=3.0576, sigma_A=3.7983, epsilon_k_K=236.77 * 1.01)]
        ),
    )
    kwargs = dict(temperature_K=300.0, volume_m3=1.0 / 7700.0, composition=[1.0])
    a = base.residual_helmholtz_temperature_derivative(**kwargs)  # type: ignore[arg-type]
    b = moved.residual_helmholtz_temperature_derivative(**kwargs)  # type: ignore[arg-type]
    assert abs(b - a) / abs(a) > 1e-3


def test_an_associating_mixture_is_refused_not_approximated() -> None:
    eos = PCSAFTEOS(components=("Water", "n-Hexane"))
    assert eos.associates()
    with pytest.raises(ModelError, match="association"):
        eos.residual_helmholtz_temperature_derivative(
            temperature_K=300.0, volume_m3=1.0 / 20000.0, composition=[0.5, 0.5]
        )


@pytest.mark.parametrize("volume", [0.0, -1e-4, float("nan")])
def test_rejects_a_non_physical_volume(volume: float) -> None:
    eos = PCSAFTEOS(components=("n-Hexane",))
    with pytest.raises(InputRangeError):
        eos.residual_helmholtz_temperature_derivative(
            temperature_K=300.0, volume_m3=volume, composition=[1.0]
        )


# ---------------------------------------------------------------------------
# Residual properties (C3)
# ---------------------------------------------------------------------------

#: ``(components, x, T/K, P/Pa, branch)``: one liquid and one vapour root.
TP_STATES = [
    (("n-Hexane",), [1.0], 300.0, 1.0e6, "liquid"),
    (("Methane", "n-Hexane"), [0.8, 0.2], 350.0, 2.0e6, "vapor"),
    (("Methane", "n-Decane"), [0.3, 0.7], 350.0, 5.0e6, "liquid"),
]


def _root(eos: PCSAFTEOS, x: list[float], t: float, p: float, branch: str) -> float:
    roots = eos.density_roots(temperature_K=t, pressure_Pa=p, composition=x)
    return max(roots) if branch == "liquid" else min(roots)


def _g_res_tp_from_ln_phi(eos: PCSAFTEOS, x: list[float], t: float, p: float, branch: str) -> float:
    rho = _root(eos, x, t, p, branch)
    ln_phi = eos.ln_fugacity_coefficients(temperature_K=t, density_mol_m3=rho, composition=x)
    return float(sum(xi * lp for xi, lp in zip(x, ln_phi)))


@pytest.mark.parametrize(("components", "x", "temperature", "pressure", "branch"), TP_STATES)
def test_gibbs_helmholtz_through_the_fugacity_route(
    components: tuple[str, ...], x: list[float], temperature: float, pressure: float, branch: str
) -> None:
    """H^res/RT = -T d(G^res/RT)/dT at fixed P, with G^res/RT = sum x ln phi.

    The right-hand side never calls the temperature derivative: it re-solves
    the density at each temperature and goes through the composition
    derivatives behind ``ln phi`` - an independent route (ledger Case P-19).
    """
    eos = PCSAFTEOS(components=components)
    rho = _root(eos, x, temperature, pressure, branch)
    props = eos.residual_properties(temperature_K=temperature, density_mol_m3=rho, composition=x)

    assert props["g_res_tp"] == pytest.approx(
        _g_res_tp_from_ln_phi(eos, x, temperature, pressure, branch), rel=1e-12, abs=1e-13
    )
    h = 1e-3 * temperature

    def g(t: float) -> float:
        return _g_res_tp_from_ln_phi(eos, x, t, pressure, branch)

    dg_dt = (-g(temperature + 2 * h) + 8 * g(temperature + h) - 8 * g(temperature - h)
             + g(temperature - 2 * h)) / (12 * h)  # fmt: skip
    assert props["h_res"] == pytest.approx(-temperature * dg_dt, rel=1e-7, abs=1e-9)


def test_the_two_references_differ_by_ln_z_and_agree_on_energies() -> None:
    import math

    eos = PCSAFTEOS(components=("n-Hexane",))
    props = eos.residual_properties(temperature_K=300.0, density_mol_m3=7700.0, composition=[1.0])
    ln_z = math.log(props["z"])
    assert props["s_res_tp"] - props["s_res_tv"] == pytest.approx(ln_z, rel=1e-14)
    assert props["g_res_tv"] - props["g_res_tp"] == pytest.approx(ln_z, rel=1e-14)
    # G = H - T S holds in both references.
    assert props["g_res_tv"] == pytest.approx(props["h_res"] - props["s_res_tv"], rel=1e-13)
    assert props["g_res_tp"] == pytest.approx(props["h_res"] - props["s_res_tp"], rel=1e-13)
    # A liquid: residual enthalpy and entropy are large and negative.
    assert props["h_res"] < -5.0 and props["s_res_tp"] < 0.0


def test_a_dilute_gas_has_vanishing_residual_properties() -> None:
    eos = PCSAFTEOS(components=("Methane",))
    props = eos.residual_properties(temperature_K=400.0, density_mol_m3=1e-3, composition=[1.0])
    for key in ("a_res", "u_res", "h_res", "s_res_tv", "s_res_tp", "g_res_tv", "g_res_tp"):
        assert abs(props[key]) < 1e-6, key


def test_residual_properties_refuse_inside_the_spinodal_and_for_association() -> None:
    with pytest.raises(ModelError, match="spinodal"):
        PCSAFTEOS(components=("Methane", "n-Hexane")).residual_properties(
            temperature_K=300.0, density_mol_m3=8000.0, composition=[0.2, 0.8]
        )
    with pytest.raises(ModelError, match="association"):
        PCSAFTEOS(components=("Water",)).residual_properties(
            temperature_K=300.0, density_mol_m3=100.0, composition=[1.0]
        )


def test_enthalpy_and_entropy_of_vaporization_satisfy_clausius_at_saturation() -> None:
    """At a saturation point dG_vap = 0, so dH_vap = T dS_vap exactly.

    The pressure is n-hexane's PC-SAFT saturation pressure at 300 K from
    teqp's own `pure_VLE_T` (ledger Case P-2), so neither side is tuned to
    make this pass; the residual G of the two roots only equals because the
    fugacities do.
    """
    eos = PCSAFTEOS(components=("n-Hexane",))
    vapor_rho, liquid_rho = eos.density_roots(
        temperature_K=300.0, pressure_Pa=21858.084278856164, composition=[1.0]
    )
    vapor = eos.residual_properties(
        temperature_K=300.0, density_mol_m3=vapor_rho, composition=[1.0]
    )
    liquid = eos.residual_properties(
        temperature_K=300.0, density_mol_m3=liquid_rho, composition=[1.0]
    )
    dh = vapor["h_res"] - liquid["h_res"]
    ds = vapor["s_res_tp"] - liquid["s_res_tp"]
    assert dh > 12.0  # about 31.5 kJ/mol in units of RT at 300 K
    assert dh == pytest.approx(ds, rel=1e-9)
