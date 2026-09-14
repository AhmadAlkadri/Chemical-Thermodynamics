"""PC-SAFT cross-check against teqp (NIST, MIT), validation Cases P-1 and P-2.

teqp is an independent implementation of the same Gross & Sadowski (2001)
model that obtains **every derivative by automatic differentiation** of one
hand-written ``alphar``. chemthermo writes all of its derivatives analytically.
The two therefore share the model and the 42 universal constants but share no
derivative code at all, which is what makes this a real cross-check of the
``Z`` and ``ln phi`` routes and not only of ``A^res``.

Skipped when ``teqp`` is not installed (``pip install -e ".[validation]"``).
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

import numpy as np
import pytest

from chemthermo.eos.pcsaft import PCSAFTEOS, R_J_PER_MOL_K

teqp = pytest.importorskip("teqp")

# Gross & Sadowski (2001) Table 1, as packaged in
# ``src/chemthermo/parameters/data/eos/pcsaft.json``. Repeated here so the
# reference model is built from values written in this file, not read out of
# the code under test.
_PARAMETERS: dict[str, tuple[float, float, float]] = {
    "Methane": (1.0000, 3.7039, 150.03),
    "Nitrogen": (1.2053, 3.3130, 90.96),
    "Carbon dioxide": (2.0729, 2.7852, 169.21),
    "n-Hexane": (3.0576, 3.7983, 236.77),
    "n-Decane": (4.6627, 3.8384, 243.87),
}


def _teqp_model(components: Sequence[str], kij_matrix: np.ndarray):
    coefficients = [
        {
            "name": name,
            "m": _PARAMETERS[name][0],
            "sigma_Angstrom": _PARAMETERS[name][1],
            "epsilon_over_k": _PARAMETERS[name][2],
            "BibTeXKey": "Gross-IECR-2001",
        }
        for name in components
    ]
    return teqp.make_model(
        {
            "kind": "PCSAFT",
            "model": {"coeffs": coefficients, "kmat": kij_matrix.tolist()},
        }
    )


#: ``(label, components, kij, composition, temperature_K, density_mol_m3)``.
#: Four pure n-hexane states plus ten mixture states: gas-like and liquid-like
#: densities, a nonzero per-pair ``kij``, a ternary, and one state inside the
#: spinodal where ``Z < 0`` and only ``A^res`` and ``Z`` can be compared.
_STATES: list[tuple[str, tuple[str, ...], object, list[float], float, float]] = [
    ("hexane 300 K / 100", ("n-Hexane",), 0.0, [1.0], 300.0, 100.0),
    ("hexane 300 K / 7700", ("n-Hexane",), 0.0, [1.0], 300.0, 7700.0),
    ("hexane 400 K / 6800", ("n-Hexane",), 0.0, [1.0], 400.0, 6800.0),
    ("hexane 500 K / 3000", ("n-Hexane",), 0.0, [1.0], 500.0, 3000.0),
    (
        "C1/C6 0.5/0.5 300 K / 200 (gas-like)",
        ("Methane", "n-Hexane"),
        0.0,
        [0.5, 0.5],
        300.0,
        200.0,
    ),
    (
        "C1/C6 0.5/0.5 300 K / 11000 (liquid-like)",
        ("Methane", "n-Hexane"),
        0.0,
        [0.5, 0.5],
        300.0,
        11000.0,
    ),
    (
        "C1/C6 0.2/0.8 450 K / 8000 (liquid-like)",
        ("Methane", "n-Hexane"),
        0.0,
        [0.2, 0.8],
        450.0,
        8000.0,
    ),
    (
        "C1/C6 0.2/0.8 300 K / 8000 (inside the spinodal, Z < 0)",
        ("Methane", "n-Hexane"),
        0.0,
        [0.2, 0.8],
        300.0,
        8000.0,
    ),
    (
        "C1/C10 kij=0.03 0.3/0.7 350 K / 100",
        ("Methane", "n-Decane"),
        {("Methane", "n-Decane"): 0.03},
        [0.3, 0.7],
        350.0,
        100.0,
    ),
    (
        "C1/C10 kij=0.03 0.3/0.7 350 K / 6500 (liquid-like)",
        ("Methane", "n-Decane"),
        {("Methane", "n-Decane"): 0.03},
        [0.3, 0.7],
        350.0,
        6500.0,
    ),
    (
        "N2/C1 0.4/0.6 150 K / 20000 (liquid-like)",
        ("Nitrogen", "Methane"),
        0.0,
        [0.4, 0.6],
        150.0,
        20000.0,
    ),
    (
        "CO2/C10 0.4/0.6 320 K / 8000 (liquid-like)",
        ("Carbon dioxide", "n-Decane"),
        0.0,
        [0.4, 0.6],
        320.0,
        8000.0,
    ),
    (
        "C1/C6/N2 0.3/0.4/0.3 250 K / 500 (ternary)",
        ("Methane", "n-Hexane", "Nitrogen"),
        0.0,
        [0.3, 0.4, 0.3],
        250.0,
        500.0,
    ),
    (
        "C1/C6/N2 0.1/0.7/0.2 250 K / 10000 (ternary, liquid-like)",
        ("Methane", "n-Hexane", "Nitrogen"),
        0.0,
        [0.1, 0.7, 0.2],
        250.0,
        10000.0,
    ),
]

_TOLERANCE = 1e-10


@pytest.mark.parametrize(
    ("label", "components", "kij", "x", "temperature", "density"),
    _STATES,
    ids=[state[0] for state in _STATES],
)
def test_residual_properties_match_teqp(
    label: str,
    components: tuple[str, ...],
    kij: object,
    x: list[float],
    temperature: float,
    density: float,
) -> None:
    """Case P-1: A^res/RT, Z, P and ln phi against teqp at 14 states."""
    eos = PCSAFTEOS(components=components, kij=kij)  # type: ignore[arg-type]
    model = _teqp_model(components, eos.kij_matrix())
    z_array = np.array(x, dtype=float)

    reference_a_res = float(model.get_Ar00(temperature, density, z_array))
    reference_z = 1.0 + float(model.get_Ar01(temperature, density, z_array))

    a_res = eos.residual_helmholtz(
        temperature_K=temperature, volume_m3=1.0 / density, composition=x
    )
    z_factor = eos.compressibility_factor(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    pressure = eos.pressure_Pa(temperature_K=temperature, density_mol_m3=density, composition=x)

    assert a_res == pytest.approx(reference_a_res, rel=_TOLERANCE, abs=_TOLERANCE)
    assert z_factor == pytest.approx(reference_z, rel=_TOLERANCE, abs=_TOLERANCE)
    assert pressure == pytest.approx(
        reference_z * density * R_J_PER_MOL_K * temperature, rel=_TOLERANCE
    )

    if reference_z <= 0.0:
        # ln phi does not exist here; chemthermo refuses rather than returning
        # a nan, and teqp's fugacity coefficients are non-positive.
        from chemthermo.exceptions import ModelError

        with pytest.raises(ModelError):
            eos.ln_fugacity_coefficients(
                temperature_K=temperature, density_mol_m3=density, composition=x
            )
        return

    reference_ln_phi = np.log(
        np.asarray(model.get_fugacity_coefficients(temperature, density * z_array))
    )
    ln_phi = eos.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=density, composition=x
    )
    np.testing.assert_allclose(ln_phi, reference_ln_phi, rtol=_TOLERANCE, atol=_TOLERANCE)


def test_the_reference_and_the_model_disagree_when_a_parameter_is_changed() -> None:
    """Negative control: the agreement above is not vacuous.

    If ``_teqp_model`` ignored the parameters it is handed, every assertion in
    this file would pass for the wrong reason. Perturbing one sigma by 1 % must
    move ``A^res`` by many orders of magnitude more than the 1e-10 tolerance.
    """
    components = ("n-Hexane",)
    eos = PCSAFTEOS(components=components)
    perturbed = {
        "kind": "PCSAFT",
        "model": {
            "coeffs": [
                {
                    "name": "n-Hexane",
                    "m": _PARAMETERS["n-Hexane"][0],
                    "sigma_Angstrom": _PARAMETERS["n-Hexane"][1] * 1.01,
                    "epsilon_over_k": _PARAMETERS["n-Hexane"][2],
                    "BibTeXKey": "perturbed",
                }
            ],
            "kmat": [[0.0]],
        },
    }
    model = teqp.make_model(perturbed)
    reference = float(model.get_Ar00(300.0, 7700.0, np.array([1.0])))
    actual = eos.residual_helmholtz(temperature_K=300.0, volume_m3=1.0 / 7700.0, composition=[1.0])
    assert abs(actual - reference) > 1e-3


# ---------------------------------------------------------------------------
# Case P-2: pure-component saturation, against teqp's own VLE solver
# ---------------------------------------------------------------------------


def _bisect(f: Callable[[float], float], lo: float, hi: float, iterations: int = 200) -> float:
    """Plain bisection; test-only, deliberately unsophisticated."""
    f_lo, f_hi = f(lo), f(hi)
    if f_lo == 0.0:
        return lo
    if f_hi == 0.0:
        return hi
    assert f_lo * f_hi < 0.0, f"bracket does not straddle a root: {f_lo!r}, {f_hi!r}"
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if f_mid == 0.0:
            return mid
        if f_lo * f_mid < 0.0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)


def _pure_saturation(
    eos: PCSAFTEOS, temperature: float, density_max: float
) -> tuple[float, float, float]:
    """Return ``(Psat, rho_liquid, rho_vapour)`` for a pure fluid.

    A deliberately simple **test-only** density root finder; it is not part of
    the public package (the density solver is the next slice). It scans the
    isotherm for the two spinodal extrema, then bisects twice on ``P(rho) = P``
    and once on ``ln phi(liquid) - ln phi(vapour) = 0``.
    """

    def pressure(rho: float) -> float:
        return eos.pressure_Pa(temperature_K=temperature, density_mol_m3=rho, composition=[1.0])

    def ln_phi(rho: float) -> float:
        return eos.ln_fugacity_coefficients(
            temperature_K=temperature, density_mol_m3=rho, composition=[1.0]
        )[0]

    grid = np.geomspace(1e-2, density_max, 6000)
    curve = np.array([pressure(float(rho)) for rho in grid])
    top = int(np.argmax(np.where(grid < 0.6 * density_max, curve, -np.inf)))
    bottom = top + int(np.argmin(curve[top:]))
    rho_vapour_spinodal = float(grid[top])
    rho_liquid_spinodal = float(grid[bottom])

    def roots(p: float) -> tuple[float, float]:
        rho_v = _bisect(lambda rho: pressure(rho) - p, 1e-8, rho_vapour_spinodal)
        rho_l = _bisect(lambda rho: pressure(rho) - p, rho_liquid_spinodal, density_max)
        return rho_l, rho_v

    def residual(p: float) -> float:
        rho_l, rho_v = roots(p)
        return ln_phi(rho_l) - ln_phi(rho_v)

    p_sat = _bisect(residual, max(float(curve[bottom]), 1e-3), 0.999999 * float(curve[top]))
    rho_l, rho_v = roots(p_sat)
    return p_sat, rho_l, rho_v


@pytest.mark.parametrize(
    "temperature",
    # 400 K is `slow`: a further temperature on the same saturation curve
    # 300 K already checks against teqp's own `pure_VLE_T`; 300 K also has its
    # own coarser sanity check below (`..._is_in_the_right_range_at_300_K`).
    [300.0, pytest.param(400.0, marks=pytest.mark.slow)],
)
def test_pure_hexane_saturation_matches_teqp_pure_vle(temperature: float) -> None:
    """Case P-2: equal fugacity on the two density roots reproduces teqp's VLE.

    ``teqp.pure_VLE_T`` solves the saturation condition with its own Newton
    iteration on the two densities; the test solves it by bisection on the
    pressure. Only the model is shared.
    """
    eos = PCSAFTEOS(components=("n-Hexane",))
    model = _teqp_model(("n-Hexane",), eos.kij_matrix())

    guess_liquid, guess_vapour = (7700.0, 10.0) if temperature < 350.0 else (6800.0, 150.0)
    rho_l_ref, rho_v_ref = model.pure_VLE_T(temperature, guess_liquid, guess_vapour, 200)
    p_ref = (
        rho_l_ref
        * R_J_PER_MOL_K
        * temperature
        * (1.0 + model.get_Ar01(temperature, rho_l_ref, np.array([1.0])))
    )

    p_sat, rho_l, rho_v = _pure_saturation(eos, temperature, 9500.0)

    assert p_sat == pytest.approx(p_ref, rel=1e-6)
    assert rho_l == pytest.approx(rho_l_ref, rel=1e-6)
    assert rho_v == pytest.approx(rho_v_ref, rel=1e-6)

    # The saturation condition itself, restated on chemthermo's own numbers.
    ln_phi_l = eos.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=rho_l, composition=[1.0]
    )[0]
    ln_phi_v = eos.ln_fugacity_coefficients(
        temperature_K=temperature, density_mol_m3=rho_v, composition=[1.0]
    )[0]
    assert ln_phi_l == pytest.approx(ln_phi_v, abs=1e-9)
    assert eos.pressure_Pa(
        temperature_K=temperature, density_mol_m3=rho_l, composition=[1.0]
    ) == pytest.approx(p_sat, rel=1e-9)


@pytest.mark.slow  # the coarser half of the check `..._matches_teqp_pure_vle[300.0]`
# makes by default, on the same saturation state
def test_hexane_saturation_pressure_is_in_the_right_range_at_300_K() -> None:
    """Model-versus-experiment sanity remark, not a validation of the model.

    n-hexane's vapour pressure at 300 K is commonly tabulated near 21.7 kPa
    (unverified against a primary source here). PC-SAFT with the packaged
    parameters gives about 21.86 kPa. This asserts only the order of magnitude,
    so it can never turn a modelling deviation into a red test.
    """
    eos = PCSAFTEOS(components=("n-Hexane",))
    p_sat, _rho_l, _rho_v = _pure_saturation(eos, 300.0, 9500.0)
    assert 1.5e4 < p_sat < 3.0e4
    assert math.isfinite(p_sat)
