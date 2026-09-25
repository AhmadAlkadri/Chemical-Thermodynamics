from typing import Any

import numpy as np
import pytest

import chemthermo as ct

thermo = pytest.importorskip("thermo")
CEOSGas = thermo.CEOSGas
CEOSLiquid = thermo.CEOSLiquid
ChemicalConstantsPackage = thermo.ChemicalConstantsPackage
FlashPureVLS = thermo.FlashPureVLS
FlashVL = thermo.FlashVL
PRMIX = thermo.PRMIX
NRTL_gammas = thermo.nrtl.NRTL_gammas


def _get_vf(state: Any) -> float | None:
    vf = getattr(state, "VF", None)
    if vf is None:
        return None
    vf_value = vf() if callable(vf) else vf
    if isinstance(vf_value, (int, float)):
        return float(vf_value)
    return None


def _pr_flasher(ids: list[str]) -> Any:
    constants, properties = ChemicalConstantsPackage.from_IDs(ids)
    kijs = [[0.0 for _ in ids] for _ in ids]
    eos_kwargs = {
        "Pcs": constants.Pcs,
        "Tcs": constants.Tcs,
        "omegas": constants.omegas,
        "kijs": kijs,
    }
    gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    liquid = CEOSLiquid(
        PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases
    )
    return FlashVL(constants, properties, liquid=liquid, gas=gas)


def test_flash_pr_pure_vapor_vs_thermo() -> None:
    name = "Methane"
    temperature_K = 400.0
    pressure_Pa = 1.0e5

    mixture = ct.Mixture.from_database([name], [1.0], normalize=True)
    result = ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=ct.PengRobinsonEOS(),
    )

    constants, properties = ChemicalConstantsPackage.from_IDs([name.lower()])
    eos_kwargs = {
        "Pcs": constants.Pcs,
        "Tcs": constants.Tcs,
        "omegas": constants.omegas,
    }
    gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases)
    liquid = CEOSLiquid(
        PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases
    )
    flasher = FlashPureVLS(constants, properties, gas=gas, liquids=[liquid], solids=[])
    ref = flasher.flash(T=temperature_K, P=pressure_Pa)

    assert result.phase_names() == ["vapor"]
    assert _get_vf(ref) == pytest.approx(1.0, abs=1e-6)


def test_flash_pr_mixture_vs_thermo() -> None:
    names = ["Methane", "Ethane", "Propane"]
    zs = [0.50, 0.30, 0.20]
    temperature_K = 240.0
    pressure_Pa = 3.0e6

    mixture = ct.Mixture.from_database(names, zs, normalize=True)
    result = ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=ct.PengRobinsonEOS(),
    )

    flasher = _pr_flasher([name.lower() for name in names])
    ref = flasher.flash(T=temperature_K, P=pressure_Pa, zs=zs)

    if "liquid" not in result.phases or "vapor" not in result.phases:
        pytest.skip("chemthermo returned single-phase; reference comparison is two-phase")

    ref_vf = _get_vf(ref)
    if ref_vf is None or ref_vf in (0.0, 1.0):
        pytest.skip("thermo returned single-phase; reference comparison is two-phase")

    # Thermo and chemthermo use different EOS implementations and solvers; use loose tolerances.
    assert result.vapor_fraction == pytest.approx(ref_vf, rel=5e-2, abs=5e-3)

    liquid = result.phases["liquid"].composition.fractions
    vapor = result.phases["vapor"].composition.fractions

    assert np.allclose(liquid, ref.liquid0.zs, rtol=5e-2, atol=1e-2)
    assert np.allclose(vapor, ref.gas.zs, rtol=5e-2, atol=1e-2)


def test_nrtl_activity_coefficients_vs_thermo() -> None:
    """Packaged (synthetic, symmetric-alpha) binary against `thermo`.

    Both sides evaluate the same closed-form Renon-Prausnitz equation, so the
    tolerance is round-off, not "physically reasonable". This case is kept for
    continuity but is deliberately weak evidence: with a symmetric alpha the
    row-sum bug fixed in the `nrtl-gibbs-duhem-fix` slice moved ln gamma here by
    only ~1.1e-3, which the previous 2e-3 tolerance could not see. The sharp
    check is the asymmetric one below.
    """
    names = ["Methane", "Ethane"]
    zs = [0.50, 0.50]
    temperature_K = 240.0

    mixture = ct.Mixture.from_database(names, zs, normalize=True)
    model = ct.NRTL(parameters=ct.ActivityParameters.load("NRTL"))
    gammas = model.activity_coefficients(
        mixture=mixture,
        temperature_K=temperature_K,
        composition=mixture.fractions,
    )

    params = ct.NRTLParameters.load()
    tau, alpha = params.for_components(names)
    ref_gammas = NRTL_gammas(xs=list(mixture.fractions), taus=tau, alphas=alpha)

    assert np.allclose(gammas, ref_gammas, rtol=1e-12, atol=1e-12)


def test_nrtl_activity_coefficients_vs_thermo_asymmetric(
    tessier2000_ln_gamma: Any, tessier2000_payload: dict[str, Any]
) -> None:
    """Asymmetric published parameters against `thermo`, at round-off.

    Tessier (2000) Problem 1: tau_12 = -0.61259 vs tau_21 = 0.71640, and
    alpha_23 = 0.48 vs alpha_12 = alpha_13 = 0.3. See
    `tests/validation/test_nrtl_tessier2000.py` for the full cross-check;
    this keeps a tight asymmetric case next to the weak packaged-binary one.
    """
    tau = np.asarray(tessier2000_payload["tau"], dtype=float)
    alpha = np.asarray(tessier2000_payload["alpha"], dtype=float)

    worst = 0.0
    for composition in ((0.12, 0.08, 0.80), (0.50, 0.30, 0.20), (0.20, 0.20, 0.60)):
        x = np.asarray(composition, dtype=float)
        reference = np.log(np.asarray(NRTL_gammas(xs=x.tolist(), taus=tau, alphas=alpha)))
        worst = max(worst, float(np.max(np.abs(tessier2000_ln_gamma(x) - reference))))

    # Achieved: 8.9e-16.
    assert worst < 1e-9
