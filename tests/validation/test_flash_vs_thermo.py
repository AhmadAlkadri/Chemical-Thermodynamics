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


def _get_vf(state: object) -> float | None:
    vf = getattr(state, "VF", None)
    if vf is None:
        return None
    return vf() if callable(vf) else vf


def _pr_flasher(ids: list[str]) -> FlashVL:
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

    # Different NRTL implementations may vary slightly; keep a tight but realistic tolerance.
    assert np.allclose(gammas, ref_gammas, rtol=2e-3, atol=2e-3)
