"""Compressibility-based phase identity (ADR-0017).

Covers ``EquationOfState.phase_identity`` for Peng-Robinson and PC-SAFT: the
kappa criterion's separation between liquid and vapor roots, the two
motivating single-root states named "liquid" instead of the pre-ADR-0017
"vapor" tie-break, a superheated state that still says "vapor", the
higher-density-is-"liquid" invariant on converged splits, and the documented
default fallback for a model that does not implement ``phase_identity``.

The full-grid measurements this module's smaller tables are drawn from (all
144 Peng-Robinson states of ``tests/test_flash_phase_detection.py``'s grid,
all 188 PC-SAFT states of ``tests/validation/test_flash_split_robustness_pcsaft.py``'s
Case F-4 grid) are reported in ADR-0017 and validation Case F-5
(``.agents/brain/validation-cases.md``); re-scanning both grids here on every
``pytest -q`` would reintroduce the runtime regression this slice fixes (see
``tests/validation/test_flash_split_robustness_pcsaft_subset.py`` and
``examples/validation/15_flash_split_robustness.py --full`` for the full-grid,
opt-in routes). This module's own tables use the full 144-state Peng-Robinson
grid (cheap: no density-root scan) and a small, fixed PC-SAFT set.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np
import pytest

import chemthermo as ct
from chemthermo.core import Mixture
from chemthermo.models.base import KAPPA_LIQUID_THRESHOLD

PR = ct.PengRobinsonEOS()
#: Matches ``chemthermo.models.peng_robinson.R_J_PER_MOL_K`` (CODATA, rounded
#: the same way); duplicated here so the cross-check below reads no private
#: constant from the package under test.
_R = 8.314462618


def _mixture(names: Sequence[str], z: Sequence[float]) -> ct.Mixture:
    return ct.Mixture.from_database(list(names), list(z), normalize=True)


# ---------------------------------------------------------------------------
# An independent (public-data-only) Peng-Robinson P(V), for the analytic vs.
# finite-difference cross-check and the kappa table below. Built from
# `Component.tc_k` / `pc_pa` / `omega` alone - never from the package's
# private `PengRobinsonEOS._mixture_parameters` - so agreement with
# `phase_identity` is a genuine check of the shipped formula, not a
# tautology.
# ---------------------------------------------------------------------------


def _pr_ab(mixture: Mixture, temperature_K: float) -> tuple[np.ndarray, np.ndarray]:
    a_i: list[float] = []
    b_i: list[float] = []
    for component in mixture.components:
        tr = temperature_K / component.tc_k
        kappa = 0.37464 + 1.54226 * component.omega - 0.26992 * component.omega**2
        alpha = (1.0 + kappa * (1.0 - math.sqrt(tr))) ** 2
        a_i.append(0.45724 * _R**2 * component.tc_k**2 / component.pc_pa * alpha)
        b_i.append(0.07780 * _R * component.tc_k / component.pc_pa)
    return np.array(a_i), np.array(b_i)


def _pr_mixture_ab(
    mixture: Mixture, temperature_K: float, composition: Sequence[float]
) -> tuple[float, float]:
    a_i, b_i = _pr_ab(mixture, temperature_K)
    y = np.asarray(composition, dtype=float)
    a_mix = float(y @ np.sqrt(np.outer(a_i, a_i)) @ y)
    b_mix = float(y @ b_i)
    return a_mix, b_mix


def _pr_pressure(temperature_K: float, volume: float, a_mix: float, b_mix: float) -> float:
    return _R * temperature_K / (volume - b_mix) - a_mix / (
        volume**2 + 2.0 * b_mix * volume - b_mix**2
    )


def _pr_volume(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    composition: Sequence[float],
    phase: str,
) -> float:
    Z = PR.compressibility_factor(
        mixture=mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        composition=composition,
        phase=phase,
    )
    return Z * _R * temperature_K / pressure_Pa


def _pr_kappa_finite_difference(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    composition: Sequence[float],
    phase: str,
) -> float:
    volume = _pr_volume(mixture, temperature_K, pressure_Pa, composition, phase)
    a_mix, b_mix = _pr_mixture_ab(mixture, temperature_K, composition)
    step = volume * 1e-6
    dP_dV = (
        _pr_pressure(temperature_K, volume + step, a_mix, b_mix)
        - _pr_pressure(temperature_K, volume - step, a_mix, b_mix)
    ) / (2.0 * step)
    return -pressure_Pa / (volume * dP_dV)


def _pr_kappa_analytic(
    mixture: Mixture,
    temperature_K: float,
    pressure_Pa: float,
    composition: Sequence[float],
    phase: str,
) -> float:
    volume = _pr_volume(mixture, temperature_K, pressure_Pa, composition, phase)
    a_mix, b_mix = _pr_mixture_ab(mixture, temperature_K, composition)
    denominator = volume**2 + 2.0 * b_mix * volume - b_mix**2
    dP_dV = (
        -_R * temperature_K / (volume - b_mix) ** 2
        + a_mix * (2.0 * volume + 2.0 * b_mix) / denominator**2
    )
    return -pressure_Pa / (volume * dP_dV)


#: The Peng-Robinson grid of ``tests/test_flash_phase_detection.py``,
#: duplicated (not imported) so this module stays self-contained - the
#: convention already used by ``tests/test_flash_refactor_bit_identity.py``.
PR_GRID_MIXTURES: tuple[tuple[tuple[str, ...], tuple[float, ...]], ...] = (
    (("Methane", "Ethane"), (0.5, 0.5)),
    (("Methane", "Propane"), (0.7, 0.3)),
    (("Ethane", "n-Heptane"), (0.7, 0.3)),
    (("Methane", "n-Pentane"), (0.6, 0.4)),
    (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2)),
    (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3)),
)
PR_GRID_T_K = (170.0, 200.0, 240.0, 280.0, 320.0, 360.0)
PR_GRID_P_PA = (2.0e5, 1.0e6, 3.0e6, 8.0e6)


def _pr_grid_roots() -> list[tuple[str, float]]:
    """``(class, kappa)`` for every liquid/vapor root PR reports over the grid.

    ``class`` is ``"single-liquid"`` / ``"single-vapor"`` / ``"two-liquid"`` /
    ``"two-vapor"`` - the phi-phi verdict and, for a single-phase result,
    which side of the threshold it landed on.
    """
    rows: list[tuple[str, float]] = []
    for names, z in PR_GRID_MIXTURES:
        mixture = _mixture(names, z)
        for temperature_K in PR_GRID_T_K:
            for pressure_Pa in PR_GRID_P_PA:
                try:
                    result = ct.flash_tp(
                        mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=PR
                    )
                except ct.ConvergenceError:
                    continue
                if len(result.phase_names()) == 2:
                    x = result.phases["liquid"].composition.fractions
                    y = result.phases["vapor"].composition.fractions
                    rows.append(
                        (
                            "two-liquid",
                            _pr_kappa_finite_difference(
                                mixture, temperature_K, pressure_Pa, x, "liquid"
                            ),
                        )
                    )
                    rows.append(
                        (
                            "two-vapor",
                            _pr_kappa_finite_difference(
                                mixture, temperature_K, pressure_Pa, y, "vapor"
                            ),
                        )
                    )
                else:
                    name = result.phase_names()[0]
                    kappa = _pr_kappa_finite_difference(
                        mixture, temperature_K, pressure_Pa, mixture.fractions, name
                    )
                    rows.append((f"single-{name}", kappa))
    return rows


# ---------------------------------------------------------------------------
# 1) The analytic dP/dV `PengRobinsonEOS.phase_identity` uses agrees with an
#    independent finite difference of the same (independently rebuilt) P(V).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("names", "z", "temperature_K", "pressure_Pa", "phase"),
    [
        (("Methane", "Ethane"), (0.5, 0.5), 240.0, 3.0e6, "liquid"),
        (("Methane", "Ethane"), (0.5, 0.5), 240.0, 3.0e6, "vapor"),
        (("Methane", "Ethane"), (0.5, 0.5), 170.0, 3.0e6, "liquid"),
        (("Methane", "Ethane"), (0.5, 0.5), 450.0, 1.0e5, "vapor"),
        (("Ethane", "n-Heptane"), (0.7, 0.3), 360.0, 8.0e6, "liquid"),
        (("Ethane", "n-Heptane"), (0.7, 0.3), 200.0, 1.0e6, "vapor"),
        (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2), 240.0, 3.0e6, "liquid"),
        (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2), 240.0, 3.0e6, "vapor"),
        (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3), 320.0, 8.0e6, "liquid"),
    ],
)
def test_pr_analytic_dP_dV_matches_finite_difference(
    names: tuple[str, ...],
    z: tuple[float, ...],
    temperature_K: float,
    pressure_Pa: float,
    phase: str,
) -> None:
    mixture = _mixture(names, z)
    analytic = _pr_kappa_analytic(mixture, temperature_K, pressure_Pa, z, phase)
    finite_difference = _pr_kappa_finite_difference(mixture, temperature_K, pressure_Pa, z, phase)
    relative = abs(analytic - finite_difference) / max(abs(finite_difference), 1e-300)
    assert relative < 1e-8, (analytic, finite_difference, relative)

    # And the identity `phase_identity` actually returns agrees with what
    # this independently-derived kappa implies.
    identity = PR.phase_identity(
        mixture=mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        composition=z,
        phase=phase,
    )
    expected = "liquid" if finite_difference < KAPPA_LIQUID_THRESHOLD else "vapor"
    assert identity == expected


# ---------------------------------------------------------------------------
# 2) Kappa separation over the in-repo Peng-Robinson grid.
# ---------------------------------------------------------------------------


def test_pr_grid_kappa_separates_liquid_and_vapor_roots() -> None:
    """Every liquid root's kappa is below the threshold, every vapor root's above.

    Measured over this module's 144-state grid (identical to
    ``tests/test_flash_phase_detection.py``'s): liquid kappa in
    [1.1e-04, 2.3e-01], vapor kappa in [7.4e-01, 1.8e+00] - a wide margin on
    both sides of ``KAPPA_LIQUID_THRESHOLD = 0.5``. See ADR-0017 for the full
    table (identical numbers, computed once and quoted there).
    """
    rows = _pr_grid_roots()
    liquid = [kappa for cls, kappa in rows if cls.endswith("liquid")]
    vapor = [kappa for cls, kappa in rows if cls.endswith("vapor")]
    assert len(liquid) >= 60, len(liquid)
    assert len(vapor) >= 60, len(vapor)
    assert max(liquid) < KAPPA_LIQUID_THRESHOLD, max(liquid)
    assert min(vapor) >= KAPPA_LIQUID_THRESHOLD, min(vapor)
    # A comfortable margin, not a coincidence at the boundary.
    assert max(liquid) < 0.3, max(liquid)
    assert min(vapor) > 0.7, min(vapor)


# ---------------------------------------------------------------------------
# 3) Kappa separation over a fixed PC-SAFT set (Case F-4 grid; see the module
#    docstring for why this is not the full 188-state scan).
# ---------------------------------------------------------------------------

#: (components, z1, T_K, P_Pa) states from the Case F-4 grid
#: (``tests/validation/test_flash_split_robustness_pcsaft.py``), a mix of
#: single-phase (dense, one density root) and two-phase states.
_PCSAFT_STATES: tuple[tuple[tuple[str, str], float, float, float], ...] = (
    (("Carbon dioxide", "n-Decane"), 0.9, 230.0, 2.5e6),  # motivating state 1
    (("Methane", "n-Hexane"), 0.5, 170.0, 2.0e6),  # motivating state 2
    (("Carbon dioxide", "n-Decane"), 0.9, 240.0, 1.0e6),  # the F-4 reference split
    (("Methane", "n-Hexane"), 0.5, 200.0, 3.5e6),
    (("Carbon dioxide", "n-Decane"), 0.6, 260.0, 1.0e6),
)


def _pcsaft_kappa_finite_difference(
    eos: ct.PCSAFTEOS, temperature_K: float, density_mol_m3: float, composition: Sequence[float]
) -> float:
    """Finite difference of the public ``pressure_Pa(T, rho, x)`` - no private access."""
    step = density_mol_m3 * 1e-6
    p_plus = eos.pressure_Pa(
        temperature_K=temperature_K, density_mol_m3=density_mol_m3 + step, composition=composition
    )
    p_minus = eos.pressure_Pa(
        temperature_K=temperature_K, density_mol_m3=density_mol_m3 - step, composition=composition
    )
    p0 = eos.pressure_Pa(
        temperature_K=temperature_K, density_mol_m3=density_mol_m3, composition=composition
    )
    dP_drho = (p_plus - p_minus) / (2.0 * step)
    return p0 / (density_mol_m3 * dP_drho)


@pytest.fixture(scope="module")
def _pcsaft_results() -> list[
    tuple[tuple[str, str], float, float, float, Mixture, ct.PCSAFTEOS, ct.FlashResult]
]:
    """Flash every ``_PCSAFT_STATES`` state once; shared by the two tests below.

    ``flash_tp`` (the stability search plus the seeded split) is what is
    expensive here, not a ``density_roots`` call on an already-known
    composition, so sharing this - not the individual kappa/density values -
    is what keeps this module cheap enough for the default ``pytest -q`` run.
    """
    results = []
    for components, z1, temperature_K, pressure_Pa in _PCSAFT_STATES:
        mixture = ct.Mixture.from_database(list(components), [z1, 1.0 - z1], normalize=True)
        eos = ct.PCSAFTEOS(components=tuple(components))
        result = ct.flash_tp(
            mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=ct.PCSAFTEOS()
        )
        results.append((components, z1, temperature_K, pressure_Pa, mixture, eos, result))
    return results


def test_pcsaft_kappa_separates_liquid_and_vapor_roots(_pcsaft_results) -> None:
    """As above, for PC-SAFT: measured over the 188-state grid liquid kappa in
    [4.9e-04, 8.5e-03], vapor kappa in [1.04, 1.87] (ADR-0017 quotes the full
    table); this test re-checks a fixed sample with an independent
    finite-difference route.
    """
    liquid: list[float] = []
    vapor: list[float] = []
    for _components, _z1, temperature_K, pressure_Pa, mixture, eos, result in _pcsaft_results:
        if len(result.phase_names()) == 2:
            x = result.phases["liquid"].composition.fractions
            y = result.phases["vapor"].composition.fractions
            rho_x = eos.density_roots(
                temperature_K=temperature_K, pressure_Pa=pressure_Pa, composition=x
            )[-1]
            rho_y = eos.density_roots(
                temperature_K=temperature_K, pressure_Pa=pressure_Pa, composition=y
            )[0]
            liquid.append(_pcsaft_kappa_finite_difference(eos, temperature_K, rho_x, x))
            vapor.append(_pcsaft_kappa_finite_difference(eos, temperature_K, rho_y, y))
        else:
            name = result.phase_names()[0]
            z = mixture.fractions
            roots = eos.density_roots(
                temperature_K=temperature_K, pressure_Pa=pressure_Pa, composition=z
            )
            rho = roots[-1] if name == "liquid" else roots[0]
            kappa = _pcsaft_kappa_finite_difference(eos, temperature_K, rho, z)
            (liquid if name == "liquid" else vapor).append(kappa)

    assert liquid, "expected at least one liquid root in _PCSAFT_STATES"
    assert max(liquid) < KAPPA_LIQUID_THRESHOLD, max(liquid)
    if vapor:
        assert min(vapor) >= KAPPA_LIQUID_THRESHOLD, min(vapor)


# ---------------------------------------------------------------------------
# 4) The two motivating states (slice declaration) and the superheated state
#    that must stay "vapor".
# ---------------------------------------------------------------------------


def test_dense_co2_n_decane_is_now_a_liquid() -> None:
    """CO2/n-decane z=(0.9, 0.1) at 230 K, 2.5 MPa: one density root, dense."""
    mixture = _mixture(("Carbon dioxide", "n-Decane"), (0.9, 0.1))
    result = ct.flash_tp(mixture, temperature_K=230.0, pressure_Pa=2.5e6, eos=ct.PCSAFTEOS())
    assert result.phase_names() == ["liquid"]
    assert result.vapor_fraction == 0.0
    assert result.diagnostics["phase_label_method"] == "compressibility"


def test_dense_methane_n_hexane_is_now_a_liquid() -> None:
    """Methane/n-hexane z=(0.5, 0.5) at 170 K, 2 MPa: one density root, dense."""
    mixture = _mixture(("Methane", "n-Hexane"), (0.5, 0.5))
    result = ct.flash_tp(mixture, temperature_K=170.0, pressure_Pa=2.0e6, eos=ct.PCSAFTEOS())
    assert result.phase_names() == ["liquid"]
    assert result.vapor_fraction == 0.0
    assert result.diagnostics["phase_label_method"] == "compressibility"


def test_superheated_methane_ethane_stays_vapor() -> None:
    """450 K, 1 bar methane/ethane: dilute, kappa near 1, correctly "vapor"."""
    mixture = _mixture(("Methane", "Ethane"), (0.5, 0.5))
    result = ct.flash_tp(mixture, temperature_K=450.0, pressure_Pa=1.0e5, eos=PR)
    assert result.phase_names() == ["vapor"]
    assert result.vapor_fraction == 1.0
    assert result.diagnostics["phase_label_method"] == "compressibility"


# ---------------------------------------------------------------------------
# 5) Splits: the phase called "liquid" always has the higher density.
# ---------------------------------------------------------------------------


def test_pr_two_phase_liquid_has_the_higher_density() -> None:
    checked = 0
    for names, z in PR_GRID_MIXTURES:
        mixture = _mixture(names, z)
        for temperature_K in PR_GRID_T_K:
            for pressure_Pa in PR_GRID_P_PA:
                try:
                    result = ct.flash_tp(
                        mixture, temperature_K=temperature_K, pressure_Pa=pressure_Pa, eos=PR
                    )
                except ct.ConvergenceError:
                    continue
                if len(result.phase_names()) != 2:
                    continue
                checked += 1
                x = result.phases["liquid"].composition.fractions
                y = result.phases["vapor"].composition.fractions
                z_liquid = PR.compressibility_factor(
                    mixture=mixture,
                    temperature_K=temperature_K,
                    pressure_Pa=pressure_Pa,
                    composition=x,
                    phase="liquid",
                )
                z_vapor = PR.compressibility_factor(
                    mixture=mixture,
                    temperature_K=temperature_K,
                    pressure_Pa=pressure_Pa,
                    composition=y,
                    phase="vapor",
                )
                # Density is proportional to 1/Z at fixed (T, P): lower Z is
                # higher density.
                assert z_liquid < z_vapor, (names, temperature_K, pressure_Pa, z_liquid, z_vapor)
    assert checked >= 30, checked


def test_pcsaft_two_phase_liquid_has_the_higher_density(_pcsaft_results) -> None:
    checked = 0
    for components, z1, temperature_K, pressure_Pa, _mixture, eos, result in _pcsaft_results:
        if len(result.phase_names()) != 2:
            continue
        checked += 1
        x = result.phases["liquid"].composition.fractions
        y = result.phases["vapor"].composition.fractions
        rho_liquid = eos.density_roots(
            temperature_K=temperature_K, pressure_Pa=pressure_Pa, composition=x
        )[-1]
        rho_vapor = eos.density_roots(
            temperature_K=temperature_K, pressure_Pa=pressure_Pa, composition=y
        )[0]
        assert rho_liquid > rho_vapor, (components, z1, temperature_K, pressure_Pa)
    assert checked >= 1, "expected at least one two-phase state in _PCSAFT_STATES"


# ---------------------------------------------------------------------------
# 6) The documented default fallback for a model without `phase_identity`.
# ---------------------------------------------------------------------------


class _NoIdentityEOS(ct.EquationOfState):
    """Wraps Peng-Robinson but does not override `phase_identity` - the
    default (`None`, "not implemented") must make `flash_tp` behave exactly
    as it did before ADR-0017."""

    name = "no-identity-eos"

    def fugacity_coefficients(self, **kwargs):  # type: ignore[override]
        return PR.fugacity_coefficients(**kwargs)


def test_a_model_without_phase_identity_keeps_the_pre_adr_0017_behavior() -> None:
    eos = _NoIdentityEOS()
    mixture = _mixture(("Methane", "Ethane"), (0.5, 0.5))

    assert (
        eos.phase_identity(
            mixture=mixture,
            temperature_K=170.0,
            pressure_Pa=3.0e6,
            composition=(0.5, 0.5),
            phase="liquid",
        )
        is None
    )

    # The same dense, single-root state that PengRobinsonEOS itself now names
    # "liquid" (see the grid above) stays the pre-ADR-0017 "vapor" tie-break
    # when the model cannot measure an identity.
    result = ct.flash_tp(mixture, temperature_K=170.0, pressure_Pa=3.0e6, eos=eos)
    assert result.phase_names() == ["vapor"]
    assert result.vapor_fraction == 1.0
    assert result.diagnostics["phase_label_method"] == "tie-break"
