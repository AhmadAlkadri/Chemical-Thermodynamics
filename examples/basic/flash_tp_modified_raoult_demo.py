"""Low-pressure TP flash from one tangent plane: modified Raoult's law.

What the model is
-----------------
At low pressure the vapor is close to an ideal gas and the pure-liquid
reference fugacity is just the saturation pressure:

    f_i^V = y_i P                 (phi_i^V = 1)
    f_i^L = x_i gamma_i(x) Psat_i(T)   (phi_i^sat = 1, Poynting = 1)

so equilibrium is **modified Raoult's law**, ``y_i P = x_i gamma_i Psat_i``.
Dividing both by ``x_i P`` puts the two phases on one Gibbs surface with the
same reference ``ln( f_i / (x_i P) )``:

    liquid candidate: ln gamma_i(w) + ln( Psat_i(T) / P )
    vapor candidate:  0

``flash_tp(..., flash_mode="modified-raoult")`` hands both candidates to
Michelsen's tangent-plane test and keeps whichever has the lower Gibbs energy
at each composition - the same rule that picks the minimum-Gibbs root of a
cubic equation of state. One test therefore answers *three* questions at once:
is the feed one phase or two, and if two, is the second phase a vapor or a
second liquid?

What this script shows
----------------------
1. **1-propanol / water** (fully miscible): a subcooled liquid, a vapor-liquid
   split, and a superheated vapor - the same call, three answers. The split is
   checked here against ``y_i P = x_i gamma_i Psat_i`` recomputed from the
   printed numbers.
2. **n-butanol / water** (partially miscible): at 330 K the same call returns a
   *liquid-liquid* tie-line, because the vapor candidate is nowhere the lower
   one. The phases are named ``liquid1`` / ``liquid2`` and ``vapor_fraction``
   is ``None``, exactly as in ``flash_mode="gamma-gamma"``.

Parameters and their provenance
-------------------------------
NRTL binary parameters are the 1-3 and 2-3 pairs of Table 1 of

    S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable phase
    stability analysis for excess Gibbs energy models", Chemical Engineering
    Science 55 (2000) 1785-1796,

where Table 1 attributes them to C. M. McDonald and C. A. Floudas, AIChE
Journal 41 (1995) 1798-1814. Table 1 prints ``G_ij`` rather than ``alpha_ij``;
the alphas below are implied by ``G_ij = exp(-alpha_ij tau_ij)``. They are
written inline so this script runs from a bare install, and they are **not**
the packaged NRTL defaults, which are synthetic placeholders.

**They are fitted to liquid-liquid data and are temperature independent**, so
the vapor-liquid numbers here are an illustration of the *method*, not a
validated correlation of these systems. Any agreement with a handbook boiling
point is partly fortuitous.

Antoine coefficients come from the packaged databank
(``ln( P^sat / bar ) = A - B / (T/K + C)``, Koretsky 2012). Their validity
ranges are enforced: a temperature outside them raises ``InputRangeError``
rather than being extrapolated.

No optional dependency is required. All inputs are SI: temperature in K,
pressure in Pa.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence

import chemthermo as ct

PRESSURE_PA = 101325.0

# Tessier et al. (2000) Table 1, pair 1-3 (n-propanol / water).
PROPANOL_WATER = ("1-Propanol", "Water")
PROPANOL_WATER_TAU = (-0.07149, 2.7425)
PROPANOL_WATER_ALPHA = 0.3

# Tessier et al. (2000) Table 1, pair 2-3 (n-butanol / water).
BUTANOL_WATER = ("n-Butanol", "Water")
BUTANOL_WATER_TAU = (0.90047, 3.51307)
BUTANOL_WATER_ALPHA = 0.48


def _model(names: Sequence[str], tau: tuple[float, float], alpha: float) -> ct.NRTL:
    parameters = ct.NRTLParameters.from_pairs([(names[0], names[1], tau[0], tau[1], alpha, alpha)])
    return ct.NRTL(parameters=parameters)


def _format_vector(label: str, names: Sequence[str], values: Iterable[float]) -> str:
    lines = [label]
    for name, value in zip(names, values):
        lines.append(f"  {name:<12} {value: .8f}")
    return "\n".join(lines)


def _saturation_pressures(names: Sequence[str], temperature_K: float) -> list[float]:
    """Psat from the databank Antoine record, recomputed here from A, B, C."""
    values = []
    for name in names:
        antoine = ct.Component.from_database(name).antoine
        assert antoine is not None
        values.append(math.exp(antoine.A - antoine.B / (temperature_K + antoine.C)) * 1.0e5)
    return values


def _report(
    names: Sequence[str],
    model: ct.NRTL,
    z: Sequence[float],
    temperature_K: float,
    note: str,
) -> ct.FlashResult:
    mixture = ct.Mixture.from_database(list(names), list(z), normalize=True)
    result = ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        flash_mode="modified-raoult",
    )
    diagnostics = result.diagnostics

    print("=" * 72)
    print("Modified-Raoult TP flash (NRTL liquid + ideal vapor, Antoine reference)")
    print(f"System: {names[0]} / {names[1]}   -   {note}")
    print(f"T [K]: {temperature_K:.4f}   P [Pa]: {PRESSURE_PA:.5g}")
    print(_format_vector("Feed z (mole):", names, z))
    psat = _saturation_pressures(names, temperature_K)
    print(_format_vector("Psat [Pa]:", names, psat))
    print(
        f"Antoine validity window for this mixture: "
        f"[{diagnostics['antoine_valid_Tmin_K']}, {diagnostics['antoine_valid_Tmax_K']}] K"
    )
    print(f"Feed candidate (lowest Gibbs at z): {diagnostics['feed_branch']}")
    print(
        f"Feed stability: {diagnostics['stability_status']}, "
        f"tpd_min = {float(diagnostics['tpd_min']):+.8e}"
    )
    print(f"Phase regime: {diagnostics['phase_regime']}")
    print(f"Phases returned: {result.phase_names()}")
    print(f"vapor_fraction: {result.vapor_fraction}")

    if len(result.phase_names()) == 1:
        print("Single phase: the feed was found stable, so no split was attempted.")
        return result

    print(f"Incipient phase candidate: {diagnostics['incipient_phase']}")
    for name in result.phase_names():
        print(
            _format_vector(
                f"{name} (phase fraction {result.phase_fractions[name]:.8f}):",
                names,
                result.phases[name].composition.fractions,
            )
        )
    print(
        "Stages: "
        f"{diagnostics['ssi_iterations']} successive substitution + "
        f"{diagnostics['second_order_iterations']} second order "
        f"(converged in: {diagnostics['converged_stage']})"
    )
    print("Verification:")
    print(f"  mass_balance_residual   = {float(diagnostics['mass_balance_residual']):.3e}")
    print(f"  equilibrium_residual    = {float(diagnostics['equilibrium_residual']):.3e}")
    print(
        f"  delta_g_split_rt        = {float(diagnostics['delta_g_split_rt']):+.6e}  (must be < 0)"
    )
    print("Post-split stability (each phase re-tested against BOTH candidates):")
    for name in result.phase_names():
        print(
            f"  {name:<8} {diagnostics[f'phase_stability_{name}']:<10}"
            f" tpd_min = {float(diagnostics[f'phase_stability_tpd_min_{name}']):+.3e}"
        )
    print(f"  post_split_status = {diagnostics['post_split_status']}")
    return result


def main() -> None:
    propanol = _model(PROPANOL_WATER, PROPANOL_WATER_TAU, PROPANOL_WATER_ALPHA)

    _report(PROPANOL_WATER, propanol, (0.50, 0.50), 330.0, "subcooled: one liquid")
    vapor_liquid = _report(
        PROPANOL_WATER, propanol, (0.50, 0.50), 361.0, "inside the two-phase band"
    )
    _report(PROPANOL_WATER, propanol, (0.50, 0.50), 380.0, "superheated: one vapor")

    print("=" * 72)
    print("Modified Raoult's law, recomputed from the printed vapor-liquid split")
    x = vapor_liquid.phases["liquid"].composition.fractions
    y = vapor_liquid.phases["vapor"].composition.fractions
    gamma = propanol.activity_coefficients(
        mixture=ct.Mixture.from_database(list(PROPANOL_WATER), list(x)),
        temperature_K=361.0,
        composition=list(x),
    )
    psat = _saturation_pressures(PROPANOL_WATER, 361.0)
    worst = 0.0
    for name, xi, yi, gi, pi in zip(PROPANOL_WATER, x, y, gamma, psat):
        left = yi * PRESSURE_PA
        right = xi * gi * pi
        worst = max(worst, abs(left - right))
        print(
            f"  {name:<12} y_i P = {left:14.6f} Pa    x_i gamma_i Psat_i = {right:14.6f} Pa"
            f"    gamma_i = {gi:.6f}"
        )
    print(f"  worst |y_i P - x_i gamma_i Psat_i| = {worst:.3e} Pa")

    butanol = _model(BUTANOL_WATER, BUTANOL_WATER_TAU, BUTANOL_WATER_ALPHA)
    _report(
        BUTANOL_WATER,
        butanol,
        (0.20, 0.80),
        330.0,
        "partially miscible: the same call returns a liquid-liquid tie-line",
    )

    print("=" * 72)
    print(
        "Notes and limits\n"
        "  * Ideal vapor: no phi^V, so this mode is for low pressure only.\n"
        "  * No Poynting correction and no phi^sat in the reference fugacity.\n"
        "  * The Antoine validity range is enforced; outside it the call raises\n"
        "    InputRangeError instead of extrapolating a vapor-pressure fit.\n"
        "  * 'stable' means no negative tangent-plane distance was found from\n"
        "    the deterministic trial set, not a global proof.\n"
        "  * This mode returns up to FlashSettings.max_phases (default 3) phases:\n"
        "    a phase set whose post-split test fails has the incipient phase\n"
        "    added and is re-solved (ADR-0011). See\n"
        "    examples/basic/flash_tp_vlle_demo.py for a three-phase answer and\n"
        "    examples/validation/10_modified_raoult_water_butanol.py for the\n"
        "    window below a three-phase temperature, where a phase is removed.\n"
        "  * flash_mode='gamma-phi' is deprecated in favour of this mode: it\n"
        "    multiplies gamma by an EOS liquid phi (double-counting the liquid's\n"
        "    nonideality) and carries no pure-liquid reference fugacity at all."
    )


if __name__ == "__main__":
    main()
