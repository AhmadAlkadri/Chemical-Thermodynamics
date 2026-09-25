"""TP flash with PC-SAFT (ADR-0015): methane / n-hexane at 300 K.

What this shows
---------------
``PCSAFTEOS`` now implements the same ``EquationOfState`` interface that
Peng-Robinson does, so it can be handed straight to ``stability_tp`` and
``flash_tp``::

    ct.flash_tp(mixture, temperature_K=..., pressure_Pa=..., eos=PCSAFTEOS())

Nothing in the stability or flash solvers changed for this. The only new piece
is the density root solver: PC-SAFT is written in ``(T, rho, x)``, so a call
that names a *pressure* and a *phase* has to solve
``P_model(T, rho, x) = P`` first. ``PCSAFTEOS.density_roots`` exposes that root
set, and ``phase="vapor"`` / ``"liquid"`` pick its lowest / highest density -
the Peng-Robinson convention verbatim.

``PCSAFTEOS`` also implements ``phase_identity`` (ADR-0017): a dense state with
only one density root is still named "liquid" or "vapor" from a measured
compressibility ratio (``kappa = P / (rho dP/drho)``, 1 for an ideal gas and
well below 1 for a liquid), not from a vapor-first Gibbs tie-break. Section 3
below is exactly that case.

The script prints, at 300 K:

1. the density roots: the feed at three pressures, where the model has only
   one and the two phase labels therefore mean the same thing, and pure
   n-hexane at its saturation pressure, where it has two and equal fugacity on
   them is what "saturation" means;
2. a two-phase flash at 3 MPa with its phase densities and the three residuals
   that make a split an answer rather than a fixed point;
3. the same feed at 8 MPa, where the tangent-plane test finds the compressed
   liquid stable and one phase comes back;
4. the same two-phase state under Peng-Robinson, for contrast. **The two
   models disagree, and neither line is evidence about the other**: PC-SAFT is
   validated against teqp in ``examples/validation/14_pcsaft_flash_vs_teqp.py``
   and Peng-Robinson against ``thermo`` elsewhere; this is a difference between
   two models, not an error in one.

Limits worth knowing before using this: the model is **non-associating** (no
alcohols, water, acids), ``phi-phi`` flash stops at **two phases**, and any
``kij`` you pass is yours to justify - the packaged parameter set carries pure
components only.

Needs no optional dependency.
"""

from __future__ import annotations

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS

COMPONENTS = ["Methane", "n-Hexane"]
FEED = [0.30, 0.70]
TEMPERATURE_K = 300.0
#: n-hexane saturation pressure at 300 K from this model, cross-checked
#: against teqp (validation Cases P-2 / P-3).
HEXANE_PSAT_PA = 21858.084278856164


def _mixture() -> ct.Mixture:
    return ct.Mixture.from_database(COMPONENTS, FEED)


def _print_roots(eos: PCSAFTEOS, mixture: ct.Mixture, pressure: float) -> None:
    roots = eos.density_roots(
        temperature_K=TEMPERATURE_K,
        pressure_Pa=pressure,
        composition=list(mixture.fractions),
        mixture=mixture,
    )
    formatted = ", ".join(f"{value:,.4f}" for value in roots)
    note = "one root: 'vapor' and 'liquid' name the same state" if len(roots) == 1 else ""
    print(f"  P = {pressure / 1e6:>5.2f} MPa -> {len(roots)} root(s): {formatted} mol/m^3  {note}")


def _print_phases(result: ct.FlashResult, eos: PCSAFTEOS, mixture: ct.Mixture) -> None:
    for name, phase in result.phases.items():
        fractions = list(phase.composition.fractions)
        amounts = ", ".join(
            f"{component:<10} {fraction: .6f}"
            for component, fraction in zip(mixture.component_names, fractions)
        )
        roots = eos.density_roots(
            temperature_K=result.temperature_K,
            pressure_Pa=result.pressure_Pa,
            composition=fractions,
            mixture=mixture,
        )
        density = roots[0] if name == "vapor" else roots[-1]
        print(
            f"    {name:<7} fraction {result.phase_fractions[name]: .6f}   "
            f"rho = {density:>10,.3f} mol/m^3"
        )
        print(f"            {amounts}")


_DIAGNOSTIC_KEYS = (
    "phase_detection",
    "stability_status",
    "tpd_min",
    "feed_branch",
    "stability_trials",
    "k_seed",
    "incipient_phase",
    "iterations",
    "converged",
    "termination_reason",
    "mass_balance_residual",
    "fugacity_residual",
    "delta_g_split_rt",
    "post_split_checked",
    "post_split_stable",
    "post_split_tpd_min",
    "phase_label_method",
)


def _print_diagnostics(result: ct.FlashResult) -> None:
    for key in _DIAGNOSTIC_KEYS:
        if key in result.diagnostics:
            value = result.diagnostics[key]
            rendered = f"{value: .6e}" if isinstance(value, float) else str(value)
            print(f"    {key:<24} {rendered}")


def main() -> None:
    mixture = _mixture()
    eos = PCSAFTEOS()

    print("=" * 78)
    print("PC-SAFT TP flash: methane / n-hexane, z = (0.30, 0.70), T = 300 K")
    print("=" * 78)

    print("\n1) Density roots")
    print("  a) at the feed composition, which is a dense fluid at 300 K:")
    for pressure in (1.0e6, 3.0e6, 8.0e6):
        _print_roots(eos, mixture, pressure)
    print("  b) pure n-hexane at its 300 K saturation pressure, where there are two:")
    hexane = ct.Mixture.from_database(["n-Hexane"], [1.0])
    roots = eos.density_roots(
        temperature_K=TEMPERATURE_K,
        pressure_Pa=HEXANE_PSAT_PA,
        composition=[1.0],
        mixture=hexane,
    )
    print(
        f"  P = {HEXANE_PSAT_PA / 1e3:>5.2f} kPa -> {len(roots)} root(s): "
        + ", ".join(f"{value:,.6f}" for value in roots)
        + " mol/m^3"
    )
    phi = [
        eos.fugacity_coefficients(
            mixture=hexane,
            temperature_K=TEMPERATURE_K,
            pressure_Pa=HEXANE_PSAT_PA,
            composition=[1.0],
            phase=phase,
        )[0]
        for phase in ("vapor", "liquid")
    ]
    print(
        f"            phi(vapor) = {phi[0]:.12f}   phi(liquid) = {phi[1]:.12f}   "
        f"|diff| = {abs(phi[0] - phi[1]):.2e}"
    )
    print(
        "            Equal fugacity on the two roots *is* the saturation condition,\n"
        "            and it is what the tangent-plane test uses to choose between them."
    )

    print("\n2) Two-phase flash at 3 MPa")
    two_phase = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=3.0e6, eos=eos)
    print(f"  phases: {len(two_phase.phases)}   vapor fraction {two_phase.vapor_fraction: .6f}")
    _print_phases(two_phase, eos, mixture)
    print("  diagnostics:")
    _print_diagnostics(two_phase)

    print("\n3) The same feed at 8 MPa: a compressed liquid, found stable")
    single = ct.flash_tp(mixture, temperature_K=TEMPERATURE_K, pressure_Pa=8.0e6, eos=eos)
    print(f"  phases: {len(single.phases)} -> {list(single.phases)}")
    _print_diagnostics(single)
    print(
        "  Only one density root exists here, so 'vapor' and 'liquid' would call the\n"
        "  same state, but the name is still measured: kappa = P / (rho dP/drho) at\n"
        "  that root is well below 1 (dense fluid), so phase_identity reports\n"
        "  'liquid' and phase_label_method records 'compressibility' - not a\n"
        "  vapor-first tie-break (ADR-0017)."
    )

    print("\n4) Peng-Robinson on the same state, for contrast (not a validation)")
    contrast = ct.flash_tp(
        mixture, temperature_K=TEMPERATURE_K, pressure_Pa=3.0e6, eos=ct.PengRobinsonEOS()
    )
    pcsaft_liquid = two_phase.phases["liquid"].composition.fractions
    pr_liquid = contrast.phases["liquid"].composition.fractions
    print(f"  PC-SAFT        x(methane) = {pcsaft_liquid[0]: .6f}")
    print(f"  Peng-Robinson  x(methane) = {pr_liquid[0]: .6f}")
    print(
        "  Two models, two answers. PC-SAFT's is cross-checked against teqp in\n"
        "  examples/validation/14_pcsaft_flash_vs_teqp.py; nothing here adjudicates."
    )

    print("\nLimits: non-associating only; phi-phi flash stops at two phases; any kij")
    print("you supply is your responsibility (the packaged set is pure components).")


if __name__ == "__main__":
    main()
