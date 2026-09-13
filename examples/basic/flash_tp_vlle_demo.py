"""Three-phase (vapor-liquid-liquid) TP flash, discovered rather than assumed.

System
------
1-propanol(1) / n-butanol(2) / water(3) at 364.0 K and 101325 Pa, with
``flash_mode="modified-raoult"``: an NRTL liquid carrying the pure-liquid
reference fugacity ``f_i^0 = Psat_i(T)`` from the databank Antoine coefficients,
against an ideal-gas vapor (ADR-0010).

The NRTL parameters are Table 1 of

    S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable phase
    stability analysis for excess Gibbs energy models", Chemical Engineering
    Science 55 (2000) 1785-1796,

where Table 1 attributes them to C. M. McDonald and C. A. Floudas, AIChE
Journal 41 (1995) 1798-1814. Table 1 prints ``G_ij``, not ``alpha_ij``; the
alphas below (0.3, 0.3, 0.48) are implied by ``G_ij = exp(-alpha_ij tau_ij)``.
They are written inline so this script runs from a bare install, and they are
deliberately **not** the packaged NRTL defaults, which are synthetic
placeholders.

**These parameters were fitted to liquid-liquid data and are temperature
independent.** Nothing printed below is a comparison against experiment; it is
what this model says, checked against the equations this model is made of.

What it shows
-------------
Nobody tells `flash_tp` there are three phases. The search is

    stability of the feed  ->  a two-phase split  ->  stability of each phase
                           ->  add the phase that was found  ->  re-solve

and it stops when every phase is stable and every phase fraction is positive.
``diagnostics["phase_set_history"]`` records the sets it passed through, and
``FlashSettings.max_phases`` (default 3) bounds it. A fourth feed, outside the
tie-triangle, comes back as two phases from the same call - the phase count is
an output, not a setting.

The two liquids are named ``liquid1`` / ``liquid2``. Those are **roles, not
identities**: which one gets which number follows the order in which the search
created them, so compare the phase *set*. ``"vapor"`` is the one name with a
model-level meaning - it is the phase the ideal-gas candidate describes.

No optional dependency is required (in particular, not `thermo`).
All inputs are SI: temperature in K, pressure in Pa.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import chemthermo as ct

NAMES = ("1-Propanol", "n-Butanol", "Water")

# Tessier et al. (2000) Table 1, Problem 1. tau[i][j] is tau_ij.
TAU = (
    (0.0, -0.61259, -0.07149),
    (0.7164, 0.0, 0.90047),
    (2.7425, 3.51307, 0.0),
)
ALPHA = (
    (0.0, 0.3, 0.3),
    (0.3, 0.0, 0.48),
    (0.3, 0.48, 0.0),
)

TEMPERATURE_K = 364.0
PRESSURE_PA = 101325.0

#: Inside the tie-triangle: its centroid, so all three phases are present in
#: comparable amounts.
FEED_INSIDE = (0.13418838, 0.08427618, 0.78153544)
#: Another interior feed, weighted towards the vapor vertex.
FEED_VAPOR_RICH = (0.15771141, 0.09296719, 0.74932140)
#: Outside the triangle, past the liquid-liquid edge: two phases, no vapor.
FEED_OUTSIDE = (0.08084578, 0.07510435, 0.84404987)


def _model() -> ct.NRTL:
    pairs = [
        (NAMES[i], NAMES[j], TAU[i][j], TAU[j][i], ALPHA[i][j], ALPHA[j][i])
        for i in range(len(NAMES))
        for j in range(i + 1, len(NAMES))
    ]
    return ct.NRTL(parameters=ct.NRTLParameters.from_pairs(pairs))


def _format_vector(label: str, names: Sequence[str], values: Iterable[float]) -> str:
    lines = [label]
    for name, value in zip(names, values):
        lines.append(f"  {name:<12} {value: .8f}")
    return "\n".join(lines)


def _report(model: ct.NRTL, z: Sequence[float], note: str) -> ct.FlashResult:
    mixture = ct.Mixture.from_database(list(NAMES), list(z), normalize=True)
    result = ct.flash_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        flash_mode="modified-raoult",
    )
    diagnostics = result.diagnostics

    print("=" * 74)
    print("TP flash, modified Raoult (NRTL liquid + Antoine reference, ideal vapor)")
    print(f"T [K]: {TEMPERATURE_K:.2f}    P [Pa]: {PRESSURE_PA:.6g}")
    print(_format_vector("Feed z (mole):", NAMES, z))
    print(f"Context: {note}")
    print(
        f"Feed stability: {diagnostics['stability_status']}, "
        f"tpd_min = {float(diagnostics['tpd_min']):+.8e}, "
        f"feed candidate = {diagnostics.get('feed_branch')}"
    )
    print(
        f"Phases returned: {result.phase_names()}   "
        f"(phase_count = {diagnostics['phase_count']}, "
        f"regime = {diagnostics['phase_regime']})"
    )
    print(f"vapor_fraction: {result.vapor_fraction}")

    if len(result.phase_names()) == 1:
        print("Single phase: the feed was found stable, so no split was attempted.")
        return result

    for name in result.phase_names():
        print(
            _format_vector(
                f"{name} (phase fraction {result.phase_fractions[name]:.8f}):",
                NAMES,
                result.phases[name].composition.fractions,
            )
        )

    history = diagnostics.get("phase_set_history")
    if history is None:
        print(
            "Phase search: not entered - the two-phase answer passed its post-split\n"
            "              stability test, so there was nothing to add."
        )
    else:
        print(f"Phase search: {history}")
        print(
            f"              L = liquid, V = vapor; "
            f"{diagnostics['phases_added']} phase(s) added, "
            f"{diagnostics['phases_removed']} removed."
        )
        print(
            "Stages: "
            f"{diagnostics['ssi_iterations']} successive substitution + "
            f"{diagnostics['second_order_iterations']} second order + "
            f"{diagnostics['rachford_rice_iterations']} Rachford-Rice Newton "
            f"(converged in: {diagnostics['converged_stage']})"
        )

    print("Verification:")
    print(f"  mass_balance_residual   = {float(diagnostics['mass_balance_residual']):.3e}")
    print(f"  equilibrium_residual    = {float(diagnostics['equilibrium_residual']):.3e}")
    print(
        f"  delta_g_split_rt        = {float(diagnostics['delta_g_split_rt']):+.6e}"
        "  (must be < 0: better than one phase)"
    )
    if "delta_g_vs_two_phase_rt" in diagnostics:
        print(
            f"  delta_g_vs_two_phase_rt = {float(diagnostics['delta_g_vs_two_phase_rt']):+.6e}"
            "  (must be < 0: better than the two-phase candidate)"
        )
    print("Post-split stability of the converged phases:")
    for name in result.phase_names():
        print(
            f"  {name:<8} {diagnostics[f'phase_stability_{name}']:<10}"
            f" tpd_min = {float(diagnostics[f'phase_stability_tpd_min_{name}']):+.3e}"
        )
    return result


def main() -> None:
    model = _model()

    first = _report(model, FEED_INSIDE, "inside the tie-triangle (its centroid)")
    second = _report(model, FEED_VAPOR_RICH, "inside the tie-triangle, nearer the vapor vertex")
    _report(model, FEED_OUTSIDE, "outside the triangle, past the liquid-liquid edge")

    print("=" * 74)
    print("One tie-triangle, two feeds: the same three phases, different amounts")
    first_set = sorted(
        tuple(first.phases[name].composition.fractions) for name in first.phase_names()
    )
    second_set = sorted(
        tuple(second.phases[name].composition.fractions) for name in second.phase_names()
    )
    worst = max(
        abs(a - b) for pair_a, pair_b in zip(first_set, second_set) for a, b in zip(pair_a, pair_b)
    )
    print(f"  max |composition difference| between the two answers: {worst:.3e}")
    print("  phase fractions:")
    for label, result in (("centroid feed", first), ("vapor-rich feed", second)):
        fractions = ", ".join(
            f"{name}={result.phase_fractions[name]:.6f}" for name in result.phase_names()
        )
        print(f"    {label:<16} {fractions}")

    print(
        "\nNotes\n"
        "  * The phase count came from the tangent-plane stability test, not from\n"
        "    an assumption. 'stable' means no negative tangent-plane distance was\n"
        "    found from the deterministic trial set; it is not a global proof, and\n"
        "    a thin three-phase region can hide from it (validation Case V-2).\n"
        "  * liquid1 / liquid2 are roles, not identities. Compare the phase set.\n"
        "  * FlashSettings(max_phases=2) turns the search off and restores the\n"
        "    pre-ADR-0011 behavior, which raises when a third phase is needed."
    )


if __name__ == "__main__":
    main()
