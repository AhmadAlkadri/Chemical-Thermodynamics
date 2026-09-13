"""Liquid-liquid TP flash with an activity model only (no equation of state).

System
------
n-butanol(1) / water(2), a partially miscible binary. The NRTL parameters are
the 2-3 pair of Table 1 of

    S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable phase
    stability analysis for excess Gibbs energy models", Chemical Engineering
    Science 55 (2000) 1785-1796,

where Table 1 attributes them to C. M. McDonald and C. A. Floudas, AIChE
Journal 41 (1995) 1798-1814. Table 1 prints ``G_ij`` rather than ``alpha_ij``;
alpha = 0.48 for this pair is implied by ``G_ij = exp(-alpha_ij tau_ij)``. tau
is dimensionless as printed, so no temperature is needed (the paper states
none); 298.15 K is passed only because the API requires a temperature.

They are written inline here so the script runs from a bare install, and they
are deliberately NOT the packaged NRTL defaults, which are synthetic
placeholders.

What it shows
-------------
``flash_tp(mixture, temperature_K=..., pressure_Pa=..., activity_model=...)``
with **no** ``eos``: the mode is inferred as ``"gamma-gamma"`` and both phases
are liquids described by the same activity model.

  * z1 = 0.10 and z1 = 0.20 -- inside the miscibility gap. The phase count is an
    *output*: the feed is found unstable by the tangent-plane test and a split
    is solved. Both feeds return the **same tie-line** (the same pair of
    conjugate compositions) and differ only in how much of each phase there is;
    the script checks that against the lever rule.
  * z1 = 0.45 -- outside the gap: one liquid, because the feed was found stable.

Every two-phase answer is verified before it is returned: material balance,
equal activities, a negative Gibbs-energy change against the single-phase feed,
and a stability test of each converged phase (the "post-split" check, which is
what would report that a third phase is needed).

The phase names ``liquid1`` / ``liquid2`` are **roles, not identities**:
``liquid1`` is the phase the split started from as feed-like. Two feeds on one
tie-line can return the same two compositions with the labels swapped, which is
why this script compares the phase *set*.

No optional dependency is required (in particular, not `thermo`).
All inputs are SI: temperature in K, pressure in Pa. Pressure is validated but
does not affect an activity-model result.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import chemthermo as ct

NAMES = ("n-Butanol", "Water")

# Tessier et al. (2000) Table 1, pair 2-3 (n-butanol / water).
TAU_12 = 0.90047
TAU_21 = 3.51307
ALPHA = 0.48

TEMPERATURE_K = 298.15
PRESSURE_PA = 101325.0


def _model() -> ct.NRTL:
    parameters = ct.NRTLParameters.from_pairs([(NAMES[0], NAMES[1], TAU_12, TAU_21, ALPHA, ALPHA)])
    return ct.NRTL(parameters=parameters)


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
    )
    diagnostics = result.diagnostics

    print("=" * 72)
    print("Liquid-liquid TP flash (NRTL activity model, no EOS)")
    print(f"T [K]: {TEMPERATURE_K:.2f}")
    print(f"P [Pa]: {PRESSURE_PA:.5g}   (validated; not used by an activity model)")
    print(_format_vector("Feed z (mole):", NAMES, z))
    print(f"Context: {note}")
    print(f"flash_mode: {diagnostics['flash_mode']}   (inferred from the models supplied)")
    print(f"Phase detection: {diagnostics['phase_detection']}")
    print(
        f"Feed stability: {diagnostics['stability_status']}, tpd_min = {diagnostics['tpd_min']:+.8e}"
    )
    print(f"Phases returned: {result.phase_names()}")
    print(f"vapor_fraction: {result.vapor_fraction}   (None: neither phase is a vapor)")

    if len(result.phase_names()) == 1:
        print("Single liquid: the feed was found stable, so no split was attempted.")
        return result

    for name in result.phase_names():
        print(
            _format_vector(
                f"{name} (phase fraction {result.phase_fractions[name]:.8f}):",
                NAMES,
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
    print(f"  mass_balance_residual   = {diagnostics['mass_balance_residual']:.3e}")
    print(f"  equilibrium_residual    = {diagnostics['equilibrium_residual']:.3e}")
    print(f"  delta_g_split_rt        = {diagnostics['delta_g_split_rt']:+.6e}  (must be < 0)")
    print("Post-split stability of the converged phases:")
    for name in result.phase_names():
        print(
            f"  {name:<8} {diagnostics[f'phase_stability_{name}']:<10}"
            f" tpd_min = {diagnostics[f'phase_stability_tpd_min_{name}']:+.3e}"
        )
    print(
        f"  post_split_status = {diagnostics['post_split_status']}"
        f" (a failure here means a third phase is needed and raises)"
    )
    return result


def main() -> None:
    model = _model()

    first = _report(model, (0.10, 0.90), "inside the miscibility gap")
    second = _report(model, (0.20, 0.80), "also inside the gap, on the same tie-line")
    _report(model, (0.45, 0.55), "butanol-rich, outside the gap")

    print("=" * 72)
    print("Same tie-line, different amounts (lever rule)")
    first_set = sorted(
        tuple(first.phases[name].composition.fractions) for name in first.phase_names()
    )
    second_set = sorted(
        tuple(second.phases[name].composition.fractions) for name in second.phase_names()
    )
    worst = max(
        abs(a - b) for pair_a, pair_b in zip(first_set, second_set) for a, b in zip(pair_a, pair_b)
    )
    print(f"  max |composition difference| between the two tie-lines: {worst:.3e}")

    # Lever rule, checked against the phase set rather than the labels: pick the
    # butanol-rich phase by composition, not by which label it happened to get.
    low, high = first_set[0][0], first_set[1][0]
    for label, feed, result in (("z1=0.10", 0.10, first), ("z1=0.20", 0.20, second)):
        rich_name = max(
            result.phase_names(), key=lambda name: result.phases[name].composition.fractions[0]
        )
        rich = result.phase_fractions[rich_name]
        expected = (feed - low) / (high - low)
        print(
            f"  {label}: butanol-rich phase is '{rich_name}', fraction = {rich:.8f}, "
            f"lever rule = {expected:.8f}, |difference| = {abs(rich - expected):.3e}"
        )

    print(
        "\nNote: 'stable' means no negative tangent-plane distance was found from\n"
        "      the deterministic trial set (one pure-component-dominant estimate\n"
        "      per component). It is not a global proof. This gamma-gamma path\n"
        "      returns at most two phases; a state needing a third is reported,\n"
        "      not solved (ADR-0011 wires phase addition to modified-raoult only,\n"
        "      because no state in this repository exercises a third liquid here).\n"
        "      For a three-phase answer see examples/basic/flash_tp_vlle_demo.py."
    )


if __name__ == "__main__":
    main()
