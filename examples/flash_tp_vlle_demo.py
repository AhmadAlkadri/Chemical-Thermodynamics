"""TP VLLE scaffold demo using NRTL (liquids) + Peng-Robinson (vapor).

This is an initial, deterministic VLLE scaffold that stages a gamma-phi VLE
solve and then splits the liquid into two phases via a fixed-point LLE update.
"""

from __future__ import annotations

from typing import Iterable, Mapping

import chemthermo as ct


def _format_composition(
    label: str, components: Iterable[ct.Component], fractions: Iterable[float]
) -> str:
    lines = [f"{label}"]
    for component, fraction in zip(components, fractions):
        lines.append(f"  {component.name:<12} {fraction: .6f}")
    return "\n".join(lines)


def _format_diagnostics(diag: Mapping[str, float | int | str | bool]) -> str:
    keys = [
        "phase_regime",
        "phase_state",
        "vle_iterations",
        "lle_iterations",
        "converged",
        "termination_reason",
        "lle_max_delta_x",
        "lle_composition_residual",
        "lle_tol",
    ]
    lines = ["Diagnostics:"]
    for key in keys:
        if key not in diag:
            continue
        value = diag[key]
        if isinstance(value, float):
            formatted = f"{value:.3e}"
        else:
            formatted = str(value)
        lines.append(f"  {key}: {formatted}")
    if len(lines) == 1:
        lines.append("  (none)")
    return "\n".join(lines)


def main() -> None:
    temperature_K = 240.0
    pressure_Pa = 3.0e6

    component_names = ("Methane", "Ethane")
    components = tuple(ct.Component.from_database(name) for name in component_names)
    z = (0.50, 0.50)

    mixture = ct.Mixture(
        components=components,
        composition=ct.Composition(fractions=z, basis="mole", normalize=False),
    )

    eos = ct.PengRobinsonEOS()
    activity_model = ct.NRTL(parameters=ct.ActivityParameters.load("NRTL"))

    result = ct.flash_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=eos,
        activity_model=activity_model,
        flash_mode="vlle",
    )

    print("TP flash (VLLE scaffold)")
    print(f"T [K]: {temperature_K:.2f}")
    print(f"P [Pa]: {pressure_Pa:.3e}")
    print(_format_composition("Feed z (mole):", components, mixture.fractions))
    print("Phase fractions:")
    for name, fraction in result.phase_fractions.items():
        print(f"  {name}: {fraction:.6f}")

    liquid1 = result.phases.get("liquid1")
    liquid2 = result.phases.get("liquid2")
    vapor = result.phases.get("vapor")

    if liquid1 is not None:
        print(_format_composition("x (liquid1):", components, liquid1.composition.fractions))
    if liquid2 is not None:
        print(_format_composition("x (liquid2):", components, liquid2.composition.fractions))
    if vapor is not None:
        print(_format_composition("y (vapor):", components, vapor.composition.fractions))

    print(_format_diagnostics(result.diagnostics))


if __name__ == "__main__":
    main()
