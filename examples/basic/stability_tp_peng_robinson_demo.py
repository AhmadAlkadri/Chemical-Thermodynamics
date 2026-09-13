"""Phase stability demo: Michelsen tangent-plane analysis with Peng-Robinson.

Runs the canonical README mixture at two states:
  - 240 K, 3 MPa, where `flash_tp` splits the feed into two phases;
  - 450 K, 1 bar, a clearly single-phase state.

All inputs are SI: temperature in K, pressure in Pa. Compositions are mole
fractions and the tangent-plane distance is dimensionless (units of RT).
"""

from __future__ import annotations

from typing import Iterable, Sequence

import chemthermo as ct


def _format_vector(label: str, names: Sequence[str], values: Iterable[float]) -> str:
    lines = [label]
    for name, value in zip(names, values):
        lines.append(f"  {name:<12} {value: .6f}")
    return "\n".join(lines)


def _report(
    names: Sequence[str], z: Sequence[float], temperature_K: float, pressure_Pa: float
) -> None:
    components = tuple(ct.Component.from_database(name) for name in names)
    composition = ct.Composition(fractions=tuple(z), basis="mole", normalize=False)
    mixture = ct.Mixture(components=components, composition=composition)

    result = ct.stability_tp(
        mixture,
        temperature_K=temperature_K,
        pressure_Pa=pressure_Pa,
        eos=ct.PengRobinsonEOS(),
    )

    print("=" * 62)
    print("Phase stability by tangent-plane distance (Peng-Robinson EOS)")
    print(f"T [K]: {temperature_K:.2f}")
    print(f"P [Pa]: {pressure_Pa:.3e}")
    print(_format_vector("Feed z (mole):", names, z))
    print(f"Verdict: {result.status.upper()} (stable={result.stable})")
    print(f"tpd_min [-]: {result.tpd_min: .8e}")
    print(f"Feed min-Gibbs root branch: {result.feed_branch}")

    if result.trial_composition is None or result.k_values is None:
        print("Minimizing trial: none (every trial collapsed onto the feed)")
    else:
        print(f"Minimizing trial: {result.diagnostics['minimizing_trial']}")
        print(f"Incipient-phase min-Gibbs root branch: {result.phase_branch}")
        print(_format_vector("Trial composition w (mole):", names, result.trial_composition))
        print(_format_vector("Incipient K-values (w_i / z_i):", names, result.k_values))
        print(f"sum(W) [-]: {float(result.diagnostics['sum_W']):.8f}")

    print("Trials:")
    header = f"  {'label':<16}{'converged':>10}{'iters':>7}{'tpd':>16}{'sum_W':>14}{'trivial':>9}"
    print(header)
    for trial in result.trials:
        print(
            f"  {trial.label:<16}{str(trial.converged):>10}{trial.iterations:>7}"
            f"{trial.tpd:>16.6e}{trial.sum_W:>14.8f}{str(trial.trivial):>9}"
        )

    print(
        "Note: 'stable' means no negative tangent-plane distance was found from\n"
        "      this deterministic trial set. It is not a global proof."
    )


def main() -> None:
    names = ("Methane", "Ethane", "Propane")
    z = (0.50, 0.30, 0.20)

    _report(names, z, temperature_K=240.0, pressure_Pa=3.0e6)
    _report(names, z, temperature_K=450.0, pressure_Pa=1.0e5)


if __name__ == "__main__":
    main()
