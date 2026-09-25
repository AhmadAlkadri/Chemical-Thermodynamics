"""Liquid-liquid phase stability with an activity model (NRTL) instead of an EOS.

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
Three feeds along the n-butanol / water composition axis:

  * z1 = 0.10 -- inside the miscibility gap: `stability_tp` reports UNSTABLE and
    the minimizing trial composition is the incipient butanol-rich phase;
  * z1 = 0.70 -- outside the gap: STABLE;
  * z1 = 0.02 -- the water-rich conjugate phase of the split, which is
    marginally stable: its tangent-plane minimum sits at zero, reached at the
    *other* equilibrium phase. Two coexisting liquids share one tangent plane.

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
        lines.append(f"  {name:<12} {value: .6f}")
    return "\n".join(lines)


def _report(model: ct.NRTL, z: Sequence[float], note: str) -> None:
    mixture = ct.Mixture.from_database(list(NAMES), list(z), normalize=True)
    result = ct.stability_tp(
        mixture,
        temperature_K=TEMPERATURE_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
    )

    print("=" * 68)
    print("Liquid-liquid stability by tangent-plane distance (NRTL activity model)")
    print(f"T [K]: {TEMPERATURE_K:.2f}")
    print(f"P [Pa]: {PRESSURE_PA:.5g}   (validated; not used by an activity model)")
    print(_format_vector("Feed z (mole):", NAMES, z))
    print(f"Context: {note}")
    print(f"Verdict: {result.status.upper()} (stable={result.stable})")
    print(f"tpd_min [-]: {result.tpd_min: .8e}")
    print(f"Model family: {result.diagnostics['model_family']}")
    print(f"Pressure dependent: {result.diagnostics['pressure_dependent']}")

    if result.trial_composition is None or result.k_values is None:
        print("Minimizing trial: none (every trial collapsed onto the feed)")
    else:
        print(f"Minimizing trial: {result.diagnostics['minimizing_trial']}")
        print(f"Converged in stage: {result.diagnostics['minimizing_trial_stage']}")
        print(_format_vector("Trial composition w (mole):", NAMES, result.trial_composition))
        print(_format_vector("Incipient K-values (w_i / z_i):", NAMES, result.k_values))
        print(f"sum(W) [-]: {float(result.diagnostics['sum_W']):.8f}")

    print("Trials:")
    header = (
        f"  {'label':<18}{'conv':>6}{'ssi':>5}{'newton':>7}{'tpd':>16}{'sum_W':>14}{'trivial':>9}"
    )
    print(header)
    for trial in result.trials:
        print(
            f"  {trial.label:<18}{str(trial.converged):>6}{trial.ssi_iterations:>5}"
            f"{trial.second_order_iterations:>7}{trial.tpd:>16.6e}"
            f"{trial.sum_W:>14.8f}{str(trial.trivial):>9}"
        )

    print(
        "Note: 'stable' means no negative tangent-plane distance was found from\n"
        "      this deterministic trial set (one pure-component-dominant estimate\n"
        "      per component). It is not a global proof."
    )


def main() -> None:
    model = _model()
    _report(model, (0.10, 0.90), "inside the miscibility gap")
    _report(model, (0.70, 0.30), "butanol-rich, outside the gap")
    _report(
        model,
        (0.019998419466946984, 0.980001580533053),
        "the water-rich conjugate phase of the split (marginally stable)",
    )


if __name__ == "__main__":
    main()
