"""Fixed density-root surfaces in the EOS stability trials, with PASS/FAIL on every check.

What is being validated
-----------------------
ADR-0021 pins every tangent-plane trial of an equation of state to one density
root: the Wilson vapour-like start iterates on the vapour root, the Wilson
liquid-like and the pure-component-dominant starts on the liquid root. Before
it, each trial re-selected the lowest-Gibbs root at *every* iterate, which is
what ADR-0012 had already found wrong for the modified-Raoult candidate pair
and had deliberately left in place here.

The consequence measured in this script is validation Case P-9 (iv). At 335 K,
above the water / n-hexane three-phase temperature ``T3 = 334.8078 K``, a
water-rich feed (``z_water = 0.7``) used to come back as **two liquids** that
are metastable by 2.495e-03 RT: from the hexane-rich liquid every trial slid
onto the liquid root and stopped at the trivial solution or at its partner
liquid, while a vapour stationary point with ``tpd = -6.52e-03`` sat unvisited.
With the trials pinned, the vapour-like trial stays on the vapour root, finds
that stationary point, and ``flash_tp`` returns the vapour-liquid pair - the
lower-Gibbs answer.

Three independent routes are used, so no single one carries the result:

1. **A 4-equation Newton written in this script** locates ``T3`` from the public
   ``fugacity_coefficients`` interface alone - no ``flash_tp``, no stability
   test.
2. **A 2-equation Newton on the liquid branch** gives the conjugate-liquid pair
   above ``T3`` independently of the flash, so the **reduced Gibbs energies** of
   the two candidate answers can be compared directly: which pair is the
   equilibrium is decided by arithmetic, not by which one the solver happened to
   return.
3. **FeOs** (feos-org/feos, MIT OR Apache-2.0; the same Gross & Sadowski model
   in Rust, every derivative by automatic differentiation) supplies chemical
   potentials **at chemthermo's converged phases and densities** - the reference
   model's own equilibrium condition, at chemthermo's answer. As in Cases P-6 to
   P-9, one input is deliberately not shared (the 42 universal constants of the
   2001 dispersion term: chemthermo packages the ten figures as printed, FeOs
   hard-codes fourteen), so that route is reported twice, as shipped and with
   FeOs's constants substituted - the substituted half under ``--full``, since
   it has to re-run the flash on FeOs's table for the comparison to mean
   anything.

``--full`` adds ``T3 + 0.05 K`` and ``T3 + 0.5 K`` (each a slower
phase-addition search than 335 K, because they sit closer to ``T3``), the
matched-universal-constants chemical potentials asserted at 1e-08, the 41-point
temperature scans at both feeds across ``[T3 - 1 K, T3 + 1 K]``, the bisected
verdict boundary to 1e-06 K, and the per-surface trial statistics over the
144-state Peng-Robinson grid. It takes several minutes; the default takes
about five seconds.

Requires the optional ``feos`` dependency for route 3::

    pip install -e ".[validation]"

Routes 1 and 2 run without it; the script says so and still exits 0.
"""

from __future__ import annotations

import argparse
import json
import math
from typing import Sequence

import numpy as np

try:
    import si_units as si
    from feos import Contributions, EquationOfState, Parameters, PureRecord, State
except ImportError:  # pragma: no cover - exercised only without the extra
    si = None  # type: ignore[assignment]

import chemthermo as ct
from chemthermo.eos import PCSAFTEOS
from chemthermo.eos import pcsaft as pcsaft_module
from chemthermo.eos.pcsaft import R_J_PER_MOL_K

BINARY = ("Water", "n-Hexane")
ATMOSPHERE_PA = 101325.0
#: The two offsets above `T3` that the slice report recomputes; `--full` only,
#: because each flash there is a phase-addition search and the closer ones are
#: the slower.
OFFSETS_K = (0.05, 0.5)
#: The temperature at which validation Case P-9 (iv) recorded the miss.
LEDGER_T_K = 335.0
Z_WATER = 0.7

#: Gross & Sadowski (2002) Table 1 for water, (2001) Table 1 for n-hexane,
#: written out here so the reference model is built from this file.
PARAMETERS: dict[str, tuple[float, float, float, float, float | None, float]] = {
    "Water": (18.015, 1.0656, 3.0007, 366.51, 0.034868, 2500.7),
    "n-Hexane": (86.177, 3.0576, 3.7983, 236.77, None, 0.0),
}

#: The Peng-Robinson grid of `tests/test_flash_phase_detection.py`, reused for
#: the trial statistics under `--full`.
PR_MIXTURES: tuple[tuple[tuple[str, ...], tuple[float, ...]], ...] = (
    (("Methane", "Ethane"), (0.5, 0.5)),
    (("Methane", "Propane"), (0.7, 0.3)),
    (("Ethane", "n-Heptane"), (0.7, 0.3)),
    (("Methane", "n-Pentane"), (0.6, 0.4)),
    (("Methane", "Ethane", "Propane"), (0.5, 0.3, 0.2)),
    (("Propane", "n-Butane", "n-Pentane"), (0.4, 0.3, 0.3)),
)
PR_T_K = (170.0, 200.0, 240.0, 280.0, 320.0, 360.0)
PR_P_PA = (2.0e5, 1.0e6, 3.0e6, 8.0e6)

POTENTIAL_TOL = 1e-8

failures: list[str] = []


def record(label: str, ok: bool) -> None:
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        failures.append(label)


# ---------------------------------------------------------------------------
# The chemthermo side, through the public interface only
# ---------------------------------------------------------------------------


def mixture(z: Sequence[float]) -> ct.Mixture:
    return ct.Mixture.from_database(list(BINARY), list(z), normalize=True)


#: One model and one mixture, reused. `PCSAFTEOS()` reads the packaged
#: parameter records on construction and `Mixture.from_database` the component
#: records, so rebuilding either inside the Newton loops below would cost more
#: than the thermodynamics they are calling.
EOS = PCSAFTEOS()
BINARY_MIXTURE = mixture((0.5, 0.5))


def ln_f(temperature_K: float, composition: Sequence[float], branch: str) -> np.ndarray:
    """``ln(x_i phi_i)`` on one named density branch."""
    x = np.asarray(composition, dtype=float)
    x = x / float(np.sum(x))
    phi = np.asarray(
        EOS.fugacity_coefficients(
            mixture=BINARY_MIXTURE,
            temperature_K=temperature_K,
            pressure_Pa=ATMOSPHERE_PA,
            composition=x.tolist(),
            phase=branch,
        ),
        dtype=float,
    )
    return np.log(x) + np.log(phi)


def reduced_g(temperature_K: float, composition: Sequence[float], branch: str) -> float:
    x = np.asarray(composition, dtype=float)
    x = x / float(np.sum(x))
    return float(np.sum(x * ln_f(temperature_K, x, branch)))


def damped_newton(residual, start: np.ndarray, steps: Sequence[float], admissible, tol: float):
    u = np.array(start, dtype=float)
    for _iteration in range(60):
        f = residual(u)
        worst = float(np.max(np.abs(f)))
        if worst < tol:
            break
        jacobian = np.zeros((f.size, u.size))
        for column in range(u.size):
            shifted = u.copy()
            shifted[column] += steps[column]
            jacobian[:, column] = (residual(shifted) - f) / steps[column]
        direction = np.linalg.solve(jacobian, -f)
        scale = 1.0
        while scale > 1e-10:
            candidate = u + scale * direction
            if admissible(candidate) and float(np.max(np.abs(residual(candidate)))) < worst:
                break
            scale *= 0.5
        u = u + scale * direction
    return u, float(np.max(np.abs(residual(u))))


def three_phase_temperature() -> tuple[float, tuple[float, float]]:
    print("\n1) The three-phase temperature, from a 4-equation Newton (no flash_tp)")
    print("-" * 78)

    def residual(u: np.ndarray) -> np.ndarray:
        first = ln_f(u[3], [u[0], 1.0 - u[0]], "liquid")
        second = ln_f(u[3], [u[1], 1.0 - u[1]], "liquid")
        vapor = ln_f(u[3], [u[2], 1.0 - u[2]], "vapor")
        return np.concatenate([first - second, first - vapor])

    u, worst = damped_newton(
        residual,
        np.array([0.9999, 0.02, 0.20, 334.5]),
        (1e-8, 1e-8, 1e-8, 3.345e-5),
        lambda v: bool(np.all(v[:3] > 0.0) and np.all(v[:3] < 1.0) and 200.0 < v[3] < 600.0),
        1e-11,
    )
    print(f"    residual    {worst:.3e}")
    print(f"    T3          {u[3]:.9f} K")
    print(f"    x_water(I)  {u[0]:.12f}")
    print(f"    x_water(II) {u[1]:.12f}")
    print(f"    y_water     {u[2]:.12f}")
    record("the 4-equation Newton converged", worst < 1e-10)
    return float(u[3]), (float(u[0]), float(u[1]))


def conjugate_liquids(temperature_K: float, start: Sequence[float]) -> tuple[np.ndarray, float]:
    """The two-liquid tie line on the liquid branch, from a 2-equation Newton.

    ``start`` is the ``T3`` tie line's two liquid compositions: the binodal
    moves slowly with temperature, so a start taken from ``T3`` converges in a
    few steps where a fixed guess needs tens of them.
    """

    def residual(u: np.ndarray) -> np.ndarray:
        return ln_f(temperature_K, [u[0], 1.0 - u[0]], "liquid") - ln_f(
            temperature_K, [u[1], 1.0 - u[1]], "liquid"
        )

    u, worst = damped_newton(
        residual,
        np.array(start, dtype=float),
        (1e-8, 1e-8),
        lambda v: bool(np.all(v > 0.0) and np.all(v < 1.0)),
        1e-12,
    )
    return u, worst


def the_flipped_verdict(
    liquids_at_t3: Sequence[float], temperatures: Sequence[tuple[str, float]]
) -> dict[str, ct.FlashResult]:
    print("\n2) Case P-9 (iv): z_water = 0.7 above T3, and which answer is lower in Gibbs")
    print("-" * 78)
    print(
        "    Before ADR-0021 the temperatures below returned two liquids. The pair is a\n"
        "    genuine stationary state of the model - it is just not the equilibrium one."
    )
    results: dict[str, ct.FlashResult] = {}
    for label, temperature in temperatures:
        result = ct.flash_tp(
            mixture((Z_WATER, 1.0 - Z_WATER)),
            temperature_K=temperature,
            pressure_Pa=ATMOSPHERE_PA,
            eos=EOS,
        )
        results[label] = result
        print(f"\n    T = {label} = {temperature:.6f} K")
        print(f"      phases           {sorted(result.phases)}")
        print(f"      regime           {result.diagnostics['phase_regime']}")
        for name in sorted(result.phases):
            fractions = result.phases[name].composition.fractions
            print(
                f"      {name:<8} x_water = {fractions[0]:.12f}"
                f"   beta = {result.phase_fractions[name]:.12f}"
            )
        record(
            f"{label} returns a vapour-liquid pair",
            sorted(result.phases) == ["liquid", "vapor"],
        )

        # The two-liquid answer, found independently of the flash.
        u, worst = conjugate_liquids(temperature, liquids_at_t3)
        record(f"the conjugate-liquid Newton converged at {label}", worst < 1e-11)
        first = np.array([u[0], 1.0 - u[0]])
        second = np.array([u[1], 1.0 - u[1]])
        # Lever rule on water.
        beta = (Z_WATER - second[0]) / (first[0] - second[0])
        g_liquids = beta * reduced_g(temperature, first, "liquid") + (1.0 - beta) * reduced_g(
            temperature, second, "liquid"
        )
        g_returned = sum(
            result.phase_fractions[name]
            * reduced_g(
                temperature,
                result.phases[name].composition.fractions,
                "vapor" if name == "vapor" else "liquid",
            )
            for name in result.phases
        )
        print(f"      G(returned)/RT   {g_returned:.12f}")
        print(f"      G(two liquids)/RT{g_liquids:.12f}")
        print(f"      gap              {g_liquids - g_returned:.6e} RT")
        record(
            f"the returned pair is lower in Gibbs than the two liquids at {label}",
            g_returned < g_liquids,
        )
    return results


def the_trial_table(liquids_at_t3: Sequence[float]) -> None:
    print("\n3) The stability trials from the hexane-rich liquid, surface by surface")
    print("-" * 78)
    print(
        "    The state validation Case P-9 (iv) pinned: 335 K, and the hexane-rich vertex\n"
        "    of the two-liquid pair that used to be returned there. Every trial reached\n"
        "    the trivial solution or the partner liquid; the vapour stationary point below\n"
        "    was never visited, and the two-liquid answer that followed is metastable by\n"
        "    2.495e-03 RT."
    )
    temperature = LEDGER_T_K
    u, worst = conjugate_liquids(temperature, liquids_at_t3)
    hexane_rich = [float(u[1]), 1.0 - float(u[1])]
    print(f"\n    T = {temperature:.6f} K, w = ({hexane_rich[0]:.12f}, {hexane_rich[1]:.12f})")
    print(f"    (the hexane-rich vertex of the metastable pair; Newton residual {worst:.2e})")
    stability = ct.stability_tp(
        mixture(hexane_rich),
        temperature_K=temperature,
        pressure_Pa=ATMOSPHERE_PA,
        eos=EOS,
    )
    print(f"\n    verdict {stability.status}   tpd_min {stability.tpd_min:.12e}")
    print(f"    trial surfaces   {stability.diagnostics.get('trial_surfaces')}")
    print(
        f"    fallbacks        {stability.diagnostics.get('surface_fallback_trial_count')}"
        f" trial(s), {stability.diagnostics.get('surface_fallback_evaluation_count')} evaluation(s)"
    )
    print(f"\n    {'trial':<16}{'surface':<9}{'branch':<9}{'iters':>6}  {'tpd':>16}  trivial")
    for trial in stability.trials:
        print(
            f"    {trial.label:<16}{str(trial.surface):<9}{str(trial.phase_branch):<9}"
            f"{trial.iterations:>6}  {trial.tpd:>16.9e}  {trial.trivial}"
        )
    record("the hexane-rich liquid is unstable", stability.status == "unstable")
    record(
        "the incipient phase it names is a vapour",
        stability.phase_branch == "vapor",
    )
    record(
        "the minimizing trial ran on the vapour root",
        stability.diagnostics.get("minimizing_trial_surface") == "vapor",
    )


# ---------------------------------------------------------------------------
# FeOs
# ---------------------------------------------------------------------------


def feos_constants() -> tuple[np.ndarray, np.ndarray]:
    """FeOs's fourteen-figure universal constants, from the Case P-6 script."""
    import importlib.util
    from pathlib import Path

    path = Path(__file__).with_name("16_pcsaft_association_vs_feos.py")
    spec = importlib.util.spec_from_file_location("_case_p6", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.FEOS_A, module.FEOS_B


def _pure_record(name: str):
    mw, m, sigma, epsilon, kappa, epsilon_ab = PARAMETERS[name]
    payload: dict[str, object] = {
        "identifier": {"name": name},
        "molarweight": mw,
        "m": m,
        "sigma": sigma,
        "epsilon_k": epsilon,
    }
    if kappa is not None:
        payload["association_sites"] = [
            {"kappa_ab": kappa, "epsilon_k_ab": epsilon_ab, "na": 1.0, "nb": 1.0}
        ]
    return PureRecord.from_json_str(json.dumps(payload))


def _feos_reduced_potentials(
    temperature_K: float, density: float, x: Sequence[float] | np.ndarray
) -> np.ndarray:
    """``mu_i / RT`` from FeOs at chemthermo's ``(T, rho, x)``, up to a constant."""
    eos = EquationOfState.pcsaft(Parameters.from_records([_pure_record(n) for n in BINARY]))
    values = np.asarray(x, dtype=float)
    state = State(
        eos,
        temperature=temperature_K * si.KELVIN,
        density=density * (si.MOL / si.METER**3),
        composition=values,
    )
    factor = R_J_PER_MOL_K * temperature_K
    z_factor = (state.pressure() / si.PASCAL) / (density * factor)
    mu_residual = np.atleast_1d(
        np.asarray(state.chemical_potential(Contributions.Residual) / si.JOULE * si.MOL)
    )
    return mu_residual / factor - math.log(z_factor) + np.log(values)


def _worst_potential_gap(result: ct.FlashResult) -> float:
    potentials = []
    for name in sorted(result.phases):
        fractions = np.asarray(result.phases[name].composition.fractions, dtype=float)
        roots = EOS.density_roots(
            mixture=BINARY_MIXTURE,
            temperature_K=result.temperature_K,
            pressure_Pa=ATMOSPHERE_PA,
            composition=fractions.tolist(),
        )
        density = roots[0] if name == "vapor" else roots[-1]
        potentials.append(_feos_reduced_potentials(result.temperature_K, density, fractions))
    return float(np.max(np.abs(potentials[0] - potentials[1])))


def the_chemical_potentials(results: dict[str, ct.FlashResult], full: bool) -> None:
    print("\n4) FeOs's chemical potentials at the new vapour-liquid phases")
    print("-" * 78)
    print(
        "    As shipped, the 42 universal constants of the 2001 dispersion term differ\n"
        "    between the two packages (ten figures as printed against FeOs's fourteen),\n"
        "    which floors the comparison near 1e-06. The asserted number is the one with\n"
        "    FeOs's own table substituted and the flash re-run on it."
    )
    for label, result in results.items():
        print(
            f"    {label}, as shipped:  max_i |mu_i^L - mu_i^V| / RT ="
            f" {_worst_potential_gap(result):.3e}"
        )

    if not full:
        print("    (--full re-runs the flash on FeOs's table and asserts the 1e-08 number)")
        return

    shipped_a, shipped_b = pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL
    matched_a, matched_b = feos_constants()
    try:
        pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = matched_a, matched_b
        for label, result in results.items():
            rerun = ct.flash_tp(
                mixture((Z_WATER, 1.0 - Z_WATER)),
                temperature_K=result.temperature_K,
                pressure_Pa=ATMOSPHERE_PA,
                eos=EOS,
            )
            if sorted(rerun.phases) != ["liquid", "vapor"]:
                record(f"the matched-constants re-run at {label} is vapour-liquid", False)
                continue
            worst = _worst_potential_gap(rerun)
            print(f"    {label}, matched:    max_i |mu_i^L - mu_i^V| / RT = {worst:.3e}")
            record(
                f"FeOs calls the {label} vapour-liquid pair an equilibrium",
                worst < POTENTIAL_TOL,
            )
    finally:
        pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = shipped_a, shipped_b


# ---------------------------------------------------------------------------
# --full
# ---------------------------------------------------------------------------


def the_scans(t3: float) -> None:
    print("\n5) The 41-point scans, and the bisected verdict boundary")
    print("-" * 78)
    for z_water in (0.3, 0.7):
        temperatures = np.linspace(t3 - 1.0, t3 + 1.0, 41)
        verdicts: list[str] = []
        errors = 0
        for temperature in temperatures:
            try:
                result = ct.flash_tp(
                    mixture((z_water, 1.0 - z_water)),
                    temperature_K=float(temperature),
                    pressure_Pa=ATMOSPHERE_PA,
                    eos=EOS,
                )
                verdicts.append(str(result.diagnostics["phase_regime"]))
            except ct.ConvergenceError:
                errors += 1
                verdicts.append("ERROR")
        switches = [i for i in range(1, 41) if verdicts[i] != verdicts[i - 1]]
        shape = "".join({"LLE": "L", "VLE": "V"}.get(v, "?") for v in verdicts)
        print(f"\n    z_water = {z_water}:  {shape}")
        print(
            f"      {len(set(verdicts))} distinct verdict(s), {len(switches)} switch(es),"
            f" {errors} ConvergenceError(s)"
        )
        record(f"z_water = {z_water}: no ConvergenceError over the window", errors == 0)
        record(f"z_water = {z_water}: exactly one verdict switch", len(switches) == 1)

        # A 2e-03 K bracket, not the whole window: the scan above has already
        # shown there is exactly one switch in it, and each bisection step here
        # costs a three-phase search.
        low, high = t3 - 1.0e-3, t3 + 1.0e-3
        below = verdicts[0]

        def regime(temperature: float) -> str:
            return str(
                ct.flash_tp(
                    mixture((z_water, 1.0 - z_water)),
                    temperature_K=float(temperature),
                    pressure_Pa=ATMOSPHERE_PA,
                    eos=EOS,
                ).diagnostics["phase_regime"]
            )

        while high - low > 1.0e-6:
            middle = 0.5 * (low + high)
            if regime(middle) == below:
                low = middle
            else:
                high = middle
        boundary = 0.5 * (low + high)
        print(f"      boundary {boundary:.9f} K, T3 - boundary = {t3 - boundary:.3e} K")
        record(f"z_water = {z_water}: the switch is at T3 to 1e-06 K", abs(boundary - t3) < 1.0e-6)


def the_trial_statistics() -> None:
    print("\n6) Which surface finds the minimizer, over the 144-state Peng-Robinson grid")
    print("-" * 78)
    eos = ct.PengRobinsonEOS()
    by_surface: dict[str, int] = {}
    by_label: dict[str, int] = {}
    iterations = 0
    trials = 0
    fallback_trials = 0
    verdicts: dict[str, int] = {}
    for names, z in PR_MIXTURES:
        for temperature in PR_T_K:
            for pressure in PR_P_PA:
                result = ct.stability_tp(
                    ct.Mixture.from_database(list(names), list(z), normalize=True),
                    temperature_K=temperature,
                    pressure_Pa=pressure,
                    eos=eos,
                )
                verdicts[result.status] = verdicts.get(result.status, 0) + 1
                for trial in result.trials:
                    trials += 1
                    iterations += trial.iterations
                    fallback_trials += int(trial.surface_fallback)
                surface = result.diagnostics.get("minimizing_trial_surface")
                if surface is not None:
                    by_surface[str(surface)] = by_surface.get(str(surface), 0) + 1
                    label = str(result.diagnostics["minimizing_trial"])
                    by_label[label] = by_label.get(label, 0) + 1
    print(f"    verdicts                     {verdicts}")
    print(f"    trials / total iterations    {trials} / {iterations}")
    print(
        f"    trials that fell back        {fallback_trials}"
        f"  (the model had one root where the trial stopped)"
    )
    print(f"    minimizing surface           {by_surface}")
    print(f"    minimizing trial             {by_label}")
    record("every state reached a verdict", verdicts.get("inconclusive", 0) == 0)
    record("both root surfaces find minimizers", len(by_surface) == 2)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="also run the 41-point scans, the bisected boundary and the grid statistics",
    )
    # `parse_known_args`, not `parse_args`: `tests/test_examples.py` runs this
    # script via `runpy.run_path` with pytest's own `sys.argv` still in place.
    args, _unknown = parser.parse_known_args()

    print("=" * 78)
    print("Fixed density-root surfaces in the EOS stability trials (Case P-11; ADR-0021)")
    print("=" * 78)
    t3, liquids_at_t3 = three_phase_temperature()
    # 335 K is the state validation Case P-9 (iv) pinned, and it is also the
    # cheapest of the three: the closer the feed sits to `T3`, the longer the
    # phase-addition search that gets to the vapour-liquid answer takes.
    temperatures: list[tuple[str, float]] = [("335 K", LEDGER_T_K)]
    if args.full:
        temperatures += [(f"T3 + {offset} K", t3 + offset) for offset in OFFSETS_K]
    results = the_flipped_verdict(liquids_at_t3, temperatures)
    the_trial_table(liquids_at_t3)
    if si is None:
        print("\nfeos is not installed; route 3 is skipped.")
        print("Install with: pip install -e '.[validation]'")
    else:
        the_chemical_potentials(results, args.full)
    if args.full:
        the_scans(t3)
        the_trial_statistics()
    else:
        print(
            "\n  (pass --full for T3 + 0.05 K and T3 + 0.5 K, the matched-constants"
            "\n   comparison, the 41-point scans, the bisected boundary and the"
            "\n   Peng-Robinson grid)"
        )

    print("\n" + "=" * 78)
    print(
        "k_ij = 0 throughout; nothing above is compared against measurement. The two\n"
        "liquids this slice replaced are a real stationary state of the model - the\n"
        "Gibbs comparison in section 2, not the solver's preference, is what makes the\n"
        "vapour-liquid pair the answer."
    )
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
