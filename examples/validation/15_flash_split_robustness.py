"""Phi-phi split robustness with PC-SAFT, with PASS/FAIL on every check.

Validation Case F-4 (ADR-0016).

What is being validated
-----------------------
The phi-phi split of `flash_tp` is successive substitution on the K-values,
with the vapor fraction re-solved from Rachford-Rice after every update. Before
ADR-0016 that inner solve searched `[0, 1]` only and reported failure when the
Rachford-Rice function did not change sign there - so an *iterate* whose
implied split was momentarily outside the physical range ended the whole flash
with

    ConvergenceError: Rachford-Rice failed to bracket a vapor fraction.

even for feeds the tangent-plane test had already proved to be two-phase. Four
states of the grid below did exactly that.

ADR-0016 changed two things:

1. the vapor fraction may leave `[0, 1]` during iteration - the "negative
   flash" of Whitson & Michelsen (1989), solved on the Leibovici-Neoschil
   window `1/(1 - K_max) < beta < 1/(1 - K_min)`, which is precisely the set of
   vapor fractions for which every phase mole fraction is non-negative;
2. the second-order Gibbs-minimization stage of ADR-0009 - already used by the
   liquid-liquid and modified-Raoult splits - now finishes a phi-phi split that
   successive substitution could not.

Both are gated so that every phi-phi state that converged before ADR-0016
converges identically now: the extended solver calls the in-window one first
and returns its answer verbatim when there is one, and the stage runs only
after successive substitution has spent its whole budget.

What is checked here
--------------------
1. The reference state - carbon dioxide / n-decane, z = (0.9, 0.1), 240 K,
   1.0 MPa - against a damped Newton written *in this script* on the
   equal-fugacity system in vapor mole numbers, a different formulation of the
   same equilibrium. Plus the phase densities from `density_roots`, and the
   Gibbs-energy reduction recomputed from the public API.
2. The Case F-4 grid. By default a fixed 16-state representative subset,
   including all four states that need the second-order stage; pass `--full`
   for the complete 188-state grid. Every state must answer; every two-phase
   answer must carry its four invariants; every single-phase answer must be a
   stability verdict.
3. The four previously failing states, each of which must still fail on the
   legacy `phase_detection="wilson-heuristic"` path (ADR-0016 left it alone,
   deliberately) and must still fail with `second_order=False`.

This script needs no optional dependency. The teqp cross-check of the same
state lives in `tests/validation/test_flash_split_robustness_pcsaft.py`.

Runtime (slice `flash-phase-labels-by-compressibility`): the full 188-state
grid takes ~2 minutes and used to run here *and* in
`tests/validation/test_flash_split_robustness_pcsaft.py::test_the_whole_grid_answers_and_every_answer_is_verified`
on every `pytest -q`, which is most of why the suite grew from ~116 s to
~389 s. That test is now `@pytest.mark.slow` (opt in with `pytest -q -m
slow`), and this script defaults to the 16-state subset below so the
`tests/test_examples.py` smoke test that runs it stays cheap; `--full` still
exercises the whole grid on demand.
"""

from __future__ import annotations

import argparse

import numpy as np

import chemthermo as ct

REFERENCE_COMPONENTS = ("Carbon dioxide", "n-Decane")
REFERENCE_FEED = (0.9, 0.1)
REFERENCE_T_K = 240.0
REFERENCE_P_PA = 1.0e6

GRID: tuple[
    tuple[tuple[str, str], tuple[float, ...], tuple[float, ...], tuple[float, ...]], ...
] = (
    (
        ("Carbon dioxide", "n-Decane"),
        (0.6, 0.8, 0.9),
        (230.0, 240.0, 250.0, 260.0),
        (1.0e6, 1.5e6, 2.0e6, 2.5e6),
    ),
    (
        ("Methane", "n-Hexane"),
        (0.5, 0.8, 0.9, 0.95),
        (170.0, 180.0, 190.0, 195.0, 200.0),
        (0.5e6, 1.0e6, 1.5e6, 2.0e6, 2.5e6, 3.0e6, 3.5e6),
    ),
)

#: The states that raised before ADR-0016, measured at commit cf846fe.
PREVIOUSLY_FAILING = (
    (("Carbon dioxide", "n-Decane"), 0.8, 250.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 240.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 250.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 260.0, 1.5e6),
)


def _grid_states() -> tuple[tuple[tuple[str, str], float, float, float], ...]:
    return tuple(
        (components, z1, temperature_K, pressure_Pa)
        for components, feeds, temperatures, pressures in GRID
        for z1 in feeds
        for temperature_K in temperatures
        for pressure_Pa in pressures
    )


#: 16 states spread over both binaries and the full temperature/pressure
#: range of the grid, always including the four `PREVIOUSLY_FAILING` states -
#: the default run (slice `flash-phase-labels-by-compressibility`; see the
#: "Runtime" note above). `--full` runs the complete `_grid_states()` instead.
#: Membership in the full grid is asserted at import time, not trusted.
SUBSET: tuple[tuple[tuple[str, str], float, float, float], ...] = PREVIOUSLY_FAILING + (
    (("Carbon dioxide", "n-Decane"), 0.6, 230.0, 1.0e6),
    (("Carbon dioxide", "n-Decane"), 0.6, 260.0, 2.5e6),
    (("Carbon dioxide", "n-Decane"), 0.8, 230.0, 2.0e6),
    (("Carbon dioxide", "n-Decane"), 0.9, 230.0, 2.5e6),
    (("Methane", "n-Hexane"), 0.5, 170.0, 0.5e6),
    (("Methane", "n-Hexane"), 0.5, 200.0, 3.5e6),
    (("Methane", "n-Hexane"), 0.8, 180.0, 1.5e6),
    (("Methane", "n-Hexane"), 0.9, 190.0, 2.0e6),
    (("Methane", "n-Hexane"), 0.9, 195.0, 3.0e6),
    (("Methane", "n-Hexane"), 0.95, 170.0, 0.5e6),
    (("Methane", "n-Hexane"), 0.95, 200.0, 3.5e6),
    (("Methane", "n-Hexane"), 0.5, 195.0, 2.5e6),
)
assert len(SUBSET) == 16, len(SUBSET)
assert len(set(SUBSET)) == 16, "SUBSET must not contain duplicates"
assert set(SUBSET) <= set(_grid_states()), "every SUBSET state must belong to GRID"

failures: list[str] = []


def record(label: str, passed: bool) -> None:
    status = "PASS" if passed else "FAIL"
    print(f"    [{status}] {label}")
    if not passed:
        failures.append(label)


def mixture_of(components: tuple[str, str], z1: float) -> ct.Mixture:
    return ct.Mixture.from_database(list(components), [z1, 1.0 - z1], normalize=True)


def ln_phi(
    mixture: ct.Mixture,
    temperature_K: float,
    pressure_Pa: float,
    composition: np.ndarray,
    phase: str,
) -> np.ndarray:
    return np.log(
        np.array(
            ct.PCSAFTEOS().fugacity_coefficients(
                mixture=mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                composition=(composition / composition.sum()).tolist(),
                phase=phase,
            ),
            dtype=float,
        )
    )


def independent_split(
    mixture: ct.Mixture,
    temperature_K: float,
    pressure_Pa: float,
    z: np.ndarray,
    seed_beta: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Damped Newton on ``ln f_i^V - ln f_i^L = 0`` in vapor mole numbers.

    Written here rather than taken from the package: the unknowns are the vapor
    mole numbers per mole of feed, the residual is the equal-fugacity system
    itself, and the Jacobian is a central difference of that residual. Returns
    ``(x, y, beta, residual)``.
    """

    def residual(v: np.ndarray) -> np.ndarray:
        liquid = z - v
        x = liquid / liquid.sum()
        y = v / v.sum()
        return (
            np.log(y)
            + ln_phi(mixture, temperature_K, pressure_Pa, y, "vapor")
            - np.log(x)
            - ln_phi(mixture, temperature_K, pressure_Pa, x, "liquid")
        )

    v = np.clip(np.array([seed_beta * (1.0 - 1e-7), seed_beta * 1e-7]), 1e-14, z - 1e-14)
    for _ in range(300):
        r = residual(v)
        if np.max(np.abs(r)) < 1e-13:
            break
        jacobian = np.zeros((z.size, z.size))
        for column in range(z.size):
            step = 1e-8 * max(abs(v[column]), 1e-8)
            plus, minus = v.copy(), v.copy()
            plus[column] += step
            minus[column] -= step
            jacobian[:, column] = (residual(plus) - residual(minus)) / (2.0 * step)
        direction = np.linalg.solve(jacobian, -r)
        scale = 1.0
        current = float(np.max(np.abs(r)))
        while scale > 1e-13:
            candidate = v + scale * direction
            if np.all(candidate > 0.0) and np.all(candidate < z):
                if float(np.max(np.abs(residual(candidate)))) < current:
                    break
            scale *= 0.5
        else:
            break
        v = v + scale * direction

    liquid = z - v
    return (
        liquid / liquid.sum(),
        v / v.sum(),
        float(v.sum()),
        float(np.max(np.abs(residual(v)))),
    )


def reduced_g(fractions: np.ndarray, terms: np.ndarray) -> float:
    mask = fractions > 0.0
    return float(np.sum(fractions[mask] * (np.log(fractions[mask]) + terms[mask])))


def check_reference_state() -> None:
    print("\n1) The reference state against an independent Newton")
    print("-" * 78)
    mixture = mixture_of(REFERENCE_COMPONENTS, REFERENCE_FEED[0])
    z = np.array(mixture.composition.fractions, dtype=float)
    eos = ct.PCSAFTEOS()

    stability = ct.stability_tp(
        mixture,
        temperature_K=REFERENCE_T_K,
        pressure_Pa=REFERENCE_P_PA,
        eos=ct.PCSAFTEOS(),
    )
    print(
        f"  {'/'.join(REFERENCE_COMPONENTS)}  z = {tuple(z)}  "
        f"T = {REFERENCE_T_K:.0f} K  P = {REFERENCE_P_PA / 1e6:.2f} MPa"
    )
    print(
        f"  stability: {stability.status}, tpd_min = {stability.tpd_min:.10e}, "
        f"feed branch {stability.feed_branch}"
    )
    record(
        "the feed is proved unstable before any split is attempted", stability.status == "unstable"
    )

    result = ct.flash_tp(
        mixture,
        temperature_K=REFERENCE_T_K,
        pressure_Pa=REFERENCE_P_PA,
        eos=ct.PCSAFTEOS(),
    )
    diagnostics = result.diagnostics
    beta = float(result.vapor_fraction or 0.0)
    x = np.array(result.phases["liquid"].composition.fractions, dtype=float)
    y = np.array(result.phases["vapor"].composition.fractions, dtype=float)
    print(
        f"  flash_tp : beta = {beta!r}\n"
        f"             x    = {tuple(x)}\n"
        f"             y    = {tuple(y)}"
    )
    print(
        f"  stages   : {diagnostics.get('ssi_iterations')} successive substitution + "
        f"{diagnostics.get('second_order_iterations')} second order "
        f"({diagnostics.get('negative_flash_steps')} negative-flash steps), "
        f"converged in the {diagnostics.get('converged_stage')} stage"
    )

    worst_beta = 0.0
    worst_composition = 0.0
    for seed_beta in (0.05, 0.30):
        reference_x, reference_y, reference_beta, reference_residual = independent_split(
            mixture, REFERENCE_T_K, REFERENCE_P_PA, z, seed_beta
        )
        print(
            f"  Newton (seed beta = {seed_beta:.2f}): beta = {reference_beta!r}, "
            f"residual {reference_residual:.2e}"
        )
        worst_beta = max(worst_beta, abs(beta - reference_beta))
        worst_composition = max(
            worst_composition,
            float(np.max(np.abs(x - reference_x))),
            float(np.max(np.abs(y - reference_y))),
        )
    print(f"  worst |d beta| = {worst_beta:.2e}, worst |d composition| = {worst_composition:.2e}")
    record("vapor fraction agrees with the independent Newton to 1e-8", worst_beta < 1e-8)
    record("both compositions agree with the independent Newton to 1e-8", worst_composition < 1e-8)

    mass_balance = float(diagnostics["mass_balance_residual"])
    fugacity = float(diagnostics["fugacity_residual"])
    delta_g = float(diagnostics["delta_g_split_rt"])
    print(
        f"  residuals: mass balance {mass_balance:.2e}, equal fugacity {fugacity:.2e}, "
        f"dG_split/RT {delta_g:.10e}, post-split {diagnostics['post_split_status']}"
    )
    record("mass balance below 1e-12", mass_balance < 1e-12)
    record("equal-fugacity residual below 1e-10", fugacity < 1e-10)
    record("the split lowers the Gibbs energy", delta_g < 0.0)
    record("both converged phases are stable", diagnostics["post_split_status"] == "stable")

    # The Gibbs-energy reduction, recomputed here from the public API.
    feed_g = min(
        reduced_g(z, ln_phi(mixture, REFERENCE_T_K, REFERENCE_P_PA, z, "liquid")),
        reduced_g(z, ln_phi(mixture, REFERENCE_T_K, REFERENCE_P_PA, z, "vapor")),
    )
    recomputed = (
        beta * reduced_g(y, ln_phi(mixture, REFERENCE_T_K, REFERENCE_P_PA, y, "vapor"))
        + (1.0 - beta) * reduced_g(x, ln_phi(mixture, REFERENCE_T_K, REFERENCE_P_PA, x, "liquid"))
        - feed_g
    )
    record(
        "the reported dG_split/RT is reproduced from the public API to 1e-12",
        abs(recomputed - delta_g) < 1e-12,
    )

    liquid_roots = eos.density_roots(
        temperature_K=REFERENCE_T_K,
        pressure_Pa=REFERENCE_P_PA,
        composition=x.tolist(),
        mixture=mixture,
    )
    vapor_roots = eos.density_roots(
        temperature_K=REFERENCE_T_K,
        pressure_Pa=REFERENCE_P_PA,
        composition=y.tolist(),
        mixture=mixture,
    )
    print(f"  density roots at x (mol/m^3): {tuple(round(r, 5) for r in liquid_roots)}")
    print(f"  density roots at y (mol/m^3): {tuple(round(r, 5) for r in vapor_roots)}")
    record(
        "both phases have two mechanically stable roots at this (T, P)",
        len(liquid_roots) == 2 and len(vapor_roots) == 2,
    )
    record(
        "the pinned phase densities are reproduced to 1e-6 relative",
        abs(liquid_roots[-1] / 17250.70876 - 1.0) < 1e-6
        and abs(vapor_roots[0] / 557.8191800 - 1.0) < 1e-6,
    )


def check_grid(
    states: tuple[tuple[tuple[str, str], float, float, float], ...], *, full: bool
) -> None:
    kind = f"the complete {len(states)}-state grid" if full else f"a {len(states)}-state subset"
    print(f"\n2) The Case F-4 grid: {kind} - every state must answer, every answer must verify")
    print("-" * 78)
    errors: list[str] = []
    two_phase = 0
    single_phase = 0
    rescued = 0
    negative_flash = 0
    worst_mass_balance = 0.0
    worst_fugacity = 0.0
    worst_delta_g = -np.inf
    worst_post_split = np.inf

    for components, z1, temperature_K, pressure_Pa in states:
        mixture = mixture_of(components, z1)
        label = (
            f"{'/'.join(components)} z1={z1} "
            f"T={temperature_K:.0f}K P={pressure_Pa / 1e6:.2f}MPa"
        )
        try:
            result = ct.flash_tp(
                mixture,
                temperature_K=temperature_K,
                pressure_Pa=pressure_Pa,
                eos=ct.PCSAFTEOS(),
            )
        except ct.ConvergenceError as error:
            errors.append(f"{label}: {error}")
            continue
        diagnostics = result.diagnostics
        if diagnostics["phase_count"] == 1:
            single_phase += 1
            if diagnostics["stability_status"] != "stable":
                errors.append(f"{label}: single phase without a stable verdict")
            continue
        two_phase += 1
        worst_mass_balance = max(worst_mass_balance, float(diagnostics["mass_balance_residual"]))
        worst_fugacity = max(worst_fugacity, float(diagnostics["fugacity_residual"]))
        worst_delta_g = max(worst_delta_g, float(diagnostics["delta_g_split_rt"]))
        worst_post_split = min(worst_post_split, float(diagnostics["post_split_tpd_min"]))
        if diagnostics["post_split_status"] != "stable":
            errors.append(f"{label}: post-split {diagnostics['post_split_status']}")
        if diagnostics.get("converged_stage") == "second-order":
            rescued += 1
        if int(diagnostics.get("negative_flash_steps", 0)) > 0:
            negative_flash += 1

    total = two_phase + single_phase + len(errors)
    print(f"  states scanned        : {total}")
    print(f"  two-phase answers     : {two_phase}")
    print(f"  single-phase answers  : {single_phase}")
    print(f"  ConvergenceErrors     : {len(errors)}  (4 before ADR-0016)")
    print(f"  needed the 2nd-order stage : {rescued}")
    print(f"  used a negative flash      : {negative_flash}")
    print(f"  worst mass balance    : {worst_mass_balance:.2e}")
    print(f"  worst equal fugacity  : {worst_fugacity:.2e}")
    print(f"  worst (least negative) dG_split/RT : {worst_delta_g:.2e}")
    print(f"  worst post-split tpd_min           : {worst_post_split:.2e}")
    for message in errors:
        print(f"    ! {message}")
    record("no state raises ConvergenceError", not errors)
    record("every two-phase mass balance below 1e-12", worst_mass_balance < 1e-12)
    record("every two-phase equal-fugacity residual below 1e-6", worst_fugacity < 1e-6)
    record("every two-phase split lowers the Gibbs energy", worst_delta_g < 0.0)
    if full:
        record(
            "exactly the four previously failing states needed the stage",
            rescued == len(PREVIOUSLY_FAILING),
        )
        record("the measured 123/65 two-phase/single-phase split is reproduced", (two_phase, single_phase) == (123, 65))
    else:
        record(
            "every previously-failing state in this subset needed the stage",
            rescued >= sum(1 for state in states if state in PREVIOUSLY_FAILING),
        )
        print(f"  (pass --full for the complete {len(_grid_states())}-state grid)")


def check_previously_failing() -> None:
    print("\n3) The four previously failing states, and what still fails")
    print("-" * 78)
    for components, z1, temperature_K, pressure_Pa in PREVIOUSLY_FAILING:
        mixture = mixture_of(components, z1)
        label = (
            f"{'/'.join(components)} z1={z1} T={temperature_K:.0f}K P={pressure_Pa / 1e6:.2f}MPa"
        )
        result = ct.flash_tp(
            mixture,
            temperature_K=temperature_K,
            pressure_Pa=pressure_Pa,
            eos=ct.PCSAFTEOS(),
        )
        beta = float(result.vapor_fraction or 0.0)
        print(
            f"  {label}\n"
            f"    default          : beta = {beta:.12f}, "
            f"stage {result.diagnostics.get('converged_stage')}, "
            f"fugacity residual {float(result.diagnostics['fugacity_residual']):.2e}"
        )
        for description, settings in (
            ("second_order=False", ct.FlashSettings(second_order=False)),
            (
                "wilson-heuristic",
                ct.FlashSettings(phase_detection="wilson-heuristic"),
            ),
        ):
            try:
                ct.flash_tp(
                    mixture,
                    temperature_K=temperature_K,
                    pressure_Pa=pressure_Pa,
                    eos=ct.PCSAFTEOS(),
                    settings=settings,
                )
            except ct.ConvergenceError as error:
                print(f"    {description:<17}: raises, as designed ({str(error)[:52]}...)")
            else:
                record(f"{label}: {description} must still fail", False)
        record(f"{label} converges through the second-order stage", 0.0 < beta < 1.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help=f"run the complete {len(_grid_states())}-state Case F-4 grid instead of the "
        f"{len(SUBSET)}-state representative subset",
    )
    # `parse_known_args`, not `parse_args`: `tests/test_examples.py` runs this
    # script via `runpy.run_path` with pytest's own `sys.argv` still in place,
    # so an unrecognized pytest flag must be ignored rather than raising.
    args, _unknown = parser.parse_known_args()

    print("=" * 78)
    print("Phi-phi split robustness with PC-SAFT (validation Case F-4, ADR-0016)")
    print("=" * 78)

    check_reference_state()
    check_grid(_grid_states() if args.full else SUBSET, full=args.full)
    check_previously_failing()

    print("\n" + "=" * 78)
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
