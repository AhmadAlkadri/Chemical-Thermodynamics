"""The three-phase neighbourhood of water / n-hexane, with PASS/FAIL on every check.

What is being validated
-----------------------
ADR-0020 wires the ADR-0011 phase addition / removal search to the
equation-of-state (phi-phi) path, with each phase pinned to its own density
root (ADR-0019). The consequence on this binary at 1 atm is the band of
temperatures just below the three-phase temperature ``T3``: the tangent-plane
test finds a **vapour** as the deepest minimum from a 50/50 feed, the
vapour-liquid pair that follows is unstable towards a second liquid, and the
three-phase solve then drives the vapour amount negative. The answer is two
conjugate liquids, reached by adding a phase and removing one. Before this
slice those temperatures raised ``ConvergenceError``.

Four independent routes are used, so no single one carries the result:

1. **A 4-equation Newton written in this script** locates ``T3`` and the three
   coexisting compositions from the public ``fugacity_coefficients`` interface
   alone - no ``flash_tp``.
2. **Two-equation Newtons** on each side of ``T3`` give the liquid-liquid and
   the vapour-liquid tie lines independently, and their **reduced Gibbs
   energies** say which pair is the equilibrium at each temperature.
3. **FeOs** (feos-org/feos, MIT OR Apache-2.0; the same Gross & Sadowski model
   in Rust, every derivative by automatic differentiation) supplies chemical
   potentials **at chemthermo's converged phases and densities**. That check
   does not depend on either flash converging.
4. **FeOs's own two-phase flash**, which at this state converges on the
   *vapour-liquid* pair - the one chemthermo's post-split stability test
   rejects. Printed and compared rather than hidden: both are stationary states
   of the same model, and the Gibbs comparison in (2) says which is the
   equilibrium.

As in Cases P-6, P-7 and P-8, one input is deliberately not shared: the 42
universal constants of the 2001 dispersion term (chemthermo packages the ten
figures as printed, FeOs hard-codes fourteen). Route 3 is therefore reported
twice, as shipped and with FeOs's constants substituted.

``--full`` adds the 41-point scan across ``[T3 - 1 K, T3 + 1 K]`` at two feeds
and the ternary water / ethanol / n-hexane tie triangle (validation Case P-10).

Requires the optional ``feos`` dependency for routes 3 and 4::

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
TERNARY = ("Water", "Ethanol", "n-Hexane")
ATMOSPHERE_PA = 101325.0
OFFSET_K = 0.05
TERNARY_T_K = 333.0

#: Gross & Sadowski (2002) Table 1 for water and ethanol, (2001) Table 1 for
#: n-hexane, written out here so the reference model is built from this file.
PARAMETERS: dict[str, tuple[float, float, float, float, float | None, float]] = {
    "Water": (18.015, 1.0656, 3.0007, 366.51, 0.034868, 2500.7),
    "Ethanol": (46.069, 2.3827, 3.1771, 198.24, 0.032384, 2653.4),
    "n-Hexane": (86.177, 3.0576, 3.7983, 236.77, None, 0.0),
}

failures: list[str] = []


def record(label: str, ok: bool) -> None:
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        failures.append(label)


# ---------------------------------------------------------------------------
# The chemthermo side, through the public interface only
# ---------------------------------------------------------------------------


def ln_f(
    names: Sequence[str],
    temperature_K: float,
    composition: Sequence[float] | np.ndarray,
    branch: str,
) -> np.ndarray:
    """``ln(x_i phi_i)`` on one named density branch."""
    mixture = ct.Mixture.from_database(list(names), [1.0 / len(names)] * len(names))
    x = np.asarray(composition, dtype=float)
    x = x / float(np.sum(x))
    phi = np.asarray(
        PCSAFTEOS().fugacity_coefficients(
            mixture=mixture,
            temperature_K=temperature_K,
            pressure_Pa=ATMOSPHERE_PA,
            composition=x.tolist(),
            phase=branch,
        ),
        dtype=float,
    )
    return np.log(x) + np.log(phi)


def reduced_g(
    names: Sequence[str],
    temperature_K: float,
    composition: Sequence[float] | np.ndarray,
    branch: str,
) -> float:
    x = np.asarray(composition, dtype=float)
    return float(np.sum(x * ln_f(names, temperature_K, x, branch)))


def damped_newton(residual, start: np.ndarray, steps: Sequence[float], admissible, tol: float):
    u = np.array(start, dtype=float)
    worst = float(np.max(np.abs(residual(u))))
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


def three_phase_point() -> tuple[float, float, float, float]:
    print("\n1) The three-phase point, from a 4-equation Newton (no flash_tp)")
    print("-" * 78)

    def residual(u: np.ndarray) -> np.ndarray:
        first = ln_f(BINARY, u[3], [u[0], 1.0 - u[0]], "liquid")
        second = ln_f(BINARY, u[3], [u[1], 1.0 - u[1]], "liquid")
        vapor = ln_f(BINARY, u[3], [u[2], 1.0 - u[2]], "vapor")
        return np.concatenate([first - second, first - vapor])

    u, worst = damped_newton(
        residual,
        np.array([0.9999, 0.02, 0.20, 334.5]),
        (1e-8, 1e-8, 1e-8, 3.345e-5),
        lambda v: bool(np.all(v[:3] > 0.0) and np.all(v[:3] < 1.0) and 200.0 < v[3] < 600.0),
        1e-11,
    )
    print(f"    residual    {worst:.3e}")
    print(f"    T3          {u[3]:.9f} K   ({u[3] - 273.15:.4f} C)")
    print(f"    x_water(I)  {u[0]:.12f}")
    print(f"    x_water(II) {u[1]:.12f}")
    print(f"    y_water     {u[2]:.12f}")
    record("the 4-equation Newton converged", worst < 1e-10)

    eos = PCSAFTEOS()
    mixture = ct.Mixture.from_database(list(BINARY), [0.5, 0.5], normalize=True)
    both = True
    for label, value in (("I", u[0]), ("II", u[1]), ("V", u[2])):
        roots = eos.density_roots(
            mixture=mixture,
            temperature_K=float(u[3]),
            pressure_Pa=ATMOSPHERE_PA,
            composition=[float(value), 1.0 - float(value)],
        )
        print(f"    density roots of {label:<2} {', '.join(f'{r:,.2f}' for r in roots)} mol/m^3")
        both = both and len(roots) == 2
    record("each of the three compositions has both a vapour and a liquid root", both)
    print(
        "    (which is exactly why each phase must be pinned to its own root, ADR-0019;\n"
        "     a single liquid/vapour branch assignment cannot hold two liquids and a vapour)"
    )
    print(
        "    Model versus experiment, a remark and not an assertion: the water / n-hexane\n"
        "    heteroazeotrope at 1 atm is commonly tabulated near 61.6 C with y_water ~ 0.21.\n"
        "    No primary source was verified for those figures here."
    )
    return float(u[3]), float(u[0]), float(u[1]), float(u[2])


def binary_pair(temperature: float, branches: tuple[str, str], start: tuple[float, float]):
    def residual(u: np.ndarray) -> np.ndarray:
        return ln_f(BINARY, temperature, [u[0], 1.0 - u[0]], branches[0]) - ln_f(
            BINARY, temperature, [u[1], 1.0 - u[1]], branches[1]
        )

    return damped_newton(
        residual,
        np.array(start, dtype=float),
        (1e-8, 1e-8),
        lambda v: bool(np.all(v > 0.0) and np.all(v < 1.0)),
        1e-11,
    )


def pair_energy(temperature: float, pair, branches: tuple[str, str]) -> tuple[float, float]:
    beta = (0.5 - pair[0]) / (pair[1] - pair[0])
    energy = (1.0 - beta) * reduced_g(
        BINARY, temperature, [pair[0], 1.0 - pair[0]], branches[0]
    ) + beta * reduced_g(BINARY, temperature, [pair[1], 1.0 - pair[1]], branches[1])
    return beta, energy


def either_side(t3: float) -> dict[float, ct.FlashResult]:
    results: dict[float, ct.FlashResult] = {}
    for section, offset in ((2, -OFFSET_K), (3, +OFFSET_K)):
        temperature = t3 + offset
        side = "below" if offset < 0 else "above"
        print(f"\n{section}) T = T3 {offset:+.2f} K = {temperature:.6f} K ({side} T3)")
        print("-" * 78)
        result = ct.flash_tp(
            ct.Mixture.from_database(list(BINARY), [0.5, 0.5], normalize=True),
            temperature_K=temperature,
            pressure_Pa=ATMOSPHERE_PA,
            eos=PCSAFTEOS(),
        )
        results[offset] = result
        names = sorted(result.phases)
        print(f"    flash_tp -> {names}, regime {result.diagnostics['phase_regime']}")
        print(
            "    phase_set_history     "
            f"{result.diagnostics.get('phase_set_history', '(search not entered)')}"
        )
        for name in names:
            fractions = result.phases[name].composition.fractions
            print(
                f"      {name:<8} beta = {result.phase_fractions[name]:.10f}  "
                f"x_water = {fractions[0]:.12f}"
            )

        liquids, liquid_residual = binary_pair(temperature, ("liquid", "liquid"), (0.9999, 0.0226))
        vapor_liquid, vl_residual = binary_pair(temperature, ("liquid", "vapor"), (0.9999, 0.213))
        record(
            "both independent two-phase Newtons converged",
            max(liquid_residual, vl_residual) < 1e-10,
        )
        beta_ll, g_ll = pair_energy(temperature, liquids, ("liquid", "liquid"))
        beta_vl, g_vl = pair_energy(temperature, vapor_liquid, ("liquid", "vapor"))
        g_feed = min(
            reduced_g(BINARY, temperature, [0.5, 0.5], branch) for branch in ("liquid", "vapor")
        )
        print(f"    independent LL tie line  x_water = {liquids[0]:.12f} / {liquids[1]:.12f}")
        print(
            f"    independent VL tie line  x_water = {vapor_liquid[0]:.12f} / {vapor_liquid[1]:.12f}"
        )
        print(f"    G/RT  two liquids {g_ll:.12f}   vapour-liquid {g_vl:.12f}   feed {g_feed:.12f}")

        if offset < 0:
            record("the two liquids have the lower Gibbs energy below T3", g_ll < g_vl < g_feed)
            record("flash_tp returned the two liquids", names == ["liquid1", "liquid2"])
            record(
                "the search route was V -> LV -> LLV -> LL",
                result.diagnostics.get("phase_set_history") == "V -> LV -> LLV -> LL",
            )
            ours = (
                float(result.phases["liquid1"].composition.fractions[0]),
                float(result.phases["liquid2"].composition.fractions[0]),
            )
            worst = max(abs(ours[0] - liquids[0]), abs(ours[1] - liquids[1]))
            print(f"    |dx| against the independent LL Newton: {worst:.3e}")
            record("the tie line matches the independent Newton to 1e-9", worst < 1e-9)
            lever = abs(float(result.phase_fractions["liquid2"]) - beta_ll)
            print(f"    |d beta| against the lever rule:        {lever:.3e}")
            record("the amounts match the lever rule to 1e-9", lever < 1e-9)
            gap = abs(float(result.diagnostics["delta_g_vs_two_phase_rt"]) - (g_ll - g_vl))
            record("delta_g_vs_two_phase_rt is G(LL) - G(VL) to 1e-9", gap < 1e-9)
        else:
            record("the vapour-liquid pair has the lower Gibbs energy above T3", g_vl < g_ll)
            record("flash_tp returned a vapour-liquid pair", names == ["liquid", "vapor"])
            record("the search was not entered", "phase_set_history" not in result.diagnostics)
            ours = (
                float(result.phases["liquid"].composition.fractions[0]),
                float(result.phases["vapor"].composition.fractions[0]),
            )
            worst = max(abs(ours[0] - vapor_liquid[0]), abs(ours[1] - vapor_liquid[1]))
            print(f"    |dx| against the independent VL Newton: {worst:.3e}")
            record("the tie line matches the independent Newton to 1e-9", worst < 1e-9)
            print(f"    (beta from the lever rule: {beta_vl:.10f})")
        record(
            "every phase is post-split stable", result.diagnostics["post_split_status"] == "stable"
        )
    return results


# ---------------------------------------------------------------------------
# The FeOs side
# ---------------------------------------------------------------------------


def _pure_record(name: str) -> "PureRecord":
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


def feos_eos(names: Sequence[str]) -> "EquationOfState":
    return EquationOfState.pcsaft(Parameters.from_records([_pure_record(n) for n in names]))


def feos_reduced_potentials(
    names: Sequence[str], temperature_K: float, density: float, x: Sequence[float]
) -> np.ndarray:
    values = np.asarray(x, dtype=float)
    state = State(
        feos_eos(names),
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


def our_densities(
    names: Sequence[str], temperature_K: float, result: ct.FlashResult
) -> list[tuple[np.ndarray, float]]:
    mixture = ct.Mixture.from_database(list(names), [1.0 / len(names)] * len(names))
    eos = PCSAFTEOS()
    out = []
    for name in sorted(result.phases):
        x = np.asarray(result.phases[name].composition.fractions, dtype=float)
        roots = eos.density_roots(
            mixture=mixture,
            temperature_K=temperature_K,
            pressure_Pa=ATMOSPHERE_PA,
            composition=x.tolist(),
        )
        out.append((x, roots[0] if name == "vapor" else roots[-1]))
    return out


def worst_potential_gap(names: Sequence[str], temperature_K: float, phases) -> float:
    potentials = [feos_reduced_potentials(names, temperature_K, rho, x) for x, rho in phases]
    return max(
        float(np.max(np.abs(potentials[i] - potentials[j])))
        for i in range(len(potentials))
        for j in range(i + 1, len(potentials))
    )


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


def the_chemical_potentials(results: dict[float, ct.FlashResult], t3: float, full: bool) -> None:
    print("\n4) FeOs's chemical potentials at chemthermo's phases")
    print("-" * 78)
    shipped_a, shipped_b = pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL
    matched_a, matched_b = feos_constants()
    # The matched-constants number needs the flash re-run with FeOs's table, so
    # by default only the below-T3 state - the one this slice made reachable -
    # pays for it; `--full` does both sides.
    matched_offsets = sorted(results) if full else [min(results)]
    for offset, result in sorted(results.items()):
        temperature = t3 + offset
        phases = our_densities(BINARY, temperature, result)
        as_shipped = worst_potential_gap(BINARY, temperature, phases)
        label = f"T3 {offset:+.2f} K ({sorted(result.phases)})"
        print(f"    {label}")
        print(f"      max_i |mu_i^a - mu_i^b| / RT   as shipped {as_shipped:.3e}")
        if offset not in matched_offsets:
            print("                                     matched    (pass --full)")
            continue
        pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = matched_a, matched_b
        try:
            matched_result = ct.flash_tp(
                ct.Mixture.from_database(list(BINARY), [0.5, 0.5], normalize=True),
                temperature_K=temperature,
                pressure_Pa=ATMOSPHERE_PA,
                eos=PCSAFTEOS(),
            )
            matched = worst_potential_gap(
                BINARY, temperature, our_densities(BINARY, temperature, matched_result)
            )
        finally:
            pcsaft_module.A_UNIVERSAL, pcsaft_module.B_UNIVERSAL = shipped_a, shipped_b
        print(f"                                     matched    {matched:.3e}")
        record(f"{label}: FeOs's potentials are equal to 1e-8 (matched constants)", matched < 1e-8)
    print(
        "    The 'as shipped' floor is the universal-constants table, not either solver:\n"
        "    chemthermo packages ten figures as printed, FeOs hard-codes fourteen."
    )


def the_reference_flash(t3: float) -> None:
    print("\n5) FeOs's own two-phase flash below T3")
    print("-" * 78)
    temperature = t3 - OFFSET_K
    state = State(
        feos_eos(BINARY),
        temperature=temperature * si.KELVIN,
        pressure=ATMOSPHERE_PA * si.PASCAL,
        composition=np.array([0.5, 0.5]),
    )
    equilibrium = state.tp_flash()
    mol_per_m3 = si.MOL / si.METER**3
    dense = np.atleast_1d(np.asarray(equilibrium.liquid.molefracs, dtype=float))
    light = np.atleast_1d(np.asarray(equilibrium.vapor.molefracs, dtype=float))
    rho_dense = float(equilibrium.liquid.density / mol_per_m3)
    rho_light = float(equilibrium.vapor.density / mol_per_m3)
    print(f"    FeOs liquid  x_water = {dense[0]:.12f}  rho = {rho_dense:,.4f} mol/m^3")
    print(f"    FeOs vapor   x_water = {light[0]:.12f}  rho = {rho_light:,.4f} mol/m^3")
    print(f"    FeOs vapour-phase fraction {float(equilibrium.vapor_phase_fraction):.10f}")
    record("FeOs's own flash returns a vapour-liquid pair here", rho_light < 200.0)

    unpoliced = ct.flash_tp(
        ct.Mixture.from_database(list(BINARY), [0.5, 0.5], normalize=True),
        temperature_K=temperature,
        pressure_Pa=ATMOSPHERE_PA,
        eos=PCSAFTEOS(),
        settings=ct.FlashSettings(post_split_stability=False),
    )
    worst = max(
        abs(float(unpoliced.phases["liquid"].composition.fractions[0]) - dense[0]),
        abs(float(unpoliced.phases["vapor"].composition.fractions[0]) - light[0]),
    )
    print(f"    chemthermo with post_split_stability=False finds the same pair to {worst:.3e}")
    record("chemthermo converges the same vapour-liquid pair first", worst < 1e-8)
    record(
        "and its post-split stability test refuses it",
        unpoliced.diagnostics["post_split_stable"] is False,
    )
    print(
        "    Not a disagreement about the model: both pairs are stationary states of it,\n"
        "    and section 2's Gibbs comparison says the two liquids are the equilibrium."
    )


# ---------------------------------------------------------------------------
# --full: the scan and the ternary triangle
# ---------------------------------------------------------------------------


def the_scan(t3: float) -> None:
    print("\n6) 41 temperatures across [T3 - 1 K, T3 + 1 K], two feeds")
    print("-" * 78)
    for z_water in (0.3, 0.7):
        verdicts = []
        raised = 0
        for temperature in np.linspace(t3 - 1.0, t3 + 1.0, 41):
            try:
                result = ct.flash_tp(
                    ct.Mixture.from_database(
                        list(BINARY), [z_water, 1.0 - z_water], normalize=True
                    ),
                    temperature_K=float(temperature),
                    pressure_Pa=ATMOSPHERE_PA,
                    eos=PCSAFTEOS(),
                )
                verdicts.append(str(result.diagnostics["phase_regime"]))
            except ct.ConvergenceError:
                raised += 1
                verdicts.append("RAISE")
        counts = {value: verdicts.count(value) for value in sorted(set(verdicts))}
        print(f"    z_water = {z_water}: {counts}")
        record(f"z_water = {z_water}: no ConvergenceError in 41 points", raised == 0)
        if z_water == 0.3:
            switches = [i for i in range(1, 41) if verdicts[i] != verdicts[i - 1]]
            boundary = float(np.linspace(t3 - 1.0, t3 + 1.0, 41)[switches[0]]) if switches else 0.0
            print(f"    one switch LLE -> VLE at {boundary:.6f} K (T3 = {t3:.6f} K)")
            record("z_water = 0.3: exactly one verdict switch, at T3", len(switches) == 1)
        else:
            print(
                "    z_water = 0.7 stays LLE on both sides. Above T3 that pair is\n"
                "    metastable: the deterministic stability trial set does not find the\n"
                "    vapour stationary point from the hexane-rich liquid. A limitation of\n"
                "    the stability test, not of the search - recorded, not worked around."
            )


def the_ternary(full_feos: bool) -> None:
    print("\n7) The ternary vapour-liquid-liquid tie triangle (validation Case P-10)")
    print("-" * 78)
    print(f"    water / ethanol / n-hexane, {TERNARY_T_K} K, 1 atm, k_ij = 0")
    vertices = None
    for feed in ((0.4, 0.3, 0.3), (0.5, 0.2, 0.3)):
        result = ct.flash_tp(
            ct.Mixture.from_database(list(TERNARY), list(feed), normalize=True),
            temperature_K=TERNARY_T_K,
            pressure_Pa=ATMOSPHERE_PA,
            eos=PCSAFTEOS(),
        )
        names = sorted(result.phases)
        print(f"\n    z = {feed} -> {names}, regime {result.diagnostics['phase_regime']}")
        for name in names:
            x = result.phases[name].composition.fractions
            print(
                f"      {name:<8} beta = {result.phase_fractions[name]:.10f}  "
                f"x = ({x[0]:.10f}, {x[1]:.10f}, {x[2]:.10f})"
            )
        record(f"z = {feed} returns three phases", len(names) == 3)
        record(f"z = {feed} is VLLE", result.diagnostics["phase_regime"] == "VLLE")
        current = {name: np.asarray(result.phases[name].composition.fractions) for name in names}
        if vertices is None:
            vertices = current
        else:
            moved = max(float(np.max(np.abs(current[n] - vertices[n]))) for n in current)
            print(f"      same triangle as the first feed to {moved:.3e}")
            record("the triangle does not move with the feed", moved < 1e-8)
        if full_feos and si is not None:
            gap = worst_potential_gap(
                TERNARY, TERNARY_T_K, our_densities(TERNARY, TERNARY_T_K, result)
            )
            print(f"      FeOs's chemical potentials across the three phases: {gap:.3e}")
            record("FeOs's potentials agree across the triangle (as shipped)", gap < 1e-5)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full",
        action="store_true",
        help="also run the 41-point scan at two feeds and the ternary tie triangle "
        "(several minutes)",
    )
    # `parse_known_args`, not `parse_args`: `tests/test_examples.py` runs this
    # script via `runpy.run_path` with pytest's own `sys.argv` still in place.
    args, _unknown = parser.parse_known_args()

    print("=" * 78)
    print("The three-phase neighbourhood of water / n-hexane (Cases P-9, P-10; ADR-0020)")
    print("=" * 78)
    t3, _x_first, _x_second, _y = three_phase_point()
    results = either_side(t3)
    if si is None:
        print("\nfeos is not installed; routes 3 and 4 are skipped.")
        print("Install with: pip install -e '.[validation]'")
    else:
        the_chemical_potentials(results, t3, args.full)
        the_reference_flash(t3)
    if args.full:
        the_scan(t3)
        the_ternary(si is not None)
    else:
        print("\n  (pass --full for the 41-point scan and the ternary tie triangle)")

    print("\n" + "=" * 78)
    print(
        "k_ij = 0 throughout; nothing above is compared against measurement. A binary at\n"
        "a fixed pressure has no three-phase region - Gibbs' phase rule leaves one degree\n"
        "of freedom, so three phases meet at a single temperature and their amounts are\n"
        "not fixed by the mass balance there. The finite three-phase region is the\n"
        "ternary of section 7."
    )
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
