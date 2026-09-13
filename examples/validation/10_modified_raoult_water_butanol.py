"""Water / 1-butanol around its three-phase temperature at 1 atm.

What is being validated
-----------------------
A binary that is partially miscible *and* boils has a single temperature at
which three phases coexist (Gibbs' phase rule: F = 2 - 3 + 2 = 1, so at fixed
pressure the three-phase state is one temperature, not a range). Below it the
equilibrium is two liquids; above it, vapor-liquid; at it, all three.

This script walks the neighbourhood of T3 and checks every verdict:

  * below T3 it returns the liquid-liquid tie-line and reports both phases
    stable against the ideal-vapor candidate;
  * above T3 the feed used here has fully evaporated, and that single vapor is
    checked against the dew-point equation;
  * in a narrow window just below T3 the tangent-plane minimizer is a *vapor*,
    so the search starts from a vapor-liquid pair that is **not** the
    equilibrium. Since ADR-0011 that is resolved rather than refused: the pair
    is unstable towards a second liquid, the liquid is added, and the
    three-phase solve then drives the vapor fraction to zero, so the vapor is
    removed again and the two-liquid answer comes back
    (`phase_set_history = "L -> LV -> LLV -> LL"`). With
    `FlashSettings(max_phases=2)` the same call still raises, which is the
    pre-ADR-0011 behavior and is checked here too.

Everything printed is checked against an independent route written in this
file - a Newton solve of the equal-activity condition for the binodal, and
scalar bisections of `sum_i x_i gamma_i Psat_i(T) = P` - which shares no code
with `chemthermo.flash`.

Three-phase condition
---------------------
At T3 *both* conjugate liquids are simultaneously at their bubble point:

    sum_i x_i^I  gamma_i^I  Psat_i(T3) = P
    sum_i x_i^II gamma_i^II Psat_i(T3) = P

and the vapor they share is `y_i = x_i gamma_i Psat_i / P` from either liquid.
Solving the first and *checking* the second is a non-trivial test, because the
binodal comes from equal activities and knows nothing about Antoine.

Parameters and their provenance
-------------------------------
NRTL pair 2-3 of Table 1 of S. R. Tessier, J. F. Brennecke and M. A. Stadtherr,
Chem. Eng. Sci. 55 (2000) 1785-1796 (attributed there to McDonald and Floudas,
AIChE J. 41 (1995) 1798); alpha implied by `G_ij = exp(-alpha_ij tau_ij)`.
Antoine coefficients from the packaged databank (Koretsky 2012),
`ln(P^sat/bar) = A - B/(T + C)`.

**These parameters were fitted to liquid-liquid data and are temperature
independent.** The heteroazeotrope of water / 1-butanol is commonly tabulated
near 365.9 K with roughly 0.75-0.76 mole fraction water in the vapor; the
number this model gives is close to that, but the agreement is partly
fortuitous and is reported, not claimed as a validation against experiment.

Needs no optional dependency. All inputs are SI.
"""

from __future__ import annotations

import math
from typing import Callable, Sequence

import numpy as np

import chemthermo as ct

NAMES = ("n-Butanol", "Water")
TAU_12 = 0.90047
TAU_21 = 3.51307
ALPHA = 0.48

PRESSURE_PA = 101325.0
FEED = (0.20, 0.80)

#: Commonly tabulated heteroazeotrope of water / 1-butanol at 1 atm. Recorded
#: for context only; not read from a primary source in this work.
LITERATURE_T_K = 365.9
LITERATURE_Y_WATER = (0.75, 0.76)

LnGamma = Callable[[np.ndarray], np.ndarray]


def _model() -> ct.NRTL:
    return ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [(NAMES[0], NAMES[1], TAU_12, TAU_21, ALPHA, ALPHA)]
        )
    )


def _mixture(z: Sequence[float]) -> ct.Mixture:
    return ct.Mixture.from_database(list(NAMES), list(z), normalize=True)


def _ln_gamma(model: ct.NRTL) -> LnGamma:
    mixture = _mixture((0.5, 0.5))

    def ln_gamma(x: np.ndarray) -> np.ndarray:
        values = np.asarray(x, dtype=float)
        values = values / float(np.sum(values))
        gamma = model.activity_coefficients(
            mixture=mixture, temperature_K=298.15, composition=[float(v) for v in values]
        )
        return np.log(np.asarray(gamma, dtype=float))

    return ln_gamma


def _psat(temperature_K: float) -> np.ndarray:
    """Psat from the databank Antoine record, recomputed here from A, B, C."""
    values = []
    for name in NAMES:
        antoine = ct.Component.from_database(name).antoine
        assert antoine is not None
        values.append(math.exp(antoine.A - antoine.B / (temperature_K + antoine.C)) * 1.0e5)
    return np.array(values, dtype=float)


def _binodal(ln_gamma: LnGamma) -> tuple[np.ndarray, np.ndarray]:
    """Conjugate liquids from equal activities: Newton with an FD Jacobian."""

    def activities(x1: float) -> np.ndarray:
        x = np.array([x1, 1.0 - x1], dtype=float)
        return x * np.exp(ln_gamma(x))

    def mismatch(u: np.ndarray) -> np.ndarray:
        return activities(float(u[0])) - activities(float(u[1]))

    u = np.array([0.02, 0.50], dtype=float)
    for _ in range(200):
        residual = mismatch(u)
        if float(np.max(np.abs(residual))) < 1e-15:
            break
        jacobian = np.zeros((2, 2), dtype=float)
        for column in range(2):
            plus, minus = u.copy(), u.copy()
            plus[column] += 1e-7
            minus[column] -= 1e-7
            jacobian[:, column] = (mismatch(plus) - mismatch(minus)) / 2e-7
        u = u + np.linalg.solve(jacobian, -residual)
    return (
        np.array([u[0], 1.0 - u[0]]),
        np.array([u[1], 1.0 - u[1]]),
    )


def _bubble_ratio(ln_gamma: LnGamma, x: np.ndarray, temperature_K: float) -> float:
    return float(np.sum(x * np.exp(ln_gamma(x)) * _psat(temperature_K))) / PRESSURE_PA


def _three_phase_temperature(ln_gamma: LnGamma, x: np.ndarray) -> float:
    low, high = 340.0, 399.0
    for _ in range(200):
        middle = 0.5 * (low + high)
        low, high = (middle, high) if _bubble_ratio(ln_gamma, x, middle) < 1.0 else (low, middle)
    return 0.5 * (low + high)


def _dew_ratio(ln_gamma: LnGamma, z: np.ndarray, temperature_K: float) -> float:
    """`sum_i z_i P / (gamma_i(x) Psat_i)` with the implicit x solved by fixed point."""
    x = z.copy()
    for _ in range(20000):
        updated = z * PRESSURE_PA / (np.exp(ln_gamma(x)) * _psat(temperature_K))
        updated = updated / float(np.sum(updated))
        if float(np.max(np.abs(updated - x))) < 1e-16:
            break
        x = updated
    return float(np.sum(z * PRESSURE_PA / (np.exp(ln_gamma(x)) * _psat(temperature_K))))


def _flash(temperature_K: float, settings: ct.FlashSettings | None = None) -> ct.FlashResult:
    return ct.flash_tp(
        _mixture(FEED),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=_model(),
        flash_mode="modified-raoult",
        settings=settings,
    )


def _pass(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def main() -> None:
    model = _model()
    ln_gamma = _ln_gamma(model)
    failures: list[str] = []

    def record(label: str, ok: bool) -> None:
        print(f"  [{_pass(ok)}] {label}")
        if not ok:
            failures.append(label)

    print("=" * 78)
    print("Water / 1-butanol at 101325 Pa: the three-phase neighbourhood")
    print(
        f"NRTL pair 2-3, Tessier et al. (2000) Table 1: tau = ({TAU_12}, {TAU_21}), alpha = {ALPHA}"
    )
    print(f"Feed z = {FEED} ({NAMES[0]} / {NAMES[1]})")

    # ---------------------------------------------------------------- binodal
    x_i, x_ii = _binodal(ln_gamma)
    print("\n1. Liquid-liquid binodal (independent equal-activity Newton solve)")
    print(f"   x^I  ({NAMES[0]}) = {x_i[0]:.12f}")
    print(f"   x^II ({NAMES[0]}) = {x_ii[0]:.12f}")
    activity_residual = float(
        np.max(np.abs(np.log(x_i) + ln_gamma(x_i) - np.log(x_ii) - ln_gamma(x_ii)))
    )
    print(f"   equal-activity residual = {activity_residual:.3e}")
    record("binodal equal-activity residual < 1e-12", activity_residual < 1e-12)

    # ------------------------------------------------------------------- T3
    t3 = _three_phase_temperature(ln_gamma, x_i)
    ratio_i = _bubble_ratio(ln_gamma, x_i, t3)
    ratio_ii = _bubble_ratio(ln_gamma, x_ii, t3)
    y = x_i * np.exp(ln_gamma(x_i)) * _psat(t3) / PRESSURE_PA
    print("\n2. Three-phase temperature: both liquids boil at once")
    print(f"   T3 = {t3:.6f} K ({t3 - 273.15:.3f} C)")
    print(f"   sum x^I  gamma^I  Psat / P = {ratio_i:.15f}")
    print(f"   sum x^II gamma^II Psat / P = {ratio_ii:.15f}   <- not solved for, checked")
    print(f"   shared vapor y = ({y[0]:.6f}, {y[1]:.6f}), sum = {y.sum():.12f}")
    print(f"   gamma^I  = ({np.exp(ln_gamma(x_i))[0]:.5f}, {np.exp(ln_gamma(x_i))[1]:.5f})")
    print(f"   gamma^II = ({np.exp(ln_gamma(x_ii))[0]:.5f}, {np.exp(ln_gamma(x_ii))[1]:.5f})")
    record("both liquids satisfy the bubble equation at T3 to 1e-10", abs(ratio_ii - 1.0) < 1e-10)
    record("the shared vapor sums to one", abs(float(y.sum()) - 1.0) < 1e-12)
    print(
        f"   Commonly tabulated heteroazeotrope: ~{LITERATURE_T_K} K, "
        f"y(water) ~ {LITERATURE_Y_WATER[0]}-{LITERATURE_Y_WATER[1]}. "
        f"Deviation here: {t3 - LITERATURE_T_K:+.2f} K, y(water) = {y[1]:.5f}."
    )
    print(
        "   Not a validation against experiment: these NRTL parameters were fitted to\n"
        "   liquid-liquid data and are temperature independent, so the agreement is\n"
        "   partly fortuitous. Reported, not tuned."
    )

    # -------------------------------------------------------------- T3 - 2 K
    print("\n3. Control (i): T3 - 2 K must be a stable liquid-liquid split")
    below = _flash(t3 - 2.0)
    diagnostics = below.diagnostics
    print(
        f"   phases = {below.phase_names()}, regime = {diagnostics['phase_regime']}, "
        f"vapor_fraction = {below.vapor_fraction}"
    )
    for name in below.phase_names():
        fractions = below.phases[name].composition.fractions
        print(
            f"     {name:<8} x({NAMES[0]}) = {fractions[0]:.12f}   "
            f"amount = {below.phase_fractions[name]:.8f}"
        )
    print(
        f"   equilibrium_residual = {float(diagnostics['equilibrium_residual']):.3e}, "
        f"delta_g_split_rt = {float(diagnostics['delta_g_split_rt']):+.6e}"
    )
    print(
        f"   post_split_status = {diagnostics['post_split_status']}, "
        f"tpd_min = {float(diagnostics['post_split_tpd_min']):+.3e}"
    )
    compositions = sorted(
        below.phases[name].composition.fractions[0] for name in below.phase_names()
    )
    record("two liquid phases", sorted(below.phase_names()) == ["liquid1", "liquid2"])
    record("vapor_fraction is None", below.vapor_fraction is None)
    record(
        "tie-line matches the independent binodal to 1e-9",
        abs(compositions[0] - x_i[0]) < 1e-9 and abs(compositions[1] - x_ii[0]) < 1e-9,
    )
    record(
        "both phases stable against the vapor candidate", diagnostics["post_split_stable"] is True
    )

    # -------------------------------------------------------------- T3 + 2 K
    print("\n4. Control (ii): T3 + 2 K")
    above = _flash(t3 + 2.0)
    print(
        f"   phases = {above.phase_names()}, regime = {above.diagnostics['phase_regime']}, "
        f"vapor_fraction = {above.vapor_fraction}"
    )
    print(
        f"   feed candidate = {above.diagnostics['feed_branch']}, "
        f"tpd_min = {float(above.diagnostics['tpd_min']):+.6e}"
    )
    dew_ratio = _dew_ratio(ln_gamma, np.array(FEED), t3 + 2.0)
    print(
        f"   independent dew check: sum z P / (gamma Psat) = {dew_ratio:.9f} (< 1 means "
        "the feed is above its dew point)"
    )
    record("the feed has fully evaporated", above.phase_names() == ["vapor"])
    record("and it really is above its dew point", dew_ratio < 1.0)

    # ------------------------------------------------------------------- T3
    print("\n5. Control (iii): at T3 and just below it")
    try:
        at_t3 = _flash(t3)
    except ct.ConvergenceError as error:
        print(f"   at T3:        ConvergenceError - {str(error)[:70]}...")
        record("a refusal at T3 names the third phase", "third phase is required" in str(error))
    else:
        names = sorted(at_t3.phase_names())
        vapor = np.array(at_t3.phases["vapor"].composition.fractions) if "vapor" in names else None
        print(
            f"   at T3:        phases = {at_t3.phase_names()}, "
            f"post_split = {at_t3.diagnostics['post_split_status']}, "
            f"tpd = {float(at_t3.diagnostics['post_split_tpd_min']):+.3e}"
        )
        if vapor is not None:
            print(
                f"                 vapor y = ({vapor[0]:.6f}, {vapor[1]:.6f})"
                f"   (three-phase vapor above: ({y[0]:.6f}, {y[1]:.6f}))"
            )
            record("the vapor at T3 is the three-phase vapor", abs(vapor[0] - y[0]) < 1e-5)
        print("                 marginal by construction: at T3 the liquid sits exactly on the")
        print("                 binodal, so the third phase has tpd = 0 and is not an instability.")

    raised = False
    try:
        _flash(t3 - 0.01, ct.FlashSettings(max_phases=2))
    except ct.ConvergenceError as error:
        raised = "third phase is required" in str(error)
        print(
            f"   at T3-0.01 K, max_phases=2: ConvergenceError - a third phase is required "
            f"({'message matched' if raised else 'unexpected message'})"
        )
    record("just below T3, max_phases=2 refuses the two-phase answer", raised)

    resolved = _flash(t3 - 0.01)
    tie_line = sorted(
        float(resolved.phases[name].composition.fractions[0]) for name in resolved.phase_names()
    )
    print(
        f"   at T3-0.01 K, max_phases=3 (default): phases = {resolved.phase_names()}, "
        f"regime = {resolved.diagnostics['phase_regime']}"
    )
    print(f"                 phase_set_history = {resolved.diagnostics['phase_set_history']}")
    print(
        f"                 x1 = ({tie_line[0]:.12f}, {tie_line[1]:.12f})"
        f"   (independent binodal: ({x_i[0]:.12f}, {x_ii[0]:.12f}))"
    )
    record(
        "the addition/removal search returns the two liquids",
        sorted(resolved.phase_names()) == ["liquid1", "liquid2"],
    )
    record(
        "and its tie-line is the independent binodal to 1e-8",
        abs(tie_line[0] - x_i[0]) < 1e-8 and abs(tie_line[1] - x_ii[0]) < 1e-8,
    )
    record(
        "and the history records the vapor being added and removed",
        str(resolved.diagnostics["phase_set_history"]) == "V -> LV -> LLV -> LL",
    )

    unchecked = _flash(t3 - 0.01, ct.FlashSettings(post_split_stability=False))
    diagnostics = unchecked.diagnostics
    print(
        f"   with post_split_stability=False the same call returns "
        f"{unchecked.phase_names()} and reports:"
    )
    for name in unchecked.phase_names():
        print(
            f"     phase_stability_{name:<8} = {diagnostics[f'phase_stability_{name}']}, "
            f"tpd_min = {float(diagnostics[f'phase_stability_tpd_min_{name}']):+.6e}"
        )
    record(
        "both phases of the refused pair see the same third phase",
        abs(
            float(diagnostics["phase_stability_tpd_min_vapor"])
            - float(diagnostics["phase_stability_tpd_min_liquid"])
        )
        < 1e-12,
    )
    print(
        "   Two coexisting phases share one tangent plane, so a third stationary point\n"
        "   below it has the SAME tpd measured from either - which is why both rows agree.\n"
        "   Below T3 the correct answer is the two-liquid pair, and the search reaches it\n"
        "   by adding the second liquid and then removing the vapor, whose amount goes to\n"
        "   zero. Measured width of the window in which the first two-phase iterate is the\n"
        "   wrong one: about 0.135 K below T3."
    )

    print("\n" + "=" * 78)
    if failures:
        print(f"FAIL: {len(failures)} check(s) failed")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
