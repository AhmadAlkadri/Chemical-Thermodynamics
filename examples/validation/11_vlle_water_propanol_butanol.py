"""The vapor-liquid-liquid tie-triangle of 1-propanol / n-butanol / water.

What is being validated
-----------------------
A ternary at fixed ``(T, P)`` can hold three phases over a whole *region* of
feed compositions, not just at one point: Gibbs' phase rule gives
``F = 3 - 3 + 2 = 2``, and fixing T and P still leaves the region - the
tie-triangle - whose three corners are the two conjugate liquids and the vapor
they share. Every feed inside that triangle returns those same three
compositions and differs only in how much of each there is.

This script checks that `flash_tp` finds the triangle without being told it is
there, and that it is the *right* triangle:

1. the three vertices are reproduced by an **independent six-equation Newton
   solve** written in this file, which shares no code with `chemthermo.flash`;
2. feeds inside the triangle (chosen by barycentric weights, so "inside" is a
   fact about the reference triangle, not about the flash) return three phases
   whose compositions and phase fractions match that reference;
3. the answer is the one with the lowest Gibbs energy:
   ``G(3 phases) < G(2-phase candidate) < G(1 phase)``, all three computed here;
4. feeds outside the triangle come back as two phases or one, each verdict
   checked by its own route (a stability test, a bubble-point sum, a
   dew-point sum);
5. the binary water / n-butanol "refusal window" of validation Case R-3 - the
   ~0.135 K band below the three-phase temperature where the first two-phase
   iterate is the wrong one - now resolves to the two-liquid pair by *removing*
   a phase, and still raises under ``FlashSettings(max_phases=2)``.

The independent tie-triangle equations
--------------------------------------
Six unknowns (two liquid compositions) and six equations:

    ln x_i^I + ln gamma_i(x^I) = ln x_i^II + ln gamma_i(x^II)          (3 eq)
    sum_i x_i^I = 1,     sum_i x_i^II = 1                              (2 eq)
    sum_i y_i = 1,       y_i = x_i^I gamma_i(x^I) Psat_i(T) / P        (1 eq)

Only liquid I is required to be at its bubble point. That liquid II is
simultaneously at *its* bubble point, and that both give the same vapor, are
consequences and are checked - which is the non-trivial part, because the
equal-activity equations know nothing about Antoine.

Parameters and their provenance
-------------------------------
NRTL Table 1 of S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, Chem. Eng.
Sci. 55 (2000) 1785-1796 (attributed there to McDonald and Floudas, AIChE J. 41
(1995) 1798); ``alpha`` implied by ``G_ij = exp(-alpha_ij tau_ij)``. Antoine
coefficients from the packaged databank (Koretsky 2012),
``ln(P^sat/bar) = A - B/(T + C)``.

**These parameters were fitted to liquid-liquid data and are temperature
independent, and no experimental ternary VLLE data is used anywhere in this
script.** Nothing printed is a comparison against measurement. `thermo` 0.6.0
cannot split two liquids over one excess-Gibbs model (validation Case L-2), so
there is no external three-phase reference for this system either, and none is
claimed.

Needs no optional dependency. All inputs are SI.
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np

import chemthermo as ct

NAMES = ("1-Propanol", "n-Butanol", "Water")
PRESSURE_PA = 101325.0

# Tessier et al. (2000) Table 1, Problem 1. tau[i][j] is tau_ij.
TAU = np.array(
    [
        [0.0, -0.61259, -0.07149],
        [0.7164, 0.0, 0.90047],
        [2.7425, 3.51307, 0.0],
    ]
)
ALPHA = np.array(
    [
        [0.0, 0.3, 0.3],
        [0.3, 0.0, 0.48],
        [0.3, 0.48, 0.0],
    ]
)
CAPITAL_G = np.exp(-ALPHA * TAU)

#: Newton starting points for the tie-triangle, one per temperature.
SEEDS = {
    365.0: (0.023616, 0.024407, 0.951977, 0.098737, 0.201955, 0.699308),
    364.0: (0.052144, 0.029172, 0.918684, 0.139639, 0.123478, 0.736883),
    363.0: (0.102828, 0.035390, 0.861782, 0.156303, 0.064227, 0.779470),
}
#: Barycentric weights (all positive) of the feeds tested inside each triangle.
INSIDE_WEIGHTS = ((1 / 3, 1 / 3, 1 / 3), (0.2, 0.3, 0.5))

# Binary control: n-Butanol / Water, validation Case R-3.
BINARY_NAMES = ("n-Butanol", "Water")
BINARY_FEED = (0.20, 0.80)
BINARY_T3_K = 366.213774
BINARY_BINODAL_X1 = (0.019998419467, 0.359999661508)

_FAILURES: list[str] = []


def record(label: str, passed: bool) -> None:
    print(f"  [{'PASS' if passed else 'FAIL'}] {label}")
    if not passed:
        _FAILURES.append(label)


# --------------------------------------------------------------------------
# The model, written out here
# --------------------------------------------------------------------------


def ln_gamma(x: np.ndarray) -> np.ndarray:
    """Renon-Prausnitz NRTL with column sums."""
    values = np.asarray(x, dtype=float)
    values = values / float(np.sum(values))
    s = CAPITAL_G.T @ values
    c = (TAU * CAPITAL_G).T @ values
    return np.array(
        [
            c[i] / s[i] + float(np.sum(values * CAPITAL_G[i, :] / s * (TAU[i, :] - c / s)))
            for i in range(len(NAMES))
        ]
    )


def _antoine(names: Sequence[str]) -> Callable[[float], np.ndarray]:
    antoine = []
    for name in names:
        record = ct.Component.from_database(name).antoine
        assert record is not None, name
        antoine.append(record)

    def psat(temperature_K: float) -> np.ndarray:
        return np.array([np.exp(a.A - a.B / (temperature_K + a.C)) * 1.0e5 for a in antoine])

    return psat


PSAT = _antoine(NAMES)


def _model() -> ct.NRTL:
    pairs = [
        (
            NAMES[i],
            NAMES[j],
            float(TAU[i][j]),
            float(TAU[j][i]),
            float(ALPHA[i][j]),
            float(ALPHA[j][i]),
        )
        for i in range(len(NAMES))
        for j in range(i + 1, len(NAMES))
    ]
    return ct.NRTL(parameters=ct.NRTLParameters.from_pairs(pairs))


# --------------------------------------------------------------------------
# The independent tie-triangle
# --------------------------------------------------------------------------


def _residual(u: np.ndarray, temperature_K: float) -> np.ndarray:
    x_i, x_ii = u[:3], u[3:]
    activity_i = np.log(np.abs(x_i)) + ln_gamma(x_i)
    activity_ii = np.log(np.abs(x_ii)) + ln_gamma(x_ii)
    y = np.exp(activity_i) * PSAT(temperature_K) / PRESSURE_PA
    return np.concatenate(
        [
            activity_i - activity_ii,
            [float(np.sum(x_i)) - 1.0, float(np.sum(x_ii)) - 1.0, float(np.sum(y)) - 1.0],
        ]
    )


def tie_triangle(temperature_K: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Damped Newton with a finite-difference Jacobian on the six equations."""
    u = np.array(SEEDS[temperature_K], dtype=float)
    residual = float("inf")
    for _ in range(200):
        f = _residual(u, temperature_K)
        residual = float(np.max(np.abs(f)))
        if residual < 1e-14:
            break
        jacobian = np.zeros((6, 6))
        step = 1e-8
        for column in range(6):
            forward = u.copy()
            forward[column] += step
            jacobian[:, column] = (_residual(forward, temperature_K) - f) / step
        direction = np.linalg.solve(jacobian, -f)
        scale = 1.0
        while np.any(u + scale * direction <= 0.0):
            scale *= 0.5
        u = u + scale * direction
    x_i, x_ii = u[:3], u[3:]
    y = np.exp(np.log(x_i) + ln_gamma(x_i)) * PSAT(temperature_K) / PRESSURE_PA
    return x_i, x_ii, y, residual


def reduced_g(x: np.ndarray, is_vapor: bool, temperature_K: float) -> float:
    """``sum_i x_i (ln x_i + t_i)``: the part of ``G/RT`` the split changes."""
    values = np.asarray(x, dtype=float)
    terms = (
        np.zeros_like(values)
        if is_vapor
        else ln_gamma(values) + np.log(PSAT(temperature_K) / PRESSURE_PA)
    )
    mask = values > 0.0
    return float(np.sum(values[mask] * (np.log(values[mask]) + terms[mask])))


def flash(
    names: Sequence[str],
    model: ct.NRTL,
    z: Sequence[float],
    temperature_K: float,
    settings: ct.FlashSettings | None = None,
) -> ct.FlashResult:
    return ct.flash_tp(
        ct.Mixture.from_database(list(names), [float(v) for v in z], normalize=True),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
        flash_mode="modified-raoult",
        settings=settings,
    )


def result_energy(result: ct.FlashResult, temperature_K: float) -> float:
    return float(
        sum(
            result.phase_fractions[name]
            * reduced_g(
                np.array(result.phases[name].composition.fractions), name == "vapor", temperature_K
            )
            for name in result.phase_names()
        )
    )


# --------------------------------------------------------------------------
# Sections
# --------------------------------------------------------------------------


def section_triangle(model: ct.NRTL) -> None:
    print("\n1. The tie-triangle from an independent Newton solve, and what flash_tp finds")
    for temperature_K in sorted(SEEDS, reverse=True):
        x_i, x_ii, y, residual = tie_triangle(temperature_K)
        print(f"\n   T = {temperature_K:.2f} K   (Newton residual {residual:.2e})")
        for label, phase in (("x^I ", x_i), ("x^II", x_ii), ("y   ", y)):
            print(f"     {label} = ({phase[0]:.8f}, {phase[1]:.8f}, {phase[2]:.8f})")
        bubble_ii = float(np.sum(x_ii * np.exp(ln_gamma(x_ii)) * PSAT(temperature_K) / PRESSURE_PA))
        y_from_ii = x_ii * np.exp(ln_gamma(x_ii)) * PSAT(temperature_K) / PRESSURE_PA
        print(f"     liquid II bubble sum = {bubble_ii:.12f}  (not solved for: a consequence)")
        record(
            f"T={temperature_K:.0f} K: both liquids boil at once",
            abs(bubble_ii - 1.0) < 1e-10,
        )
        record(
            f"T={temperature_K:.0f} K: they share one vapor",
            float(np.max(np.abs(y_from_ii - y))) < 1e-10,
        )

        vertices = np.column_stack([x_i, x_ii, y])
        for weights in INSIDE_WEIGHTS:
            tag = f"T={temperature_K:.0f} K, beta={tuple(round(w, 3) for w in weights)}"
            z = vertices @ np.array(weights)
            z = z / float(np.sum(z))
            exact = np.linalg.solve(vertices, z)
            result = flash(NAMES, model, z, temperature_K)
            names = sorted(result.phase_names())
            computed = {name: np.array(result.phases[name].composition.fractions) for name in names}
            worst_x = 0.0
            worst_beta = 0.0
            for index, reference in enumerate((x_i, x_ii, y)):
                nearest = min(
                    computed, key=lambda n: float(np.max(np.abs(computed[n] - reference)))
                )
                worst_x = max(worst_x, float(np.max(np.abs(computed[nearest] - reference))))
                worst_beta = max(
                    worst_beta, abs(result.phase_fractions[nearest] - float(exact[index]))
                )
            print(
                f"     feed at barycentric {tuple(round(w, 4) for w in weights)}: "
                f"{names}, history = {result.diagnostics['phase_set_history']}"
            )
            print(
                f"       worst |dx| = {worst_x:.2e}, worst |dbeta| = {worst_beta:.2e}, "
                f"equilibrium residual = {float(result.diagnostics['equilibrium_residual']):.2e}, "
                f"mass balance = {float(result.diagnostics['mass_balance_residual']):.2e}"
            )
            record(f"{tag}: three phases", len(names) == 3)
            record(f"{tag}: compositions to 1e-6", worst_x < 1e-6)
            record(f"{tag}: fractions to 1e-6", worst_beta < 1e-6)
            record(
                f"{tag}: equilibrium residual < 1e-10",
                float(result.diagnostics["equilibrium_residual"]) < 1e-10,
            )
            record(
                f"{tag}: mass balance < 1e-12",
                float(result.diagnostics["mass_balance_residual"]) < 1e-12,
            )
            record(
                f"{tag}: every phase post-split stable",
                bool(result.diagnostics["post_split_stable"]),
            )

            two = flash(
                NAMES, model, z, temperature_K, ct.FlashSettings(post_split_stability=False)
            )
            g3 = result_energy(result, temperature_K)
            g2 = result_energy(two, temperature_K)
            g1 = reduced_g(z, False, temperature_K)
            print(f"       G3/RT = {g3:.9f} < G2/RT = {g2:.9f} < G1/RT = {g1:.9f}")
            record(f"{tag}: G3 < G2 < G1", g3 < g2 < g1)


def section_negative_controls(model: ct.NRTL) -> None:
    print("\n2. Negative controls at 364.0 K: feeds outside the tie-triangle")
    temperature_K = 364.0
    x_i, x_ii, y, _ = tie_triangle(temperature_K)
    vertices = np.column_stack([x_i, x_ii, y])

    # (a) Vapor-liquid region: past the x^II - y edge, away from x^I.
    z = vertices @ np.array([-0.25, 0.45, 0.80])
    z = z / float(np.sum(z))
    result = flash(NAMES, model, z, temperature_K)
    liquid = np.array(result.phases["liquid"].composition.fractions)
    liquid_stability = ct.stability_tp(
        ct.Mixture.from_database(list(NAMES), list(liquid), normalize=True),
        temperature_K=temperature_K,
        pressure_Pa=PRESSURE_PA,
        activity_model=model,
    )
    print(f"\n   (a) VL region, z = ({z[0]:.6f}, {z[1]:.6f}, {z[2]:.6f})")
    print(
        f"       barycentric weights in the triangle: {np.round(np.linalg.solve(vertices, z), 6)}"
    )
    print(f"       phases = {result.phase_names()}, vapor_fraction = {result.vapor_fraction:.8f}")
    print(
        f"       the liquid tested on its own against a second liquid: "
        f"{liquid_stability.status}, tpd_min = {float(liquid_stability.tpd_min):+.3e}"
    )
    record(
        "(a) two phases, vapor and one liquid", sorted(result.phase_names()) == ["liquid", "vapor"]
    )
    record("(a) the liquid is stable against a second liquid", liquid_stability.status == "stable")

    # (b) Liquid-liquid region: past the x^I - x^II edge, away from y.
    z = vertices @ np.array([0.55, 0.60, -0.15])
    z = z / float(np.sum(z))
    result = flash(NAMES, model, z, temperature_K)
    bubbles = [
        float(
            np.sum(
                np.array(result.phases[name].composition.fractions)
                * np.exp(ln_gamma(np.array(result.phases[name].composition.fractions)))
                * PSAT(temperature_K)
                / PRESSURE_PA
            )
        )
        for name in result.phase_names()
    ]
    print(f"\n   (b) LL region, z = ({z[0]:.6f}, {z[1]:.6f}, {z[2]:.6f})")
    print(
        f"       barycentric weights in the triangle: {np.round(np.linalg.solve(vertices, z), 6)}"
    )
    print(f"       phases = {result.phase_names()}, vapor_fraction = {result.vapor_fraction}")
    print(
        "       both liquids below their bubble point: "
        + ", ".join(f"sum x gamma Psat / P = {value:.9f}" for value in bubbles)
    )
    for name in result.phase_names():
        print(
            f"       phase_stability_{name:<8} = {result.diagnostics[f'phase_stability_{name}']}, "
            f"tpd_min = {float(result.diagnostics[f'phase_stability_tpd_min_{name}']):+.3e}"
        )
    record("(b) two liquids", sorted(result.phase_names()) == ["liquid1", "liquid2"])
    record("(b) vapor_fraction is None", result.vapor_fraction is None)
    record("(b) neither liquid can boil", all(value < 1.0 for value in bubbles))

    # (c) Water-rich corner: a single liquid.
    z = np.array([0.01, 0.005, 0.985])
    result = flash(NAMES, model, z, temperature_K)
    bubble = float(np.sum(z * np.exp(ln_gamma(z)) * PSAT(temperature_K) / PRESSURE_PA))
    print(f"\n   (c) water-rich corner, z = ({z[0]:.6f}, {z[1]:.6f}, {z[2]:.6f})")
    print(
        f"       phases = {result.phase_names()}, feed candidate = "
        f"{result.diagnostics['feed_branch']}, sum z gamma Psat / P = {bubble:.9f}"
    )
    record("(c) a single liquid", result.phase_names() == ["liquid"])
    record("(c) and it is below its bubble point", bubble < 1.0)

    # (d) Superheated: a single vapor at 380 K.
    z = np.array([0.20, 0.15, 0.65])
    result = flash(NAMES, model, z, 380.0)
    x = z.copy()
    for _ in range(500):
        x = z * PRESSURE_PA / (np.exp(ln_gamma(x)) * PSAT(380.0))
        x = x / float(np.sum(x))
    ratio = float(np.sum(z * PRESSURE_PA / (np.exp(ln_gamma(x)) * PSAT(380.0))))
    print(f"\n   (d) superheated at 380.0 K, z = ({z[0]:.6f}, {z[1]:.6f}, {z[2]:.6f})")
    print(
        f"       phases = {result.phase_names()}, feed candidate = "
        f"{result.diagnostics['feed_branch']}, dew sum z P / (gamma Psat) = {ratio:.9f}"
    )
    record("(d) a single vapor", result.phase_names() == ["vapor"])
    record("(d) and it is above its dew point", ratio < 1.0)


def section_refusal_window() -> None:
    print("\n3. The binary refusal window (validation Case R-3 -> V-3)")
    model = ct.NRTL(
        parameters=ct.NRTLParameters.from_pairs(
            [(BINARY_NAMES[0], BINARY_NAMES[1], 0.90047, 3.51307, 0.48, 0.48)]
        )
    )
    print(
        f"   n-Butanol / Water, z = {BINARY_FEED}, T3 = {BINARY_T3_K:.6f} K.\n"
        "   A binary can hold three phases only at T3 (F = 2 - 3 + 2 = 1), so just\n"
        "   below it the three-phase Rachford-Rice has no finite solution and the\n"
        "   vapor's amount runs negative - which is what removes it."
    )
    for offset in (-0.05, -0.10):
        result = flash(BINARY_NAMES, model, BINARY_FEED, BINARY_T3_K + offset)
        tie_line = sorted(
            float(result.phases[name].composition.fractions[0]) for name in result.phase_names()
        )
        deviation = max(
            abs(tie_line[0] - BINARY_BINODAL_X1[0]), abs(tie_line[1] - BINARY_BINODAL_X1[1])
        )
        print(
            f"\n   T3 {offset:+.2f} K: {result.phase_names()}, "
            f"history = {result.diagnostics['phase_set_history']}"
        )
        print(
            f"       x1 = ({tie_line[0]:.12f}, {tie_line[1]:.12f}), "
            f"worst deviation from the binodal = {deviation:.2e}"
        )
        record(
            f"T3{offset:+.2f} K: two liquids",
            sorted(result.phase_names()) == ["liquid1", "liquid2"],
        )
        record(f"T3{offset:+.2f} K: binodal to 1e-8", deviation < 1e-8)
        record(
            f"T3{offset:+.2f} K: the history shows the add then the removal",
            str(result.diagnostics["phase_set_history"]) == "V -> LV -> LLV -> LL",
        )

        raised = False
        try:
            flash(
                BINARY_NAMES,
                model,
                BINARY_FEED,
                BINARY_T3_K + offset,
                ct.FlashSettings(max_phases=2),
            )
        except ct.ConvergenceError as error:
            raised = "third phase is required" in str(error)
        record(f"T3{offset:+.2f} K: max_phases=2 still raises the documented error", raised)

    above = flash(BINARY_NAMES, model, BINARY_FEED, BINARY_T3_K + 0.05)
    binary_psat = _antoine(BINARY_NAMES)(BINARY_T3_K + 0.05)

    def binary_ln_gamma(x: np.ndarray) -> np.ndarray:
        mixture = ct.Mixture.from_database(list(BINARY_NAMES), [0.5, 0.5], normalize=True)
        return np.log(
            np.array(
                model.activity_coefficients(
                    mixture=mixture,
                    temperature_K=298.15,
                    composition=[float(v) for v in np.asarray(x) / float(np.sum(x))],
                )
            )
        )

    x = np.array(above.phases["liquid"].composition.fractions)
    y = np.array(above.phases["vapor"].composition.fractions)
    bubble = float(np.sum(x * np.exp(binary_ln_gamma(x)) * binary_psat / PRESSURE_PA))
    raoult = float(
        np.max(np.abs(y * PRESSURE_PA - x * np.exp(binary_ln_gamma(x)) * binary_psat)) / PRESSURE_PA
    )
    print(f"\n   T3 +0.05 K: {above.phase_names()}, vapor_fraction = {above.vapor_fraction:.8f}")
    print(
        f"       liquid at its bubble point: sum x gamma Psat / P = {bubble:.12f}; "
        f"modified Raoult residual = {raoult:.2e}"
    )
    record("T3+0.05 K: a vapor-liquid pair", sorted(above.phase_names()) == ["liquid", "vapor"])
    record("T3+0.05 K: the liquid is exactly at its bubble point", abs(bubble - 1.0) < 1e-10)
    record("T3+0.05 K: modified Raoult's law holds", raoult < 1e-10)


def main() -> None:
    print("=" * 78)
    print("Ternary VLLE: 1-propanol / n-butanol / water at 101325 Pa")
    print("Modified Raoult (NRTL liquid + Antoine reference fugacity, ideal vapor)")
    print("Parameters: Tessier et al. (2000) Table 1; LLE-fitted, temperature independent.")
    print("No experimental data is used; nothing here is a comparison against measurement.")
    print("=" * 78)

    model = _model()
    section_triangle(model)
    section_negative_controls(model)
    section_refusal_window()

    print("\n" + "=" * 78)
    if _FAILURES:
        print(f"FAIL: {len(_FAILURES)} check(s) failed")
        for failure in _FAILURES:
            print(f"  - {failure}")
        raise SystemExit(1)
    print("PASS: every check above passed")


if __name__ == "__main__":
    main()
