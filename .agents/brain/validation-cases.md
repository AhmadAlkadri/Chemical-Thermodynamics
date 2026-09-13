# Validation Cases Ledger ("thermodynamics exam")

Purpose: one entry per case that a thermodynamics reviewer could re-run by hand.
Every entry states the source, where in the source, the assumptions, the
components and units, the parameters and their provenance, the expected outcome,
the tolerance, the independent route used, and the test path.

Rules:
- Never record an expected value that was not read from a source or produced by
  an independent route. No fabricated reference numbers.
- Record the tolerance actually achieved, not the tolerance asserted, when they
  differ.
- An entry whose "independent route" is an internal invariant must say so.

---

## Case S-1: Michelsen tangent-plane identity at the feed

- **Source:** M. L. Michelsen, "The isothermal flash problem. Part I.
  Stability", Fluid Phase Equilibria 9 (1982) 1-19.
- **Location:** Definition of the tangent plane distance; restated as equations
  (1)-(3) in the module docstring of `src/chemthermo/stability/tp.py`.
- **Assumptions:** Fixed T, P; single-phase reference tangent plane at the feed.
- **Components / units:** Methane / Ethane / Propane, z = 0.5 / 0.3 / 0.2,
  T = 240 K, P = 3 MPa. tpd and tm are dimensionless (units of RT).
- **Parameters and provenance:** Peng-Robinson with Tc, Pc, omega from the
  packaged chemthermo databank (`src/chemthermo/data/components.json`,
  Koretsky 2012); kij = 0.
- **Expected outcome:** tpd(w = z) = 0 and tm(W = z) = 0 exactly (analytic
  identity, not a fitted number).
- **Tolerance:** asserted 1e-12; achieved exactly 0.0 in double precision.
- **Independent route:** internal invariant (the test recomputes tpd and tm from
  the definitions rather than calling the solver).
- **Test path:** `tests/test_stability_tp.py::test_tpd_and_tm_vanish_at_the_feed_composition`

---

## Case S-2: Stationarity and the tm / tpd / sum(W) relations

- **Source:** Michelsen (1982), op. cit.; M. L. Michelsen and J. M. Mollerup,
  "Thermodynamic Models: Fundamentals and Computational Aspects", chapter on
  stability analysis.
- **Location:** Derivation of the modified tangent-plane function tm(W) and its
  stationarity condition; restated as equations (4)-(7) in
  `src/chemthermo/stability/tp.py`.
- **Assumptions:** Gibbs-Duhem holds at fixed T, P, so
  `sum_i w_i d ln phi_i = 0`; this is what collapses the gradient of tm.
- **Components / units:** as Case S-1.
- **Parameters and provenance:** as Case S-1.
- **Expected outcome:** at each converged non-trivial trial,
  `max_i |ln W_i + ln phi_i(w) - d_i| = 0`;
  the central finite-difference gradient of tm w.r.t. W is zero;
  `tm* = 1 - sum_i W_i`; `tpd = -ln(sum_i W_i)`; `tm* = 1 - exp(-tpd)`.
- **Tolerance:** asserted residual < 1e-9, |d tm / d W_k| < 1e-7 at step 1e-6,
  relations to 1e-10. Achieved: residual <= 1.8e-11, max |d tm / d W_k| = 6.7e-10,
  relations agree to ~2e-16.
- **Independent route:** internal invariant plus numerical differentiation of an
  independently coded tm.
- **Test path:** `tests/test_stability_tp.py::test_converged_non_trivial_trials_are_stationary_points_of_tm`

---

## Case S-3: Marginal stability of converged equilibrium phases

- **Source:** Michelsen (1982), op. cit. (equilibrium phases are points of
  tangency of a common hyperplane).
- **Location:** Tangent-plane criterion; see the test docstring for the
  argument.
- **Assumptions:** `flash_tp` has converged to a genuine equal-fugacity
  solution at 240 K, 3 MPa.
- **Components / units:** Methane / Ethane / Propane, feed z = 0.5 / 0.3 / 0.2;
  the stability feeds are the converged equilibrium liquid x and vapor y.
- **Parameters and provenance:** as Case S-1.
- **Expected outcome:** with x as feed, tpd_min >= 0 to numerical precision and
  the minimizing trial composition equals y; symmetrically with y as feed.
- **Tolerance:** asserted tpd_min >= -1e-6 and composition match atol 1e-6.
  Achieved: tpd_min = +5.41e-10 (feed = x) and -3.60e-11 (feed = y);
  composition match to ~2e-10.
- **Independent route:** cross-solver invariant (stability solver vs flash
  solver; they share no iteration code).
- **Test path:** `tests/test_stability_tp.py::test_equilibrium_phases_are_marginally_stable`

---

## Case S-4: Cross-check against `thermo` 0.6.0 Michelsen stability test

- **Source:** Caleb Bell, `thermo` 0.6.0 (open source). Functions
  `FlashVL.stability_test_Michelsen` and
  `thermo.flash.flash_utils.stability_iteration_Michelsen`.
- **Location:** `.venv/lib/python3.11/site-packages/thermo/flash/flash_vl.py`
  (~line 526) and `.../flash_utils.py` (~line 3917).
- **Assumptions:** The `thermo` PRMIX mixture is built from EXACTLY the Tc, Pc
  and omega values carried by chemthermo's `Component` objects (not thermo's own
  databank) with all kijs = 0, wrapped in `CEOSGas` / `CEOSLiquid` and
  `FlashVL`. Both codes use lowest-Gibbs fugacities. The reference stationary
  point is taken as the largest `sum(W)` over the same deterministic trial set.
- **Components / units:** seven states, SI (K, Pa), mole fractions:
  1. C1/C2/C3 0.5/0.3/0.2 @ 240 K, 3.0e6 Pa  -> unstable
  2. C1/C2 0.5/0.5 @ 450 K, 1.0e5 Pa         -> stable
  3. C1/C2/C3 0.5/0.3/0.2 @ 300 K, 5.0e7 Pa  -> stable
  4. C1/C2/C3 0.5/0.3/0.2 @ 200 K, 1.0e8 Pa  -> stable
  5. C1/C2/C3 0.5/0.3/0.2 @ 220 K, 2.0e6 Pa  -> unstable
  6. C1/C2 0.5/0.5 @ 200 K, 2.0e6 Pa         -> unstable
  7. C1/C3 0.7/0.3 @ 250 K, 5.0e6 Pa         -> unstable
- **Parameters and provenance:** Tc, Pc, omega from
  `src/chemthermo/data/components.json` (Koretsky 2012), passed verbatim to
  thermo. kij = 0 on both sides.
- **Expected outcome:** identical stable/unstable verdict on all seven states;
  matching `sum(W)`, tpd and stationary composition on the four unstable ones.
- **Tolerance:** asserted |d tpd| <= 1e-3 and max |d w| <= 1e-3.
  Achieved: verdicts identical on 7/7; max |d tpd| = 5.689e-5 (case 7);
  max |d w| = 2.300e-5 (case 7); |d sum(W)| <= 6.7e-5.
- **Known residual difference (explains the whole gap):** chemthermo uses the
  rounded Peng-Robinson constants 0.45724 / 0.07780, while `thermo` uses the
  exact roots 0.4572355289213822 / 0.0777960739038885. At these states that
  alone shifts ln(phi) by up to ~2.3e-4 and Z by ~4e-5, which dominates the
  differences above. This is an EOS-implementation difference, not a stability
  algorithm difference, and is out of scope for the `stability-tpd-pr` slice.
- **Independent route:** external reference implementation.
- **Test path:** `tests/validation/test_stability_vs_thermo.py`
- **Script:** `examples/validation/06_stability_vs_thermo.py`

---

## Case S-5 (NOT FOUND): published PR-based tangent-plane example

- **Status: not recorded. No usable open source was located.**
- Searched (September 2026) for a published Peng-Robinson tangent-plane
  stability example that states all of: component Tc / Pc / omega, the feed
  composition, T and P, kij = 0, and a numerical expected outcome (a TPD value
  or a stationary composition). Candidates considered:
  - Michelsen (1982) - uses SRK, so not directly usable for a PR check.
  - Hoteit and Firoozabadi (2006), Nichita et al., Li and Firoozabadi (2012) -
    the numerical tables were behind paywalls in the sources reachable from
    here, and the free secondary sources found gave mutually inconsistent
    accounts of the benchmark mixtures (for example two different compositions
    were reported for the "Y8" mixture).
- No expected value was invented. Until a primary source is read directly, the
  exam for this slice rests on Cases S-1 to S-4.
- **Follow-up:** if the original Hoteit and Firoozabadi (2006) or Nichita et al.
  tables become available, add the case here with the exact table values and the
  achieved tolerance.
