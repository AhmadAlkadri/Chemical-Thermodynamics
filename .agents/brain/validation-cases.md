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

---

## Case K-1: Per-pair `kij` cross-check against `thermo` 0.6.0 PRMIX

- **Source:** Caleb Bell, `thermo` 0.6.0 (open source). `thermo.PRMIX` /
  `CEOSGas` / `CEOSLiquid` / `FlashVL` / `FlashVL.stability_test_Michelsen`.
- **Location:** `tests/validation/test_pr_kij_vs_thermo.py`.
- **Assumptions:** The `thermo` PRMIX mixture is built from EXACTLY the Tc,
  Pc, omega values carried by chemthermo's `Component` objects (not thermo's
  own databank), with the SAME nonzero kij matrix passed to both codes, so
  the comparison isolates the `pr-kij-matrix` mixing-rule fix rather than
  component-data or kij-value differences.
- **Components / units:** SI (K, Pa), mole fractions.
  - Binary: Methane / n-Decane, kij(Methane, n-Decane) = 0.0411. States
    (320 K, 2.0 MPa) and (380 K, 4.0 MPa); compositions
    (0.3, 0.7) / (0.5, 0.5) / (0.7, 0.3); both vapor and liquid roots.
  - Three-component (synthetic matrix, exercises the general n x n path):
    Methane / Ethane / n-Decane with kij(Methane, Ethane) = -0.0026,
    kij(Methane, n-Decane) = 0.0411, kij(Ethane, n-Decane) = 0.0170. Same two
    states; compositions (0.5, 0.3, 0.2) / (0.3, 0.3, 0.4) / (0.2, 0.5, 0.3).
  - End-to-end: Methane / n-Decane, z = 0.5 / 0.5, kij = 0.0411,
    T = 350 K, P = 3.0 MPa (independently confirmed two-phase on both sides).
- **Parameters and provenance:** Tc, Pc, omega from
  `src/chemthermo/data/components.json` (Koretsky 2012), passed verbatim to
  `thermo`. The kij value 0.0411 for Methane/n-Decane is an **illustrative,
  literature-order-of-magnitude value for a light-heavy alkane pair, not
  sourced from a specific publication read in this session** -- it must not
  be read as a validated physical parameter. The 3x3 matrix's other two
  entries are entirely synthetic (chosen only to give three distinct
  off-diagonal values and exercise the matrix-building code); they carry no
  physical claim at all.
- **Expected outcome:** ln(phi) and Z agree with `thermo` up to the known
  rounded-PR-constants gap (see below); `flash_tp` and `stability_tp` verdicts
  and quantities agree with `thermo`'s `FlashVL` / `stability_test_Michelsen`
  using the same kij.
- **Tolerance:** asserted max \|d ln phi\| <= 1e-3, max \|dZ\| <= 2e-4 (phi/Z
  cross-checks); asserted \|d beta\| <= 5e-4, max \|dx\|/\|dy\| <= 2e-3
  (end-to-end flash). Achieved:
  - Binary (2 states x 3 compositions x 2 roots): max \|d ln phi\| = 5.783e-4,
    max \|dZ\| = 6.468e-5.
  - Three-component: max \|d ln phi\| = 5.823e-4, max \|dZ\| = 1.833e-5.
  - Pure-component limit within a kij != 0 binary (y = [1, 0]):
    \|d ln phi\| well under 1e-3 -- this is exactly the check the pre-fix
    diagonal bug would have failed (see "Known residual difference" below).
  - End-to-end flash (T=350 K, P=3.0 MPa): \|d beta\| = 2.94e-6,
    max \|dx\| = 5.27e-6, max \|dy\| = 8.32e-7.
  - End-to-end stability: verdict "unstable" on both sides (matches).
- **Known residual difference (explains the phi/Z gap):** as in Case S-4,
  chemthermo uses the rounded Peng-Robinson constants 0.45724 / 0.07780,
  while `thermo` uses the exact roots 0.4572355289213822 / 0.0777960739038885.
  This is an EOS-implementation difference, not a kij-handling difference,
  and ADR-0006 keeps the rounded constants unchanged.
- **Pre-fix vs post-fix evidence (not asserted in a persisted test; measured
  during development and reported alongside this slice):** at T=350 K,
  P=2.0 MPa, Methane/n-Decane z=0.5/0.5, kij=0.0411, the pre-fix formula
  (`aij = sqrt(outer(a_i,a_i)) * (1 - kij)` applied to the diagonal too) gave
  max \|d ln phi\| = 0.365 against `thermo` (single-real-root state, vapor and
  liquid branches identical); the post-fix formula gives max
  \|d ln phi\| = 3.40e-4 at the same state -- roughly a 1000x reduction,
  and the remainder is explained entirely by the rounded-constants gap above.
- **Permutation invariance:** a positional-matrix alternative was rejected in
  ADR-0006 partly because it would not be permutation-invariant. Verified
  directly (no `thermo` dependency, since `thermo`'s own CEOSLiquid root
  solver was observed to be numerically order-sensitive at some
  near-critical-locus states tried during development, which is a property of
  the reference solver, not of chemthermo):
  `tests/test_pr_eos.py::test_pr_eos_kij_permutation_invariance` and
  `tests/validation/test_pr_kij_vs_thermo.py::test_pr_eos_kij_permutation_invariance_with_matrix_kij`.
- **Independent route:** external reference implementation.
- **Test path:** `tests/validation/test_pr_kij_vs_thermo.py`;
  regression/unit coverage in `tests/test_pr_eos.py`.
- **Script:** `examples/basic/tp_flash_pr_kij_demo.py` (golden path; contrasts
  kij=0.0 against the mapping form on the same feed, no `thermo` dependency).

---

## Case N-1: NRTL Gibbs-Duhem consistency and binary reduction

- **Source:** H. Renon and J. M. Prausnitz, "Local compositions in
  thermodynamic excess functions for liquid mixtures", AIChE Journal 14 (1968)
  135-144 (the NRTL equation and its two-component form); the Gibbs-Duhem
  relation at fixed T, P.
- **Location:** Working equation restated as equations (1)-(4) in the module
  docstring of `src/chemthermo/models/nrtl.py`; the binary expressions are
  written out in the test docstring.
- **Assumptions:** Fixed T, P. `ln gamma` must be the composition derivative of
  a single reduced excess Gibbs energy `g^E`, hence
  `sum_i x_i d ln gamma_i = 0` along any direction in the composition simplex.
- **Components / units:** 1-Propanol / n-Butanol / Water (dimensionless
  compositions and activity coefficients). Four interior compositions
  (0.12, 0.08, 0.80), (0.50, 0.30, 0.20), (0.20, 0.20, 0.60),
  (0.05, 0.05, 0.90) and four simplex directions (1,-1,0), (1,0,-1), (0,1,-1),
  (1,1,-2). Binary reduction uses a synthetic asymmetric pair
  (tau_12 = 0.6, tau_21 = -0.35, alpha = 0.3) at six compositions.
- **Parameters and provenance:** Tessier (2000) Table 1 via
  `tests/fixtures/nrtl/tessier2000_problem1.json` (see Case N-3 for
  provenance). The binary-reduction parameters are synthetic and carry no
  physical claim; only the algebraic identity is under test.
- **Expected outcome:** Gibbs-Duhem residual 0; multicomponent code equals the
  textbook two-component NRTL formulas exactly; reordering components permutes
  `ln gamma` exactly; `tau = 0` gives `gamma = 1`.
- **Tolerance:** asserted |residual| < 1e-7 (and < 1e-8 for the worst case),
  binary abs 1e-14, permutation rtol 1e-12. Achieved: worst Gibbs-Duhem
  residual 2.07e-10 at step 1e-6; binary agreement 1.1e-16; permutation
  agreement at round-off.
- **Pre-fix vs post-fix evidence (measured during the
  `nrtl-gibbs-duhem-fix` slice, not persisted as a test):** the previous
  implementation summed `G` along rows and used per-term denominators in the
  first term. Gibbs-Duhem residuals with the same parameters were
  +1.2549e-01 at (0.12, 0.08, 0.80) along (1,-1,0), +5.8681e-01 along
  (1,0,-1), and -4.5359e-02 at (0.50, 0.30, 0.20) along (1,-1,0). `ln gamma`
  differed from the standard equation by up to 0.113 at (0.12, 0.08, 0.80),
  0.649 at (0.50, 0.30, 0.20) and 0.889 at (0.0597, 0.0282, 0.9120).
  A regression guard pins `ln gamma` at (0.12, 0.08, 0.80) to
  (0.9120183964, 1.1593117705, 0.1849764954) within 1e-9; the pre-fix code
  returned (0.857845, 1.045984, 0.174959).
- **Independent route:** internal invariant (Gibbs-Duhem by central
  differences) plus an independently written two-component formula. External
  cross-check is Case N-2.
- **Test path:** `tests/test_activity_nrtl.py`
  (`test_nrtl_satisfies_gibbs_duhem_for_asymmetric_parameters`,
  `test_nrtl_binary_reduces_to_textbook_two_component_formulas`,
  `test_nrtl_is_permutation_invariant`,
  `test_nrtl_regression_guard_against_row_sum_bug`,
  `test_nrtl_gamma_unity_with_zero_tau`,
  `test_nrtl_fully_symmetric_ternary_is_permutation_symmetric`).

---

## Case N-2: NRTL against an independent implementation (`thermo` 0.6.0)

- **Source:** `thermo` 0.6.0 (`thermo.NRTL` and `thermo.nrtl.NRTL_gammas`) as
  an external implementation of the same closed-form Renon-Prausnitz equation.
- **Location:** `tests/validation/test_nrtl_tessier2000.py::test_nrtl_matches_thermo_for_asymmetric_parameters`.
- **Assumptions:** Both sides evaluate the same closed form at the same
  dimensionless tau and alpha, so agreement is expected at round-off, not at a
  "physically reasonable" tolerance. `thermo` receives the parameters as
  temperature-independent coefficients (`tau_coeffs` / `alpha_coeffs` with only
  the constant term non-zero).
- **Components / units:** 1-Propanol / n-Butanol / Water at six compositions:
  (0.12, 0.08, 0.80), (0.50, 0.30, 0.20), (0.0597449, 0.0282358, 0.9120193),
  (0.20, 0.20, 0.60), (0.80, 0.10, 0.10), (0.05, 0.90, 0.05). `ln gamma` is
  dimensionless.
- **Parameters and provenance:** as Case N-3 (Tessier 2000 Table 1 fixture).
  Strongly asymmetric: tau_12 = -0.61259 vs tau_21 = 0.71640, and
  alpha_23 = 0.48 vs alpha_12 = alpha_13 = 0.3.
- **Expected outcome:** identical `ln gamma` from both implementations.
- **Tolerance:** asserted max |d ln gamma| < 1e-9 per composition and < 1e-12
  overall. Achieved: 8.88e-16 against both `thermo.NRTL` and
  `thermo.nrtl.NRTL_gammas`.
- **Known weakness this replaces:** the pre-existing cross-check
  (`tests/validation/test_flash_vs_thermo.py::test_nrtl_activity_coefficients_vs_thermo`)
  used a 50/50 Methane/Ethane binary with tau = 0.2/0.1 and symmetric
  alpha = 0.3 at rtol/atol 2e-3. The row-sum bug moved `ln gamma` by only
  ~1.1e-3 there, so that test passed both before and after the fix. A
  symmetric-`G` binary cannot detect a row/column mix-up.
- **Independent route:** external reference implementation.
- **Test path:** `tests/validation/test_nrtl_tessier2000.py`.

---

## Case N-3: Tessier, Brennecke & Stadtherr (2000) Table 2 stationary points

- **Source:** S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable
  phase stability analysis for excess Gibbs energy models", Chemical
  Engineering Science 55 (2000) 1785-1796. Author copy of the accepted
  manuscript: https://academicweb.nd.edu/~markst/srt2000.pdf
- **Location:** Section 2.2 (NRTL form), Table 1 (Problem 1 parameters),
  Table 2 (stationary points of the tangent-plane distance D). Table 1
  attributes the parameters to C. M. McDonald and C. A. Floudas, AIChE Journal
  41 (1995) 1798-1814.
- **Assumptions:** The paper's `g^E` is the standard Renon-Prausnitz form
  (verified against section 2.2 equations 5-7), so `D` is the activity-based
  tangent-plane distance
  `D(x) = sum_i x_i [ln x_i + ln gamma_i(x) - ln z_i - ln gamma_i(z)]`.
  tau is dimensionless and temperature-independent as printed, so no
  temperature is needed (the paper states none for Problem 1).
- **Components / units:** n-propanol(1) / n-butanol(2) / water(3), mapped to
  the chemthermo databank names `1-Propanol`, `n-Butanol`, `Water`. Four feeds:
  (0.148, 0.052, 0.80), (0.12, 0.08, 0.80), (0.13, 0.07, 0.80),
  (0.12, 0.05, 0.83). D is dimensionless.
- **Parameters and provenance:** `tests/fixtures/nrtl/tessier2000_problem1.json`
  (deliberately NOT in the packaged default data). Table 1 prints `G_ij` and
  `tau_ij` but not `alpha_ij`; alpha is implied by
  `alpha_ij = -ln(G_ij)/tau_ij`, which evaluates to 0.3 for pairs 1-2 and 1-3
  and 0.48 for pair 2-3 (recovered values 0.29999945 to 0.30000003 and
  0.48000005). Re-exponentiating the rounded alpha reproduces the printed G to
  max |dG| = 3.99e-08 (asserted < 1e-7).
- **Expected outcome:** every printed stationary composition is a stationary
  point of D (residual `max_i |ln w_i + ln gamma_i(w) - d_i - k| = 0` with
  `D = k`), the trivial point `w = z` gives `D = 0` exactly, and the recomputed
  D matches the printed D.
- **Tolerance:** asserted stationarity residual < 1e-10, composition match
  atol 1e-3 (the paper prints three digits), `|D/D_printed - 1| < 2e-5` for
  the undisputed points. Achieved: residuals 1.6e-16 to 3.7e-15; composition
  agreement 1.5e-05 to 3.8e-04; `D = k` to ~1e-16.
- **Recomputed vs printed D (all eight non-trivial points):**

  | feed | refined w | printed D | recomputed D | rel. diff |
  | --- | --- | --- | --- | --- |
  | (0.148, 0.052, 0.80) | (0.143614, 0.049882, 0.806503) | +4.5711e-08 | +4.571050e-08 | 1.09e-05 |
  | (0.148, 0.052, 0.80) | (0.114336, 0.035993, 0.849671) | -9.9851e-06 | -9.851037e-06 | **1.34e-02** |
  | (0.12, 0.08, 0.80) | (0.129742, 0.089090, 0.781167) | -3.0693e-06 | -3.069309e-06 | 3.03e-06 |
  | (0.12, 0.08, 0.80) | (0.059745, 0.028236, 0.912019) | -7.4818e-04 | -7.481797e-04 | 4.15e-07 |
  | (0.13, 0.07, 0.80) | (0.137539, 0.075645, 0.786817) | -8.6268e-07 | -8.626898e-07 | 1.14e-05 |
  | (0.13, 0.07, 0.80) | (0.073787, 0.030312, 0.895901) | -3.2762e-04 | -3.276225e-04 | 7.76e-06 |
  | (0.12, 0.05, 0.83) | (0.157573, 0.072897, 0.769530) | -5.7360e-05 | -5.735988e-05 | 2.04e-06 |
  | (0.12, 0.05, 0.83) | (0.093970, 0.034854, 0.871177) | -3.0088e-05 | -3.088768e-05 | **2.66e-02** |

  The four trivial points (`w = z`) reproduce `D = 0` to 0.0 or 1.8e-17.
- **Two printed D values are treated as typographical errors, not as model
  disagreement (recorded honestly rather than choosing whichever number
  passes):**
  1. Feed (0.148, 0.052, 0.80), root near (0.114, 0.036, 0.850): printed
     -9.9851e-06, recomputed -9.851037e-06. The paper's own section 4 text
     reports the verified interval
     `x1 = [0.11433639929296194, 0.11433639934254627]` for this root; feeding
     the printed `G` matrix verbatim into the same stationarity solve returns
     `x1 = 0.11433639931776886`, **inside** that interval, with
     D = -9.850999e-06. So the composition is confirmed to 17 digits while the
     printed D digit string is not reproducible; it looks like a duplicated
     leading digit of 9.8510.
  2. Feed (0.12, 0.05, 0.83), root near (0.094, 0.0349, 0.871): printed
     -3.0088e-05, recomputed -3.088768e-05 (-3.088768e-05 with the printed G
     matrix as well). Looks like 3.0888 printed as 3.0088.
  Both mismatches are two to three orders of magnitude larger than the
  round-off of a 5-digit printed value and than the spread of the six
  reproducing points, and both are insensitive to whether the rounded alpha or
  the printed G matrix is used. The tests assert against the **recomputed**
  values and additionally assert that the mismatch exceeds 1e-3 relative, so a
  future change that accidentally "fixed" them would fail.
- **Solver-behaviour note (asserted, not incidental):** successive substitution
  `ln W_i = d_i - ln gamma_i(w)` converges to six of the seven non-trivial
  roots but drifts to the trivial solution `w = z` for the near-plait-point
  root of the z1 = 0.148 feed (printed D = +4.5711e-08, a saddle rather than a
  minimum). That is precisely the initialization dependence the paper sets out
  to eliminate. The reported stationary points are therefore obtained with a
  damped Newton solve of the constrained stationarity system, whose Jacobian is
  built by central differences and shares no derivative code with the model.
- **Independent route:** published reference values plus a from-scratch
  tangent-plane distance and stationary-point solve written in the test module
  (no `thermo` dependency); cross-checked against `thermo` via Case N-2.
- **Test path:** `tests/validation/test_nrtl_tessier2000.py`
  (`test_published_G_matrix_is_reproduced_by_the_implied_alpha`,
  `test_tessier2000_table2_stationary_points_are_reproduced`,
  `test_successive_substitution_reaches_the_printed_minima`).
- **Script:** `examples/validation/07_nrtl_tessier_stationary_points.py`
  (golden path; prints printed vs recomputed points and D values with a
  PASS / KNOWN-TYPO / FAIL line per point; no optional dependency).
- **Data provenance caution (packaged defaults):** the pairs shipped in
  `src/chemthermo/parameters/data/activity/nrtl.json` (Methane/Ethane,
  Benzene/Water) are synthetic illustrative placeholders added ad hoc in
  commits e5ccd8f and ecd476e with no source. They are now labelled as such in
  the file's `provenance` block and per-pair `source` field, and in README /
  `examples/README.md`. They must never be cited as physical parameters.
