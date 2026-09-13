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

---

## Case S-6: `stability_tp` with NRTL reproduces the Tessier (2000) Problem 1 global minima

- **Source:** S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable
  phase stability analysis for excess Gibbs energy models", Chemical
  Engineering Science 55 (2000) 1785-1796. Author copy:
  https://academicweb.nd.edu/~markst/srt2000.pdf
- **Location:** Table 1 (NRTL parameters), Table 2 (stationary points and D).
  Where Case N-3 tests the *model* (is the published stationary point a
  stationary point of our `ln gamma`?), this case tests the *solver*: given
  only the feed and the model, does `stability_tp` find the global minimum and
  report the right verdict?
- **Assumptions:** Both phases are liquids with the same pure-liquid reference
  state, so `tpd(w) = sum_i w_i [ln w_i + ln gamma_i(w) - ln z_i -
  ln gamma_i(z)]` is exactly the paper's `D(x)`. tau is dimensionless and
  temperature independent as printed; T = 298.15 K and P = 101325 Pa are passed
  only because the API requires them (`pressure_dependent` is recorded as
  `False`).
- **Components / units:** n-propanol(1) / n-butanol(2) / water(3), databank
  names `1-Propanol`, `n-Butanol`, `Water`. `tpd` is dimensionless.
- **Parameters and provenance:** `tests/fixtures/nrtl/tessier2000_problem1.json`
  (see Case N-3). Not packaged runtime data.
- **Expected outcome:** all four feeds UNSTABLE, `tpd_min` equal to the feed's
  lowest stationary D (recomputed, Case N-3), minimizing `w` equal to that
  stationary point.
- **Results (defaults: `tol = 1e-10`, `tpd_tol = 1e-8`, `ssi_iterations = 50`):**

  | feed | verdict | `tpd_min` | reference D (refined) | rel. diff | minimizing `w` | max abs `dw` | winning trial | ssi + newton |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | (0.148, 0.052, 0.80) | unstable | -9.8510373256e-06 | -9.8510373257e-06 | 1.9e-11 | (0.11433639, 0.03599266, 0.84967095) | 1.5e-10 | `pure-Water` | 50 + 5 |
  | (0.12, 0.08, 0.80) | unstable | -7.4817968990e-04 | -7.4817968990e-04 | 1.6e-14 | (0.05974494, 0.02823583, 0.91201923) | 6.6e-17 | `pure-Water` | 50 + 3 |
  | (0.13, 0.07, 0.80) | unstable | -3.2762254313e-04 | -3.2762254313e-04 | 8.9e-14 | (0.07378748, 0.03031163, 0.89590089) | 2.3e-12 | `pure-Water` | 50 + 3 |
  | (0.12, 0.05, 0.83) | unstable | -5.7359882819e-05 | -5.7359882819e-05 | 1.7e-12 | (0.15757285, 0.07289671, 0.76953044) | 1.2e-10 | `pure-n-Butanol` | 50 + 4 |

  Every trial of every feed converged (`converged_stage = "second-order"` for
  all of them). Against the *printed* five-digit D the relative differences are
  1.34e-02 (the Case N-3 typo), 4.1e-07, 7.8e-06 and 2.0e-06.
- **Tolerance:** asserted `tpd_min` to rel 1e-5 against the recomputed minima,
  `w` to abs 1e-6 against an independent damped-Newton refinement written in
  the test module. Achieved: 1.9e-11 relative and 1.5e-10 absolute (worst).
- **`tpd_tol` choice:** the default 1e-8. The smallest |D| here is 9.85e-06,
  three orders of magnitude above it, so the verdicts do not depend on the
  tolerance. A `tpd_tol` of 1e-5 or looser would wrongly call the z1 = 0.148
  feed stable.
- **Second-order stage is load-bearing here (asserted, not incidental):** with
  `second_order=False` and `max_iter=1000`, *no* trial at feed
  (0.148, 0.052, 0.80) reaches `tol = 1e-10` - all three end with
  `termination_reason = "max_iter"` and the result is `"inconclusive"`. The
  fixed-point map's contraction ratio is too close to one near the plait point.
- **Stationary points reached / not reached (honest record):** seven of the
  eight non-trivial printed Table 2 points are reached by some trial. The one
  that is not is the near-plait-point *saddle* of the z1 = 0.148 feed (printed
  D = +4.5711e-08): trials heading toward it are pulled into the trivial
  solution. It has positive D, so it cannot change a verdict.
- **Independent route:** a from-scratch damped-Newton refinement of the
  constrained stationarity system written in the test module (no shared code
  with `chemthermo.stability`), plus `thermo` 0.6.0 `ln gamma` (Case S-7's
  cross-check covers both problems).
- **Test path:** `tests/validation/test_stability_nrtl_tessier2000.py`
  (`test_problem1_stability_tp_finds_the_table2_global_minima`,
  `test_problem1_reaches_every_printed_table2_point_except_the_plait_saddle`),
  `tests/test_stability_activity.py::test_successive_substitution_alone_cannot_solve_the_near_plait_feed`.
- **Script:** `examples/validation/08_stability_nrtl_tessier2000.py`.

---

## Case S-7: `stability_tp` with NRTL reproduces the Tessier (2000) Problem 2 global minima

- **Source:** Tessier, Brennecke and Stadtherr (2000), op. cit.
- **Location:** Section 4.2, Table 4 (NRTL parameters), Table 5 (stationary
  points and D). The paper states "All but the second feed listed are
  unstable", which is the verdict under test.
- **Assumptions:** as Case S-6.
- **Components / units:** n-propanol(1) / n-butanol(2) / benzene(3) / water(4),
  databank names `1-Propanol`, `n-Butanol`, `Benzene`, `Water`.
- **Parameters and provenance:**
  `tests/fixtures/nrtl/tessier2000_problem2.json`. Table 4 prints `G_ij` and
  `tau_ij` but not `alpha_ij`; alpha is implied by
  `alpha_ij = -ln(G_ij)/tau_ij`. The two directions of each pair agree to
  max |alpha_ij - alpha_ji| = **2.597e-05**, i.e. alpha is symmetric to the
  precision the five-decimal printed G supports; the fixture stores the
  symmetric average rounded to three decimals
  (0.494, 0.286, 0.282, 0.297, 0.344, 0.281). Re-exponentiating it reproduces
  the printed G to max |dG| = **4.413e-06**, which moves the recomputed D by at
  most 9.0e-05 relative (measured against using the printed G matrix verbatim).
  **Redistribution caution:** Table 4 attributes these values to Gmehling et
  al., DECHEMA Chemistry Data Series (1977-1990). They are fitted third-party
  data reproduced here solely as a cited test fixture; they are NOT packaged
  runtime data and must never be loaded as defaults.
- **Expected outcome:** feed 2 STABLE (printed D values 0.03079, 0.06532,
  0.00000, all >= 0); the other four UNSTABLE with the printed global minima.
- **Results (defaults):**

  | feed | verdict | `tpd_min` | reference D (refined) | rel. diff vs refined | printed D | rel. diff vs printed | winning trial | ssi + newton |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | (0.148, 0.052, 0.600, 0.200) | unstable | -3.3982528336e-01 | -3.3982528336e-01 | 6.5e-16 | -0.33982 | 1.56e-05 | `pure-Water` | 18 + 0 |
  | (0.25, 0.25, 0.25, 0.25) | **stable** | +3.0793112777e-02 | +3.0793112777e-02 | 1.9e-15 | +0.03079 | 1.01e-04 | `pure-Water` | 35 + 0 |
  | (0.148, 0.052, 0.700, 0.100) | unstable | -3.1097303625e-01 | -3.1097303625e-01 | 8.9e-16 | -0.31097 | 9.76e-06 | `pure-n-Butanol` | 30 + 0 |
  | (0.25, 0.15, 0.40, 0.20) | unstable | -3.8665151099e-02 | -3.8665151099e-02 | 5.4e-15 | -0.03867 | 1.25e-04 | `pure-Water` | 31 + 0 |
  | (0.25, 0.15, 0.35, 0.25) | unstable | -7.3625780741e-02 | -7.3625780741e-02 | 9.4e-16 | -0.07363 | 5.73e-05 | `pure-Water` | 28 + 0 |

  The minimizing composition matches the independent refinement to
  max |dw| = 3.7e-12 in every case, and every trial of every feed converged.
- **Tolerance:** asserted `tpd_min` to rel 1e-6 against the refined values and
  `w` to abs 1e-6; asserted the refined D against the printed D at rel 2.5e-04
  (five printed digits alone justify ~5e-05; the alpha recovery adds up to
  9.0e-05). Achieved: 5.4e-15 relative against the refinement, 1.25e-04 against
  the printed digits.
- **Stationary points reached / not reached (honest record):** of the eight
  non-trivial printed Table 5 points, five are the global minima above and are
  all reached. Of the remaining three non-global points, two are reached - the
  benzene-rich points of the z3-rich feeds, printed D = -0.03365 (recomputed
  -3.3651657e-02) and -3.1279e-03 (recomputed -3.1280532e-03) - and **three are
  not**: the positive-D points +0.06532 (feed 2), +0.02268 (feed 4) and
  +0.01066 (feed 5). Trials heading toward those collapse onto the trivial
  solution. All three have D > 0, so none can change a verdict, but the trial
  set is genuinely not exhaustive and that is what "stable" is bounded by.
- **One printed D is not reproduced and is recorded, not accommodated:** feed
  (0.25, 0.15, 0.40, 0.20), stationary point near
  (0.195, 7.86e-2, 0.114, 0.613). Refining the printed composition to a
  stationarity residual of 6.7e-16 gives
  w = (0.194545, 0.078562, 0.113961, 0.612933) - within 3.4e-04 of the printed
  three-digit composition - but D = **+2.66799e-02**, not the printed
  +2.26800e-02 (17.6% apart). The recomputed value is +2.667985e-02 with the
  rounded alpha and +2.667994e-02 with the printed G matrix verbatim, so it is
  not a parameter-rounding artefact, and every other Table 5 point reproduces
  to 2.5e-04 relative or better. Treated as a typographical error
  (2.6680 printed as 2.2680). The test asserts the mismatch exceeds 1e-2
  relative, so a future change that "fixed" it would fail.
- **Independent route:** the same from-scratch damped-Newton refinement as
  Case S-6, plus `thermo` 0.6.0. Recomputing the tangent-plane distance at
  `stability_tp`'s reported minimizer with `thermo.NRTL`'s `ln gamma` agrees to
  max |dD| = **4.72e-16** over the nine feeds of Problems 1 and 2 (asserted
  1e-10 per feed, 1e-12 overall).
- **Not attempted:** confirming the two-liquid split of feed 1 with `thermo`'s
  `FlashVLN`. That needs a full `thermo` property package (pure-component
  `HeatCapacityGas`, `VaporPressure`, volume correlations) for the four
  components before an LLE flash can run, and the paper's parameters are
  dimensionless taus with no temperature attached, so the flash would be
  answering a different question than the table does. The `ln gamma`
  cross-check above is the part that is actually comparable, and it is exact.
- **Test path:** `tests/validation/test_stability_nrtl_tessier2000.py`
  (`test_problem2_alpha_is_symmetric_and_reproduces_the_printed_G_matrix`,
  `test_problem2_components_map_to_the_databank`,
  `test_problem2_stability_tp_reproduces_the_table5_verdicts_and_minima`,
  `test_problem2_records_which_other_printed_points_are_reached`,
  `test_problem2_disputed_printed_D_is_recorded_not_accommodated`,
  `test_stability_minima_agree_with_thermo_ln_gamma`).
- **Script:** `examples/validation/08_stability_nrtl_tessier2000.py`.

---

## Case S-8: n-butanol / water LLE control, and one tangent plane for two liquids

- **Source:** Internal invariant (the definition of liquid-liquid equilibrium),
  applied to the cited n-butanol / water NRTL pair. This is the activity-model
  analogue of Case S-3.
- **Location:** `tests/test_stability_activity.py::test_butanol_water_binary_splits_and_both_phases_share_one_tangent_plane`.
- **Assumptions:** At an LLE split the two liquid phases are the two points
  where one common hyperplane touches the Gibbs surface. Each of them, used as
  a feed, must therefore be marginally stable (`tpd_min = 0`) and must find the
  *other* phase as its stationary point.
- **Components / units:** n-Butanol(1) / Water(2), T = 298.15 K,
  P = 101325 Pa (inert). Mole fractions; `tpd` dimensionless.
- **Parameters and provenance:** the 2-3 pair of Tessier et al. (2000) Table 1:
  tau_12 = 0.90047, tau_21 = 3.51307, alpha = 0.48 (implied by the printed
  `G_23`, `G_32`). Same fixture provenance as Case N-3. Note these are ternary
  parameters used here as a binary sub-system; the claim under test is the
  tangent-plane geometry, not the physical n-butanol / water phase diagram.
- **Expected outcome and results:**
  - Feed z1 = 0.10 (inside the gap): UNSTABLE, `tpd_min = -2.9994488835e-02`
    at w = (0.419473294157, 0.580526705843), found by `pure-n-Butanol`
    (50 ssi + 1 Newton). The second trial finds a second negative stationary
    point at D = -5.2965300e-03. Recomputing `tpd` at the reported minimizer
    from the definition agrees to 1e-12.
  - Conjugate phases from an independent equal-activity solve
    (`x_i gamma_i` equal in both phases, Newton with a finite-difference
    Jacobian, residual 1.1e-16): x1 = 0.019998419467 and x1 = 0.359999661508.
    The feed lies strictly between them.
  - Feed = x1 = 0.0199984: STABLE, `tpd_min = +1.77e-16`, stationary point at
    x1 = 0.359999662 (the other phase, to 1e-9).
  - Feed = x1 = 0.3599997: STABLE, `tpd_min = -1.34e-16`, stationary point at
    x1 = 0.019998419 (the other phase, to 1e-9).
- **Tolerance:** asserted `tpd_min` to rel 1e-9 (unstable feed) and abs 1e-10
  (marginal feeds), conjugate composition to abs 1e-9. Achieved: 1.8e-16 on the
  marginal feeds.
- **Other negative controls in the same module:** a one-component feed and a
  component at zero mole fraction are reported STABLE with `tpd = 0` and a
  single trivial trial; an ideal solution (`tau = 0`, so `gamma = 1`) is STABLE
  at z = (0.5, 0.5), (0.1, 0.9) and (0.9, 0.1), which is the analytic answer
  since `tpd(w) = sum_i w_i ln(w_i / z_i) >= 0`.
- **Independent route:** internal invariant plus an equal-activity binodal
  solve sharing no code with `chemthermo.stability`.
- **Test path:** `tests/test_stability_activity.py`.
- **Script:** `examples/basic/stability_tp_nrtl_lle_demo.py`.

---

## Case F-1: `flash_tp` phase verdicts against `thermo` 0.6.0 `FlashVL`

- **Source:** `thermo` 0.6.0 (`FlashVL` with `PRMIX`), used as an independent
  implementation of the same Peng-Robinson model. Not a literature case: it is a
  cross-implementation check, and it is recorded as such.
- **Location:** `tests/validation/test_flash_phase_detection_vs_thermo.py::test_phase_verdicts_over_a_grid_match_thermo_more_often_than_the_heuristic`.
- **Assumptions:** Both sides use **the same** `Tc`, `Pc` and `omega`, read from
  the chemthermo databank and handed to `thermo`'s `ChemicalConstantsPackage`,
  and `kij = 0` everywhere. `thermo`'s ideal-gas heat capacities come from its
  own databank and do not enter an isothermal-isobaric VLE verdict. A verdict
  difference is therefore a *solver* difference.
- **Components / units:** 5 systems - Methane/Ethane (0.5, 0.5),
  Methane/n-Pentane (0.6, 0.4), Ethane/n-Heptane (0.7, 0.3),
  Methane/Ethane/Propane (0.5, 0.3, 0.2), Propane/n-Butane/n-Pentane
  (0.4, 0.3, 0.3). T in {170, 175, 200, 240, 280, 320, 360} K, P in
  {2e5, 1e6, 1.778e6, 3e6, 8e6} Pa. 175 states. SI units; `tpd` dimensionless.
- **Parameters and provenance:** packaged chemthermo databank
  (`src/chemthermo/data/components.json`, Koretsky 2012); `kij = 0`.
- **Expected outcome:** the tangent-plane path's 1-vs-2 phase verdict should
  match `FlashVL` at least as often as the legacy Wilson heuristic does.
- **Results:**

  | path | verdicts matching `thermo` | misses |
  | --- | --- | --- |
  | `phase_detection="tangent-plane"` (default) | **175 / 175** | 0 |
  | `phase_detection="wilson-heuristic"` (legacy) | 166 / 175 | 8 `ConvergenceError` + 1 wrong single-phase |

  All 8 legacy non-convergences are states where `thermo` says single phase and
  the tangent-plane path now returns a single phase. The 1 wrong single-phase is
  Case F-2 below. Over the 56 states that are two-phase on both sides, the worst
  `|beta - beta_thermo|` is **1.601e-04**.
- **Wider scan (recorded, not asserted):** the same comparison over 1144 states
  (8 databank mixtures, T 150-450 K in 11 steps, P 1e5-3.2e7 Pa in 13 geometric
  steps) gives: tangent-plane 1138/1144 verdicts matching `thermo`, legacy
  1061/1144. 75 states move from `ConvergenceError` to a single-phase answer
  (`thermo` agrees single-phase on 75/75); 2 move from single-phase to
  two-phase; **0** move from two-phase to single-phase. Of the 6 remaining
  tangent-plane disagreements, 5 are Methane/Propane/n-Decane (0.7, 0.2, 0.1) at
  1.98e7-1.22e7 Pa where chemthermo splits (`tpd_min` from -5.2e-2 to -4.0e-3,
  `delta_g_split_rt < 0`) and `thermo` returns `VF = 0`; those 5 disagree with
  `thermo` **before and after** this slice, so they are not caused by it. They
  are now **adjudicated** (see below) rather than left open. The 6th is the
  near-critical Methane/n-Pentane state below.
- **Adjudication of the 5 Methane/Propane/n-Decane disagreements:** at
  (T, P) = (270 K, 1.22e7 Pa), (300 K, 1.5e7 Pa), (330 K, 1.98e7 Pa),
  (360 K, 1.5e7 Pa) and (450 K, 1.5e7 Pa), `thermo`'s own fugacity model says
  the feed is unstable too, so the disagreement is a `thermo` stability-*search*
  miss, not a chemthermo error. Method: build a `thermo` `PRMIX` `CEOSLiquid`
  with chemthermo's own `Tc`, `Pc` and `omega` (`kij = 0`, the same constants
  `flash_tp` used), take chemthermo's `stability_tp` minimizing trial
  composition `w`, and evaluate `tpd(w) = sum_i w_i [ln w_i + ln phi_i(w) -
  ln z_i - ln phi_i(z)]` using *thermo's* `lnphis_at_zs(zs, most_stable=True)`
  (the minimum-Gibbs root) at both `w` and the feed `z`. Results (chemthermo
  `tpd_min` vs. thermo-evaluated tpd at `w`, both negative; thermo `FlashVL`
  `VF` at each state):

  | T (K) | P (Pa) | chemthermo `tpd_min` | thermo tpd at `w` | thermo `VF` |
  | --- | --- | --- | --- | --- |
  | 270 | 1.22e7 | -5.278114e-2 | -5.278752e-2 | 0.0 |
  | 300 | 1.5e7  | -3.909584e-2 | -3.910660e-2 | 0.0 |
  | 330 | 1.98e7 | -6.908276e-3 | -6.914199e-3 | 0.0 |
  | 360 | 1.5e7  | -4.712275e-2 | -4.713469e-2 | 0.0 |
  | 450 | 1.5e7  | -3.590357e-2 | -3.593922e-2 | 0.0 |

  All 5 pairs agree to within 6.6e-5 absolute (the known rounded-PR-constants
  gap between the two implementations); `thermo`'s own model puts a negative
  tpd at chemthermo's stationary point in every case, and chemthermo's
  `flash_tp` split lowers the Gibbs energy (`delta_g_split_rt < 0`) at all 5.
  Asserted (not just recorded, with a 5e-5 absolute tolerance on the tpd match)
  in `tests/validation/test_flash_thermo_disagreements_adjudicated.py`, which
  is skipped if `thermo` is not installed and records `thermo`'s `VF` without
  asserting its value, so the test stays valid (and prints a note instead of
  failing) if a future `thermo` release finds the split.
- **Regression:** over the 47 states of the in-repo grid
  (`tests/test_flash_phase_detection.py`) where both paths find two phases, the
  worst relative vapor-fraction difference is **8.60e-07** (asserted < 1e-6).
  The legacy path reproduces the pre-slice values bit-identically
  (`0.6745181801306899` binary, `0.46829043780053325` ternary).
- **Known failure kept honest:** Methane/n-Pentane (0.6, 0.4) at 390 K,
  1.2236e7 Pa is weakly unstable (`tpd_min = -1.17e-3`, near-critical) and the
  successive-substitution split hits the iteration limit in both paths. This
  slice adds no acceleration to the *flash*, so it still raises
  `ConvergenceError`.
- **Tolerance:** asserted 175/175 verdict match, `tangent >= legacy`, and
  `|beta - beta_thermo| < 1e-3`. Achieved 1.601e-04.
- **Independent route:** `thermo` 0.6.0 `FlashVL`, a separate implementation of
  the same EOS with its own stability test and phase-split solver.
- **Test path:** `tests/validation/test_flash_phase_detection_vs_thermo.py`,
  `tests/test_flash_phase_detection.py::test_both_paths_agree_on_every_state_where_both_find_two_phases`,
  `tests/test_flash_phase_detection.py::test_legacy_path_reproduces_the_pre_slice_numbers_exactly`,
  `tests/validation/test_flash_thermo_disagreements_adjudicated.py` (the 5
  Methane/Propane/n-Decane adjudication above).
- **Script:** `examples/basic/flash_tp_auto_phase_demo.py`,
  `examples/validation/00_reference_case.py`.

---

## Case F-2: a state where the Wilson heuristic and the tangent plane disagree

- **Source:** found by scanning the databank's Peng-Robinson states (see Case
  F-1); adjudicated by `thermo` 0.6.0 `FlashVL` and by the Gibbs-energy
  criterion. Not a literature case.
- **Location:** `tests/validation/test_flash_phase_detection_vs_thermo.py::test_the_disagreement_state_is_two_phase_in_thermo_and_lowers_the_gibbs_energy`
  and `tests/test_flash_phase_detection.py::test_heuristic_calls_a_two_phase_feed_single_phase_and_tangent_plane_does_not`.
- **Assumptions:** Peng-Robinson, `kij = 0`, chemthermo databank constants on
  both sides. n-Pentane's normal melting point is ~143 K, so 175 K is a
  physically liquid-pentane state; Peng-Robinson has no solid phase, so the
  question under test is the fluid-phase answer.
- **Components / units:** Methane(0.6) / n-Pentane(0.4), T = 175.0 K,
  P = 1.778e6 Pa. Mole fractions; `tpd` and `delta_g/RT` dimensionless.
- **Parameters and provenance:** packaged chemthermo databank. Methane
  Tc = 190.6 K, Pc = 4.600e6 Pa, omega = 0.008; n-Pentane Tc = 469.6 K,
  Pc = 3.374e6 Pa, omega = 0.251.
- **Expected outcome:** two phases.
- **Results:**
  - **Legacy heuristic:** single `liquid`, `termination_reason = "rr_no_root"`.
    The Wilson K-values straddle 1 (`k_min = 2.3121e-05`, `k_max = 1.5964`), so
    the K-bound test does not fire, but Rachford-Rice cannot bracket a root for
    them: `f(0) = -4.213e-02` and `f(1) = -1.730e+04` have the same sign.
  - **Tangent-plane path:** `stability_status = "unstable"`,
    `tpd_min = -4.265697e-02`, incipient phase vapor-like
    (`w = (0.99998, 1.824e-05)`), converged split
    `beta = 0.0792360`, `x = (0.565580, 0.434420)`,
    `y = (0.9999798, 2.0201e-05)`, `mass_balance_residual = 7.2e-15`,
    `fugacity_residual = 6.2e-09`, `delta_g_split_rt = -1.76661e-03`.
  - **`thermo` `FlashVL` (same constants, kij = 0):** `VF = 0.07912756`,
    `x = (0.5656310, 0.4343690)`, `y = (0.9999798, 2.0188e-05)`.
    `|d beta| = 1.08e-04`, `max |dx| = 6.8e-05`, `max |dy| = 1.9e-07`.
  - **Gibbs-energy criterion, evaluated at `thermo`'s split with chemthermo's
    own fugacity coefficients on the minimum-Gibbs root:**
    `G_feed/RT = -5.0263909`, `G_split/RT = -5.0281575`, so
    `dG/RT = -1.7666e-03 < 0`. The split is the lower-Gibbs state, independently
    of which solver produced it.
- **Two further states of the same kind** (same system, recorded not asserted):
  150 K / 6.8399e5 Pa (`tpd_min = -4.5567e-02`, chemthermo `beta = 0.080184`,
  `thermo` 0.080020, `dG/RT = -1.9074e-03`) and 210 K / 4.6784e6 Pa
  (`tpd_min = -2.2673e-02`, chemthermo `beta = 0.049352`, `thermo` 0.049299,
  `dG/RT = -5.7644e-04`).
- **Tolerance:** asserted `beta` to abs 5e-4 against `thermo`, compositions to
  abs 1e-3, `dG/RT` to rel 1e-3 against -1.7666e-03, and `dG/RT < 0`.
- **Independent route:** two, and they agree - `thermo` 0.6.0 `FlashVL`, and a
  Gibbs-energy comparison computed from the definition (not from the flash
  solver) at `thermo`'s compositions.
- **Test path:** the two tests named above.
- **Script:** `examples/basic/flash_tp_auto_phase_demo.py` (case 3).

---

## Case F-3: verification invariants of every converged two-phase flash

- **Source:** internal invariants (material balance, equal fugacities, the
  Gibbs-energy criterion for a phase split, Michelsen's instability proof).
  This entry's "independent route" is an internal invariant, stated as required
  by the ledger rules.
- **Location:** `tests/test_flash_phase_detection.py::test_grid_invariants`.
- **Assumptions:** A converged isothermal-isobaric two-phase solution must
  satisfy `z = beta y + (1 - beta) x`, `x_i phi_i^L = y_i phi_i^V` for every
  component present in both phases, `0 < beta < 1`, and must have a lower molar
  Gibbs energy than the single-phase feed. A feed that splits must have been
  found unstable (`tpd_min < 0`), and a feed reported as one phase must have
  been found stable.
- **Components / units:** 6 binary and ternary hydrocarbon systems
  (Methane/Ethane, Methane/Propane, Ethane/n-Heptane, Methane/n-Pentane,
  Methane/Ethane/Propane, Propane/n-Butane/n-Pentane) over T in
  {170, 200, 240, 280, 320, 360} K and P in {2e5, 1e6, 3e6, 8e6} Pa - 144
  states, of which 47 are two-phase and 97 single-phase (0 refused). Mole
  fractions; residuals dimensionless.
- **Parameters and provenance:** packaged chemthermo databank; `kij = 0`.
- **Expected outcome and results (worst values over the 47 two-phase states):**

  | invariant | required | achieved (worst) |
  | --- | --- | --- |
  | `max_i |z_i - (beta y_i + (1-beta) x_i)|` | < 1e-10 | **2.014e-13** |
  | `max_i |ln(x_i phi_i^L) - ln(y_i phi_i^V)|` | < 1e-6 | **8.246e-09** |
  | `delta_g_split_rt` | < 0 | **-1.4575e-04** (least negative) |
  | `beta` | in (0, 1) | all |
  | feed `tpd_min` | < 0 | all |
  | single-phase `stability_status` | `"stable"` | all 97 |

  The material balance is re-checked in the test from the returned phase
  compositions, not from the residual the solver wrote into diagnostics.
- **Permutation invariance and determinism:** over 4 states (two-phase and
  single-phase, binary and ternary) and every component permutation, the phase
  names are identical, `|d beta| = 0.0` and the worst composition difference
  after undoing the permutation is **2.22e-16**. Repeated calls are
  bit-identical including the full diagnostics mapping.
- **Failure semantics:** `StabilitySettings(max_iter=1, second_order=False)`
  makes the stability analysis `"inconclusive"`; `flash_tp` then raises
  `ConvergenceError` in tangent-plane mode and still returns
  `0.46829043780053325` through `phase_detection="wilson-heuristic"`.
  `FlashSettings(max_iter=1, tol=1e-12)` still raises `ConvergenceError` for
  the split, as before this slice.
- **Tolerance:** as tabulated above; also asserted `>= 30` two-phase states so
  the sweep cannot silently shrink.
- **Independent route:** internal invariant (the test recomputes the material
  balance from the returned compositions; the fugacity and Gibbs residuals come
  from the definitions in `chemthermo.flash.tp._verify_split`, which the solver
  does not use to iterate).
- **Test path:** `tests/test_flash_phase_detection.py`.
- **Script:** `examples/basic/flash_tp_auto_phase_demo.py`.

---

## Case L-1: Liquid-liquid tie-lines of Tessier (2000) Problem 1

- **Source:** S. R. Tessier, J. F. Brennecke and M. A. Stadtherr, "Reliable
  phase stability analysis for excess Gibbs energy models", Chemical
  Engineering Science 55 (2000) 1785-1796, for the **parameters and the feeds**
  (Table 1, Table 2). The paper publishes *stationary points of the
  tangent-plane distance*, **not tie-lines**, so no printed tie-line exists to
  compare against and none is claimed. The published verdicts for these feeds
  (all unstable) are checked in Case S-6.
- **Location:** `tests/validation/test_flash_lle_tessier2000.py::test_problem1_tie_lines`.
- **Assumptions:** Two liquid phases at one temperature share a pure-liquid
  reference state, so it cancels and equilibrium is `x_i^I gamma_i^I =
  x_i^II gamma_i^II`. The tie-line is a property of (T, P, model, feed) that any
  correct solver must reproduce; the feed must lie on it.
- **Components / units:** 1-Propanol(1) / n-Butanol(2) / Water(3), T = 298.15 K,
  P = 101325 Pa (validated, inert). Mole fractions; residuals dimensionless.
- **Parameters and provenance:** `tests/fixtures/nrtl/tessier2000_problem1.json`
  (Table 1; alpha implied by the printed G matrix). Same fixture and provenance
  as Cases N-3 and S-6. Not packaged runtime data.
- **Expected outcome and results** (phases ordered by the 1-Propanol fraction;
  `beta` is the fraction of the second one):

  | feed z | phase A | phase B | fraction of B | `delta_g_split_rt` |
  | --- | --- | --- | --- | --- |
  | (0.12, 0.08, 0.80) | (0.063945, 0.030845, 0.905210) | (0.147546, 0.104156, 0.748298) | 0.670502 | -1.952259e-04 |
  | (0.13, 0.07, 0.80) | (0.077843, 0.032599, 0.889559) | (0.153205, 0.086640, 0.760155) | 0.692089 | -7.868224e-05 |
  | (0.12, 0.05, 0.83) | (0.093662, 0.034385, 0.871953) | (0.156101, 0.071402, 0.772496) | 0.421822 | -4.073533e-05 |
  | (0.148, 0.052, 0.80) | (0.116727, 0.037073, 0.846200) | (0.154474, 0.055090, 0.790436) | 0.828496 | -1.065209e-06 |

  Residuals (worst over the four feeds): equal-activity
  `max_i |ln(x_i^I gamma_i^I) - ln(x_i^II gamma_i^II)| = 4.441e-16`, material
  balance `2.220e-16`, post-split `tpd_min` `-1.414e-16` (both phases reported
  `"stable"` on every feed, so none of these states needs a third phase).
  Iterations per stage: 50 successive substitutions (the budget) plus 5, 6, 7
  and 5 second-order steps respectively.
- **Tolerance:** asserted composition and phase fraction to abs 1e-6 against the
  independent route, equal-activity residual < 1e-10, material balance < 1e-12,
  `delta_g_split_rt < 0`. Achieved against the independent route: worst
  difference **1.081e-12** (the near-plait feed); the other three are <= 9.7e-14.
- **Independent route:** a from-scratch solve written in the same test module,
  sharing no code with `chemthermo.flash`: its own successive-substitution loop
  (with its own Rachford-Rice bisection) followed by a damped Newton solve of
  the full system in `(x^I, x^II, beta)` - `n` equal-activity equations, `n - 1`
  material balances and both normalization constraints - with a
  central-difference Jacobian. Final residual <= 1.8e-15. It needs 536, 794,
  1315 and 3922 substitutions before the Newton stage, which is the measurement
  behind the claim that a second-order stage is required (ADR-0009 decision 3).
- **Test path:** `tests/validation/test_flash_lle_tessier2000.py`.
- **Script:** `examples/validation/09_lle_tessier2000_tie_lines.py`.

---

## Case L-2: Tessier (2000) Problem 2 tie-lines, stable control, and post-split stability

- **Source:** Tessier, Brennecke and Stadtherr (2000), section 4.2, Table 4
  (parameters) and Table 5 (feeds and stability verdicts). As in Case L-1, the
  paper prints no tie-lines; the verdicts are checked in Case S-7.
- **Location:** `tests/validation/test_flash_lle_tessier2000.py::test_problem2_tie_lines_and_stable_control`.
- **Assumptions:** As Case L-1, plus: the paper's four unstable feeds must split
  into a benzene-rich and a water-rich liquid, its stable feed must come back as
  a single liquid, and **both** phases of every split must themselves be stable
  (otherwise a third liquid exists and the two-phase answer is wrong).
- **Components / units:** 1-Propanol(1) / n-Butanol(2) / Benzene(3) / Water(4),
  T = 298.15 K, P = 101325 Pa (inert). Mole fractions.
- **Parameters and provenance:** `tests/fixtures/nrtl/tessier2000_problem2.json`
  (Table 4, regressed from the DECHEMA Chemistry Data Series; test fixture only,
  never packaged runtime data). Same provenance note as Case S-7.
- **Expected outcome and results:**

  | feed z | benzene-rich phase | water-rich phase | water-rich fraction |
  | --- | --- | --- | --- |
  | (0.148, 0.052, 0.600, 0.200) | (0.166960, 0.061916, 0.717646, 0.053478) | (0.052871, 0.002251, 0.009743, 0.935136) | 0.166189 |
  | (0.148, 0.052, 0.700, 0.100) | (0.154123, 0.055232, 0.744766, 0.045880) | (0.053605, 0.002170, 0.009833, 0.934392) | 0.060911 |
  | (0.25, 0.15, 0.40, 0.20) | (0.255929, 0.154162, 0.411130, 0.178779) | (0.041250, 0.003483, 0.008157, 0.947110) | 0.027619 |
  | (0.25, 0.15, 0.35, 0.25) | (0.265823, 0.161123, 0.375985, 0.197069) | (0.041786, 0.003625, 0.008064, 0.946525) | 0.070626 |

  `delta_g_split_rt` = -4.358433e-02, -1.171946e-02, -5.581938e-04 and
  -2.884102e-03. Worst equal-activity residual **1.821e-14**, worst material
  balance **1.110e-16**. Iterations per stage: 38+1, 31+1, 35+1 and 38+1
  (successive substitution + second order); these feeds converge the first stage
  inside its budget and the second-order stage only polishes them from ~6e-09 to
  ~1e-14.
  **Post-split:** all eight phases report `"stable"`; worst post-split `tpd_min`
  is **-7.116e-16**. No feed here needs a third phase.
  **Stable control:** z = (0.25, 0.25, 0.25, 0.25) returns a single phase named
  `"liquid"` with `tpd_min = +3.079311e-02` and `vapor_fraction is None`.
- **Tolerance:** as Case L-1. Achieved against the independent route: worst
  difference **4.441e-16**.
- **Independent route:** the same from-scratch solve as Case L-1 (its own
  successive substitution plus a damped Newton on the full system), which needs
  50-56 substitutions on these feeds.
- **Test path:** `tests/validation/test_flash_lle_tessier2000.py`.
- **Script:** `examples/validation/09_lle_tessier2000_tie_lines.py`.

---

## Case L-3: n-butanol / water binodal, lever rule, and a `thermo` cross-check

- **Source:** Internal invariants (a tie-line does not depend on the feed; the
  lever rule) applied to the cited n-butanol / water NRTL pair, plus `thermo`
  0.6.0 as an external route. This is the flash-level companion of Case S-8.
- **Location:** `tests/test_flash_lle.py::test_binary_binodal_is_independent_of_the_feed_and_obeys_the_lever_rule`
  and `tests/validation/test_flash_lle_tessier2000.py::test_thermo_agrees_on_the_binary_binodal_but_its_flash_will_not_split`.
- **Assumptions:** For a binary at fixed T and P the miscibility gap has exactly
  one tie-line: every feed strictly inside it must return the *same* two
  conjugate compositions, with only the phase amounts changing, and those
  amounts must follow the lever rule exactly. A feed outside the gap must return
  one phase.
- **Components / units:** n-Butanol(1) / Water(2), T = 298.15 K, P = 101325 Pa
  (inert). Mole fractions.
- **Parameters and provenance:** the 2-3 pair of Tessier et al. (2000) Table 1:
  tau_12 = 0.90047, tau_21 = 3.51307, alpha = 0.48 (implied by the printed
  `G_23`, `G_32`). Ternary parameters used as a binary sub-system; the claim
  under test is the equilibrium arithmetic, not the physical n-butanol / water
  phase diagram.
- **Expected outcome and results:**
  - Binodal from an equal-activity solve written in the test (Newton, FD
    Jacobian, residual < 1e-14): x1 = **0.019998419467** and
    **0.359999661508** - the same pair as Case S-8.
  - Feeds 0.05, 0.10, 0.20 and 0.30 all split and all return that pair;
    worst composition deviation over the four feeds **1.7e-12**.
  - Phase fractions of the butanol-rich phase: 0.088240, 0.235298, 0.470586 and
    0.176469; lever-rule error **<= 2.8e-12**.
  - `delta_g_split_rt` = -4.105485e-03, -1.053007e-02, -8.580513e-03,
    -1.622934e-03; worst equal-activity residual 7.654e-13; worst material
    balance 1.110e-16; all four post-split checks `"stable"` with worst
    `tpd_min` -1.358e-13.
  - Feed 0.45 (above the upper branch) returns a single `"liquid"` with
    `tpd_min = +4.078314e-02`.
  - **Label independence:** at z1 = 0.10 the butanol-rich phase comes back as
    `liquid2`; at z1 = 0.20 it comes back as `liquid1`. The labels are roles
    assigned by the seed (ADR-0009 decision 2), so every assertion compares the
    phase *set*.
  - **Correction to the slice brief:** the brief listed z1 = 0.30 as a
    single-phase control ("outside the gap"). It is not: 0.30 < 0.3599997, so it
    is inside, and it splits. 0.45 is used as the single-phase control instead.
    `thermo` independently calls 0.30 unstable (below).
- **Independent route / external check:** `thermo` 0.6.0 with the same taus and
  alphas.
  - `thermo.FlashVLN` built from two `thermo.GibbsExcessLiquid` phases over one
    `thermo.NRTL` model plus a `CEOSGas` **will not return a liquid-liquid
    split**: it reports `unique_liquid_count == 1` (the two liquid phases are
    deduplicated because they share an excess-Gibbs model) and returns a single
    phase for all four in-gap feeds. No phase-fraction comparison is therefore
    possible and none is claimed. This is recorded, not worked around; the
    assertion `flasher.unique_liquid_count == 1` will fail if a future `thermo`
    changes it.
  - `thermo`'s own `stability_test_Michelsen` on the same model **does** find
    the split, calls all four feeds unstable, and returns the conjugate pair
    x1 = 0.0199994 and 0.3600248. Against this package's 0.0199984 / 0.3599997
    that is **1.0e-06** and **2.5e-05**, both inside the asserted 1e-04.
    `thermo`'s own stationarity residual there is ~1.9e-06, which is the
    precision limit of the comparison.
- **Tolerance:** binodal asserted to abs 1e-6 (achieved 1.7e-12), lever rule to
  abs 1e-6 (achieved 2.8e-12), `thermo` to abs 1e-4 (achieved 2.5e-05).
- **Test path:** `tests/test_flash_lle.py`,
  `tests/validation/test_flash_lle_tessier2000.py`.
- **Script:** `examples/basic/flash_tp_nrtl_lle_demo.py`.

---

## Case L-4: post-split stability of every two-phase phi-phi flash

- **Source:** Internal invariant (two coexisting phases share one tangent plane;
  a converged phase set is only an answer if each phase is itself stable). This
  entry's "independent route" is an internal invariant, stated as required by
  the ledger rules. It extends Case F-3, which verified the split but not the
  phases.
- **Location:** `tests/test_flash_phase_detection.py::test_every_two_phase_grid_state_passes_the_post_split_stability_check`,
  and the guard itself in `tests/test_flash_lle.py::test_a_converged_phase_that_finds_its_partner_is_marginal_not_unstable`
  and `::test_post_split_failure_raises_and_post_split_stability_false_returns`.
- **Assumptions:** Feeding a converged equilibrium phase back into
  `stability_tp` must find its partner phase as the tangent-plane minimizer with
  `tpd = 0` (Case S-3). Any *other* negative stationary point means a third
  phase exists and the two-phase answer is wrong.
- **Components / units:** the 144-state Peng-Robinson grid of Case F-3 (6 binary
  and ternary hydrocarbon systems, T 170-360 K, P 2e5-8e6 Pa, `kij = 0`), of
  which 47 states are two-phase. Mole fractions; `tpd` dimensionless.
- **Parameters and provenance:** packaged chemthermo databank; `kij = 0`.
- **Expected outcome and results:**
  - All 47 two-phase states pass: 94 phases, every one reported `"stable"`,
    `post_split_status == "stable"`, nothing raised. **No state in this
    repository genuinely needs a third phase.**
  - Most negative post-split `tpd_min` anywhere on the grid: **-7.0055e-09**, at
    Ethane/n-Heptane (0.7, 0.3), 360 K, 1 MPa, on the vapor phase. That is
    inside the default `tpd_tol = 1e-8` by a factor of only **1.4**, which is
    the reason the converged-onto-the-partner ("marginal") rule exists.
  - At the default `tol = 1e-8` the minimizer *is* the partner phase to
    `max_i |w_i - x_i^partner| <= 7.0e-09` and
    `sum_i ln(w_i / x_i^partner)^2 <= 4.0e-16` - twelve orders below
    `trivial_tol = 1e-8` - so the rule never had to fire on this grid
    (`marginal` count 0).
  - The rule is exercised directly by loosening the split tolerance:
    at `FlashSettings(tol=1e-4)` the canonical ternary's liquid phase reports a
    raw `tpd_min` of **-2.0e-05** (far past `tpd_tol`) at a point that *is* the
    partner, is classified `"marginal"`, and the flash still returns.
  - **Synthetic failure control:** at `FlashSettings(tol=1e-3)` the
    Ethane/n-Heptane state's phases are only accurate to ~1e-3, the minimizer is
    no longer recognisable as the partner, and `flash_tp` raises
    `ConvergenceError` ("not a stable phase set"). With
    `post_split_stability=False` the same call returns the two-phase result with
    `post_split_stable = False`, `post_split_status = "unstable"` and
    `post_split_tpd_min < -1e-8`. This is **not** a genuine three-phase state -
    it is an under-converged split - and the test says so; it exists because
    that is the only way to reach the failure path in this repository today.
- **Tolerance:** asserted `post_split_tpd_min > -1e-8` on every grid state and
  the worst value pinned to rel 1e-3.
- **Not covered:** `phase_detection="wilson-heuristic"` and `gamma-phi` results
  are not post-split checked (no stability test is available for gamma-phi,
  ADR-0007, and the legacy path exists to reproduce old behavior unchanged).
  Both report `post_split_checked = False` with a reason; asserted in
  `tests/test_flash_lle.py::test_gamma_phi_and_legacy_paths_declare_that_they_were_not_checked`.
- **Independent route:** internal invariant.
- **Test path:** `tests/test_flash_phase_detection.py`, `tests/test_flash_lle.py`.
- **Script:** `examples/basic/flash_tp_nrtl_lle_demo.py` (prints the post-split
  block for a liquid-liquid split).

---

## Case R-1: Modified-Raoult liquid-liquid agrees with the `gamma-gamma` path

- **Source:** Internal invariant (the pure-liquid reference fugacity cancels
  between two liquid phases), applied to the cited n-butanol / water NRTL pair.
  This entry's "independent route" is an internal invariant, stated as required
  by the ledger rules. It is the modified-Raoult companion of Case L-3.
- **Location:** `tests/test_flash_modified_raoult.py::test_liquid_liquid_matches_the_gamma_gamma_path_exactly`,
  `::test_a_feed_outside_the_gap_is_one_liquid`,
  `::test_both_liquid_phases_are_stable_against_the_vapor_candidate`, and
  `tests/test_stability_candidates.py::test_a_liquid_liquid_tangent_plane_is_unchanged_by_the_reference_offset`.
- **Assumptions:** In `flash_mode="modified-raoult"` the liquid candidate's term
  is `ln gamma_i + ln(Psat_i/P)`. When both the feed reference and the trial are
  the liquid candidate, the offset `ln(Psat_i/P)` appears on both sides of the
  tangent-plane distance and cancels identically, and in the split it cancels
  from `K_i = exp(t_i(x) - t_i(y)) = gamma_i(x)/gamma_i(y)`. So at a temperature
  where the vapor candidate is never the lower-Gibbs one, the two modes must
  return the *same* tie-line - not merely a similar one.
- **Components / units:** n-Butanol(1) / Water(2), T = 330 K, P = 101325 Pa.
  Mole fractions. 330 K is inside both Antoine windows (n-Butanol [288, 404] K,
  Water [284, 441] K).
- **Parameters and provenance:** Tessier et al. (2000) Table 1 pair 2-3:
  tau_12 = 0.90047, tau_21 = 3.51307, alpha = 0.48 (implied by the printed
  `G_23`, `G_32`). Antoine from the packaged databank (Koretsky 2012).
- **Expected outcome and results:**
  - Feeds z1 = 0.05, 0.10, 0.20, 0.30 all return `phase_regime == "LLE"`,
    phases `liquid1`/`liquid2`, `vapor_fraction is None`.
  - Worst composition deviation from the `gamma-gamma` answer over the four
    feeds: **1.7e-15** (asserted 1e-10). Worst phase-fraction deviation
    **2.7e-15**.
  - Tie-line: x1 = **0.019998419467** / **0.359999661508**, the Case L-3 /
    Case S-8 binodal.
  - z1 = 0.45 returns a single `"liquid"` with `feed_branch == "liquid"` and
    `tpd_min = +4.078314e-02`.
  - Post-split at z1 = 0.20: both liquids `"stable"` against the vapor
    candidate, `post_split_tpd_min = -1.354e-13`.
  - Stability-level check: at the same state `tpd_min` with `vapor="ideal"`
    equals `tpd_min` with `vapor="none"` to **< 1e-14**, and the minimizing
    compositions agree to **5.1e-12** (the two evaluators reach the same
    stationary point from different trial sets, so the agreement is the
    solver's `tol = 1e-10`, not round-off).
- **Tolerance:** asserted 1e-10 on compositions and phase fractions (achieved
  1.7e-15); 1e-14 on the tangent-plane distance (achieved < 1e-14).
- **Independent route:** internal invariant plus the pre-existing `gamma-gamma`
  path, which Case L-3 validated against `thermo` and an independent
  equal-activity solve.
- **Test path:** `tests/test_flash_modified_raoult.py`,
  `tests/test_stability_candidates.py`.
- **Script:** `examples/basic/flash_tp_modified_raoult_demo.py`.

---

## Case R-2: Modified-Raoult VLE - bubble, dew, split invariants and the azeotrope

- **Source:** The modified-Raoult equilibrium condition itself
  (`y_i P = x_i gamma_i Psat_i`, e.g. Koretsky 2012 ch. 8; Smith, Van Ness &
  Abbott ch. 10), evaluated independently in the test file. Literature value for
  the azeotrope is **commonly tabulated and was not read from a primary source
  in this work**; it is recorded as context, not as a validation anchor.
- **Location:** `tests/test_flash_modified_raoult.py::test_bubble_temperature_from_flash_verdicts_satisfies_the_scalar_equation`,
  `::test_dew_temperature_from_flash_verdicts_satisfies_the_scalar_equation`,
  `::test_the_verdict_boundary_sits_exactly_at_the_stability_tolerance`,
  `::test_vapor_liquid_split_satisfies_modified_raoults_law`,
  `::test_the_azeotrope_matches_an_independent_solve_and_is_a_vanishing_band`.
- **Assumptions:** `flash_tp` is treated as a black box returning a phase
  *count*. The temperature at which the count changes from 1 to 2 is the bubble
  point and from 2 to 1 (upwards) the dew point. What is asserted is the scalar
  identity at that temperature, computed by code written in the test file.
- **Components / units:** 1-Propanol(1) / Water(2), P = 101325 Pa, mole
  fractions. Antoine windows: 1-Propanol [285, 400] K, Water [284, 441] K, so
  the mixture window is [285, 400] K and is reported in diagnostics.
- **Parameters and provenance:** Tessier et al. (2000) Table 1 pair 1-3:
  tau_13 = -0.07149, tau_31 = 2.7425, alpha = 0.3. **LLE-fitted and temperature
  independent**; with these parameters the binary is fully miscible (checked:
  `stability_tp(..., vapor="none")` is `"stable"` at z1 = 0.1 ... 0.9, 330 K).
- **Expected outcome and results:**
  - **Verdict boundary identity.** At a stationary point `tpd = -ln(sum_i W_i)`
    and the first vapor substitution from a liquid feed gives
    `W_i = x_i gamma_i Psat_i / P`, so a feed is called unstable exactly when
    `sum_i x_i gamma_i Psat_i / P > exp(tpd_tol)`. Measured at the default
    `tpd_tol = 1e-8`: the bisected boundary gives a ratio of
    **1.000000010000000** (i.e. `exp(1e-8)` to 1e-11). The boundary tests
    therefore use `StabilitySettings(tpd_tol=1e-14)`.
  - **Bubble** (tpd_tol = 1e-14): x1 = 0.3 -> T = **360.976515439 K**, residual
    `sum x gamma Psat / P - 1 = 4.7e-15`; x1 = 0.5 -> T = **360.988638079 K**,
    residual 3.3e-15.
  - **Dew** (implicit x solved in the test by fixed point): y1 = 0.3 ->
    T = **364.328227469 K**, x1 = 0.04181356, residual
    `sum y P/(gamma Psat) - 1 = 2.6e-14`; y1 = 0.5 -> T = **361.611862978 K**,
    x1 = 0.63295515, residual 1.4e-14.
  - **Split invariants** at 361 K, z = (0.5, 0.5): x = (0.50549118, 0.49450882),
    y = (0.44232093, 0.55767907), beta = **0.08692668740410228**.
    `max_i |y_i P - x_i gamma_i Psat_i| / P = 1.6e-15`; mass balance **0.0**;
    `delta_g_split_rt = -2.004104e-05`; lever rule
    `(z1 - x1)/(y1 - x1) - beta = < 1e-12`; post-split `"stable"` on both
    phases. Stages: 32 successive substitutions + 1 second-order step.
  - **Incipient vapor equals the equilibrium vapor.** At `T_bubble + 1e-8 K`
    the stability minimizer reproduces `y(x)` from the independent bubble solve
    to **< 1e-8** (x1 = 0.3: 0.40599129; x1 = 0.4: 0.41626188).
  - **Azeotrope** (independent bisection of `y_1(x_1) - x_1` on the bubble
    curve, written in the test): x1 = **0.419874**, T = **360.917975 K =
    87.768 C**. Commonly tabulated (unverified here): ~87.7 C, x1 ~ 0.43.
    Deviation **+0.07 K** and **-0.012 in x1**. Reported, not tuned; the
    parameters are LLE-fitted and temperature independent, so agreement this
    close is partly fortuitous.
  - **Flash-level azeotrope signature:** at x1 = xa the two-phase band has zero
    width - `flash_tp` returns one `"liquid"` 0.01 K below and one `"vapor"`
    0.01 K above. The incipient vapor is richer in propanol at xa - 0.05 and
    leaner at xa + 0.05, bracketing the azeotrope from the package side.
- **Tolerance:** 1e-8 relative on both scalar equations (achieved 4.7e-15 and
  2.6e-14); 1e-10 on modified Raoult's law (achieved 1.6e-15); 1e-12 on mass
  balance (achieved 0.0).
- **Independent route:** scalar bubble/dew/azeotrope solves written in the test
  file, sharing no code with `chemthermo.flash`; plus `thermo` 0.6.0 in
  Case R-4.
- **Test path:** `tests/test_flash_modified_raoult.py`.
- **Script:** `examples/basic/flash_tp_modified_raoult_demo.py`.

---

## Case R-3: The three-phase neighbourhood of water / 1-butanol (negative controls)

- **Source:** Gibbs' phase rule (a binary at fixed pressure has three phases at
  a single temperature) plus the two simultaneous bubble equations that define
  it. Literature heteroazeotrope is **commonly tabulated, not read from a
  primary source in this work**.
- **Location:** `tests/test_flash_modified_raoult.py::test_three_phase_temperature_is_where_both_liquids_boil`,
  `::test_two_kelvin_below_t3_is_a_stable_liquid_liquid_split`,
  `::test_two_kelvin_above_t3_the_feed_has_evaporated`,
  `::test_at_t3_the_two_phase_answer_is_refused`, and
  `examples/validation/10_modified_raoult_water_butanol.py`.
- **Assumptions:** At T3 both conjugate liquids are at their bubble point
  simultaneously and share one vapor. The binodal is obtained from equal
  activities (no Antoine); T3 is then solved from the *water-rich* branch only
  and **checked** on the butanol-rich branch, which is a non-trivial identity.
- **Components / units:** n-Butanol(1) / Water(2), P = 101325 Pa, feed
  z = (0.20, 0.80). Mole fractions. T3 is inside both Antoine windows.
- **Parameters and provenance:** as Case R-1.
- **Expected outcome and results:**
  - Binodal (independent Newton, residual 1.1e-16): x1 = **0.019998419467** and
    **0.359999661508**.
  - **T3 = 366.213774 K (93.064 C).** Both liquids give
    `sum x gamma Psat / P = 0.999999999999999` (checked to 1e-10 on the branch
    that was not solved for). Shared vapor
    **y = (0.234063, 0.765937)**, sum = 1.000000000000.
    gamma^I = (30.55830, 1.00658), gamma^II = (1.69755, 1.54133).
    These reproduce the independent orchestrator reference (T3 = 366.2138 K,
    y = (0.23406, 0.76594), same gammas) to 1e-5 or better.
  - Commonly tabulated heteroazeotrope ~365.9 K, y(water) ~ 0.75-0.76:
    deviation **+0.31 K**, y(water) = **0.76594**. Partly fortuitous - the NRTL
    pair is LLE-fitted and temperature independent.
  - **(i) T3 - 2 K:** liquid-liquid split, tie-line equal to the binodal to
    1e-9, `equilibrium_residual = 7.7e-13`,
    `delta_g_split_rt = -8.580513e-03`, `post_split_stable = True`,
    `post_split_tpd_min = -1.354e-13`. Both liquids are stable against the
    vapor candidate.
  - **(ii) T3 + 2 K:** the feed has **fully evaporated** - a single `"vapor"`
    with `feed_branch == "vapor"` and `tpd_min = +3.393303e-02`. Verified
    independently: `sum z P/(gamma(x) Psat) = 0.966636 < 1`, i.e. the feed is
    above its dew point. (The brief allowed either this or a vapor-liquid pair;
    the model gives this.)
  - **(iii) T3 (to 1e-7 K):** `flash_tp` returns a vapor-liquid pair with
    vapor **(0.234063, 0.765937)** - the three-phase vapor to 1e-6 - and liquid
    on the water-rich binodal, `post_split_tpd_min = -1.4e-15`, i.e. *marginal*:
    at T3 the third phase lies exactly on the tangent plane, so it is not an
    instability. Which side of that knife edge a run lands on depends on how
    exactly T3 is known, so the test accepts either the marginal two-phase
    answer or the `ConvergenceError`, and asserts the numbers in each branch.
  - **(iii, continued) T3 - 0.01 K:** `flash_tp` **raises**
    `ConvergenceError("... a third phase is required ...")`. With
    `FlashSettings(post_split_stability=False)` the same call returns
    `['vapor', 'liquid']` with `post_split_status = "unstable"` and
    `phase_stability_tpd_min_vapor == phase_stability_tpd_min_liquid =
    -6.122901e-04` - identical because two coexisting phases share one tangent
    plane, so a third stationary point below it has the same tpd measured from
    either.
  - **Measured width of the refusal window: about 0.135 K below T3** for this
    feed (LLE returned at T3 - 0.14 K, refusal from T3 - 0.13 K to T3).
- **Finding, recorded not accommodated:** inside that window the *correct*
  answer is the two-liquid pair, not three phases. The solver seeds its split
  from the deepest tangent-plane minimum, which there is the vapor; the pair it
  converges is genuinely not the equilibrium, so refusing is right, but the
  message's diagnosis ("a third phase is required") is only half the story -
  resolving it needs phase addition **and removal** (the vapor amount would go
  to zero). No retry-from-the-second-minimum heuristic was added: that is
  `flash-vlle-phase-addition`. See ADR-0010 "Known limitation, measured".
- **Tolerance:** 1e-10 on the bubble equation at T3 (achieved 1e-15); 1e-5 on
  the shared vapor against the independent reference (achieved ~3e-6); 1e-9 on
  the tie-line below T3.
- **Independent route:** binodal and T3 solved in the test file and in the
  script, sharing no code with `chemthermo.flash`; cross-checked against the
  orchestrator's separately written reference.
- **Test path:** `tests/test_flash_modified_raoult.py`.
- **Script:** `examples/validation/10_modified_raoult_water_butanol.py`
  (prints the whole neighbourhood with PASS/FAIL).

---

## Case R-4: Modified-Raoult VLE against `thermo` 0.6.0

- **Source:** `thermo` 0.6.0 as an independent implementation:
  `GibbsExcessLiquid` (`use_Poynting=False`, `use_phis_sat=False`) over
  `thermo.NRTL`, an `IdealGas` vapor, and `FlashVL`.
- **Location:** `tests/validation/test_modified_raoult_vs_thermo.py`.
- **Assumptions:** The comparison is only meaningful if both packages evaluate
  the same `Psat_i(T)`, so that is established first. `thermo`'s
  `VaporPressure.add_correlation(model="Antoine", ..., base=e)` computes
  `base**(A - B/(T+C))` in **Pa**, while the chemthermo databank stores
  `ln(P/bar)`, so the only conversion needed is `A -> A + ln(1e5)`.
- **Components / units:** 1-Propanol(1) / Water(2), P = 101325 Pa, mole
  fractions; `thermo` IDs `1-propanol`, `water`.
- **Parameters and provenance:** as Case R-2; identical taus and alphas are
  passed to both packages, and chemthermo's own Antoine records are pushed into
  `thermo`.
- **Expected outcome and results:**
  - **Psat:** worst relative difference **2.3e-15** over 300, 330, 361, 380 and
    399 K (asserted 1e-10). The shared reference is established.
  - **Bubble / dew**, each package bisected on its own verdicts:
    z1 = 0.3 bubble 360.976515439 vs 360.976515464 (**2.6e-08 K**), dew
    364.328227469 vs 364.328227441 (**2.9e-08 K**); z1 = 0.5 bubble
    360.988638079 vs 360.988638104 (**2.6e-08 K**), dew 361.611862978 vs
    361.611862978 (**3.4e-13 K**).
  - **VL flash** at 361 K, z = (0.5, 0.5): |d beta| = **4.3e-06**,
    worst |dx| = **3.0e-07**, worst |dy| = **3.4e-13** (asserted 1e-5).
    Adjudicated with the modified-Raoult residual
    `max_i |ln(x_i gamma_i Psat_i / P) - ln y_i|`: chemthermo **1.6e-15**,
    `thermo` **2.0e-07**. Both satisfy their own material balance to round-off,
    so the difference is `thermo`'s convergence tolerance, not a model
    disagreement.
- **Recorded limitation (not worked around):** `thermo`'s
  vapor-fraction-specified flash on this phase pair is unusable -
  `flash(P=101325, zs=[0.5,0.5], VF=0)` returns **T = 1.8e5 K** and `VF=1`
  returns **T = 9.2e3 K**, for a feed whose true boundaries are 360.99 K and
  361.61 K. The T,P flashes on the same objects are correct, so this is a
  solver-path problem in `thermo`. Bubble and dew points from `thermo` are
  therefore obtained by bisecting its own T,P flashes, and the failure is
  asserted (`test_thermo_vapor_fraction_specified_flash_is_unusable_here`) so a
  future `thermo` that fixes it is noticed. Same spirit as Cases S-7 and L-2.
- **Tolerance:** 1e-10 relative on Psat (achieved 2.3e-15); 1e-5 K on bubble and
  dew (achieved 2.9e-08 K); 1e-5 on beta, x and y (achieved 4.3e-06).
- **Independent route:** `thermo` 0.6.0 (optional dependency; the module skips
  without it).
- **Test path:** `tests/validation/test_modified_raoult_vs_thermo.py`.
- **Script:** none (the `thermo` comparison lives in tests only).

---

## Case V-1: The ternary vapor-liquid-liquid tie-triangle

- **Source:** Internal invariant plus an independent solve. A ternary at fixed
  `(T, P)` holds three phases over a *region* of feeds (Gibbs' phase rule:
  F = 3 - 3 + 2 = 2, and fixing T and P still leaves the tie-triangle), whose
  corners are the two conjugate liquids and the vapor they share. **No
  experimental ternary VLLE data is used anywhere in this case.** This entry's
  "independent route" is a separately written Newton solve, stated as required
  by the ledger rules.
- **Location:** `tests/validation/test_vlle_water_propanol_butanol.py`
  (`::test_the_independent_tie_triangle_is_a_genuine_three_phase_state`,
  `::test_a_feed_inside_the_tie_triangle_returns_three_verified_phases`,
  `::test_the_three_phase_state_has_the_lowest_gibbs_energy`,
  `::test_three_phase_results_are_deterministic_and_permutation_invariant`) and
  `examples/validation/11_vlle_water_propanol_butanol.py`.
- **Assumptions:** modified Raoult (`f_i^0 = Psat_i(T)`, `phi^sat = 1`,
  Poynting = 1, ideal vapor). Six unknowns - the two liquid compositions - and
  six equations: three equal activities, two normalizations, and the bubble
  condition on liquid I only. That liquid II is *also* at its bubble point, and
  that both liquids give the same vapor, are consequences and are checked.
- **Components / units:** 1-Propanol(1) / n-Butanol(2) / Water(3),
  P = 101325 Pa, T = 363.0, 364.0 and 365.0 K. Mole fractions; `tpd` and
  `G/RT` dimensionless. All three Antoine windows cover these temperatures.
- **Parameters and provenance:** NRTL Table 1 of Tessier, Brennecke and
  Stadtherr, Chem. Eng. Sci. 55 (2000) 1785-1796 (attributed there to McDonald
  and Floudas, AIChE J. 41 (1995) 1798), fixture
  `tests/fixtures/nrtl/tessier2000_problem1.json`; alpha implied by
  `G_ij = exp(-alpha_ij tau_ij)`. Antoine from the packaged databank
  (Koretsky 2012). **LLE-fitted and temperature independent.**
- **Expected outcome and results:**
  - Tie-triangle vertices (independent Newton, residual <= 6.7e-16):

    | T [K] | x^I | x^II | y |
    |---|---|---|---|
    | 365.0 | (0.02361601, 0.02440696, 0.95197703) | (0.09873744, 0.20195471, 0.69930785) | (0.12724459, 0.14989592, 0.72285949) |
    | 364.0 | (0.05214424, 0.02917157, 0.91868419) | (0.13963944, 0.12347805, 0.73688251) | (0.21078146, 0.10017893, 0.68903962) |
    | 363.0 | (0.10282779, 0.03539032, 0.86178189) | (0.15630287, 0.06422717, 0.77946996) | (0.28312929, 0.06046985, 0.65640086) |

    These reproduce the orchestrator's separately written reference (same
    vertices to the 8 digits it recorded).
  - The *unsolved* equations hold: liquid II's bubble sum is
    1.000000000000 to <= 2e-15 at all three temperatures, and the vapor
    computed from liquid II equals the vapor computed from liquid I to
    <= 1e-15.
  - Six feeds inside the triangle (barycentric weights (1/3, 1/3, 1/3) and
    (0.2, 0.3, 0.5) at each temperature) return **three phases**
    `liquid1`/`liquid2`/`vapor`. Worst deviations over the six:
    compositions **2.2e-14**, phase fractions **1.6e-13**,
    `equilibrium_residual` **1.8e-15**, `mass_balance_residual` **6.9e-18**,
    post-split `tpd_min` **-4.1e-16**, every phase `"stable"`.
  - Independently reproduced with the orchestrator's own reference feeds at
    364.0 K: z = (0.134, 0.084, 0.782) -> (LI, LII, V) =
    (0.336081159, 0.329853965, 0.334064877) against the reference
    (0.336081, 0.329854, 0.334065); z = (0.10, 0.06, 0.84) ->
    (0.621821904, 0.170607495, 0.207570601); z = (0.16, 0.09, 0.75) ->
    (0.218035553, 0.227613768, 0.554350679). `delta_g_split_rt` for the first
    is **-0.002712735**, i.e. G3/RT = -0.69315706 against a single liquid at
    -0.69044432, the reference values; the other two feeds give
    **-0.001641130** and **-0.006438596**.
  - **Gibbs ordering**, computed in the test file: G3 < G2 < G1 at every one of
    the six feeds, e.g. at 364 K, centroid: **-0.693912758 < -0.693543323 <
    -0.691218474**. `delta_g_vs_two_phase_rt` ranges from -1.2e-05 (363 K,
    centroid) to -1.58e-03 (365 K, centroid) and is negative everywhere.
  - `phase_set_history` is `"L -> LV -> LLV"` or `"L -> LL -> LLV"` depending
    on which candidate the feed's deepest tangent-plane minimum is; both routes
    reach the same triangle.
  - **Invariants:** deterministic (identical diagnostics on a repeat run);
    permutation invariant over all 6 component orderings with worst restored
    composition difference **4.3e-15** and worst phase-fraction difference
    **1.5e-14**; the `vapor` phase stays the vapor under every permutation;
    `vapor_fraction == phase_fractions["vapor"]` exactly.
- **Tolerance:** asserted 1e-6 on compositions and phase fractions (achieved
  2.2e-14 and 1.6e-13); 1e-10 on the equilibrium residual (achieved 1.8e-15);
  1e-12 on mass balance (achieved 6.9e-18); 1e-9 on permutation invariance
  (achieved 1.5e-14).
  - **External cross-check: attempted, and it does not exist.** Case L-2 had
    found `thermo` 0.6.0 collapsing two `GibbsExcessLiquid` phases into one for
    a *binary*; that was re-checked here directly for the three-phase state,
    with identical NRTL parameters and chemthermo's own Antoine records pushed
    into `thermo` so `Psat` is shared.
    `FlashVLN(..., liquids=[liquid])` reports `unique_liquid_count == 1` and
    returns a **two-phase** answer at the 364 K centroid feed, which is inside
    the tie-triangle: vapor (0.21027032, 0.09978016, 0.68994952), liquid
    (0.09928882, 0.07716435, 0.82354683), betas (0.31446284, 0.68553716). That
    liquid is neither conjugate liquid. `FlashVLN(..., liquids=[liquid, liquid])`
    still reports `unique_liquid_count == 1` and raises
    `TypeError: 'NoneType' object is not subscriptable`.
    **Adjudicated by Gibbs energy**, computed in the test file: chemthermo's
    three-phase answer **-0.693912758** against `thermo`'s two-phase
    **-0.693543321**. And the *intermediate* is externally confirmed:
    `thermo`'s two-phase answer has the same Gibbs energy as the two-phase
    candidate chemthermo converged before adding the third phase, to
    **2.4e-09** - so the disagreement is about the phase count, not about the
    two-phase thermodynamics. The failure is asserted
    (`::test_thermo_cannot_hold_two_distinct_excess_gibbs_liquids`) so a future
    `thermo` that fixes it is noticed.
- **Independent route:** the six-equation Newton solve written in the test file
  and again in the script, sharing no code with `chemthermo.flash`; plus the
  Gibbs-energy ordering computed there. **No external VLLE reference is
  available for this system** (see the bullet above) and none is claimed.
- **Test path:** `tests/validation/test_vlle_water_propanol_butanol.py`,
  `tests/test_flash_vlle.py`.
- **Script:** `examples/validation/11_vlle_water_propanol_butanol.py`,
  `examples/basic/flash_tp_vlle_demo.py`.

---

## Case V-2: Feeds outside the tie-triangle, and one the stability test misses

- **Source:** Internal invariants, each verified by a different independent
  route (a liquid-liquid stability test, a bubble-point sum, a dew-point sum).
  Stated as required by the ledger rules.
- **Location:** `tests/validation/test_vlle_water_propanol_butanol.py`
  (`::test_a_feed_in_the_vapor_liquid_region_returns_two_phases`,
  `::test_a_feed_in_the_liquid_liquid_region_returns_two_liquids`,
  `::test_the_water_rich_corner_is_a_single_liquid`,
  `::test_a_superheated_feed_is_a_single_vapor`,
  `::test_a_thin_tie_triangle_can_hide_from_the_deterministic_trial_set`),
  section 2 of `examples/validation/11_vlle_water_propanol_butanol.py`.
- **Assumptions:** as Case V-1. "Outside" is decided by the barycentric weights
  of the feed in the *independent* tie-triangle, not by what the flash returns.
- **Components / units:** as Case V-1; T = 364.0 K except the superheated
  control at 380.0 K (inside every Antoine window: 1-Propanol [285, 400] K,
  n-Butanol [288, 404] K, Water [284, 441] K).
- **Parameters and provenance:** as Case V-1.
- **Expected outcome and results:**
  - **(a) Vapor-liquid region**, weights (-0.25, 0.45, 0.80),
    z = (0.21842685, 0.12841537, 0.65315777): `['vapor', 'liquid']`,
    `vapor_fraction = 0.71129282`, vapor (0.22553936, 0.10455392, 0.66990672),
    liquid (0.20090366, 0.18720324, 0.61189310). Independent check: that liquid
    on its own is `"stable"` against a second liquid (`stability_tp` with
    `vapor="none"`, `tpd_min = 0.0`), so no third phase exists there.
    Post-split `tpd_min = -2.4e-16`.
  - **(b) Liquid-liquid region**, weights (0.55, 0.60, -0.15),
    z = (0.08084578, 0.07510435, 0.84404987): `['liquid1', 'liquid2']`,
    `vapor_fraction is None`, fractions 0.41992979 / 0.58007021, compositions
    (0.13181652, 0.14030670, 0.72787677) and
    (0.04394656, 0.02790247, 0.92815097). Independent check: **both** liquids
    have `sum_i x_i gamma_i Psat_i / P = 0.990647948 < 1`, i.e. neither can
    boil. Per-phase post-split `tpd_min` +1.03e-13 and -1.27e-14, both
    `"stable"`.
  - **(c) Water-rich corner**, z = (0.01, 0.005, 0.985): a single `"liquid"`,
    `feed_branch = "liquid"`, `tpd_min = 0.0`; independently
    `sum_i z_i gamma_i Psat_i / P = 0.859278566 < 1`.
  - **(d) Superheated**, 380.0 K, z = (0.20, 0.15, 0.65): a single `"vapor"`,
    `feed_branch = "vapor"`, `vapor_fraction = 1.0`; independently the dew sum
    `sum_i z_i P / (gamma_i(x) Psat_i) = 0.574785505 < 1`, i.e. above the dew
    point.
  - None of (a)-(d) reports `phase_set_history`: the search was never entered,
    which is the correct signal that nothing was added or removed.
  - **Recorded limitation, not accommodated.** At 363.0 K the two liquid
    vertices differ by only 0.053 in x_1 (near the plait point). For the feed
    at barycentric weights (0.5, 0.3, 0.2), z = (0.15493061, 0.04905728,
    0.79601210) - which **is** inside the triangle - every trial of the
    deterministic stability set collapses onto the trivial solution,
    `tpd_min = 0.0`, and `flash_tp` returns a single `"liquid"`. That is wrong
    for this model. The failure is in the *stability* test, not in the phase
    search, which is never entered. Pinned by test so that a future improvement
    to the trial set is noticed rather than silently absorbed.
- **Tolerance:** verdicts are exact (phase names and counts); the independent
  bubble and dew sums are compared against 1 with no tolerance needed
  (0.859-0.991 and 0.575).
- **Independent route:** a liquid-liquid `stability_tp` call for (a); bubble
  and dew sums written in the test file for (b), (c) and (d).
- **Test path:** `tests/validation/test_vlle_water_propanol_butanol.py`.
- **Script:** `examples/validation/11_vlle_water_propanol_butanol.py`.

---

## Case V-3: The binary refusal window, resolved by removing a phase

- **Source:** Gibbs' phase rule (a binary at fixed pressure has three phases at
  a single temperature only, F = 2 - 3 + 2 = 1) plus the binodal and T3 of
  Case R-3. This entry continues Case R-3, whose "Finding, recorded not
  accommodated" this case discharges.
- **Location:** `tests/test_flash_vlle.py`
  (`::test_below_t3_the_window_resolves_to_the_two_liquids`,
  `::test_just_above_t3_the_answer_is_a_vapor_liquid_pair`,
  `::test_the_window_answer_is_the_same_as_a_direct_liquid_liquid_flash`,
  `::test_max_phases_two_reproduces_the_pre_adr_0011_refusal`,
  `::test_post_split_stability_false_still_returns_the_two_phase_pair`),
  `tests/test_flash_modified_raoult.py::test_at_t3_the_two_phase_answer_is_refused`,
  section 3 of `examples/validation/11_vlle_water_propanol_butanol.py`, and
  section 5 of `examples/validation/10_modified_raoult_water_butanol.py`.
- **Assumptions:** as Case R-3. The binodal of this NRTL pair is temperature
  independent (the fitted tau are), so the two-liquid answer below T3 is the
  *same* tie-line at every temperature and the lever rule fixes the amounts.
- **Components / units:** n-Butanol(1) / Water(2), P = 101325 Pa,
  z = (0.20, 0.80), T = T3 - 0.05 K, T3 - 0.10 K and T3 + 0.05 K with
  T3 = 366.213774 K.
- **Parameters and provenance:** as Case R-1 / R-3.
- **Expected outcome and results:**
  - **Below T3 (both -0.05 K and -0.10 K):** `['liquid1', 'liquid2']`,
    `phase_regime = "LLE"`, `vapor_fraction is None`, tie-line
    x1 = **0.019998419467 / 0.359999661508** - the Case L-3 / R-1 / R-3
    binodal - with worst deviation **3.3e-13**. Phase fractions
    **0.470585520652 / 0.529414479348**, equal to the lever rule to
    **0.0**. `equilibrium_residual` 8.9e-16 and 3.1e-15,
    `mass_balance_residual` 0.0, `delta_g_split_rt` -5.269806e-03 and
    -7.185443e-03, both phases post-split `"stable"`.
  - **The route:** `phase_set_history = "V -> LV -> LLV -> LL"`,
    `phases_added = 1`, `phases_removed = 1`. The feed's own lowest-Gibbs
    candidate here is the *vapor*, so the search starts from `V`; the
    vapor-liquid pair is unstable towards a second liquid, which is added; the
    three-phase Rachford-Rice for a binary away from T3 has no finite solution,
    the vapor's amount runs negative, and it is removed. **This is the case
    that needs removal and not only addition.**
  - **Same answer, different route:** the tie-line reached through the window
    equals the one a direct `gamma-gamma` flash returns at the same state to
    **1.8e-12**, which is the two routes' own convergence tolerances rather
    than round-off.
  - **Above T3 (+0.05 K):** `['vapor', 'liquid']`,
    `vapor_fraction = 0.847371275`, vapor (0.23249827, 0.76750173), liquid
    (0.01957460, 0.98042540). Verified independently in the test file: the
    liquid is exactly at its bubble point
    (`sum x gamma Psat / P = 1.000000000000`, |deviation| < 1e-12) and modified
    Raoult's law holds to **5.8e-15**; the liquid is also `"stable"` against a
    second liquid. (Case R-3(ii) recorded a *single vapor* two kelvin above T3;
    0.05 K above it the feed is still between its bubble and dew points, so a
    pair is what the physics gives.)
  - **Regression of the documented behavior:** `FlashSettings(max_phases=2)` at
    T3 - 0.05 K still raises
    `ConvergenceError("... a third phase is required ...")`, and so does
    `max_phases=1`. `FlashSettings(post_split_stability=False)` still returns
    the vapor-liquid pair with `post_split_status = "unstable"` and no
    `phase_set_history` key.
- **Tolerance:** 1e-8 on the tie-line against the independent binodal
  (achieved 3.3e-13); 1e-10 on the lever rule (achieved 0.0); 1e-10 on the
  bubble equation and on modified Raoult's law above T3 (achieved 5.8e-15);
  1e-9 between the two routes to the same tie-line (achieved 1.8e-12).
- **Independent route:** the binodal and T3 of Case R-3, solved in
  `tests/test_flash_modified_raoult.py` and in
  `examples/validation/10_modified_raoult_water_butanol.py` with no
  `chemthermo.flash` code; the bubble and Raoult checks written in
  `tests/test_flash_vlle.py`; and the `gamma-gamma` path as a second route to
  the same tie-line.
- **Test path:** `tests/test_flash_vlle.py`,
  `tests/test_flash_modified_raoult.py`.
- **Script:** `examples/validation/11_vlle_water_propanol_butanol.py`,
  `examples/validation/10_modified_raoult_water_butanol.py`.

---

## Case V-4: Multiphase Rachford-Rice against Okuno et al. (2010) Table 1

- **Source:** R. Okuno, R. T. Johns and K. Sepehrnoori, "A new algorithm for
  Rachford-Rice for multiphase compositional simulation", SPE Journal 15 (2010)
  313-325. Table 1 gives four overall compositions with constant K-values; the
  text prints the solution of Example 3 as `(beta_1, beta_2) = (0.87, 2.2e-6)`
  and gives Example 4 as three phase compositions `x_ij` for a feed outside the
  tie-triangle. Open-access PDF text used:
  `scratchpad/okuno2010_multiphase_rr.txt`, Table 1 and the "Comparisons in
  Standalone Calculations" section.
- **Location:** `tests/test_multiphase_rr.py`.
- **Assumptions:** the Rachford-Rice stage takes `(z, K)` and returns `beta`
  and never calls a thermodynamic model, so constant-K data isolates exactly
  this solver. The reference phase is the paper's phase 3, so the two K columns
  are phases 1 and 2 against it.
- **Components / units:** 7 components for Examples 1-3, 3 for Example 4; mole
  fractions, dimensionless K.
- **Parameters and provenance:** the printed Table 1 values verbatim.
- **Expected outcome and results:**
  - **Example 1:** beta = (0.68683289, 0.06019424), reference phase
    0.25297286451, 8 Newton iterations, scaled residual 3.1e-16.
  - **Example 2** (the case the paper's reference root-finder cannot solve at
    all): beta = (0.46945316, 0.47024452), reference 0.06030232018,
    7 iterations, residual 1.3e-14.
  - **Example 3** (next to a critical end point): beta =
    **(0.8701633569, 2.180303e-06)**, reference 0.1298344631, **4 iterations**,
    residual 1.7e-15. The paper prints **(0.87, 2.2e-6)**; agreement 1.9e-04
    absolute on beta_1 and 0.9% on beta_2 against the two printed digits. The
    paper reports 4 iterations for its own algorithm on this example.
  - **Example 4** (a *negative* flash near a critical end point): the printed
    K-values are the printed compositions against phase 3 to **< 1e-8**, so the
    phase fractions are the exact solution of `z = sum_j beta_j x^j`:
    **(1.2, 14.66, -14.86)**. Achieved
    **(1.2000000000, 14.6599999, -14.8599999033)** in 5 iterations, residual
    2.3e-13. One fraction is negative, which is the signal `flash_tp` uses to
    *remove* a phase.
  - **Iteration counts differ from the paper** for Examples 1 and 2 (8 and 7
    here against 4 and 4 printed). The initial estimate and the line search are
    reimplemented from the description, the stopping test here is a *scaled*
    residual, and the paper's tolerance is 1e-8 against 1e-12 here. The
    solutions agree; the counts are recorded, not matched.
  - **Convexity, checked not assumed:** the Hessian
    `A^T diag(z / t^2) A` is positive definite at all four solutions, and `F`
    does not decrease along any of eight deterministic feasible perturbations
    at each.
  - **Feasible region:** every point of `S` tested satisfies
    `t_i >= max(z_i, max_j K_i^j z_i)`, i.e. `S` contains no pole.
  - **Constructed negative flash:** three chosen phase compositions and weights
    (0.8, 0.45, -0.25) are recovered to **1e-10**.
  - **No-solution control:** the `(z, K)` the flash reaches for the binary
    three-phase set at T3 - 0.05 K (Case V-3) is reported as a recession, with
    the reference phase (the vapor) named as the one whose amount runs to minus
    infinity - not as a convergence failure.
  - **Degenerate control:** `K^j = 1` makes phase `j` the reference phase; `F`
    is then flat in `beta_j`, the Hessian is singular (smallest eigenvalue
    0 to 1e-14), and the duplicate is removed upstream by the trivial-solution
    metric rather than guessed at here.
- **Tolerance:** Rachford-Rice equations asserted < 1e-11 (achieved 1.5e-14);
  material balance < 1e-12; Example 3 against the printed digits at abs 5e-3
  and rel 5e-2; Example 4 against the exact linear solve at abs 1e-6.
- **Independent route:** published constant-K data plus, for Example 4, the
  exact 3x3 linear solve from the paper's own printed compositions. The
  Rachford-Rice equations, the objective and the Hessian are all rewritten in
  the test file rather than imported.
- **Test path:** `tests/test_multiphase_rr.py`.
- **Script:** none (a solver-level unit test; the flash-level scripts are
  Cases V-1 and V-3).
