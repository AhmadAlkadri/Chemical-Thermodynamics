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

## Table of contents: case families

One family letter per capability area; case numbers are listed in the order
they appear below, not necessarily numeric order (an amended case keeps its
original number). "Slice (ADR)" names the slice that created the family, in
delivery order; later slices that amended a case in place are named in that
case's own entry.

| Family | Cases | Capability | Slice (ADR) |
| --- | --- | --- | --- |
| **S** | S-1..S-4, S-5 (not found), S-6..S-8 | Michelsen tangent-plane stability (`stability_tp`): the feed-tangent identity, stationarity/tm/tpd relations, marginal stability, cross-checks against `thermo` and against Tessier (2000)'s published NRTL stationary points | `stability-tpd-pr` (ADR-0005/0007), `stability-tpd-nrtl` (ADR-0007) |
| **K** | K-1 | Peng-Robinson per-pair `kij` (scalar or name-keyed mapping), the fixed diagonal-corruption bug | `pr-kij-matrix` (ADR-0006) |
| **N** | N-1..N-3 | NRTL activity coefficients: Gibbs-Duhem consistency, cross-check against an independent implementation, reproduction of Tessier (2000) Table 2 stationary points | `nrtl-gibbs-duhem-fix` (no ADR; signature unchanged) |
| **F** | F-1..F-5 | `flash_tp` phi-phi phase detection and split robustness: verdicts against `thermo`'s `FlashVL`, the Wilson-vs-tangent-plane disagreement, split verification invariants, the negative-flash/second-order stage, phase identity by compressibility | `flash-auto-phase-detection` (ADR-0008), `flash-phi-phi-second-order` (ADR-0016), `flash-phase-labels-by-compressibility` (ADR-0017) |
| **L** | L-1..L-4 | `gamma-gamma` liquid-liquid `flash_tp` (activity model, no EOS): Tessier (2000) tie-lines, binodal/lever-rule cross-check, post-split stability of every two-phase phi-phi flash | `flash-lle-activity` (ADR-0009) |
| **R** | R-1..R-4 | `modified-raoult` mode (activity liquid + ideal vapor): VLE/LLE from one tangent plane, the three-phase neighbourhood refusal window, cross-check against `thermo` | `flash-modified-raoult` (ADR-0010) |
| **V** | V-1..V-5 | Three-phase (VLLE) discovery: the ternary tie-triangle, stability misses fixed by fixed candidate surfaces, the binary refusal window, multiphase Rachford-Rice against Okuno et al. (2010), the ternary verdict map | `flash-vlle-phase-addition` (ADR-0011), `stability-candidate-surfaces` (ADR-0012) |
| **P** | P-0..P-18 | PC-SAFT: residual Helmholtz/derivatives against teqp, density roots, flash/stability integration, split robustness re-measured on the PC-SAFT grid, phase identity, association against FeOs, per-phase density roots (LLE), EOS phase addition/removal (VLLE), fixed EOS stability surfaces, polymer components, the polymer liquid-liquid split, the polymer vapour-liquid split in log mole numbers, the log-space stability normalization, the curvature safeguard and the Rachford-Rice denominator that close the 0.3-3.6 MPa sweep, the seed ladder that retires the map's polymer refusals (P-17), and the reversible removal plus multiphase log-space stage that retire its nine `eos-three-phase` refusals (P-18) | `pcsaft-residual-helmholtz` (ADR-0014), `pcsaft-density-roots-flash` (ADR-0015), `pcsaft-association` (ADR-0018), `flash-eos-per-phase-roots` (ADR-0019), `flash-phase-addition-eos` (ADR-0020), `stability-eos-root-surfaces` (ADR-0021), `pcsaft-polymer-solvent` (ADR-0022), `pcsaft-polymer-vle` (ADR-0024), `stability-log-space-sums` (ADR-0025), `flash-polymer-edge-cases` (ADR-0026), `flash-polymer-robustness` (ADR-0028), `flash-eos-multiphase-robustness` (ADR-0029) |
| **B** | B-1, B-2 | The `chemthermo.bench` measurement harness and its acceptance rule (baseline + after + identical result hashes + measured ratio) (B-1), and the call-local equation-of-state solve memo plus the structured multiphase Hessian that was measured and declined (B-2) | `perf-baseline-and-root-reuse` (ADR-0023), `perf-eos-solve-cache` (ADR-0030) |
| **R-MAP** | R-MAP-1, R-MAP-2 | The `chemthermo.bench robustness` coverage map: 2110 states over the original six model families (R-MAP-1), grown to 2505 states over ten families - EOS three-phase windows, the legacy `gamma-phi` path, near-critical Peng-Robinson states and associating ternaries (R-MAP-2) - classified into a phase verdict, an invariant violation or one of eight refusal classes, with the ranked classes and their diagnoses. Both are amended in place as later slices retire what they found: R-MAP-2 by ADR-0029 (Case P-18), which leaves 5 refusals, all of them the deprecated `gamma-phi` path's | `robustness-map` (ADR-0027), `robustness-map-coverage` (ADR-0027 amendment) |

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
  - **AMENDED (slice `stability-candidate-surfaces`, ADR-0012): the recorded
    miss is resolved.** As first recorded: at 363.0 K the two liquid vertices
    differ by only 0.053 in x_1 (near the plait point), and for the feed at
    barycentric weights (0.5, 0.3, 0.2), z = (0.15493061, 0.04905728,
    0.79601210) - which **is** inside the triangle - every trial of the
    deterministic stability set collapsed onto the trivial solution,
    `tpd_min = 0.0`, and `flash_tp` returned a single `"liquid"`, wrong for
    this model.

    The cause was diagnosed as candidate switching *inside* a trial, not a
    reachability failure of the trial set. At that feed the lowest-Gibbs
    candidate is the liquid, and the reduced tangent-plane distance is
    **+8.07e-04** at x^I, **+1.33e-03** at x^II and **-9.92e-03** at the
    equilibrium vapor y, so the feed was provably unstable. On the vapor
    surface the stationary point is
    w = (0.30324249, 0.05009844, 0.64665907) with
    `tpd = -ln(sum W) = -1.1680295426e-02`, reachable in **one** successive
    substitution because the ideal-gas term is zero (`ln W_i = d_i`). The old
    iteration never got there: from the Raoult-vapor start the liquid
    candidate has the lower Gibbs energy at the intermediate compositions, so
    the update used the liquid terms and the iterate was dragged onto the
    liquid surface and onto the trivial solution.

    After ADR-0012 (one fixed candidate surface per trial, same starting
    points, same trial count): `stability_tp(..., vapor="ideal")` reports
    **unstable**, `tpd_min = -1.1680295426e-02`, `feed_branch = "liquid"`,
    `phase_branch = "vapor"`, `minimizing_trial = "raoult-vapor"`,
    `minimizing_trial_surface = "vapor"`, `trial_surfaces = "vapor:1,liquid:4"`;
    the vapor trial converges in **2** successive-substitution iterations with
    residual exactly 0.0. `flash_tp` then returns
    `['liquid1', 'liquid2', 'vapor']` with compositions equal to the tie-triangle
    above and fractions (0.5, 0.3, 0.2), both to **< 1e-8**. Pinned in
    `tests/validation/test_vlle_water_propanol_butanol.py::test_the_thin_tie_triangle_at_363_k_is_found`
    and `tests/test_flash_vlle.py::test_the_near_plait_ternary_feed_at_363_k_is_a_three_phase_state`.
  - **What remains at that feed, recorded not accommodated.** The `pure-Water`
    liquid-surface trial still does not converge: `converged = False`,
    `tpd = nan` (the documented sentinel for a failed trial, not a number that
    blew up - every iterate stays finite), `termination_reason =
    "second_order_no_progress"`, residual **5.4e-04** at
    w = (0.13762, 0.04163, 0.82075). The cause is the plait point, not
    arithmetic: the stationarity Jacobian `dg/d(ln W)` there has eigenvalues
    {**3.26e-09**, 1.106, 1.000}, condition number **8.2e+08**, so the Newton
    direction is dominated by the near-null direction and the line search
    cannot reduce the residual. Successive substitution on the same surface
    *does* reach the trivial solution (residual 0.0) but needs ~1.2e+04
    iterations, against `StabilitySettings.max_iter = 300`. The verdict does
    not depend on it: 4 of the 5 trials converge and the instability is found
    on the vapor surface. Pinned in
    `tests/test_stability_candidates.py::test_the_pure_water_trial_stalls_near_the_plait_point`.
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

---

## Case V-5: The verdict map of the ternary VLLE region

- **Source:** none external. This is a self-adjudicated case: the reference is
  an independent lowest-Gibbs classifier written in the test and example files,
  which share no code with `chemthermo.flash`. `thermo` 0.6.0 cannot hold two
  distinct excess-Gibbs liquids over one model (Case L-2, re-checked directly
  for this ternary in Case V-1), so no external three-phase reference exists
  and none is claimed.
- **Location:** `tests/validation/test_vlle_verdict_map.py` and
  `examples/validation/12_vlle_verdict_map.py`.
- **Assumptions:** modified Raoult (`f_i^0 = Psat_i(T)`, `phi^sat = 1`,
  Poynting = 1, ideal vapor), as Case V-1. The adjudication rule is that the
  equilibrium state is the admissible state of **least Gibbs energy**; equal
  fugacities alone are satisfied by more than one of the candidates below.
- **Components / units:** 1-Propanol(1) / n-Butanol(2) / Water(3),
  P = 101325 Pa, T = 363.0, 364.0 and 365.0 K. Mole fractions; `G/RT`
  dimensionless.
- **Parameters and provenance:** as Case V-1 (Tessier et al. 2000 Table 1,
  fixture `tests/fixtures/nrtl/tessier2000_problem1.json`; Antoine from the
  packaged databank). **LLE-fitted, temperature independent, no experimental
  ternary VLLE data is used.**
- **Method:** for each feed, every state the model admits is built
  independently and scored:
  1. the tie-triangle (the six-equation Newton of Case V-1), admissible when
     the feed's barycentric weights in it are all positive;
  2. a vapor-liquid state, from a four-equation Newton solve of
     `z_i - (1 - b) x_i - b K_i(x) x_i = 0`, `sum_i x_i = 1`, with
     `K_i = gamma_i(x) Psat_i / P`, from five deterministic starts;
  3. a liquid-liquid state, from a seven-equation Newton solve of the three
     equal activities, two normalizations and two mass balances, from two
     deterministic starts;
  4. the single-phase state: the feed on its lowest-Gibbs candidate.

  The lowest-Gibbs admissible state is the **expected** answer, compared feed
  by feed against `flash_tp(..., flash_mode="modified-raoult")`.
- **Grid:** per temperature, 21 feeds strictly inside the triangle (barycentric
  lattice, all weights >= 1/8) plus every point of a 1/12 mole-fraction lattice
  lying outside it by more than 1e-3 in barycentric coordinates - 55, 54 and 55
  feeds respectively, so **76 / 75 / 76** feeds in total. Deterministic.
- **Expected outcome and results:**
  - **Zero disagreements at all three temperatures.** Confusion matrices
    (rows expected, columns obtained; all off-diagonal entries are 0):

    | T [K] | (1, 1) | (2, 2) | (3, 3) | feeds |
    |---|---|---|---|---|
    | 363.0 | 48 | 7 | 21 | 76 |
    | 364.0 | 44 | 10 | 21 | 75 |
    | 365.0 | 41 | 14 | 21 | 76 |

    Every one of the 21 inside-feeds per temperature is three-phase and every
    outside-feed is one or two phases, so the geometric statement "inside the
    tie-triangle means three phases" is a *result* here, not an input.
  - **Three-phase answers are the same triangle whatever the feed.** Worst
    composition deviation from the independent tie-triangle over all 63
    three-phase feeds: **9.50e-12** (363 K), 2.88e-12 (364 K), 1.03e-12
    (365 K). Worst phase-fraction deviation from the feed's barycentric
    weights: **1.63e-11**, 7.45e-12, 3.10e-12. Phase names are always
    `liquid1`/`liquid2`/`vapor`.
  - **Two-phase answers** satisfy `ln x_i + t_i` equal across the two phases to
    **6.22e-15** (363 K), 3.78e-13 (364 K), 8.44e-15 (365 K), and mass balance
    `sum_j beta_j x^j = z` to **1.11e-16** or better.
  - **Single-phase answers** return the feed composition itself (to < 1e-12)
    on the candidate of lower Gibbs energy, with
    `diagnostics["phase_regime"] == "single-phase"`.
  - **Two feeds needed the multi-start classifier, not the flash.** At 363 K,
    z = (0.08333333, 0.16666667, 0.75) is liquid-liquid
    (`flash_tp` G/RT = -0.8091073 against -0.8084371 for a single liquid) and
    at 365 K, z = (0.16666667, 0.08333333, 0.75) is vapor-liquid
    (-0.7234793 against -0.7214637 for a single vapor). A single-seed
    independent solver missed both and would have recorded two false
    disagreements; the seeds were widened, and the flash was right both times.
- **Tolerance:** verdicts are exact (phase counts); 1e-8 asserted on
  three-phase compositions and phase fractions (achieved 9.5e-12 and 1.6e-11);
  1e-8 on the two-phase equilibrium residual (achieved 3.8e-13); 1e-10 on mass
  balance (achieved 1.1e-16); 1e-12 on the single-phase composition;
  1e-14 on the tie-triangle Newton residual (achieved 1.1e-15).
- **Independent route:** the four solvers listed under **Method**, all written
  in the test and example files with forward-difference Jacobians and their own
  line searches, plus the NRTL equations and the Antoine form rewritten there.
- **Test path:** `tests/validation/test_vlle_verdict_map.py`.
- **Script:** `examples/validation/12_vlle_verdict_map.py`.
- **Side measurement taken on the same grid (ADR-0012 diagnostics).** Over the
  227 feeds, the minimizing stability trial is `raoult-vapor` on the vapor
  surface for 193 of them, a `pure-<name>` trial on the liquid surface for 26,
  and `raoult-liquid` on the liquid surface for 8. In 99 of the 193 vapor-surface
  cases the *lowest-Gibbs* candidate at the converged point is the liquid, so
  `surface != phase_branch` there and `tpd` is correctly taken from the liquid;
  `|tpd_from_sum_W - tpd_min|` is then non-zero (worst 0.370) but only ever on a
  **stable** verdict, whose `tpd_min` is a positive number that decides nothing.
  Restricted to the unstable verdicts - the ones that seed a split - the two
  agree to **3.0e-12**, so equation (7) still holds where it is load bearing.
  The only unconverged trials anywhere on the grid are `pure-Water` on the
  liquid surface: 7 with `second_order_no_progress` and 3 with
  `second_order_max_iter`, all near the plait point, none changing a verdict.
- **Honesty note:** this case is evidence that the phase-count verdict is right
  *on this grid, for this model*. It is not a proof of global correctness, and
  it is not a comparison against measurement.

---

## Case P-0: PC-SAFT invariants, constants, and the Eq. (A.11) typo

- **Source:** J. Gross and G. Sadowski, "Perturbed-Chain SAFT: An Equation of
  State Based on a Perturbation Theory for Chain Molecules", Ind. Eng. Chem.
  Res. 40 (2001) 1244-1260 (DOI 10.1021/ie0003887). Erratum on Eq. (A.11)
  stated by NIST TRC, https://trc.nist.gov/TDE/TDE_Help/eos-PC-SAFT.htm.
- **Location:** Appendix equations (A.3)-(A.19) and (A.31)-(A.35); the paper's
  Table 1 (both the 21 + 21 universal constants and the pure-component
  parameters). Restated in the module docstring of
  `src/chemthermo/eos/pcsaft.py` with the index conventions.
- **Assumptions:** Non-associating, non-polar. Hard-chain plus dispersion only.
  Fixed `(T, molar density, x)`; no density root solving.
- **Components / units:** Methane, n-Hexane, n-Decane, Nitrogen and their
  binaries/ternaries. `sigma` in Angstrom, `epsilon/k` in K, densities in
  mol/m^3, volumes in m^3/mol, `k_B = 1.380649e-23` J/K,
  `N_A = 6.02214076e23` /mol (exact SI).
- **Parameters and provenance:** **The primary source was NOT read.** ACS
  returns HTTP 403 for `pubs.acs.org/doi/10.1021/ie0003887` from this
  environment and no open copy was located. Both tables were therefore taken
  from independent secondary sources that agree digit for digit:
  - *42 universal constants*: teqp (NIST, MIT) `src/data/PCSAFT.cpp`, namespace
    `teqp::saft::PCSAFT::PCSAFTMatrices::GrossSadowski2001`; and the table in
    the Wikipedia article "PC-SAFT", section "Dispersion Term". Checked
    entry-by-entry: **0 differences in 42 values**.
  - *11 pure-component parameter triples*: FeOs `parameters/pcsaft/gross2001.json`
    and Clapeyron.jl `database/SAFT/PCSAFT/PCSAFT_like.csv` (rows whose `source`
    column is the DOI of this paper). **0 differences.** teqp's built-in
    `PCSAFTLibrary` independently confirms Methane / Ethane / Propane with the
    BibTeXKey `Gross-IECR-2001`.
- **Expected outcome:**
  - `C1` is the **reciprocal** of the bracket of Eq. (A.11) (the printed
    equation drops the outer `-1` on its right-hand side); `C1 -> 1` as
    `eta -> 0`, and `0 < C1 < 1` at any `eta > 0` with `mbar >= 1`.
  - Ideal-gas limit: `A^res/RT -> 0`, `Z -> 1`, `ln phi_i -> 0` as `rho -> 0`.
  - Pure-component limit: the mixture code equals independently written pure
    formulas (Carnahan-Starling `a_hs`, `g = (1 - eta/2)/(1 - eta)^3`).
  - `Z = 1 + rho (d a_res/d rho)`, `mu_i^res/RT = d(n a_res)/dn_i|_{T,V}`.
  - Euler identity `sum_i x_i ln phi_i = A^res/RT + Z - 1 - ln Z`.
  - Gibbs-Duhem **at fixed (T, rho)**: `sum_i x_i d ln phi_i = (1 - 1/Z) dZ`.
    (The familiar `= 0` form holds at fixed T **and P**; these `ln phi` are
    evaluated at fixed density, where `P` moves with composition, so the right
    identity is derived in the test docstring rather than assumed.)
- **Tolerance achieved** (asserted -> measured on the states in the test file):
  - `C1 * bracket - 1`: 1e-14 -> exact to round-off.
  - `dC1/deta` vs central difference (step 1e-7): 1e-6 rel -> **8.7e-11**.
    `dC1/dmbar` likewise: 1e-6 rel -> **3.9e-08**.
  - Hard-chain and dispersion density derivatives, separately, vs central
    differences (step `rho * 1e-6`): 1e-7 rel, both pass.
  - Ideal-gas limit at `rho = 1e-6` mol/m^3: `< 1e-8` -> **1.3e-09** for
    `|A^res/RT|`, `|Z - 1|` and `max |ln phi|` alike (they coincide there).
  - Pure-component formulas: 1e-13 rel -> **2.2e-14** worst over 4 fluids
    x 3 states.
  - `Z` vs FD of `A^res` in density: 1e-8 rel -> **5.9e-09**
    (finite-difference limited).
  - `mu_i^res` vs FD of `n A^res` in mole numbers at fixed `(T, V)`:
    1e-7 rel -> **1.4e-08** (finite-difference limited).
  - Euler identity: 1e-13 -> **2.0e-15** here, **8.9e-16** worst over the
    Case P-1 state set.
  - Gibbs-Duhem at fixed `(T, rho)`: 1e-6 -> **9.4e-09** (finite-difference
    limited).
- **Independent route:** internal invariants plus independently written pure
  formulas; no external library is used in this case.
- **Test path:** `tests/test_pcsaft.py`.

---

## Case P-1: PC-SAFT residual properties against teqp 0.23.2

- **Source:** teqp (NIST, MIT licence), https://github.com/usnistgov/teqp -
  an independent implementation of Gross & Sadowski (2001) whose derivatives
  are all obtained by **automatic differentiation** of one hand-written
  `alphar`, where chemthermo writes them analytically. The two share the model
  and the universal constants and share no derivative code.
- **Location:** `teqp.make_model({'kind': 'PCSAFT', ...})`, methods
  `get_Ar00` (= `A^res/RT`), `get_Ar01` (= `rho dAr/drho` = `Z - 1`),
  `get_fugacity_coefficients(T, rhovec)`.
- **Assumptions:** Non-associating PC-SAFT; `k_ij` only (no `l_ij`); the state
  is fixed by `(T, rho, x)` on both sides, so no root selection is involved.
- **Components / units:** 14 states -
  pure n-hexane at (300 K, 100), (300 K, 7700), (400 K, 6800), (500 K, 3000)
  mol/m^3; methane/n-hexane 0.5/0.5 at (300 K, 200) and (300 K, 11000);
  methane/n-hexane 0.2/0.8 at (450 K, 8000) and at (300 K, 8000) (inside the
  spinodal, `Z = -0.775`); methane/n-decane with `k_ij = 0.03` at 0.3/0.7,
  (350 K, 100) and (350 K, 6500); nitrogen/methane 0.4/0.6 at (150 K, 20000);
  carbon dioxide/n-decane 0.4/0.6 at (320 K, 8000); and two ternaries
  methane/n-hexane/nitrogen, 0.3/0.4/0.3 at (250 K, 500) and 0.1/0.7/0.2 at
  (250 K, 10000). Nine of the ten mixture states have `Z > 0`; five of them are
  liquid-like (`eta` between 0.29 and 0.45).
- **Parameters and provenance:** as Case P-0. The teqp model in the test is
  built from parameter values written **in the test file**, not read from the
  package under test.
- **Expected outcome:** `A^res/RT`, `Z`, `P = Z rho R T` and `ln phi_i` equal
  teqp's to 1e-10 (relative and absolute). At the spinodal state teqp's `Z` is
  negative and chemthermo **refuses** to return `ln phi` (`ModelError`) rather
  than returning a `nan`.
- **Tolerance:** asserted 1e-10. **Achieved:** worst `|d A^res/RT| = 4.44e-15`,
  worst `|dZ| = 2.58e-14`, worst `max_i |d ln phi_i| = 2.66e-14` - i.e. four
  orders inside the asserted tolerance and at the level of double-precision
  summation order. Per-state numbers are printed by the golden-path script.
- **Negative control:** perturbing one `sigma` by 1 % moves `A^res/RT` by
  4.15e-03, so the agreement is not vacuous (asserted `> 1e-3`).
- **Not compared:** teqp's `get_Ar10` (the temperature derivative).
  chemthermo does not implement one; recorded here as a gap, not skipped
  silently.
- **Independent route:** teqp (external, autodiff).
- **Test path:** `tests/validation/test_pcsaft_vs_teqp.py::test_residual_properties_match_teqp`
- **Script:** `examples/validation/13_pcsaft_vs_teqp.py`.

---

## Case P-2: Pure n-hexane saturation from equal fugacity on two density roots

- **Source:** teqp's `pure_VLE_T` (its own Newton solve of the pure saturation
  condition) as the reference; the model is Gross & Sadowski (2001).
- **Location:** `model.pure_VLE_T(T, rho_liquid_guess, rho_vapour_guess, 200)`.
- **Assumptions:** Pure fluid, so equality of fugacity reduces to
  `ln phi(rho_L) = ln phi(rho_V)` at one pressure. chemthermo ships **no**
  density solver in this slice, so the test and the script write their own
  (scan the isotherm for the two spinodal extrema, bisect `P(rho) = P` on each
  branch, then bisect on the fugacity difference). That solver is test-only and
  deliberately unsophisticated.
- **Components / units:** n-hexane at 300 K and 400 K; pressures in Pa,
  densities in mol/m^3.
- **Parameters and provenance:** as Case P-0 (m = 3.0576, sigma = 3.7983 A,
  eps/k = 236.77 K).
- **Expected outcome (teqp):**
  - 300 K: `Psat = 21858.084278 Pa`, `rho_L = 7518.498734`,
    `rho_V = 8.868596` mol/m^3.
  - 400 K: `Psat = 463846.275296 Pa`, `rho_L = 6367.988628`,
    `rho_V = 158.932771` mol/m^3.
- **Tolerance:** asserted 1e-6 relative on all six numbers. **Achieved:**
  300 K - `Psat` 2.9e-12, `rho_L` 1.2e-16, `rho_V` 4.2e-12;
  400 K - `Psat` 3.0e-13, `rho_L` 1.4e-15, `rho_V` 1.0e-13.
  The saturation condition restated on chemthermo's own numbers gives
  `|ln phi_L - ln phi_V| <= 1e-9`.
- **Independent route:** teqp's Newton solve versus bisection written in the
  test; only the model is shared.
- **Model versus experiment (remark, not an assertion):** PC-SAFT with these
  parameters gives 21.858 kPa at 300 K. The n-hexane vapour pressure at 300 K
  is **commonly tabulated near 21.7 kPa; that figure was not verified against a
  primary source here**, so the ~0.7 % difference is reported and nothing
  asserts it. This case validates one implementation of PC-SAFT against
  another, not PC-SAFT against measurement.
- **Test path:** `tests/validation/test_pcsaft_vs_teqp.py::test_pure_hexane_saturation_matches_teqp_pure_vle`
- **Script:** `examples/validation/13_pcsaft_vs_teqp.py`.

---

## Case P-3: PC-SAFT density roots at a specified (T, P, x)

- **Source:** teqp 0.23.2 (`pure_VLE_T`) for the saturation state; the model is
  Gross & Sadowski (2001). The `dP/drho` check is an internal invariant
  (central difference of the *other* module's pressure).
- **Location:** `chemthermo.eos._pcsaft_density.solve_density_roots`, reached
  publicly through `PCSAFTEOS.density_roots` and
  `PCSAFTEOS.fugacity_coefficients`. ADR-0015.
- **Assumptions:** Fixed `(T, x)`; the model is a function of the packing
  fraction `eta` alone, so the roots of `P_model(T, rho, x) = P` are the sign
  changes of one scalar function on `0 < eta < 0.7405`. A root is *admissible*
  only if `dP/drho > 0`.
- **Components / units:** n-hexane at 300 K and 400 K; methane at 300 K;
  methane / n-hexane 0.3 / 0.7 at 300 K. Pressures in Pa, densities in
  mol/m^3.
- **Parameters and provenance:** as Case P-0 (packaged
  `src/chemthermo/parameters/data/eos/pcsaft.json`, Gross & Sadowski 2001
  Table 1); `kij = 0`.
- **Expected outcome and achieved values:**
  - n-hexane, 300 K, `P = Psat = 21858.084278856164 Pa` (teqp): **exactly two**
    admissible roots, `rho_V = 8.868596301913758`,
    `rho_L = 7518.498733715524`. teqp gives `8.868596301925571` and
    `7518.498733715526`; relative deviation **1.33e-12** (vapour) and
    **2.22e-16** (liquid), against an asserted 1e-8. Three brackets were found
    and one - the spinodal branch, `dP/drho < 0` - was discarded.
  - n-hexane, 400 K, `P = Psat = 463846.2753046097 Pa` (teqp): two roots,
    `158.9327710580466` and `6367.988627933045`; relative deviation
    **2.41e-13** and **1.11e-16**. Three brackets, one discarded.
  - n-hexane, 300 K, 10 MPa: **one** root, `7663.9424167785` (liquid only).
  - n-hexane, 300 K, 600 kPa: **one** root, `7527.542618141328`. 600 kPa is
    above the vapour spinodal maximum of this isotherm (524.4 kPa), which is
    where the vapour branch ceases to exist.
  - Methane, 300 K, 5 MPa (supercritical): **one** root,
    `2197.6188798617954`.
  - Methane / n-hexane 0.3 / 0.7, 300 K, 8 MPa: **one** root.
  - `dP/drho > 0` at every returned root, by construction and by assertion.
- **Deviation from the slice brief, recorded rather than hidden:** the brief
  expected "at `P = 0.5 Psat` only the vapour root". That is **wrong for this
  fluid**, and chemthermo returns **two** roots there
  (`4.407635978116232` and `7518.326944709169`). The liquid spinodal of
  PC-SAFT n-hexane at 300 K sits at **-39.99 MPa**, so the stretched
  (metastable) liquid root exists at every positive pressure below saturation;
  a cubic behaves identically. A root set is a *mechanical* statement. Which
  root is the phase is a Gibbs-energy question, answered by
  `_ln_phi_min_gibbs`, and below saturation it answers "vapour" - which is what
  the test asserts instead.
- **Residual note:** `|P_model - P| / P` at a returned root is limited by
  cancellation in `Z = 1 + eta a'(eta)`, not by the iteration. Achieved
  **2.2e-14** at 10 MPa, **1.6e-13** at 600 kPa, **6.5e-12** at 21.9 kPa (the
  dense liquid root, where `Z = 1.17e-3`) and **1.4e-11** at 0.5 Psat. The
  density is unaffected (see the 2.2e-16 above) and is what is pinned.
  `DensityRoots.max_relative_residual` reports the number.
- **Analytic `dP/drho` versus central finite difference** of
  `PCSAFTEOS.pressure_Pa` (step `1e-6 rho`): n-hexane 300 K at rho = 10,
  **4.9e-11**; at rho = 7518.5, **2.4e-10**; methane / n-hexane 0.3 / 0.7 at
  rho = 9000, **5.1e-11**; methane / n-decane 0.4 / 0.6 at 350 K, rho = 5000
  (inside the spinodal, `dP/drho = -3096.44`), **8.9e-11**.
- **The one-variable rewrite versus the original module:** `a`, `Z` and `P` from
  `_pcsaft_density` against `residual_helmholtz`, `compressibility_factor` and
  `pressure_Pa` at eight `eta` per state over five states - asserted 1e-12
  relative, achieved at round-off.
- **Fugacity interface:** `fugacity_coefficients(..., phase=)` equals
  `exp(ln_fugacity_coefficients)` at the corresponding root to 1e-14; a
  single-root state returns identical values for both labels; at `Psat` the two
  roots of pure n-hexane tie to `|phi_L - phi_V| = 5.6e-12`, and bisecting on
  that difference recovers `Psat` to 1e-10 relative.
- **Tolerance:** as stated per item above.
- **Independent route:** teqp for the saturation densities; central finite
  differences and the other in-tree derivation for the derivatives.
- **Test path:** `tests/test_pcsaft_density.py`,
  `tests/validation/test_pcsaft_flash_vs_teqp.py::test_pure_hexane_saturation_roots_match_teqp`
- **Script:** `examples/basic/flash_tp_pcsaft_demo.py`,
  `examples/validation/14_pcsaft_flash_vs_teqp.py`.

---

## Case P-4: Tangent-plane stability with PC-SAFT

- **Source:** the two-phase boundary of the PC-SAFT methane / n-hexane isotherm
  at 300 K, established against teqp in Case P-5. The identities are internal
  invariants (Michelsen 1982, equations (1)-(7) of
  `src/chemthermo/stability/tp.py`).
- **Location:** `stability_tp(mixture, temperature_K=300, pressure_Pa=...,
  eos=PCSAFTEOS())`. The stability solver is **unchanged** by ADR-0015;
  PC-SAFT enters through `EquationOfState.fugacity_coefficients`.
- **Assumptions:** Non-associating PC-SAFT, `kij = 0`; the Wilson trial
  estimates use Tc, Pc and omega from the packaged databank, which are
  *starting points* for the trials and not a model statement.
- **Components / units:** Methane / n-hexane at 300 K; `z1` is the methane mole
  fraction; pressures in Pa.
- **Expected outcome and achieved values** (four trials each: `wilson-vapor`,
  `wilson-liquid`, `pure-Methane`, `pure-n-Hexane`):

  | z1 | P | verdict | tpd_min |
  |---|---|---|---|
  | 0.30 | 3 MPa | unstable | -5.588192e-01 |
  | 0.50 | 5 MPa | unstable | -5.213477e-01 |
  | 0.95 | 1 MPa | unstable | -5.927738e-01 |
  | 0.30 | 8 MPa | stable | +1.982012e-01 |
  | 0.99 | 1 MPa | stable | +8.753476e-01 |
  | 0.02 | 0.5 MPa | stable | +2.996219e-01 |

- **Deviation from the slice brief, recorded rather than hidden:** the brief
  expected `z1 = 0.95` at 1 MPa to be "a vapour -> stable". It is **unstable**:
  the 1 MPa tie line runs from `x1 = 0.055615` to `y1 = 0.973473` (Case P-5),
  so `z1 = 0.95` is inside the two-phase region. `z1 = 0.99` is the stable
  vapour, and is what the test uses.
- **Identities:** `tpd(w = z) = 0` to 1e-12; every converged non-trivial trial
  satisfies `|ln W_i + ln phi_i(w) - d_i| < 1e-9` (achieved, asserted 1e-9) and
  `tpd = -ln(sum W)` to 1e-10, with the directly recomputed `tpd` agreeing to
  1e-12.
- **Direction:** at `z1 = 0.30`, 3 MPa the minimizing stationary point is a
  methane-rich vapour (`w1 > 0.9`) with `K_methane > 1 > K_hexane`, which is the
  direction of the tie line the flash then finds.
- **Permutation invariance:** reversing the component order reverses the
  stationary composition to 1e-10 and leaves `tpd_min` unchanged to 1e-10.
- **Determinism:** repeated calls return bit-identical `tpd_min` and trial
  composition.
- **Naming note:** `feed_branch` is `"vapor"` for the stable compressed liquid
  at 8 MPa. Only one density root exists there, so both phase labels evaluate
  the same state and the min-Gibbs tie-break keeps the first candidate. That is
  the documented convention (brain.md, "Vapor/liquid naming is a convention
  where the model cannot tell"), not a misclassification.
- **Tolerance:** as stated per item.
- **Independent route:** the verdicts are adjudicated by Case P-5's teqp tie
  lines; the identities are internal invariants recomputed from the
  definitions.
- **Test path:** `tests/test_pcsaft_flash.py` (stability section).
- **Script:** `examples/validation/14_pcsaft_flash_vs_teqp.py` section 3.

---

## Case P-5: PC-SAFT phi-phi flash against teqp's 300 K isotherm

- **Source:** teqp 0.23.2. Reference tie lines from
  `trace_VLE_isotherm_binary` (numerical continuation along the isotherm,
  `polish=True`, `integration_order=5`, `max_steps=10000`) polished at each
  target pressure by `mix_VLE_Tp`; bubble pressures from `mix_VLE_Tx` at a
  specified liquid composition. The model is Gross & Sadowski (2001).
- **Location:** `flash_tp(mixture, temperature_K=300, pressure_Pa=...,
  eos=PCSAFTEOS())`. The flash solver is **unchanged** by ADR-0015.
- **Assumptions:** Non-associating PC-SAFT; `kij = 0` for methane / n-hexane.
  teqp and chemthermo share the model and its 42 universal constants and
  nothing else: teqp differentiates by autodiff and continues along the
  isotherm, chemthermo differentiates analytically and reaches the tie line
  through Michelsen's tangent-plane test plus a Rachford-Rice /
  successive-substitution split.
- **Components / units:** Methane / n-hexane at 300 K, feeds taken as the
  midpoint of teqp's tie line at each pressure; densities in mol/m^3.
- **Expected outcome and achieved values** (`d` = chemthermo minus teqp;
  "teqp f" is the worst relative difference between the two phases' fugacities
  computed by **teqp's own** `get_fugacity_coefficients` at chemthermo's
  compositions and densities):

  | P / MPa | x1 (teqp) | dx1 | y1 (teqp) | dy1 | d rho_L | d rho_V | teqp f | dG/RT |
  |---|---|---|---|---|---|---|---|---|
  | 0.50 | 0.027547 | +5.8e-14 | 0.951840 | -3.6e-12 | +3.2e-14 | +3.8e-13 | 7.5e-11 | -0.8415 |
  | 1.00 | 0.055615 | +1.9e-11 | 0.973473 | -3.6e-12 | +1.1e-11 | +7.2e-13 | 3.4e-10 | -0.7602 |
  | 2.00 | 0.109597 | +6.7e-11 | 0.983863 | -2.0e-13 | +3.9e-11 | +7.9e-14 | 5.9e-10 | -0.5079 |
  | 3.00 | 0.160870 | +9.9e-11 | 0.986862 | -2.4e-14 | +5.9e-11 | +1.4e-14 | 5.8e-10 | -0.3746 |
  | 5.00 | 0.255990 | +1.6e-10 | 0.988103 | -9.9e-15 | +1.0e-10 | +1.1e-14 | 5.6e-10 | -0.2280 |
  | 7.00 | 0.342234 | +3.6e-10 | 0.986961 | -1.4e-14 | +2.3e-10 | +2.1e-14 | 8.8e-10 | -0.1464 |
  | 8.50 | 0.401747 | +1.97e-9 | 0.985041 | -1.7e-14 | +1.3e-9 | +3.4e-14 | 3.8e-9 | -0.1054 |

- **Tolerance:** asserted `|dx1|, |dy1| <= 1e-6`, densities `<= 1e-5` relative,
  teqp equal-fugacity `< 1e-8` relative. Achieved as tabulated: worst
  `|dx1| = 2.0e-9`, worst `|dy1| = 3.7e-12`, worst density deviation
  `1.3e-9`, worst teqp fugacity mismatch `3.8e-9` (all at 8.5 MPa, the state
  nearest the mixture critical region).
- **Split verification (every row):** `mass_balance_residual < 1e-12`,
  `fugacity_residual < 1e-8`, `delta_g_split_rt < 0` (tabulated above), and
  every converged phase passes the post-split stability test.
- **Bubble pressures from the stability *verdict* alone** (chemthermo solves no
  bubble-point equation: the pressure is where `stability_tp` flips from
  `unstable` to `stable` at a fixed feed, located by 40 bisections):

  | x1 | teqp `mix_VLE_Tx` / Pa | chemthermo verdict flip / Pa | relative |
  |---|---|---|---|
  | 0.10 | 1818388.518385 | 1818388.498552 | 1.09e-08 |
  | 0.20 | 3798775.343023 | 3798775.297263 | 1.20e-08 |
  | 0.30 | 5996498.051410 | 5996497.969490 | 1.37e-08 |

  Asserted 1e-6 relative; achieved 1.4e-8 worst.
- **One `kij` case:** methane / n-decane at 350 K, 5 MPa, feed `z1 = 0.4`, with
  `kij = 0.03`. **That value is illustrative** - it is the nonzero value used
  in the ADR-0014 cross-check states, **not** a literature-validated binary
  parameter for this pair; the check is that two independent codes agree once
  both are given it. Achieved: vapour fraction 0.26321712,
  `x1 = 0.18642106`, `y1 = 0.99783844`, `rho_L = 5548.0481`,
  `rho_V = 1805.7640` mol/m^3; teqp equal-fugacity mismatch **1.94e-10**
  relative; `delta_g_split_rt = -0.113917`; `fugacity_residual = 1.94e-10`.
  teqp's own `mix_VLE_Tp`, seeded from chemthermo's answer, does not move it
  (within 1e-6 mole fraction).
- **Negative control:** perturbing n-hexane's `sigma` by 1 % in the reference
  moves the saturation liquid density by more than 1e-3 relative, so the
  agreement above is not vacuous.
- **Two-phase limit:** phi-phi stops at two phases whatever
  `FlashSettings.max_phases` says (ADR-0009, ADR-0011), and the test asserts
  that a `max_phases=3` run returns two phases and no `phase_set_history`.
- **Amendment (ADR-0016, Case F-4):** the phi-phi split gained an extended
  (negative-flash) Rachford-Rice and a second-order stage after this case was
  recorded. **Every number in this entry is unchanged**: all seven tie lines,
  the three bubble pressures and the `kij` state converge in the first stage,
  where the extended solver returns the in-window solver's `float` bit for
  bit, so none of them enters the new code path. Re-measured, not assumed -
  `tests/validation/test_pcsaft_flash_vs_teqp.py` passes unchanged.
- **Independent route:** teqp's continuation + Newton equilibrium solver, and
  teqp's own fugacity coefficients evaluated at chemthermo's answer.
- **Test path:** `tests/validation/test_pcsaft_flash_vs_teqp.py`,
  `tests/test_pcsaft_flash.py` (flash section).
- **Script:** `examples/validation/14_pcsaft_flash_vs_teqp.py`.

---

## Case F-4: phi-phi split robustness - the negative flash and the second-order stage

- **Source:** The defect and the fix are internal (ADR-0016). The *method* is
  published: the negative flash is C. H. Whitson and M. L. Michelsen, "The
  negative flash", *Fluid Phase Equilibria* **53** (1989) 51-71, and the
  admissible window is the two-phase case of the `t_i > 0` region of
  C. F. Leibovici and J. Neoschil, *Fluid Phase Equilibria* **112** (1995)
  217-221. The reference equilibrium at the reference state is produced by an
  **independent route**, not read from a source: a damped Newton on the
  equal-fugacity system in vapor mole numbers, written in the test and in the
  script, with a finite-difference Jacobian of that residual.
- **Location:** `tests/test_flash_phi_phi_second_order.py`,
  `tests/validation/test_flash_split_robustness_pcsaft.py`,
  `tests/test_rachford_rice_extended.py`.
- **Assumptions:** PC-SAFT (Gross & Sadowski 2001) hard chain + dispersion, no
  association, `kij = 0`. Mole fractions; T in K, P in Pa, densities in
  mol/m^3; `tpd`, `beta` and `dG/RT` dimensionless.
- **Components / units:** carbon dioxide / n-decane and methane / n-hexane.
  Grid: CO2/n-decane `z1` in {0.6, 0.8, 0.9} x T in {230, 240, 250, 260} K x
  P in {1.0, 1.5, 2.0, 2.5} MPa (48 states); methane/n-hexane `z1` in
  {0.5, 0.8, 0.9, 0.95} x T in {170, 180, 190, 195, 200} K x P in
  {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5} MPa (140 states). 188 in total.
- **Parameters and provenance:** packaged PC-SAFT records (Gross & Sadowski
  2001 Table 1, `src/chemthermo/parameters/data/eos/pcsaft.json`); the
  validation test rebuilds teqp's model from the same four records written out
  in that file. No `kij`.

- **The defect (measured at commit `cf846fe`, before ADR-0016).** Four of the
  188 states raised `ConvergenceError("Rachford-Rice failed to bracket a vapor
  fraction.")` on the tangent-plane path - all four also on the legacy path:

  | binary | z1 | T / K | P / MPa |
  |---|---|---|---|
  | CO2 / n-decane | 0.8 | 250 | 1.0 |
  | CO2 / n-decane | 0.9 | 240 | 1.0 |
  | CO2 / n-decane | 0.9 | 250 | 1.0 |
  | CO2 / n-decane | 0.9 | 260 | 1.5 |

  At 240 K / 1.0 MPa the stability test is unambiguous (`tpd_min =
  -5.668015362e-02`, feed branch liquid, minimizer `w = (1 - 4.24e-08,
  4.24e-08)`), the seeded K-values bracket (`f(0) = +0.058`, `f(1) = -2.2e+05`)
  and successive substitution then oscillates: `K_CO2` = 1.176, 1.095, 1.196,
  1.073, 1.225, 1.042, 1.258, 1.009, 0.991, with the Rachford-Rice root
  negative at iterations 1, 3, 5 and 7 and **no root at all** at iteration 9.

- **The reference equilibrium at 240 K / 1.0 MPa** (independent damped Newton,
  two different starting vapor fractions, residual `<= 7.3e-14`):

  | quantity | value |
  |---|---|
  | `beta` | 0.1835824343 |
  | `x` (CO2, n-decane) | (0.87751368, 0.12248632) |
  | `y` (CO2, n-decane) | (1 - 7.0725e-08, 7.0725e-08) |
  | `dG_split/RT` | -5.5909263161e-03 |
  | `rho` at `x` / mol m^-3 | 699.97874 (vapor-like), 17250.70876 (liquid-like) |
  | `rho` at `y` / mol m^-3 | 557.81918 (vapor-like), 24161.93786 (liquid-like) |

  `flash_tp` after ADR-0016 agrees to **worst `|d beta| = 7.8e-13`** and
  **worst `|d composition| = 1.2e-13`** (asserted 1e-8), with
  `mass_balance_residual = 0.0`, `fugacity_residual = 6.7e-13`,
  `post_split_status = "stable"`, and the reported `delta_g_split_rt`
  reproduced from the public API to `< 1e-12`. The split converges in
  **9 successive-substitution + 9 second-order iterations**, using
  **4 negative-flash steps**.

- **teqp cross-check at the reference state.** teqp 0.23.2's own
  `get_fugacity_coefficients`, evaluated at chemthermo's converged
  compositions and at the densities `PCSAFTEOS.density_roots` returns for them,
  makes the two phases' fugacities equal to **5.5e-13** relative (asserted
  1e-8). No reference tie line is used. Negative control: perturbing CO2's
  `sigma` by 1 % breaks the agreement past 1e-3.

- **The grid, after ADR-0016** (tangent-plane path, default `FlashSettings`):

  | quantity | before (`cf846fe`) | after |
  |---|---|---|
  | states | 188 | 188 |
  | `ConvergenceError` | **4** | **0** |
  | two-phase answers | 119 | 123 |
  | single-phase answers | 65 | 65 |
  | states needing the second-order stage | - | 4 |
  | states using a negative flash | - | 4 |

  Worst invariants over the 123 two-phase answers: `mass_balance_residual`
  **1.86e-13**, `fugacity_residual` **3.58e-08** (a state converged by
  successive substitution at `tol = 1e-8` on the K-update; the four rescued
  states are at `<= 6.7e-13`), least negative `delta_g_split_rt`
  **-3.95e-04**, worst post-split `tpd_min` **-8.13e-09** - inside the default
  `tpd_tol = 1e-8` by a factor of 1.2, which is a *tighter* margin than the
  -7.0e-09 recorded for the Peng-Robinson grid in Case L-4 and is worth
  watching. Every single-phase answer is a `"stable"` stability verdict, never
  a fallback.

- **Bit-identity (the load-bearing claim).** All **184** states that converged
  before ADR-0016 return the same phase fractions, compositions and residuals
  afterwards, compared with `==`. Independently, over a **1248-state**
  Peng-Robinson grid (the 6 mixtures of Case F-3, T 150-400 K in 10 K steps,
  8 pressures from 0.1 to 8 MPa) **1247 states are unchanged** and the single
  changed state is the weakly unstable near-critical
  methane/ethane/propane (0.5, 0.3, 0.2) at 290 K and 8 MPa, which previously
  exhausted the iteration limit (`max_delta_k = 3.578e-04`) and now converges
  through the second-order stage to `beta = 0.44737586`,
  `fugacity_residual = 1.8e-15`, `delta_g_split_rt = -1.487e-05`, post-split
  stable. That was the open item recorded in `brain.md` as "an accelerated /
  second-order phi-phi split"; it is discharged.
- **The extended Rachford-Rice itself** is checked separately in
  `tests/test_rachford_rice_extended.py`: the window is exactly the set where
  every phase mole fraction is non-negative (scanned); `f` is strictly
  decreasing on it (scanned, 5000 points); a root below 0 and a root above 1
  are each found to a residual `< 1e-14` and, for a binary, agree with the
  closed form `beta = -(z1 a1 + z2 a2) / (a1 a2)` derived in the test to 1e-12;
  all-`K`-on-one-side gives no window and no root; and **every in-window root
  is the pre-ADR-0016 solver's `float` compared with `==`** - on 400 random
  `(z, K)` draws and on 60 K-vectors taken from real Peng-Robinson iterations.
- **Error semantics.** A converged vapor fraction outside `(0, 1)` while the
  feed is unstable raises `ConvergenceError` naming `beta` and `tpd_min`. **No
  state in this repository reaches that path**, so it is exercised by a
  **synthetic** fixture (the converged split is replaced; the feed, the model
  and the stability verdict are real), and the test says so.
- **Not covered:** the legacy `phase_detection="wilson-heuristic"` path, which
  still raises on all four states **by design** (ADR-0016 decision 8) and is
  asserted to; gamma-phi, gamma-gamma and modified-Raoult, whose Rachford-Rice
  is unchanged; and any claim that the second-order stage always succeeds -
  what is measured is that it succeeds on these four states and on the one
  Peng-Robinson state above.
- **Tolerance:** asserted `|d beta| <= 1e-8` and `|d composition| <= 1e-8`
  against the independent Newton (achieved 7.8e-13 / 1.2e-13); teqp
  equal-fugacity `< 1e-8` (achieved 5.5e-13); grid `mass_balance < 1e-12`,
  `fugacity_residual < 1e-6`, `delta_g_split_rt < 0` on every two-phase state.
- **Independent route:** a damped Newton on the equal-fugacity system written
  in the test and in the script (a different formulation from the solver's
  Gibbs minimization), plus teqp 0.23.2's own fugacity coefficients.
- **Test path:** `tests/test_flash_phi_phi_second_order.py`,
  `tests/validation/test_flash_split_robustness_pcsaft.py`,
  `tests/test_rachford_rice_extended.py`.
- **Script:** `examples/validation/15_flash_split_robustness.py`.

## Case F-5: Phase identity by compressibility, and a runtime trim

- **Source:** The criterion is internal (ADR-0017); the alternative it was
  weighed against and rejected is published - G. Venkatarathnam and
  L. R. Oellrich, "Identification of the phase of a fluid using partial
  derivatives of pressure, volume, and temperature without reference to
  saturation properties: Applications in phase equilibria calculations",
  *Fluid Phase Equilibria* **301** (2011) 200-203 (the `Pi` criterion,
  rejected for needing `(dP/dT)_v` and `Cp`, which this package's PC-SAFT
  implementation does not have). All numbers below are produced by an
  **independent route**: the analytic `kappa` the shipped code computes is
  cross-checked against a finite difference built without reading the
  package's private mixing-rule code (Peng-Robinson: rebuilt from
  `Component.tc_k/pc_pa/omega` directly; PC-SAFT: a finite difference of the
  already-public `pressure_Pa(T, rho, x)`).
- **Location:** `tests/test_phase_identity.py`,
  `tests/test_flash_refactor_bit_identity.py`,
  `tests/validation/test_flash_split_robustness_pcsaft_subset.py`.
- **Assumptions:** Peng-Robinson with `kij = 0`; PC-SAFT (Gross & Sadowski
  2001) hard chain + dispersion, no association, `kij = 0`. Mole fractions; T
  in K, P in Pa, densities in mol/m^3; `kappa` dimensionless.
- **Components / units, grids:**
  - Peng-Robinson: the 144-state grid of `tests/test_flash_phase_detection.py`
    / `tests/test_flash_refactor_bit_identity.py` - 6 mixtures (Methane/Ethane,
    Methane/Propane, Ethane/n-Heptane, Methane/n-Pentane,
    Methane/Ethane/Propane, Propane/n-Butane/n-Pentane) x T in {170, 200, 240,
    280, 320, 360} K x P in {0.2, 1.0, 3.0, 8.0} MPa.
  - PC-SAFT: the 188-state Case F-4 grid of
    `tests/validation/test_flash_split_robustness_pcsaft.py` - carbon dioxide /
    n-decane and methane / n-hexane, as in Case F-4.
- **Parameters and provenance:** Peng-Robinson `Tc`/`Pc`/`omega` from the
  packaged databank (unchanged by this slice); PC-SAFT records as in Case F-4.

### 1) The criterion and the threshold

`kappa = P / (rho (dP/drho)_T)` (equivalently `-P / (V (dP/dV)_T)`): exactly 1
for an ideal gas, well below 1 for a liquid. `KAPPA_LIQUID_THRESHOLD = 0.5`
(`chemthermo.models.base`) labels a root `"liquid"` when `kappa < 0.5`, else
`"vapor"`. Both derivatives are analytic in the shipped code (Peng-Robinson:
`dP/dV` of the cubic's own `P = RT/(V-b) - a/(V^2+2bV-b^2)`; PC-SAFT: the
`dP/drho` `PCSAFTIsotherm.pressure_and_slope` already computes for the
density-root solver, ADR-0015) - no new root-finding, no finite difference in
the shipped implementation.

**Analytic vs. finite difference (Peng-Robinson).** Over 9 representative
states (`tests/test_phase_identity.py::test_pr_analytic_dP_dV_matches_finite_difference`,
parametrized), the shipped analytic `kappa` and an independently rebuilt
central-difference `kappa` (step `1e-6` relative in `V`, `a_mix`/`b_mix` from
`Component.tc_k/pc_pa/omega` alone, never the package's private
`_mixture_parameters`) agree to **< 1e-8 relative** on every state, and
`PengRobinsonEOS.phase_identity`'s returned label matches what the
finite-difference `kappa` alone implies on every one.

### 2) Kappa separation over both grids

Full-grid measurements (script logic reproduced in
`tests/test_phase_identity.py`'s docstrings; the pytest module itself uses a
smaller, cheap fixed sample of each grid rather than re-scanning both in
full on every `pytest -q`, to protect the runtime trim in part 4):

| grid | class | n roots | kappa min | kappa max |
|---|---|---:|---:|---:|
| Peng-Robinson, 144-state | two-phase liquid root | 47 | 1.13e-04 | 2.30e-01 |
| Peng-Robinson, 144-state | two-phase vapor root | 47 | 1.01e+00 | 1.48e+00 |
| Peng-Robinson, 144-state | single-phase liquid root | 52 | 1.35e-04 | 1.95e-01 |
| Peng-Robinson, 144-state | single-phase vapor root | 45 | 7.44e-01 | 1.76e+00 |
| PC-SAFT, 188-state (Case F-4) | two-phase liquid root | 123 | 4.88e-04 | 8.50e-03 |
| PC-SAFT, 188-state | two-phase vapor root | 123 | 1.04e+00 | 1.87e+00 |
| PC-SAFT, 188-state | single-phase liquid root | 65 | 1.17e-03 | 3.20e-02 |
| PC-SAFT, 188-state | single-phase vapor root | 0 | - | - |

Over both grids combined: every liquid root has `kappa <= 0.230`; every vapor
root has `kappa >= 0.744`. `0.5` separates every one of the 599 roots across
both grids with a margin `>= 0.24` on each side - it is not a value chosen to
make one state come out right. (PC-SAFT's grid has no single-phase vapor
state: a property of its feed/T/P range, not of the criterion - the two-phase
vapor roots on that same grid are still `>= 1.04`.)

**Near-critical states**, where any label is a convention: none exist as a
*two-phase* result on either validated grid.
`phase_label_method == "compressibility"` on all 47 Peng-Robinson and all 123
PC-SAFT two-phase states measured - the Wilson-ranking fallback
(`_orient_two_phase_labels`, triggered when both converged phases land on the
same side of the threshold) is not exercised by either grid, and is recorded
as such rather than hidden.

### 3) The two motivating states, and one that must not move

| state | before (`32f693b`) | after |
|---|---|---|
| PC-SAFT, CO2/n-decane `z=(0.9,0.1)`, 230 K, 2.5 MPa (one density root, 18,676.8 mol/m^3) | `"vapor"`, `vapor_fraction=1.0` | `"liquid"`, `vapor_fraction=0.0`, `phase_label_method="compressibility"` |
| PC-SAFT, methane/n-hexane `z=(0.5,0.5)`, 170 K, 2.0 MPa (one density root, 12,871.6 mol/m^3) | `"vapor"`, `vapor_fraction=1.0` | `"liquid"`, `vapor_fraction=0.0`, `phase_label_method="compressibility"` |
| Peng-Robinson, methane/ethane `z=(0.5,0.5)`, 450 K, 1 bar (dilute, genuinely superheated) | `"vapor"`, `vapor_fraction=1.0` | unchanged: `"vapor"`, `vapor_fraction=1.0`, `phase_label_method="compressibility"` |

### 4) Splits: the phase called "liquid" always has the higher density

Checked over every two-phase result of both grids (47 Peng-Robinson, 123
PC-SAFT, both measured via the kappa tables above using the same converged
`x`/`y`): `Z(liquid) < Z(vapor)` on all 47 Peng-Robinson states (`Z(liquid)`
in `[0.0069, 0.302]`, `Z(vapor)` in `[0.615, 0.990]` - lower `Z` is higher
density at fixed `(T, P)`), and `rho(liquid) > rho(vapor)` on all 123 PC-SAFT
states. Zero violations on either grid.

### 5) The bit-identity audit: every changed state

`tests/test_flash_refactor_bit_identity.py` pins 155 states, bit for bit
(floats compared with `==`). Auditing the old fixture
(`refactor_bit_identity_v1.json`, HEAD `e927623`-through-`32f693b`, kept for
history) against the ADR-0017 code, state by state:

- **46 of the 144 Peng-Robinson phi-phi grid states changed** - every one
  single-phase, every one `feed_branch`/name `"vapor"` (tie-break) ->
  `"liquid"` (measured), `vapor_fraction` `1.0 -> 0.0`,
  `diagnostics["phase_state"]` mirrors the name change, and a new
  `diagnostics["phase_label_method"] = "compressibility"` key. **Zero**
  composition changes, **zero** fraction-*set* changes (a single-phase
  result's fraction set `{1.0}` is unchanged regardless of which name it is
  attached to), **zero** other diagnostics-number changes, checked
  programmatically for every one of the 155 states before the v2 fixture was
  written. The full table (mixture, T, P, the finite-difference `kappa` at
  that state, old label, new label):

  | mixture | T / K | P / MPa | kappa | old | new |
  |---|---:|---:|---:|---|---|
  | Ethane/n-Heptane | 170.0 | 1.00 | 8.04e-04 | vapor | liquid |
  | Ethane/n-Heptane | 170.0 | 3.00 | 2.33e-03 | vapor | liquid |
  | Ethane/n-Heptane | 170.0 | 8.00 | 5.74e-03 | vapor | liquid |
  | Ethane/n-Heptane | 200.0 | 1.00 | 1.27e-03 | vapor | liquid |
  | Ethane/n-Heptane | 200.0 | 3.00 | 3.65e-03 | vapor | liquid |
  | Ethane/n-Heptane | 200.0 | 8.00 | 8.80e-03 | vapor | liquid |
  | Ethane/n-Heptane | 240.0 | 1.00 | 2.38e-03 | vapor | liquid |
  | Ethane/n-Heptane | 240.0 | 3.00 | 6.72e-03 | vapor | liquid |
  | Ethane/n-Heptane | 240.0 | 8.00 | 1.56e-02 | vapor | liquid |
  | Ethane/n-Heptane | 280.0 | 3.00 | 1.31e-02 | vapor | liquid |
  | Ethane/n-Heptane | 280.0 | 8.00 | 2.83e-02 | vapor | liquid |
  | Ethane/n-Heptane | 320.0 | 8.00 | 5.53e-02 | vapor | liquid |
  | Ethane/n-Heptane | 360.0 | 8.00 | 1.24e-01 | vapor | liquid |
  | Methane/Ethane | 170.0 | 3.00 | 1.03e-02 | vapor | liquid |
  | Methane/Ethane | 170.0 | 8.00 | 2.31e-02 | vapor | liquid |
  | Methane/Ethane | 200.0 | 3.00 | 2.52e-02 | vapor | liquid |
  | Methane/Ethane | 200.0 | 8.00 | 4.88e-02 | vapor | liquid |
  | Methane/Ethane | 240.0 | 8.00 | 1.95e-01 | vapor | liquid |
  | Methane/Ethane/Propane | 170.0 | 3.00 | 8.00e-03 | vapor | liquid |
  | Methane/Ethane/Propane | 170.0 | 8.00 | 1.84e-02 | vapor | liquid |
  | Methane/Ethane/Propane | 200.0 | 3.00 | 1.72e-02 | vapor | liquid |
  | Methane/Ethane/Propane | 200.0 | 8.00 | 3.58e-02 | vapor | liquid |
  | Methane/Ethane/Propane | 240.0 | 8.00 | 1.09e-01 | vapor | liquid |
  | Methane/Propane | 170.0 | 3.00 | 1.05e-02 | vapor | liquid |
  | Methane/Propane | 170.0 | 8.00 | 2.35e-02 | vapor | liquid |
  | Methane/Propane | 200.0 | 8.00 | 4.88e-02 | vapor | liquid |
  | Methane/Propane | 240.0 | 8.00 | 1.81e-01 | vapor | liquid |
  | Methane/n-Pentane | 170.0 | 3.00 | 4.18e-03 | vapor | liquid |
  | Methane/n-Pentane | 170.0 | 8.00 | 1.00e-02 | vapor | liquid |
  | Methane/n-Pentane | 200.0 | 8.00 | 1.66e-02 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 170.0 | 1.00 | 6.67e-04 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 170.0 | 3.00 | 1.94e-03 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 170.0 | 8.00 | 4.81e-03 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 200.0 | 1.00 | 1.03e-03 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 200.0 | 3.00 | 2.96e-03 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 200.0 | 8.00 | 7.22e-03 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 240.0 | 1.00 | 1.83e-03 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 240.0 | 3.00 | 5.21e-03 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 240.0 | 8.00 | 1.23e-02 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 280.0 | 1.00 | 3.41e-03 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 280.0 | 3.00 | 9.47e-03 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 280.0 | 8.00 | 2.12e-02 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 320.0 | 3.00 | 1.88e-02 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 320.0 | 8.00 | 3.83e-02 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 360.0 | 3.00 | 4.64e-02 | vapor | liquid |
  | Propane/n-Butane/n-Pentane | 360.0 | 8.00 | 7.54e-02 | vapor | liquid |

  (Feed composition for each mixture is the one fixed value used throughout
  the grid: Methane/Ethane (0.5, 0.5), Methane/Propane (0.7, 0.3),
  Ethane/n-Heptane (0.7, 0.3), Methane/n-Pentane (0.6, 0.4),
  Methane/Ethane/Propane (0.5, 0.3, 0.2), Propane/n-Butane/n-Pentane
  (0.4, 0.3, 0.3).)

- **Zero of the 47 two-phase Peng-Robinson grid states changed**: kappa
  already agreed with the historical Wilson-ranking orientation on every one.
- **The other 11 fixture states are untouched**: 5 gamma-gamma binary feeds, 1
  gamma-gamma single-component feed, 1 Tessier (2000) near-plait feed, 2
  gamma-phi cases, 2 legacy `wilson-heuristic` cases - none build an
  `_EOSTangentPlane`, so `identity_label` is never reached for them.
- **The CLI fixture** (`tests/fixtures/cli/tp_flash_v1.json`,
  methane/ethane/propane 240 K/3 MPa, two-phase, well-separated) is
  unaffected, confirmed by `tests/test_cli_tp_flash.py` (unchanged, passing).

### 6) Runtime trim

The 188-state Case F-4 grid ran twice on every `pytest -q` -
`tests/validation/test_flash_split_robustness_pcsaft.py`'s own test and
`examples/validation/15_flash_split_robustness.py`'s smoke test in
`tests/test_examples.py` - which is most of why the suite grew from ~116 s to
~389 s between when Case F-4 was added and this slice. Decision recorded
here: **CI's default `pytest -q` does not fit the full grid**, so
- the full-grid test is `@pytest.mark.slow`, deselected by default
  (`pyproject.toml` `addopts = "-m 'not slow'"`) and run explicitly with
  `pytest -q -m slow`;
- a new fixed **16-state** subset (not the 24 first considered, trimmed for
  runtime; always including the four states that need the ADR-0016
  second-order stage) covers the grid in the default run
  (`tests/validation/test_flash_split_robustness_pcsaft_subset.py`);
- `examples/validation/15_flash_split_robustness.py` defaults to the same
  16-state idea (its own list, not imported, matching this repository's
  "duplicate small grids for self-containment" convention) with `--full` for
  the complete grid. **Amended by ADR-0020:** its four previously-failing
  states moved behind `--full` as well (they run in
  `tests/test_flash_phi_phi_second_order.py` on every default `pytest -q`),
  which took the example from ~26.9 s to ~16.5 s. The grid subset is
  unchanged at 16 states.

Measured on this machine: `pytest -q` **~389 s -> ~188 s** (553 passed, 1
deselected). The full grid, run explicitly (`pytest -q -m slow`), still
passes with the pre-slice counts unchanged: 188 states, 0
`ConvergenceError`s, 123 two-phase, 65 single-phase, 4 rescued by the
second-order stage, worst mass balance `1.86e-13`, worst fugacity residual
`3.58e-08`, worst `delta_g_split_rt` `-3.95e-04` - identical to Case F-4's own
numbers, since this slice changes no numeric result of that grid. The
remaining ~8 s gap against the ~180 s target is two pre-existing, unrelated
slow examples out of this slice's scope -
`examples/validation/12_vlle_verdict_map.py` (~16.5 s) and
`examples/validation/14_pcsaft_flash_vs_teqp.py` (~16.3 s), both present
before this slice and neither touching the PC-SAFT-grid-duplication
regression this slice fixes. (Both gained a `--full` flag of their own in
ADR-0020: 16.4 -> 7.0 s and 16.2 -> 4.6 s; see Case P-9's runtime note.)

- **Not covered:** whether `0.5` is optimal in any formal sense (it is a
  value that separates the two measured grids with a wide margin, not a
  fitted or globally optimal cut); any claim that the Wilson-ranking fallback
  correctly orients a near-critical split, since no such state exists in
  either validated grid to check it against.
- **Tolerance:** analytic vs. finite-difference `kappa` (Peng-Robinson)
  `< 1e-8` relative (achieved on 9 states); threshold separation with `>= 0.24`
  margin on both sides over both full grids (599 roots).
- **Independent route:** a finite-difference `kappa` built without the
  package's private mixing-rule/isotherm code (Peng-Robinson: from
  `Component.tc_k/pc_pa/omega`; PC-SAFT: from the public `pressure_Pa(T, rho,
  x)`).
- **Test path:** `tests/test_phase_identity.py`,
  `tests/test_flash_refactor_bit_identity.py`,
  `tests/validation/test_flash_split_robustness_pcsaft_subset.py`,
  `tests/validation/test_flash_split_robustness_pcsaft.py` (the `slow`-marked
  full grid).
- **Script:** `examples/basic/flash_tp_pcsaft_demo.py`,
  `examples/validation/15_flash_split_robustness.py` (`--full` for the
  complete grid).

---

## Case P-6: PC-SAFT association term by term against FeOs 0.10.1

- **Source:** FeOs (feos-org/feos), https://github.com/feos-org/feos, MIT OR
  Apache-2.0 - an independent Rust implementation of Gross & Sadowski (2001)
  **and** (2002) whose derivatives are all obtained by automatic
  differentiation (`num-dual`), where chemthermo writes them analytically.
  **teqp cannot serve here**: its `PCSAFT` kind implements no association term.
  FeOs reproduces teqp on non-associating n-hexane (`A^res/RT` at 300 K /
  7700 mol/m^3: teqp `-5.783742760059240`, chemthermo `-5.783742760059239`,
  FeOs `-5.783742760694397`; see the universal-constants note below), so the
  two external references corroborate each other where they overlap.
- **Location:** `PureRecord.from_json_str` with an `association_sites` entry,
  `Parameters.new_pure` / `Parameters.from_records` (an empty binary-record
  list, i.e. `k_ij = 0`), `EquationOfState.pcsaft`, `State(eos, temperature=,
  density=, composition=)`, then
  `State.residual_molar_helmholtz_energy_contributions()` (per-term J/mol,
  keys `Hard Sphere`, `Hard Chain`, `Dispersion`, `Association`),
  `State.pressure()` and `State.chemical_potential(Contributions.Residual)`.
- **Assumptions:** 2B association (`na = nb = 1`) for every associating
  component; `k_ij = 0` throughout; the state is fixed by `(T, rho, x)` on both
  sides, so no root selection is involved.
- **Components / units:** 18 states. Pure water at (300 K, 55000), (350 K,
  50000), (373.15 K, 100), (550 K, 41241.19) and (550 K, 1250.01) mol/m^3;
  pure ethanol at (300 K, 17000), (350 K, 15928.62), (450 K, 691.11) and
  (500 K, 9929.27); water/ethanol at 0.2/0.8, 0.5/0.5 and 0.8/0.2 (320 K,
  liquid-like), at 0.5/0.5 (400 K, 94.62, gas-like) and at 0.2/0.8 (351 K,
  18574.96); water/n-hexane at 0.3/0.7 and 0.9/0.1 (298.15 K, liquid-like) and
  0.5/0.5 (400 K, 93.93, gas-like); and the ternary methanol/water/n-hexane
  0.3/0.4/0.3 at (320 K, 16684.84). Every state has `Z > 0`.
- **Parameters and provenance:** Gross & Sadowski, *Ind. Eng. Chem. Res.* **41**
  (2002) 5510-5515 (DOI 10.1021/ie010954d), Table 1, for the associating five -
  Water 1.0656 / 3.0007 / 366.51 / 0.034868 / 2500.7; Methanol 1.5255 / 3.2300
  / 188.90 / 0.035176 / 2899.5; Ethanol 2.3827 / 3.1771 / 198.24 / 0.032384 /
  2653.4; 1-Propanol 2.9997 / 3.2522 / 233.40 / 0.015268 / 2276.8; n-Butanol
  (the paper's "1-butanol") 2.7515 / 3.6139 / 259.59 / 0.006692 / 2544.6
  (`m`, `sigma`/A, `eps/k` in K, `kappa^AB`, `eps^AB/k` in K) - and the 2001
  Table 1 record for n-hexane. **The 2002 paper was not read**:
  `https://pubs.acs.org/doi/10.1021/ie010954d` returned HTTP 403 from the
  environment that produced this entry, exactly as the 2001 paper did. The
  values were transcribed from two independent secondary sources that both
  cite that DOI and agree digit for digit: FeOs's
  `parameters/pcsaft/gross2002.json` and Clapeyron.jl's
  `database/SAFT/PCSAFT/PCSAFT_like.csv` + `PCSAFT_assoc.csv` (whose `source`
  column is the DOI, and whose `n_H = n_e = 1` confirms the 2B scheme). The
  reference model in the test is built from values written **in the test
  file**, not read from the package under test.
- **Expected outcome:** every term of `A^res/(R T)` separately, plus the total,
  `Z` and `ln phi_i`, equal FeOs's.
- **The one shared input that is not shared - and what it costs.** chemthermo
  packages the 2001 paper's 42 universal constants **as printed**, to ten
  significant figures, which is also what teqp uses (the two agree on
  `A^res/RT` to 4.4e-15, Case P-1). FeOs hard-codes the same constants to
  **fourteen** figures (`crates/feos/src/pcsaft/eos/dispersion.rs`); the two
  tables differ by up to **4.8e-09** termwise. That is a difference in an
  *input*, not in the model or the derivatives, and it floors any comparison
  of the dispersion term. Every dispersion-dependent quantity is therefore
  compared twice. **Achieved worst |difference| over all 18 states:**

  | quantity | as shipped | with FeOs's constants |
  |---|---:|---:|
  | hard chain (`a_hs + a_hc`) | 7.11e-15 | 7.11e-15 |
  | dispersion | 7.13e-10 | 8.88e-15 |
  | **association** | **3.11e-15** | **3.11e-15** |
  | `A^res/RT` | 7.13e-10 | 7.11e-15 |
  | `Z` | 3.86e-09 | 3.40e-14 |
  | `ln phi_i` | 1.74e-06 | 1.31e-11 |

  The association term does **not** depend on those constants and matches to
  round-off either way, which is the comparison this case exists for. The
  `ln phi` figure as shipped is the largest because a dilute component's
  `ln phi` amplifies the difference (worst state: water/n-hexane 0.9/0.1 at
  298 K). Substituting FeOs's own constants is done by monkeypatching
  `chemthermo.eos.pcsaft.A_UNIVERSAL` / `B_UNIVERSAL` in the test; the shipped
  package is **not** changed, because its stated provenance is the paper's
  table as printed and Cases P-0 to P-5 are pinned against it.
- **Tolerance:** asserted 1e-10 for the association and hard-chain terms with
  either table, and for *every* quantity with FeOs's table; 1e-8 for
  `A^res/RT` and `Z` as shipped; 1e-5 for `ln phi` as shipped. All achieved
  values are in the table above.
- **The `sigma^3` versus `d^3` question, settled numerically.** The
  association strength is
  `Delta = sigma_ij^3 g_ij^hs(d_ij) kappa^{AB}_ij [exp(eps^{AB}_ij/kT) - 1]`.
  Both `sigma_ij^3` and `d_ij^3` appear in the literature for this prefactor,
  and at 300 K water's `d` is 2.9915 A against `sigma = 3.0007` A, so the
  choice moves the cube by 0.93 % and `a_assoc` by 0.16 % of itself at a liquid
  density. Since the paper could not be read, the
  test recomputes the pure-water 2B term from scratch, in the test file, both
  ways, at (300 K, 55000 mol/m^3): **`sigma^3` reproduces FeOs's
  `-5.703948225068` to 0.0 (exactly, in double precision) and `d^3` is
  8.9e-03 away.** `sigma_ij^3` is what is implemented. Since the parameters
  and the convention were regressed together, the other spelling would be
  wrong with these parameters whatever the paper prints.
- **The cross-association rules, confirmed against a third source.** The
  Wolbach-Sandler rules (*IECR* **37** (1998) 2917) -
  `eps^{AB}_ij = (eps_ii + eps_jj)/2` and
  `kappa^{AB}_ij = sqrt(kappa_ii kappa_jj) [sqrt(sigma_i sigma_j) /
  (0.5(sigma_i + sigma_j))]^3` - reproduce Clapeyron.jl's **explicitly stored**
  water/ethanol cross pair (`epsilon_assoc = 2577.05` K,
  `bondvol = 0.03356196748232913`, source DOI 10.1021/ie010954d) digit for
  digit from the two pure records. The plain geometric mean
  (`sqrt(kappa_ii kappa_jj) = 0.03360305509920192`) does not, and using it
  moves `a_assoc` for a 0.5/0.5 water/ethanol liquid by 5.8e-04 - so the
  sigma-ratio factor is not decorative.
- **Derivative discipline (internal invariants, no external dependency),
  ten states** covering pure water, pure ethanol, water/ethanol,
  water/n-hexane and the ternary: `Z` against a central difference of
  `A^res/RT` in the density, `ln phi_i` against a central difference of
  `n A^res/RT` in the mole numbers at fixed `(T, V)`, the Euler identity
  `sum_i x_i ln phi_i = A^res/RT + Z - 1 - ln Z`, the packing-fraction path's
  `a'` and `a''` against central differences, `dP/drho` against a central
  difference of the public `pressure_Pa`, the mass-action residual
  `X_a (1 + rho (Delta W X)_a) - 1`, and Michelsen-Hendriks stationarity
  `dQ/dX = 0` with `Q` written out in the test file from its definition (both
  analytically and by central difference). **Achieved:** `|dZ| <= 4.13e-09`,
  `max_i |d ln phi_i| <= 7.26e-09`, Euler `<= 1.78e-15`, `a'` `<= 4.94e-09`,
  `a''` `<= 4.76e-09`, `dP/drho` `<= 3.72e-10`, mass action `<= 4.44e-16`
  (and `<= 4.44e-16` over the whole 1599-point density scan grid),
  stationarity `< 1e-12` analytically and `< 1e-6` by finite difference at a
  1e-7 relative step. Asserted: 1e-8 for `Z`, 1e-7 for `ln phi`, 1e-13 for
  Euler, 1e-7 for `a'`/`a''`, 1e-13 for mass action, 1e-12 for stationarity.
  The finite-difference figures are step-limited, not accuracy-limited.
- **Negative control:** a 1 % change in water's `kappa^AB` moves `a_assoc` by
  9.60e-03, and a 1 % change in `sigma` moves `A^res/RT` by 4.15e-03
  (Case P-1's control, still asserted), so neither agreement is vacuous.
- **Why the site solve is damped substitution + Newton and not substitution
  alone** (measured, to a 1e-14 mass-action residual): plain undamped
  substitution needs **926** steps for pure water at 300 K / 55 kmol/m^3 and
  **does not converge within 200,000 steps** for a 0.2/0.8 water/ethanol
  liquid at 300 K / 40 kmol/m^3, where the map oscillates (500 plain steps
  there leave `a_assoc` 4.2e-02 wrong). Damping to 0.5 needs 14 and 52 steps.
  The shipped solver takes 12 damped steps and then Newton: 0 further steps for
  pure water, 2 for the water/ethanol state, and over the whole 1599-point
  density scan grid the 12 damped steps already leave a worst residual of
  4.4e-16.
- **Robustness sweep (internal):** five associating systems (pure water,
  water/ethanol, water/n-hexane, methanol/water/n-hexane,
  1-propanol/n-butanol/water) at 250 / 298.15 / 350 / 500 / 800 K over 40
  geometrically spaced densities from 1e-6 to 60,000 mol/m^3 - 1000 states.
  Every state either returns finite `A^res/RT`, `Z` and site fractions in
  `(0, 1]`, or raises the pre-existing `eta` out-of-range `ModelError`. Zero
  `nan`s, zero unexpected exceptions.
- **Bit-identity of the non-associating path:** n-hexane at 300 K and
  7700 mol/m^3 still gives exactly `A^res/RT = -5.783742760059239`,
  `Z = 0.661534529144653`, `ln phi = -5.709015132378622`, and its density
  roots at `Psat(300 K)` are still exactly
  `(8.868596301913758, 7518.498733715524)` - asserted with `==`, not a
  tolerance, because the association code does not run when no component has
  sites.
- **Permutation invariance:** reordering methanol/water/n-hexane reproduces
  `A^res/RT` to 1.1e-16 and `ln phi_i` to **3.8e-13** relative. That is looser
  than the non-associating path's round-off because the site sums are
  accumulated in component order; recorded rather than asserted tightly.
- **Not compared:** any temperature derivative (FeOs has them, chemthermo has
  none - recorded as a gap, not skipped silently); any association scheme other
  than 2B; induced association (not implemented).
- **Independent route:** FeOs (external, autodiff, Rust), plus Clapeyron.jl's
  stored cross-association pair for the combining rules and finite differences
  written in the test for the derivatives.
- **Test path:** `tests/validation/test_pcsaft_association_vs_feos.py`,
  `tests/test_pcsaft_association.py`
- **Script:** `examples/validation/16_pcsaft_association_vs_feos.py`.

---

## Case P-7: Phase equilibrium with associating PC-SAFT

- **Source:** FeOs 0.10.1, as Case P-6: `PhaseEquilibrium.pure(eos, T)` (its
  own Newton solve of the pure saturation condition), `State.tp_flash()` (its
  own two-phase flash) and `State.chemical_potential(Contributions.Residual)`
  evaluated at **chemthermo's** converged phases and densities.
- **Location:** `stability_tp(mixture, ..., eos=PCSAFTEOS())` and
  `flash_tp(mixture, ..., eos=PCSAFTEOS())`. **No solver changed** in the
  `pcsaft-association` slice: `stability/`, `flash/` and `models/` are
  untouched by it.
- **Assumptions:** 2B association; `k_ij = 0`; the packaged 2002 parameters.
  The equilibrium comparisons against FeOs run with **FeOs's** universal
  constants substituted into chemthermo (see Case P-6), because FeOs evaluates
  fugacities at chemthermo's own densities and the table difference otherwise
  floors the residual; the tie lines are the same to eight decimal places
  either way, and the pure-water saturation comparison below is run on the
  **shipped** constants.
- **Components / units:** water at 373.15 K; water/ethanol at 351 K and
  80 kPa, `z = (0.7, 0.3)`; water/n-hexane at 298.15 K, `z = (0.5, 0.5)`, at
  101,325 Pa and at 1 MPa. Pressures in Pa, densities in mol/m^3.

### (i) Pure-water saturation at 373.15 K

- **Route:** bisection (90 steps) on
  `ln phi(liquid root) - ln phi(vapour root) = 0`, with chemthermo's own
  `density_roots` supplying the two branches - the Case P-2 route updated to
  use the ADR-0015 root solver. FeOs reaches the same state by a Newton
  iteration of its own, so only the model is shared.
- **Expected outcome (FeOs):** `Psat = 100,890.27301264 Pa`,
  `rho_L = 48,755.50956363`, `rho_V = 33.12715221` mol/m^3.
- **Achieved (chemthermo, shipped constants):** `100,890.27305023 Pa`,
  `48,755.50956045`, `33.12715222` - relative **3.73e-10**, **6.52e-11** and
  **3.36e-10** against an asserted 1e-6. The residual difference is the
  universal-constants table, not the solvers. The saturation condition
  restated on chemthermo's own numbers gives
  `|ln phi_L - ln phi_V| < 1e-9`.
- **Model versus experiment (remark, not an assertion):** 373.15 K *is*
  water's normal boiling point, so the experimental saturation pressure there
  is 101,325 Pa by definition. PC-SAFT with the 2002 parameters is **0.43 %
  low**. Nothing asserts that; this case validates one implementation against
  another.

### (ii) Water / ethanol vapour-liquid flash at 351 K and 80 kPa

- **Expected outcome:** `stability_tp` unstable, `flash_tp` returns a verified
  two-phase VLE result, and FeOs agrees the two phases are in equilibrium.
- **Achieved:** `stability_tp` -> `unstable`, `tpd_min = -2.582408e-01`.
  `flash_tp` -> `vapor_fraction = 0.60468054`,
  liquid `x = (0.95844408, 0.04155592)` at 45,864.1053 mol/m^3,
  vapour `y = (0.53103810, 0.46896190)` at 28.2488 mol/m^3,
  `phase_regime = "VLE"`, `phase_label_method = "compressibility"`,
  `delta_g_split_rt = -5.928e-02`, `mass_balance_residual = 8.8e-14`,
  `fugacity_residual = 6.58e-10`, `post_split_status = "stable"`.
  **FeOs's own fugacities at chemthermo's phases and densities give
  `max_i |ln(x_i phi_i)^I - ln(x_i phi_i)^II| = 6.12e-10`** against an
  asserted 1e-8. (With the shipped universal constants the same residual is
  2.53e-06 and the compositions are unchanged to eight decimal places.)
- The two-phase window at 351 K is narrow with these parameters: the same feed
  is a single stable liquid at 1 atm, and at `z = (0.5, 0.5)` the unstable
  band found on a 5 kPa scan is 85-95 kPa.

### (iii) Water / n-hexane at 298.15 K - what the phi-phi path can and cannot do

> **Amended 2026-09-13 by ADR-0019 (`flash-eos-per-phase-roots`), superseded by
> Case P-8.** The limitation recorded below is **resolved**: the split no longer
> pins one phase to the liquid root and the other to the vapour root, so at 1 atm
> `flash_tp` returns the two liquids FeOs returns, and at 1 MPa the same tie line
> comes back named `liquid1` / `liquid2` with `vapor_fraction = None` instead of
> falling through to the Wilson ranking. Everything below is kept **as written**,
> because it is the measured record of the defect and of what the pre-ADR-0019
> code actually produced; Case P-8 carries the numbers that hold now. The two
> tests named at the end of Case P-7 were renamed and rewritten accordingly.

- **At 1 atm, the documented failure (pre-ADR-0019).** `stability_tp` is **right**: the feed
  is `unstable` with `tpd_min = -9.281926e-01`. `flash_tp` then **raises
  `ConvergenceError`** ("a third phase is required"). With
  `FlashSettings(post_split_stability=False)` the converged pair is a
  water-rich liquid `(0.99992175, 0.00007826)` on the liquid root against a
  hexane-rich phase `(0.03206195, 0.96793805)` on the **vapour** root
  (43.3 mol/m^3), with `delta_g_split_rt = +0.2594` - i.e. a "split" whose
  Gibbs energy is *above* the feed's - and both phases individually unstable
  (`tpd = -1.53` and `-0.78`). The post-split test catches it. **The cause is
  structural, not numerical**: `_split._ln_phi_function` evaluates one phase
  with `phase="liquid"` and the other with `phase="vapor"`, so the phi-phi
  path cannot put both phases on the liquid density branch, and at 1 atm a
  vapour root still exists. The true answer at that state is two liquids -
  water's and n-hexane's vapour pressures sum to about 23 kPa, far below
  1 atm - and **FeOs's own `State.tp_flash()` at the same state returns
  exactly that**: `x = (0.00631, 0.99369)` at 7578.0 mol/m^3 and
  `(0.99998, 0.00002)` at 51,174.6 mol/m^3, both liquid densities. Recorded as
  a limitation of the flash path, not worked around; the brief for this slice
  forbade touching the solvers and the fix belongs to
  `flash-phase-addition-eos`.
- **At 1 MPa, the liquid-liquid split is reachable and correct.** Above about
  0.6 MPa the isotherm has a **single** density root at every composition
  checked (0.5/0.5, 0.999/0.001, 0.001/0.999), so `"vapor"` and `"liquid"`
  name the same root and the unchanged machinery expresses a genuine
  liquid-liquid split. Achieved: `stability_tp` -> `unstable`; `flash_tp` ->
  water-rich `(0.99998318, 0.00001682)` at 51,185.85 mol/m^3 and hexane-rich
  `(0.00630483, 0.99369517)` at 7,591.86 mol/m^3,
  `delta_g_split_rt = -5.079e-01`, `fugacity_residual = 2.56e-13`,
  `post_split_status = "stable"`. `PCSAFTEOS.phase_identity` measures **both**
  phases as `"liquid"`. FeOs's fugacities at those phases give an
  equal-fugacity residual of **1.22e-11** against an asserted 1e-8 (2.62e-07
  with the shipped constants). The compositions agree with FeOs's own 1 atm
  `tp_flash` to the printed five decimal places, as they should - a liquid
  tie line barely moves between 1 atm and 1 MPa.
- **The labels actually produced, recorded as fact.** Both phases are liquids
  by the ADR-0017 compressibility criterion, so the two-phase orientation rule
  falls through to its documented fallback: the phases come back named
  `"liquid"` and `"vapor"`, `diagnostics["phase_label_method"] ==
  "wilson-ranking"`, and `vapor_fraction = 0.503164` is really the hexane-rich
  **liquid**'s fraction. ADR-0017 anticipated the near-critical version of
  this; a liquid-liquid EOS split is a second, more common way in. Naming them
  `liquid1` / `liquid2` requires `_detect` to grow that vocabulary for
  phi-phi, which is a flash-module change and out of this slice's scope.
- **Model versus experiment (remark, not an assertion):** this model gives
  **1.68e-05** mole fraction hexane in the water-rich phase and **6.30e-03**
  water in the hexane-rich phase. Commonly tabulated experimental values are
  about **2e-6** and **5e-4** at 298 K - figures that were **not verified
  against a primary source here**. The model is an order of magnitude out on
  both, which is the expected behaviour of PC-SAFT with `k_ij = 0` for a
  water/hydrocarbon pair. **This case is a check of the code against another
  implementation, not of the model against measurement.**
- **Tolerance:** asserted 1e-6 relative on the three saturation numbers
  (achieved <= 3.73e-10) and 1e-8 on both equal-fugacity residuals (achieved
  6.12e-10 and 1.22e-11).
- **Not covered (at the time of writing):** a three-phase (vapour + two
  liquids) EOS state, which this system has at 1 atm near 330 K and which no
  path in this package can return - **still not covered**, and now the only
  remaining phi-phi phase-count gap (Case P-8 brackets the window); any
  liquid-liquid EOS split at a pressure where a vapour root still exists -
  **covered since ADR-0019**, see Case P-8.
- **Independent route:** FeOs's own saturation Newton solve, its own two-phase
  flash, and its own chemical potentials evaluated at chemthermo's converged
  states.
- **Test path:** `tests/validation/test_pcsaft_association_vs_feos.py`
  (`test_pure_water_saturation_matches_feos`,
  `test_water_ethanol_vapor_liquid_flash_matches_feos`,
  `test_water_hexane_at_one_atm_is_two_liquids_and_feos_agrees` - renamed from
  `test_water_hexane_is_unstable_and_the_phi_phi_split_cannot_express_two_liquids`
  by ADR-0019, `test_water_hexane_liquid_liquid_split_above_the_vapour_root`)
- **Script:** `examples/validation/16_pcsaft_association_vs_feos.py`,
  `examples/basic/pcsaft_association_demo.py`.

---

## Case P-8: Liquid-liquid equilibrium from an equation of state

- **Source:** FeOs 0.10.1 (feos-org/feos, MIT OR Apache-2.0), as Cases P-6 and
  P-7: `State.tp_flash()` (its own two-phase flash, with its own stability
  analysis) and `State.chemical_potential(Contributions.Residual)` evaluated at
  **chemthermo's** converged phases and densities.
- **Location:** `flash_tp(mixture, ..., eos=PCSAFTEOS())` and
  `flash_tp(mixture, ..., eos=PengRobinsonEOS())` on the tangent-plane path.
  The slice is `flash-eos-per-phase-roots` (ADR-0019); it changes
  `chemthermo/flash/_split.py` and `_detect.py` and nothing else in `src/`.
- **Assumptions:** 2B association for water; `k_ij = 0`; the packaged Gross &
  Sadowski (2002) / (2001) parameters. As in Case P-7, the comparisons that
  evaluate **FeOs** at **chemthermo's own densities** are run twice - with the
  shipped ten-figure universal constants and with FeOs's fourteen-figure ones -
  because that table difference, not either solver, floors the residual. Both
  numbers are recorded below.
- **Components / units:** water / n-hexane at 298.15 K, at 101,325 Pa and at
  1 MPa, `z = (0.5, 0.5)`, `(0.2, 0.8)` and `(0.8, 0.2)`. Pressures in Pa,
  densities in mol/m^3, compositions in mole fractions.

### (i) The state Case P-7(iii) could not express: 298.15 K, 1 atm, z = 0.5/0.5

- **Premise, measured:** `stability_tp` -> `unstable`,
  `tpd_min = -9.281926e-01`, `feed_branch = "liquid"`,
  `phase_branch = "liquid"`; the isotherm still has **two** density roots at
  the feed composition (43.01 and 13,419.37 mol/m^3), so the pre-ADR-0019 fixed
  liquid/vapour pairing had a vapour root to land on and did (Case P-7(iii)).
- **Expected outcome (FeOs `tp_flash`):** water-rich phase
  `x = (0.9999832572043167, 1.6742795683233053e-05)` at
  `rho = 51,174.64281417162`; hexane-rich phase
  `x = (0.0063122298802306175, 0.9936877701197695)` at
  `rho = 7,578.020786675787`; fraction of the hexane-rich phase
  `0.5031677924139043`. FeOs's container calls the hexane-rich phase `vapor`;
  **both densities are liquid densities**, which the test asserts rather than
  trusting the names.
- **Achieved (chemthermo, shipped constants):** `phases = {"liquid1",
  "liquid2"}`, `vapor_fraction = None`, `phase_regime = "LLE"`,
  `phase_label_method = "compressibility"`,
  `phase_i_branch = phase_ii_branch = "liquid"`.
  `liquid1 = (0.999983257204303, 1.6742795697016702e-05)` at
  `rho = 51,174.64281031051`, fraction `0.49683221079476425`;
  `liquid2 = (0.006312223543853401, 0.9936877764561466)` at
  `rho = 7,578.020745569911`, fraction `0.5031677892052357`.
  `delta_g_split_rt = -5.084812e-01`, `fugacity_residual = 7.09e-12`,
  `mass_balance_residual = 9.33e-14`, `post_split_status = "stable"`.
- **Agreement with FeOs's flash (shipped constants):** compositions
  **1.38e-14** (water-rich) and **6.34e-09** (hexane-rich) absolute; densities
  **7.54e-11** and **5.42e-09** relative; phase fraction **3.21e-09**. With
  FeOs's constants substituted the water-rich numbers fall to **6.66e-16** and
  **3.55e-15**; the hexane-rich ones do **not** (6.34e-09, 5.34e-09), because
  they are FeOs's own flash tolerance rather than a model difference - see the
  chemical-potential check next.
- **FeOs's chemical potentials at chemthermo's phases:**
  `max_i |mu_i^I - mu_i^II| / RT` = **5.265e-12** with matched constants and
  **2.606e-06** as shipped, against an asserted 1e-8. This is the check that
  does not depend on FeOs's flash converging: only the two compositions and
  the two densities are chemthermo's.
- **Both phases are liquids, independently:** `kappa = P / (rho dP/drho)`
  recomputed from the public `PCSAFTEOS.pressure_Pa` by central difference (not
  from `phase_identity`'s own analytic derivative) is **2.475e-05** for
  `liquid1` and **2.081e-04** for `liquid2`, against
  `KAPPA_LIQUID_THRESHOLD = 0.5` (ADR-0017).

### (ii) The same tie line from three feeds (the lever rule)

- **Expected outcome:** a tie line is a property of the state; changing the
  feed must change only the amounts.
- **Achieved:** `z = (0.2, 0.8)` and `(0.8, 0.2)` return the *same* two
  compositions as `z = (0.5, 0.5)` to better than 1e-10, with fractions
  `0.19492143 / 0.80507857` and `0.79874299 / 0.20125701`; the lever-rule
  residual `max_i |z_i - ((1-beta) x_i^I + beta x_i^II)|` is **9.3e-14**,
  **1.1e-16** and **0.0** respectively. `liquid1` is the water-rich phase at
  every feed, because ADR-0019 orders the two names by the first component's
  mole fraction rather than by seed role.
- **Observed limit of the reference, recorded not worked around:** FeOs's own
  `tp_flash` raises `"stability analysis did not converge"` at
  `z = (0.2, 0.8)` on this binary (feos 0.10.1). `z = (0.8, 0.2)` converges. So
  the three feeds are compared against the `z = 0.5/0.5` FeOs tie line, which
  all three reproduce to 1e-8.

### (iii) The 1 MPa tie line, now named `liquid1` / `liquid2`

- **Expected outcome:** the same numbers Case P-7(iii) recorded at 1 MPa, with
  the `wilson-ranking` fallback no longer firing.
- **Achieved:** `liquid1 = (0.9999831809598392, 1.68190401608668e-05)` at
  `rho = 51,185.848656244416`, fraction `0.4968359921458614`;
  `liquid2 = (0.006304831246250744, 0.9936951687537494)` at
  `rho = 7,591.861938295979`, fraction `0.5031640078541386`;
  `vapor_fraction = None`, `phase_regime = "LLE"`,
  `phase_label_method = "compressibility"` (it was `"wilson-ranking"` before),
  `delta_g_split_rt = -5.079126e-01`, `fugacity_residual = 2.56e-13`,
  `post_split_status = "stable"`. Against FeOs's flash: compositions 1.49e-14
  and 6.47e-09, densities 7.54e-11 and 5.53e-09 relative, fraction 3.28e-09.
  FeOs's chemical potentials at these phases: **1.222e-11** matched,
  **2.620e-07** as shipped. The tie line moves by less than 1e-5 between 1 atm
  and 1 MPa, as a liquid one should.

### (iv) Peng-Robinson on the same binary (reported, not forced)

- **Expected outcome:** none asserted. The question was only whether a cubic
  with `k_ij = 0` shows a miscibility gap here, and whether whatever it returns
  is a verified equilibrium **of that model**.
- **Achieved:** at 298.15 K / 1 atm, `stability_tp` -> `unstable`,
  `tpd_min = -2.481382`, both branches `"liquid"`. `flash_tp` returns
  `liquid1 = (0.99999999999851219, 1.4878377535914647e-12)` (fraction
  `0.49135398481939774`) and
  `liquid2 = (0.016998098723775595, 0.98300190127622444)` (fraction
  `0.5086460151806023`), `delta_g_split_rt = -1.121704189194578`,
  `fugacity_residual = 6.59e-13`, `post_split_status = "stable"`, both phases
  measured `"liquid"` by `phase_identity`. **Before this slice the same call
  raised** `ConvergenceError` ("a third phase is required", post-split
  `tpd = -1.523828`), which was confirmed by re-running with the pre-slice
  branch pinning forced.
- **Model versus experiment (remark, not an assertion):** Peng-Robinson with
  `k_ij = 0` puts essentially **zero** hexane in the water-rich phase
  (1.5e-12 mole fraction) and 1.70e-02 water in the hexane-rich phase. Neither
  is a useful prediction and nothing asserts them; this sub-case records that
  the machinery is model-agnostic, not that the cubic is right here.

### Bit-identity and regression

- **Expected outcome:** every previously pinned number unchanged.
- **Achieved:** `tests/test_flash_refactor_bit_identity.py` passes **unchanged**
  against `refactor_bit_identity_v2.json` (155 states, floats compared with
  `==`); **no fixture regeneration was needed and none was done**, so the audit
  policy's v3 path was not taken. Separately,
  `test_the_pinned_root_is_the_historical_branch_on_the_whole_phi_phi_grid`
  replays all 144 Peng-Robinson phi-phi states of that fixture with `_PhaseRoot`
  instrumented and asserts that on all **47** two-phase states, at **every**
  iterate, each phase's fugacity coefficients are exactly (`==`) the ones the
  pre-slice `phase="liquid"` / `phase="vapor"` assignment would have produced
  (the other 97 states are single-phase and build no `_PhaseRoot`). Cases
  P-1..P-7, F-4 and its 16-state subset, the 188-state slow grid, the teqp
  cross-checks, the modified-Raoult / VLLE / gamma-gamma numbers and the CLI
  contract fixtures are unchanged.
- **A rejected design, measured:** re-selecting the lowest-Gibbs root at every
  *iterate* (rather than pinning the branch the stability test reported) breaks
  the ADR-0016 reference state. PC-SAFT carbon dioxide / n-decane,
  `z = (0.9, 0.1)`, 240 K, 1.0 MPa: the liquid phase reaches
  `x = (0.99066, 0.00934)` at the fifteenth iterate, where the vapour root has
  the lower Gibbs energy, both phases then collapse onto one root and the split
  returns `beta = -3.247203e+09`. Recorded in ADR-0019 "Alternatives
  considered"; the shipped rule pins the branch and lets the post-split
  stability test enforce the lowest-Gibbs condition at the solution.

### Still not covered (at the time of ADR-0019; resolved by ADR-0020)

~~A three-phase (vapour + two liquids) EOS state. **Bracketed, measured with
this slice's code**, water / n-hexane at 1 atm, `z = (0.5, 0.5)`: 322, 324,
326 and 328 K return a stable two-liquid result; 335, 336, 340 and 350 K
return a stable vapour-liquid result; **330 K and 334 K raise**
`ConvergenceError` because the post-split stability test finds a converged
phase unstable. That window is the target of the next slice,
`flash-phase-addition-eos`.~~

**Resolved by ADR-0020**, and the diagnosis of the window turned out to be the
Case R-3 one rather than a three-phase region: a binary at fixed pressure has
no three-phase region at all (Gibbs' phase rule), and below `T3` the search
runs `V -> LV -> LLV -> LL`, reaching the two conjugate liquids by *removing*
a phase that addition had to add first. 330 K and 334 K now return verified
two-liquid answers. Genuine three-phase EOS states are ternary and are
recorded in Case P-10; the binary window is Case P-9.

- **Tolerance:** asserted 1e-8 absolute on compositions and phase fractions
  against FeOs's flash (achieved <= 6.47e-09), 1e-6 relative on densities
  (achieved <= 5.53e-09), 1e-8 on FeOs's chemical potentials at chemthermo's
  phases with matched constants (achieved <= 1.22e-11), 1e-9 on the
  fugacity residual (achieved <= 7.09e-12) and 1e-12 on the mass balance
  (achieved <= 9.33e-14).
- **Independent route:** FeOs's own two-phase flash and its own chemical
  potentials; the lever rule; a `kappa` rebuilt from the public
  `pressure_Pa` by finite difference; and, for Peng-Robinson, `stability_tp`
  re-run on each converged phase from the public API.
- **Negative control:** perturbing water's `epsilon^AB / k` by 1 % moves
  `x_water` in the hexane-rich phase by **4.79e-04**, so the agreement is not a
  statement about numbers that stopped depending on the model.
- **Test path:** `tests/test_flash_eos_lle.py` (no optional dependency) and
  `tests/validation/test_pcsaft_lle_vs_feos.py` (skipped without `feos`).
- **Script:** `examples/basic/flash_tp_pcsaft_lle_demo.py`,
  `examples/validation/17_pcsaft_lle_vs_feos.py`.

## Case P-9: The three-phase neighbourhood of water / n-hexane (binary)

- **Source:** an independent 4-equation Newton written in the test and in
  `examples/validation/18_pcsaft_vlle_water_hexane.py` (two separate copies,
  agreeing to 1e-09 K), independent two-equation Newtons for each two-phase
  pair, reduced Gibbs energies computed from the public
  `EquationOfState.fugacity_coefficients`, and **FeOs 0.10.1** (feos-org/feos,
  MIT OR Apache-2.0): `State.chemical_potential(Contributions.Residual)`
  evaluated at chemthermo's converged phases and densities, and `State.tp_flash`.
- **Location:** `flash_tp(mixture, ..., eos=PCSAFTEOS())` on the tangent-plane
  path. The slice is `flash-phase-addition-eos` (ADR-0020); it changes
  `chemthermo/flash/_multiphase.py` and `_detect.py` and nothing else in `src/`.
- **Assumptions:** 2B association for water; `k_ij = 0`; the packaged Gross &
  Sadowski (2002) / (2001) parameters. As in Cases P-6 to P-8, every comparison
  that evaluates **FeOs** at **chemthermo's own densities** is run twice - with
  the shipped ten-figure universal constants and with FeOs's fourteen-figure
  ones - because that table difference, not either solver, floors the residual.
- **Components / units:** water / n-hexane at 101,325 Pa, `z_water` of 0.3, 0.5
  and 0.7. Temperatures in K, densities in mol/m^3, compositions in mole
  fractions.

### The three-phase point

- **Expected outcome:** none available from a source. `T3` and the three
  coexisting compositions are *computed* by a 4-equation Newton in
  `(x^I, x^II, y, T)` - two liquids on the model's liquid branch, a vapour on
  its vapour branch - started from a coarse bracket `(0.9999, 0.02, 0.20,
  334.5)`.
- **Achieved:** residual **1.74e-12** in 7 iterations.
  `T3 = 334.807826336 K` (61.6578 C); `x_water(I) = 0.999935973994`
  (water-rich liquid), `x_water(II) = 0.022598600659` (hexane-rich liquid),
  `y_water = 0.213124406737`.
- **Why per-phase roots are needed:** at `T3` every one of the three
  compositions has **two** mechanically stable density roots - I (37.91,
  49,987.86), II (37.93, 7,289.79), V (37.59, 8,795.95) mol/m^3. A split that
  pins one phase to the liquid branch and another to the vapour branch cannot
  describe two liquids and a vapour at once (ADR-0019).
- **Weak external sanity check, not a reference value:** the water / n-hexane
  heteroazeotrope at 1 atm is commonly tabulated near 61.6 C with
  `y_water ~ 0.21`. No primary source was verified, and the test asserts only
  `330 K < T3 < 340 K` and `0.15 < y_water < 0.28`.

### (i) Below T3: the refusal window, resolved by add-then-remove

- **Premise, measured with ADR-0019's code:** 330 K and 334 K raised
  `ConvergenceError` (post-split test finds a converged phase unstable).
- **Expected outcome:** the two conjugate liquids, from an independent
  two-equation Newton on the liquid branch, with the amounts the lever rule
  gives.
- **Achieved** at `T = T3 - 0.05 K = 334.757826 K`, `z = (0.5, 0.5)`:
  `phases = {"liquid1", "liquid2"}`, `vapor_fraction = None`,
  `phase_regime = "LLE"`, `phase_label_method = "compressibility"`,
  `phase_set_history = "V -> LV -> LLV -> LL"`, `phases_added = 1`,
  `phases_removed = 1`, `post_split_status = "stable"`.
  `liquid1 = (0.999936079895, 6.392010504959e-05)` at `beta = 0.4884895575`;
  `liquid2 = (0.022563892766, 0.977436107234)` at `beta = 0.5115104425`.
  `equilibrium_residual = 8.0e-13`, `mass_balance_residual = 0.0`,
  `delta_g_split_rt = -0.196766`.
- **Against the independent Newton:** `|dx| = 2.285e-13`; against the lever
  rule `|d beta| = 1.198e-13`.
- **Gibbs ordering:** `G(LL)/RT = -0.916796623689 < G(VL)/RT = -0.915715633921
  < G(feed)/RT = -0.720030204244`, and
  `delta_g_vs_two_phase_rt = -1.080990e-03` equals `G(LL) - G(VL)` to better
  than 1e-09 - i.e. the diagnostics key *is* the comparison against the
  two-phase pair the search started from.
- **FeOs's chemical potentials at chemthermo's two phases:**
  `max_i |mu_i^I - mu_i^II| / RT` = **6.828e-11** with matched constants and
  **2.193e-06** as shipped, against an asserted 1e-8.
- **What the reference's own flash does here, recorded not worked around:**
  FeOs's `State.tp_flash` converges on the **vapour-liquid** pair
  (`x_water = 0.999935940756` at 49,989.45 mol/m^3 against
  `y_water = 0.212638848885` at 37.60 mol/m^3, vapour fraction 0.6350029054),
  which is exactly the pair chemthermo converges first and then refuses:
  `FlashSettings(post_split_stability=False)` reproduces it to **6.35e-10**.
  Both are stationary states of the same model; the Gibbs comparison above is
  what decides between them, and it prefers the two liquids by 1.081e-03 RT.
- **`max_phases = 2`** reproduces the pre-ADR-0020 refusal at this state.

### (ii) Above T3: the vapour-liquid answer, and the search is not entered

- **Achieved** at `T = T3 + 0.05 K = 334.857826 K`, `z = (0.5, 0.5)`:
  `phases = {"liquid", "vapor"}`, `phase_regime = "VLE"`, **no**
  `phase_set_history` key, `post_split_status = "stable"`.
  `liquid = (0.999936007302, 6.399269782050e-05)`,
  `vapor = (0.213610914465, 0.786389085535)`, `vapor_fraction = 0.6357879354`.
- **Against the independent VL Newton:** `|dx| = 6.612e-12`.
- **Gibbs ordering:** `G(VL)/RT = -0.914057391035 < G(LL)/RT = -0.912975942037`,
  so the vapour-liquid pair is the equilibrium here and the two-liquid pair -
  which still exists as a stationary state - is not.
- **FeOs's chemical potentials:** **4.078e-11** matched, **2.498e-06** as
  shipped.

### (iii) At T3: what is and is not claimed

- Gibbs' phase rule gives `F = 2 - 3 + 2 = 1`, so on a binary at fixed pressure
  three phases coexist at **one** temperature, and there the three phase
  *amounts* solve an underdetermined system (three unknowns, two independent
  balances). **No three-phase `FlashResult` is produced or claimed for the
  binary**, and none should be.
- **Achieved** at `T = T3`: `flash_tp` returns the vapour-liquid edge of the
  tie triangle - `liquid = 0.999935974`, `vapor = 0.213124406` - matching
  vertices I and V of the 4-equation Newton to better than **1e-06**. The third
  vertex is a **zero** of the tangent-plane distance from the returned liquid
  (`|tpd| < 1e-09`, computed in the test), which is "three phases coexist here"
  written in the stability test's own terms, and is why the post-split test
  reports the pair stable rather than unstable.
- **The verdict boundary locates T3.** Bisecting the `LLE` / `VLE` verdict of
  `flash_tp` over `T3 +/- 0.05 K` puts the switch at **334.8078261 K**, which is
  **2.4e-07 K** below the independently computed `T3` (45-step bisection,
  scratch measurement). The shipped test bisects a `+/- 1e-03 K` bracket ten
  times and asserts agreement to 1e-05 K.

### (iv) The scan across the window, and one honest miss - **retired by ADR-0021**

- **Expected outcome:** no `ConvergenceError` anywhere in `[T3 - 1 K, T3 + 1 K]`,
  and a verdict that switches once.
- **Achieved:** 41 temperatures at `z_water = 0.3` and 41 at `z_water = 0.7`,
  **82 states, zero `ConvergenceError`**. At `z_water = 0.3`: 20 `LLE` then 21
  `VLE`, exactly one switch, at `T3`. A third scan at `z_water = 0.5` (scratch,
  not shipped) gives the same shape.
- **The miss, as measured on the `flash-phase-addition-eos` slice.** At
  `z_water = 0.7` the verdict was `LLE` at all 41 temperatures, including above
  `T3` where it should be `VLE`. At 335 K the returned two-liquid pair has
  `G/RT = -1.1618107137` against the vapour-liquid pair's `-1.1643059308`: the
  answer was **metastable by 2.495e-03 RT**. The cause was the *stability*
  test, not the search: all four deterministic trials from the hexane-rich
  liquid `(0.02273, 0.97727)` converged to the trivial solution or to its
  partner (`tpd_min = -3.007e-09`, verdict `"stable"`), while the vapour
  stationary point at `y = (0.21500, 0.78500)` has `tpd = -6.530e-03`. That was
  `_EOSTangentPlane`'s trial set, which ADR-0012 deliberately left with
  per-iterate minimum-Gibbs root selection and no fixed surfaces, and it was the
  invariant "a phase count is never better than the stability test that
  produced it" made concrete again.
- **AMENDMENT (slice `stability-eos-root-surfaces`, ADR-0021): the miss is
  gone.** With each equation-of-state trial pinned to one density root, the
  `wilson-vapor` trial stays on the vapour root and reaches that stationary
  point in 7 iterations (`tpd = -6.5237494612e-03` at
  `w = (0.213559327395, 0.786440672605)`, `phase_branch = "vapor"`,
  `minimizing_trial_surface = "vapor"`), the hexane-rich liquid is reported
  `"unstable"`, and `flash_tp` returns the vapour-liquid pair. Both 41-point
  scans now switch `LLE -> VLE` **exactly once**, with **zero
  `ConvergenceError`** (`z_water = 0.3`: 20 `LLE` then 21 `VLE`;
  `z_water = 0.7`: 21 then 20). Bisecting each verdict boundary to 1e-06 K puts
  it **2.09e-07 K below** `T3` at `z_water = 0.3` and **2.68e-07 K above** it at
  `z_water = 0.7`. See **Case P-11** for the full evidence, and note that
  `tests/test_flash_vlle_eos.py::test_the_window_scan_never_raises` no longer
  pins the miss - it pins the switch, at both feeds.

### Negative controls

- `T = T3 + 20 K = 354.81 K`, `z = (0.5, 0.5)` -> a single `"vapor"`,
  `phase_regime = "single-phase"`, no search key.
- `z_water = 0.99999`, 300 K -> a single `"liquid"`. **`z_water = 0.999` is
  not a valid negative control for this model**: it puts 1e-03 mole fraction
  n-hexane into a water-rich phase whose binodal composition is 1.67e-05
  (Case P-8), so two liquids is the correct answer there and the test records
  that instead of calling it a miss.
- The 155-state bit-identity fixture carries **no** `phase_set_history`,
  `phases_added`, `phases_removed`, `rachford_rice_iterations` or
  `delta_g_vs_two_phase_rt` key, asserted directly on the JSON, so the search
  is provably not entered by any pinned state.

- **Tolerance:** asserted 1e-09 absolute on compositions and phase amounts
  against the independent Newtons (achieved <= 6.61e-12), 1e-09 on the
  equal-fugacity residual (achieved <= 8.0e-13), 1e-12 on the mass balance
  (achieved 0.0), 1e-06 K at `T3` on the verdict boundary (achieved 2.4e-07 K),
  1e-08 on FeOs's chemical potentials with matched constants (achieved
  <= 6.83e-11).
- **Independent route:** the 4-equation and 2-equation Newtons written in the
  test and in the example; the lever rule; reduced Gibbs energies from the
  public fugacity-coefficient interface; FeOs's chemical potentials and its own
  two-phase flash.
- **Negative control:** perturbing water's `epsilon^AB / k` by 1 % moves
  `x_water` in the hexane-rich phase by more than 1e-04 (slow-marked).
- **Test path:** `tests/test_flash_vlle_eos.py` (no optional dependency) and
  `tests/validation/test_pcsaft_vlle_water_hexane.py` (skipped without `feos`).
- **Script:** `examples/basic/flash_tp_pcsaft_vlle_demo.py`,
  `examples/validation/18_pcsaft_vlle_water_hexane.py`.

## Case P-10: Three-phase equation-of-state tie triangles (ternary)

- **Source:** independent 6-equation Newton solves written in the test (equal
  fugacity across three phases, each on its own named density branch,
  parameterized by `ln(x_k / x_last)` so every iterate stays inside the
  simplex), an independent mass-balance solve for the phase amounts, and
  **FeOs 0.10.1** chemical potentials at chemthermo's three phases.
- **Location:** `flash_tp(mixture, ..., eos=...)` on the tangent-plane path,
  slice `flash-phase-addition-eos` (ADR-0020).
- **Assumptions:** `k_ij = 0` throughout; 2B association for water and ethanol;
  packaged parameters. **Nothing here is compared against measurement**, and
  neither model is claimed to be right for these systems - what is checked is
  that the returned phase set is a verified equilibrium *of that model*,
  discovered rather than assumed.
- **Components / units:** water / ethanol / n-hexane at 101,325 Pa.

### (i) PC-SAFT vapour-liquid-liquid at 333 K

- **Expected outcome:** a tie triangle - the same three vertices from every
  feed inside it, with amounts fixed by the mass balance alone.
- **Achieved:** `phases = {"liquid1", "liquid2", "vapor"}`,
  `phase_regime = "VLLE"`, `phase_count = 3`,
  `phase_set_history = "L -> LL -> LLV"`, `phases_added = 1`,
  `phases_removed = 0`, `phase_label_method = "compressibility"`,
  `post_split_status = "stable"`, `vapor_fraction` = the vapour's fraction.
  Vertices:
  `liquid1 = (0.9613879781718806, 0.03842241608698442, 1.8960574113e-04)`,
  `liquid2 = (0.2553513856079896, 0.3913146548516096, 0.35333395954040087)`,
  `vapor = (0.194162215411462, 0.19904827565750266, 0.6067895089310353)`;
  densities 46,638.36 / 12,961.80 / 37.83 mol/m^3.
  `equilibrium_residual = 6.795e-12`, `mass_balance_residual = 0.0`,
  `delta_g_split_rt = -0.034026`, `delta_g_vs_two_phase_rt = -2.893e-04`
  (so `G3 < G2 < G1`).
- **Against the independent 6-equation Newton:** `|dx| <= 2.6e-11` on all
  three vertices, residual 1.7e-12.
- **Against the independent mass balance:** over nine feeds inside the
  triangle, `|d beta| <= 5.5e-11` and `|dx| <= 2.6e-11`; every feed returns the
  same three vertices.
- **FeOs's chemical potentials across the three phases:**
  `max |mu_i^a - mu_i^b| / RT` = **1.418e-11** with matched constants and
  **2.232e-06** as shipped, against an asserted 1e-8.
- **The region is finite, and it was scanned.** A 36-feed grid (mole fractions
  in steps of 0.1) at 333 K gives 9 `VLLE`, 7 `VLE`, 6 `LLE`, 13 single-phase
  and **1 `ConvergenceError`**; at 335 K, 9 `VLLE`, 13 `VLE`, 2 `LLE`, 12
  single-phase and **0** errors; at 337 K, 10 `VLLE`, 15 `VLE`, 11 single-phase
  and 0 errors; at 328 K and 331 K there is **no** three-phase region at all
  (22 `LLE`, 14 single-phase), which is consistent with the binary `T3` of
  334.81 K being the top of the two-liquid band.
- **The failing feed, recorded not worked around.** `z = (0.1, 0.1, 0.8)` at
  333 K raises `ConvergenceError` ("a two-phase set converged to a non-positive
  phase fraction"): the three-phase solve from the stability seed does not
  converge (residual 0.59 after the budget, fractions `(-0.76, -4.31, 6.07)`),
  and the two-phase set its removal leaves does not converge either. **That
  feed also raised before this slice** - with `max_phases = 2` it still does -
  so it is not a regression. It sits near the edge of the triangle where the
  water-rich liquid's amount is tiny.

### (ii) Peng-Robinson three *liquid* phases at 280 K

- **Premise:** ADR-0019 and Case P-8 recorded "no pure Peng-Robinson
  three-phase case found in the databank with `k_ij = 0`". A scan over 15
  ternaries x 6 temperatures x 3 pressures x 28 feeds found one.
- **Achieved:** water / ethanol / n-hexane, 280 K, 1 atm, `k_ij = 0` ->
  `phases = {"liquid1", "liquid2", "liquid3"}`, `vapor_fraction = None`,
  `phase_regime = "LLE"`, `phase_set_history = "L -> LL -> LLL"`,
  `phases_added = 1`, `phases_removed = 0`, `phase_label_method =
  "compressibility"`, `post_split_status = "stable"`, in **0.097 s**.
  `liquid1 = (0.9988087221058589, 1.1912778940878e-03, 5.324007536133e-14)`,
  `liquid2 = (0.10277098734252052, 0.8598745315402574, 0.037354481117222046)`,
  `liquid3 = (0.018917958864194933, 0.41169611644509935, 0.5693859246907057)`.
  `equilibrium_residual = 8.303e-09`, `mass_balance_residual = 1.11e-16`,
  `delta_g_split_rt = -0.250606`, `delta_g_vs_two_phase_rt = -2.905e-03`.
- **Against the independent 6-equation Newton:** residual 8.4e-15 and
  `|dx| <= 9.1e-09`; against the independent mass balance over four feeds,
  `|d beta| <= 1.2e-08`. Both are floored by chemthermo's own
  `equilibrium_residual` of 8.3e-09 (successive substitution stopped at
  `tol = 1e-8`; the second-order stage could not improve it, because
  `liquid1`'s n-hexane mole fraction of 5.3e-14 makes the Hessian nearly
  singular).
- **Model versus reality, stated and not asserted:** `k_ij = 0` between water
  and a hydrocarbon is not a serious parameterization and a three-liquid split
  for this ternary at 280 K is not a claim about the real system. What this
  sub-case establishes is that the search is model-agnostic and that a **cubic**
  reaches it, at a cost small enough for the default test suite.

- **Tolerance:** asserted 1e-08 absolute on the PC-SAFT vertices and amounts
  (achieved <= 5.5e-11), 1e-07 on the Peng-Robinson ones (achieved <= 1.2e-08),
  1e-09 on the PC-SAFT equal-fugacity residual (achieved 6.8e-12), 1e-08 on the
  Peng-Robinson one (achieved 8.3e-09), 1e-08 on FeOs's chemical potentials
  with matched constants (achieved 1.42e-11).
- **Independent route:** 6-equation Newton solves and a mass-balance solve
  written in the test; FeOs's chemical potentials.
- **Negative control:** at 328 K and 331 K the same ternary has **no**
  three-phase region, so the 333 K result is not something the search produces
  everywhere.
- **Test path:** `tests/test_flash_vlle_eos.py`
  (`test_peng_robinson_returns_three_liquid_phases` and
  `test_the_peng_robinson_triangle_matches_an_independent_newton` run by
  default; `test_the_pcsaft_ternary_returns_a_vapor_liquid_liquid_tie_triangle`
  is `slow`), `tests/validation/test_pcsaft_vlle_water_hexane.py`
  (`slow`, skipped without `feos`).
- **Script:** `examples/basic/flash_tp_pcsaft_vlle_demo.py --full`,
  `examples/validation/18_pcsaft_vlle_water_hexane.py --full`.

## Case P-11: Fixed density-root surfaces in the equation-of-state stability trials

- **Source:** the two independent 4-equation and 2-equation Newtons of Case
  P-9, re-run here; reduced Gibbs energies from the public
  `EquationOfState.fugacity_coefficients`; **FeOs 0.10.1** chemical potentials
  at chemthermo's converged phases; and, for the "nothing moved" half, a
  state-by-state diff of this repository against itself at HEAD `58190f5`
  (144-state Peng-Robinson grid, 188-state PC-SAFT Case F-4 grid, a 1144-state
  Peng-Robinson flash scan, and the 155-state bit-identity fixture).
- **Location:** `chemthermo/stability/_evaluator.py` and
  `chemthermo/stability/tp.py`; slice `stability-eos-root-surfaces`
  (ADR-0021). Nothing in `chemthermo/flash/`, `chemthermo/eos/` or
  `chemthermo/models/` changed.
- **Assumptions:** `k_ij = 0` throughout; packaged parameters; 2B association
  for water. Nothing here is compared against measurement.
- **Components / units:** water / n-hexane at 101,325 Pa (PC-SAFT) and the
  light-hydrocarbon grids of Cases F-1 and F-4 (Peng-Robinson, PC-SAFT).
  Temperatures in K, pressures in Pa, distances in units of RT.

### (i) The defect, and the repair

- **Expected outcome:** the hexane-rich liquid of Case P-9 (iv) is unstable
  towards a vapour, and `flash_tp` returns the lower-Gibbs pair.
- **Achieved, at the state Case P-9 (iv) pinned** (335 K, 1 atm,
  `w = (0.02273, 0.97727)`): `status = "unstable"`,
  `tpd_min = -6.5237494612e-03` at `w = (0.213559327395, 0.786440672605)`,
  `feed_branch = "liquid"`, `phase_branch = "vapor"`,
  `minimizing_trial = "wilson-vapor"`, `minimizing_trial_surface = "vapor"`,
  7 iterations. The other three trials are what the whole answer used to be:
  `wilson-liquid` and `pure-n-Hexane` trivial, `pure-Water` on a
  `tpd = +6.26e-05` liquid-surface point.
- **The flash that follows**, `z_water = 0.7`, 335 K:
  `phases = {"liquid", "vapor"}`, `phase_regime = "VLE"`,
  `liquid = (0.999936102393, 6.3897607e-05)`,
  `vapor = (0.214999486349, 0.785000513651)`,
  `beta_vapor = 0.382115060327`, `phase_set_history = "L -> LL -> LLV -> LV"`.
  Reduced Gibbs energies, from the public interface and the lever rule:
  `G(VL)/RT = -1.164305930752` against the two-liquid pair's
  `-1.161810713699`, a gap of **2.495217e-03 RT** - the ledger's own P-9 (iv)
  number, reproduced from the other side.
- **Against FeOs** at chemthermo's converged phases and densities, with FeOs's
  fourteen-figure universal constants substituted and the flash re-run on them:
  `max_i |mu_i^L - mu_i^V| / RT` = **9.04e-12** at 335 K, **3.62e-11** at
  `T3 + 0.05 K` and **2.75e-11** at `T3 + 0.5 K`, against an asserted 1e-08. As
  shipped (the ten-figure table), the same comparisons are 2.496e-06,
  2.498e-06 and 2.491e-06 - the constants difference, not either solver, as in
  Cases P-6 to P-9.
- **The two offsets the slice brief named**, both `z_water = 0.7` (quoted to
  twelve decimals from the example's own run; the last two digits move with
  the `T3` the 4-equation Newton returns, which is why they are not pinned in
  a test): at `T3 + 0.05 K` the flash returns
  `liquid = (0.999936007302, ...)`, `vapor = (0.213610914456, ...)`,
  `beta_vapor = 0.381440208422`, with
  `G(VL)/RT = -1.168271522770` against the two liquids' `-1.167622708774`
  (gap **6.488140e-04 RT**); at `T3 + 0.5 K`,
  `liquid = (0.999936310223, ...)`, `vapor = (0.218032446374, ...)`,
  `beta_vapor = 0.383597426858`, `G(VL)/RT = -1.155747388671` against
  `-1.149246614657` (gap **6.500774e-03 RT**).
- **The scans:** 41 temperatures per feed across `[T3 - 1 K, T3 + 1 K]`,
  **82 states, zero `ConvergenceError`**, **exactly one verdict switch per
  feed** (`z_water = 0.3`: 20 `LLE` then 21 `VLE`; `z_water = 0.7`: 21 then
  20). Bisecting the whole `[T3 - 1 K, T3 + 1 K]` window 25 times (resolution
  6e-08 K, scratch measurement) puts the switch **2.09e-07 K below** the
  independent `T3 = 334.807826336 K` at `z_water = 0.3` and **2.68e-07 K
  above** it at `z_water = 0.7`. The shipped example bisects a 2e-03 K bracket
  to 1e-06 K instead - 11 flashes per feed rather than 25, each one a
  phase-addition search - and therefore reports `|T3 - boundary| = 4.88e-07 K`
  at both feeds, which is that bisection's half-width and not a different
  answer.

### (ii) What did not move

- **The 144-state Peng-Robinson stability grid** (Cases F-1 / F-5 geometry):
  **0 verdict changes** (47 unstable / 97 stable, before and after), **0
  branch-label changes**, worst `|delta tpd_min|` where both runs found a
  stationary point **1.78e-15**, worst `|delta w|` **7.27e-12**.
- **The 188-state PC-SAFT grid** (Case F-4): **0 verdict changes** (123
  unstable / 65 stable), **0 branch-label changes**, worst
  `|delta tpd_min|` **1.24e-12**, worst `|delta w|` **2.22e-15**.
- **A 1144-state Peng-Robinson `flash_tp` scan** (8 databank mixtures,
  T 150-450 K in 11 steps, P 1e5-3.2e7 Pa in 13 geometric steps - the shape of
  Case F-1's wider scan): **0 verdict changes**, **0 `ConvergenceError` before
  and after**, worst `|delta composition|` **1.11e-15**, phase fractions
  **identical**.
- **The one thing that does change, and why it cannot change a verdict.** On
  18 of the 144 Peng-Robinson grid states, 38 of the 188 PC-SAFT ones and 118
  of the 1144 scan states, `tpd_min` moves from `0.0` - the value the summary
  reports when *no* non-trivial stationary point was reachable - to a
  **positive** number (smallest seen, 0.0968 on the 1144-state scan). A
  vapour-root-pinned trial now reaches a real stationary point that lies
  *above* the tangent plane. A positive `tpd_min` decides nothing: every one of
  those states is single-phase and `"stable"` before and after.
- **The bit-identity fixture, audited state by state before regeneration**
  (`refactor_bit_identity_v2.json` -> `v3.json`, 155 states): **122
  bit-identical**; **18** changed `diagnostics["tpd_min"]` from `0.0` to a
  positive number and *nothing else*; **15** changed only in the last bits of
  iterative quantities - worst composition move **1.11e-16** (one ulp), worst
  diagnostics move **5.33e-15** (on `max_delta_k`), against an audited
  tolerance of 1e-09. **No state changed its phase names, its phase set, any
  composition or any phase fraction beyond those last bits, and no state gained
  or lost a diagnostics key.** `v1` and `v2` are kept, unreferenced, for
  history.

### (iii) Trial statistics, and two decisions they settled

- **Who finds the minimizer.** Over the 144-state Peng-Robinson grid the
  minimizing trial runs on the **vapour** root on 32 states and the **liquid**
  root on 45 (by label: `wilson-vapor` 32, `wilson-liquid` 18,
  `pure-<name>` 27). Over the 188-state PC-SAFT grid: vapour 132, liquid 48
  (`wilson-vapor` 132, `pure-Methane` 29, `wilson-liquid` 13,
  `pure-n-Hexane` 6). Neither surface is decorative.
- **How often the fallback fires.** Measured once with an instrumented
  both-branches build (not what ships; see ADR-0021 decision 4): **3889
  single-root evaluations against 4378 trial iterations** on the Peng-Robinson
  grid and **14975 against 16898** on the PC-SAFT one. The two counters have
  different denominators (the second-order stage evaluates the model several
  times per iteration), so this is a ratio of about nine in ten and not an
  exact percentage. The single-root case is the common one, so The fixed surface therefore changes the
  iteration at about one evaluation in nine, and that ninth is where
  P-9 (iv) lived. What *ships* counts the same condition where the solver
  compares the branches anyway - the trial's stopping point - so on the
  Peng-Robinson grid 448 of 624 trials and on the PC-SAFT grid 589 of 752 stop
  in a one-root region.
- **Iteration counts.** Total trial iterations rose from 4080 to 4378 (+7.3 %)
  on the Peng-Robinson grid and from 16291 to 16898 (+3.7 %) on the PC-SAFT
  one - trials that used to collapse onto the feed now walk to a real
  stationary point. Wall time went the other way, because a pinned trial asks
  the model for one root instead of two: `stability_tp` over the 188-state
  PC-SAFT grid **76.6 s -> 44.3 s** (-42 %), and the default test suite
  **244.2 s -> 231.7 s** before the new golden path is added; with it,
  **247.3 s for 750 tests** against 244.2 s for 738. `pytest -q -m slow` is
  **827 s (13:47) for 22 tests**, against 855 s for 21.
- **Pure-component starts on the vapour root: measured, then not added.**
  Running the `n` pure-component-dominant estimates on the vapour root as well
  was tried over all 334 states above: **0 verdict changes** and **not one
  state whose `tpd_min` fell by more than 2.51e-13** (they are the nominal
  minimizer on 41 states only by tying to the last bits with a trial that
  already found the same point). The trial count therefore stays `n + 2`.

### Negative controls and what is *not* claimed

- The ternary water / ethanol / n-hexane feed `z = (0.1, 0.1, 0.8)` at 333 K,
  which Case P-10 recorded as raising, **still raises** the same
  `ConvergenceError` ("a two-phase set converged to a non-positive phase
  fraction") in 5.6 s. It is a multiphase-solver failure at the edge of the
  tie triangle, not a stability miss, and this slice does not touch it.
- `stability_tp`'s honesty note is unchanged. Fixing the surfaces enlarges the
  set of reachable stationary points; it does not turn a local search into a
  global proof.
- Nothing here is compared against measurement, and `k_ij = 0` between water
  and a hydrocarbon is not a serious parameterization.

- **Tolerance:** asserted 1e-12 absolute on the pinned `tpd_min` and minimizing
  composition (achieved exactly, the numbers are pinned); 1e-09 on every
  previously validated composition, fraction and `tpd_min` (achieved
  <= 7.27e-12 on compositions, <= 1.24e-12 on `tpd_min`, 0.0 on phase
  fractions over the 1144-state scan); 1e-08 on FeOs's chemical potentials with
  matched constants (achieved <= 3.62e-11); 1e-06 K on the bisected verdict
  boundary against the independent `T3` (achieved <= 2.68e-07 K).
- **Independent route:** the 4-equation and 2-equation Newtons of Case P-9,
  written again in the example; the lever rule; reduced Gibbs energies from the
  public fugacity-coefficient interface; FeOs's chemical potentials; and the
  repository at HEAD `58190f5` as its own before-state.
- **Negative control:** the 1144-state Peng-Robinson scan and the 155-state
  bit-identity fixture, which must *not* move - and do not.
- **Test path:** `tests/test_stability_eos_surfaces.py` (11 by default, 1
  `slow`), `tests/test_flash_refactor_bit_identity.py` (v3),
  `tests/test_flash_vlle_eos.py::test_the_window_scan_never_raises` (`slow`),
  `tests/test_stability_candidates.py`, `tests/test_stability_tp.py`
  (pinned Peng-Robinson numbers, unchanged).
- **Script:** `examples/validation/19_eos_stability_surfaces.py` (about 5 s;
  `--full` for both offsets, the matched-constants comparison, the two
  41-point scans, the bisected boundary and the Peng-Robinson grid).

### Cross-platform (ADR-0032, 2026-09-25)

The pins above were captured on macOS arm64. On the first Linux runs
(GitHub Actions `ubuntu-latest` run 36121300741, and the cloud-baseline host:
Linux 6.18 x86_64, CPython 3.11.15, numpy 2.4.6) three captured-value guards
failed and nothing else did:

| guard | Linux measurement | discrete fields |
| --- | --- | --- |
| `refactor_bit_identity_v3.json`, 155 states | deterministic run to run; 83 states / 688 floats move; worst 2.44e-13 absolute (`k_min`); compositions <= 1.2e-14, phase fractions <= 2.1e-13, `tpd_min` <= 3.6e-15 | all identical: phase names, key sets, statuses, stages, iteration counts |
| PC-SAFT n-hexane literals (Cases P-1/P-3) | cloud host: `Z(7700)` 32 ULP (3.5e-15), `ln phi(7700)` 2 ULP, rest exact; CI runner: `a_res(7700)` 3 ULP | n/a |
| grid `minimizing_trial_surface` count | liquid/vapour 45/32 (macOS), 46/31 (CI), 47/30 (cloud) | verdicts 47/97 identical |

The surface count is last-bit-decided on 20 of its 77 states: both surfaces
reach the same stationary point (best `tpd` equal to <= 5e-16, trial
compositions to <= 6e-13). On the other 57 the winner is decisive
(19 vapour, 38 liquid; the two decisive two-surface gaps are 1.7e-02 and 1.0).

- **Rule adopted (ADR-0032):** exact on macOS arm64; elsewhere discrete fields
  exact, floats to `atol = rtol = 1e-12` (flash fixture) and `atol = 5e-14`
  (PC-SAFT literals), and the surface statement is `{vapor: 19, liquid: 38,
  tie: 20}` with `TIE_MARGIN = 1e-10`, the 45/32 split still pinned on macOS.
- **Negative control:** `tests/test_capture_identity.py` - a 1e-11 move, any
  discrete change, a key added or removed, and a non-finite mismatch all fail
  the bounded comparison; a last-bit move passes.
- **Not claimed:** bit-identity on Linux. The in-process A/B guards
  (ADR-0023/0030) remain exact on every platform.

---

## Case P-12: PC-SAFT properties for a polymer, against FeOs

- **Source:** **FeOs 0.10.1** (feos-org/feos, MIT OR Apache-2.0), the same
  Gross & Sadowski model in Rust with every derivative by automatic
  differentiation. FeOs packages no polymer parameters but takes a segment
  number directly, so it is given the **derived** `m = (m/M) Mw` this package
  computes from the segments-per-mass record - which makes the comparison a
  check of the convention as well as of the equations.
- **Location:** `chemthermo/parameters/pcsaft.py` (`PCSAFTRecord.segments_per_g`),
  `chemthermo/core/component.py` (`Component.custom`), `chemthermo/eos/pcsaft.py`;
  slice `pcsaft-polymer-solvent` (ADR-0022).
- **Parameters and provenance.** Polyethylene: `m/M = 0.0263` mol/g,
  `sigma = 4.0217 A`, `eps/k = 247.5 K`. n-pentane: the packaged Gross &
  Sadowski (2001) record, `m = 2.6896`, `sigma = 3.7729 A`, `eps/k = 231.20 K`.
  The polymer row is **as tabulated by Martini, Cismondi, Barbosa & Brignole,
  *Sep. Sci. Technol.* 44(11) (2009) (author manuscript, CONICET open
  repository), Table 1, citing Gross & Sadowski, *IECR* 41 (2002) 1084 - and is
  NOT verified against that primary table**, which is paywalled and was not
  read. A search for a second open source printing the same three numbers found
  none (FeOs and Clapeyron.jl carry no polymer records; every open hit leads
  back to this one manuscript). The values are therefore a cited **test
  fixture**, `tests/fixtures/pcsaft/martini2009_polymers.json`, and are never
  packaged runtime data. An earlier recollection of `eps/k = 252.0` for PE and
  `m/M = 0.0205` for PS is recorded in the fixture's notes as *not* reproduced
  anywhere.
- **Assumptions:** `k_ij = 0` on both sides (FeOs cannot be given one - see
  Case P-13); the polymer is **monodisperse**, one chain length and one
  component, while the samples these parameters describe have polydispersities
  of 1.14 to 2.94. Nothing here is compared against measurement.
- **Components / units:** polyethylene (`Mw = 16400` and `53000` g/mol, so
  `m = 431.32` and `1393.9`) and n-pentane at 453.0 K. Densities in mol/m^3,
  pressures in Pa, everything compared dimensionless.

### (i) Segment number from the mass-based parameter

- **Expected outcome:** `m = (m/M) * Mw`, derived in the record rather than by
  the caller.
- **Achieved:** `431.32` for `Mw = 16400` and `1393.9` for `Mw = 53000`, equal
  to `0.0263 * Mw` to the last bit. A record giving both `m` and
  `segments_per_g`, or neither, or `segments_per_g` with no `MW_g_mol`, raises
  `PCSAFTParameterError`.

### (ii) Density roots of the pure melt

- **Expected outcome:** one mechanically stable root at every pressure over
  1-30 MPa, of a plausible melt magnitude, measuring as a liquid.
- **Achieved** (453.0 K, `Mw = 16400`): one root at every pressure, densities
  `46.1654`, `46.9073`, `47.6454`, `48.3130` mol/m^3 at 1, 10, 20 and 30 MPa,
  i.e. **0.757112, 0.769280, 0.781384, 0.792333 g/cm^3**. For `Mw = 53000` the
  same pressures give 0.757789 to 0.792889 g/cm^3. Monotone in pressure, no
  overflow anywhere in the `eta` scan, and `phase_identity` returns `"liquid"`
  at each - independently confirmed by a `kappa = P/(rho dP/drho)` recomputed
  in the test from the public `pressure_Pa`, all far below the 0.5 threshold.
- **Unverified remark, not an assertion:** commonly tabulated polyethylene melt
  densities near 450 K are around **0.77-0.80 g/cm^3**. That figure was not
  read from a source here; it is recorded to say the magnitude is not absurd,
  and the test brackets `0.75 < rho < 0.80` as a sanity range only.

### (iii) `A^res/RT`, `Z` and `ln phi` against FeOs

- **Expected outcome:** agreement to 1e-10 with matched universal constants,
  over eleven states whose densities are chemthermo's **own** liquid roots
  (pure melt at both molar masses at 1 / 10 / 30 MPa; 5, 10 and 15 wt% polymer
  in n-pentane at 10 MPa; 5 and 15 wt% at 15 MPa).
- **Achieved, matched constants:** worst `|dA^res/RT| = 9.10e-13`,
  `|dZ| = 5.94e-12`, `max |d ln phi| = 5.00e-12`.
- **Achieved, as shipped:** `2.85e-07`, `1.12e-06`, `1.41e-06`. That residual
  is the 42 universal constants of the 2001 dispersion term (chemthermo
  packages the ten printed figures, FeOs hard-codes fourteen) - an *input*
  difference, as in Cases P-6 to P-11.
- **Cross-check on the root itself:** FeOs's pressure at chemthermo's melt
  density reproduces the target pressure to better than 1e-5 relative at every
  melt state.

### (iv) The exponential's range

- **Expected outcome:** the model is finite where `exp(ln phi)` is not, and the
  ADR-0022 guard is what carries it.
- **Achieved:** polyethylene `Mw = 53000` (`m = 1393.9`) at 5 wt% in n-pentane,
  453 K, 10 MPa: `ln phi_polymer = -1690.615`, and `exp` of that is an exact
  `0.0` (the smallest positive double is `exp(-744.44)`). At 8 MPa the same
  quantity is `-1678.189`. `PCSAFTEOS.log_fugacity_coefficients` returns the
  same doubles as the `(T, rho, x)` route, `==`.
- **Before the guard:** `stability_tp` raised `ModelError("No usable
  fugacity-coefficient branch for stability analysis (vapor: non-finite or
  non-positive fugacity coefficients; liquid: ...)")` at every pressure tried
  (5, 10, 15, 20 MPa).
- **Guard instrumentation** (counted at
  `chemthermo.flash._common.eos_branch_terms`): PC-SAFT water / n-hexane at
  298.15 K / 1 atm, **0 of 274** branch evaluations in log space; the
  `Mw = 16400` polymer flash at 8 MPa, **0 of 600**; the `Mw = 53000` one,
  **1025 of 1025**. The 144-state Peng-Robinson replay test additionally
  asserts `terms.phi is not None` at every split iterate.

- **Tolerance:** asserted 1e-10 on `A^res/RT`, `Z` and `ln phi` with matched
  constants (achieved <= 5.94e-12); 1e-5 as shipped (achieved 1.41e-06); 1e-5
  relative on the melt pressure round trip.
- **Independent route:** FeOs, which shares no code and takes every derivative
  by automatic differentiation.
- **Negative control:** the packaged parameter set contains no component whose
  name matches `poly` (asserted), so nothing here reaches runtime data.
- **Test path:** `tests/validation/test_pcsaft_polymer_vs_feos.py`,
  `tests/test_pcsaft_polymer.py`.
- **Script:** `examples/validation/20_pcsaft_polymer_vs_feos.py` (about 4 s;
  `--full` for the FeOs-flash survey, the `Mw = 53000` chain and the cloud
  point).

---

## Case P-13: The polyethylene / n-pentane liquid-liquid split

- **Source:** an equal-fugacity Newton written independently in the test and in
  the example (two unknowns, carried as logarithms so the 1e-05 branch keeps
  its digits, finite-difference Jacobian, started a thousandth away from
  `flash_tp`'s answer); **FeOs 0.10.1** chemical potentials at chemthermo's
  converged phases and densities; and FeOs's own `tp_flash` at the one pressure
  where it converges on this system.
- **Location:** `chemthermo/flash/_split.py`, `chemthermo/flash/_common.py`,
  `chemthermo/stability/_evaluator.py`; slice `pcsaft-polymer-solvent`
  (ADR-0022). No solver equation changed.
- **Parameters and provenance:** as Case P-12, plus `k_ij = -0.006`, which the
  same secondary source's Table 3 gives for the `Mw = 16400` (polydispersity
  1.16) sample and which **that source fitted** to the cloud-point data of
  Kiran & Zhuang, *Polymer* 33 (1992) 5259. Those experimental data were not
  read here and nothing below is compared against them.
- **Assumptions:** monodisperse polymer; 453.0 K throughout; feeds stated as
  polymer **mass** fractions and converted with the component's molar mass.

### (i) The cloud point, by bisecting the stability verdict

- **Expected outcome:** two phases at low pressure, one phase at high pressure,
  one switch.
- **Achieved** (5 wt% polymer, `k_ij = -0.006`): `unstable` at 5.0 and 8.0 MPa,
  `stable` at 10.0, 15.0, 20.0, 25.0 and 30.0 MPa. Bisecting the verdict to
  1 Pa gives **9,751,759 Pa = 9.7518 MPa = 97.518 bar**.
- The cited source's figures put separation pressures for this system in the
  tens-to-few-hundred bar range. **The figures were not digitized**; this is a
  magnitude remark, not a comparison, and no number here is fitted to them.

### (ii) The split at 8 MPa

- **Expected outcome:** two liquid phases, all three verification residuals
  satisfied, phase fractions in `(0, 1)`, post-split stable.
- **Achieved:** `liquid1` / `liquid2`, `vapor_fraction = None`,
  `phase_regime = "LLE"`, `phase_label_method = "compressibility"`,
  `phase_i_branch = phase_ii_branch = "liquid"`.
  - polymer-rich phase: `x_polymer = 9.751377e-04`, fraction `0.2282983`,
    `rho = 5782.65` mol/m^3, `kappa = 6.805e-02`;
  - solvent-rich phase: `x_polymer = 1.1478711e-05`, fraction `0.7717017`,
    `rho = 6218.31` mol/m^3, `kappa = 1.225e-01`;
  - feed `x_polymer = 2.314804e-04`.
  - `mass_balance_residual = 2.71e-20`, `fugacity_residual = 3.41e-13`,
    `delta_g_split_rt = -1.6914873e-04`, `post_split_status = "stable"` with
    `post_split_tpd_min = -1.30e-14`, `converged_stage = "second-order"`
    (100 successive substitutions then 6 Newton steps).
  - Both `kappa` values are recomputed in the test from the public
    `pressure_Pa` by finite difference, so the liquid identities do not rest on
    `phase_identity`'s own derivative.
- **Independent Newton:** `x_polymer = 1.1478711124e-05` and `9.7513772864e-04`,
  i.e. **max `|dx_polymer| = 1.70e-15`** against `flash_tp`, with its own
  residual at 1.0e-12.
- **A tie line, not a feed property:** a 10 wt% feed at the same state returns
  the same two compositions to 1e-9 relative.
- **FeOs chemical potentials at chemthermo's phases** (8 MPa, `k_ij = 0` on
  both sides): **4.55e-13** with matched universal constants, **4.13e-08** as
  shipped. FeOs's own flash *raises* at this state, and this route does not
  need it to converge - only to evaluate.

### (iii) FeOs's own flash, and where it is not a reference

- **Measured**, 5 wt% feed, 453 K, `k_ij = 0`, feos 0.10.1:

  | P / MPa | FeOs `tp_flash` | chemthermo |
  | --- | --- | --- |
  | 3 | degenerate: both phases equal to nine figures, vapour fraction exactly 0.5 | two liquids, `dG/RT = -3.41e-03` |
  | 5 | `RuntimeError: `rachford_rice` encountered illegal values during the iteration` | two liquids, `dG/RT = -1.35e-03` |
  | 8 | same `RuntimeError` | two liquids, `dG/RT = -1.69e-04` |
  | 10 | converges | two liquids |
  | 11, 15 | `RuntimeError: No phase split according to stability analysis` | one liquid |

- **At 10 MPa**, the one usable comparison: the polymer-rich branch agrees to
  `|dx| = 5.77e-13` (matched constants; 1.88e-12 as shipped), the solvent-rich
  branch to `1.23e-09`, densities to `9.3e-08` relative and the phase fraction
  to `1.41e-06`. That is **not** the 1e-11 of Case P-12 because 10 MPa is
  within 8 % of this system's `k_ij = 0` cloud point (10.768 MPa) and both
  solvers are near a plait point: measured in chemthermo's own model,
  chemthermo's converged pair has an equal-fugacity residual of **1.06e-07**
  and FeOs's has **9.07e-06**, so chemthermo's is the tighter stationary point
  by about 86x. Recorded as conditioning, not as disagreement.
- **`k_ij` could not be given to FeOs.** `BinaryRecord.from_json_str` and
  `Parameters.from_records` both accept one, but `EquationOfState.pcsaft` then
  raises `RuntimeError: missing field `k_ij`` in feos 0.10.1 for every
  serialization tried (bare float, `{"k_ij": x}`, `{"k_ij": [x]}`, with and
  without `l_ij`, and via `Parameters.new_binary`). Every FeOs comparison in
  this case therefore runs at `k_ij = 0` on both sides.

### (iv) LCST-type behaviour and the `k_ij` trend (both qualitative)

- **At a fixed 10 MPa**, 5 wt%: `stable` at 400, 405, ..., 450 and 453 K;
  `unstable` at 455 and 460 K. Bisecting gives a switch at **454.56 K**. The
  polymer comes out of solution on **heating**, which is the direction the
  cited manuscript's Figure 1 describes (the region above its cloud-point curve
  is single phase). **The figure was not digitized**; the check is qualitative
  and is labelled as such in the test and the examples.
- **Cloud-point pressure against `k_ij`** (5 wt%, 453 K, bisected to 1 kPa):
  `k_ij = -0.006` -> **97.518 bar**, `k_ij = 0` -> **107.680 bar**,
  `k_ij = +0.02` -> **215.074 bar**. Monotone increasing, which is the
  manuscript's own statement. Again qualitative: no value of theirs is
  reproduced.

### (v) Asymmetry stress

- **Permutation invariance:** reversing the component order moves the two
  phase compositions by at most **3.47e-15** and the phase fractions by
  **9.31e-13**.
- **Determinism:** two identical calls return `==` compositions and `==` phase
  fractions.
- **Extreme feeds at 8 MPa:** 0.1 wt% (`x_polymer = 4.40e-06`) and 40 wt%
  (`x_polymer = 2.92e-03`) both return a **single stable liquid** and neither
  raises. Both are outside the binodal, which runs from `1.15e-05` to
  `9.75e-04` there - the dilute feed below the solvent-rich branch and the
  concentrated one above the polymer-rich branch.
- **Three components:** polyethylene / n-pentane / n-hexane (5 wt% polymer,
  equal solvent masses) at 5 MPa returns two liquids with
  `mass_balance_residual = 5.55e-17`, `fugacity_residual = 8.06e-08`,
  `dG/RT = -1.73e-05`, post-split stable.
- **The long chain:** `Mw = 53000` (`m = 1393.9`) at 8 MPa returns
  `x_polymer = 3.729971e-04` and `6.269597e-10`, fractions `0.192063` /
  `0.807937`, `mass_balance_residual = 0.0`,
  `fugacity_residual = 4.55e-13`, `dG/RT = -4.1379e-04`, post-split stable.
  Every one of its 1025 branch evaluations runs in log space (Case P-12 (iv)).
- **No post-split refusal was observed** on any state in this case: the worst
  `post_split_tpd_min` measured is `-1.30e-14`, inside `tpd_tol = 1e-8` by six
  orders.
- **The second-order stage's mole-number box is not reached.** The smallest
  mole number it sees on these systems is `1.2e-20` (the `Mw = 53000`
  solvent-rich phase at 5 MPa), against its existing `1e-300` floor. The guard
  the slice design anticipated there was therefore **not added**: there is no
  measured failure behind it.

### (vi) A gap in the phi-phi split, pinned and not patched - **RESOLVED by ADR-0024**

> **Amendment (slice `pcsaft-polymer-vle`, ADR-0024).** Both states below now
> converge and are pinned the other way round in **Case P-14**. The text is
> kept as written because it is the measured record of what the pre-ADR-0024
> code produced, and because the diagnosis in it - "a gap in the split, not a
> property of the polymer support" - is what the fix turned out to confirm. The
> two numbered facts that changed: `flash_tp` at 1 and 2 MPa returns a verified
> vapour-liquid split, and the ternary at 3 MPa returns a verified
> liquid-liquid one. The named later candidate
> `pcsaft-polymer-vle-ethylene` was delivered as `pcsaft-polymer-vle`.

- **Below about 3 MPa** at 453 K, n-pentane is subcritical and the mixture has
  **two** density roots, so the equilibrium in question is vapour-liquid rather
  than liquid-liquid. `stability_tp` still reports the 5 wt% feed `unstable`
  (2 roots at 1 and 2 MPa), and the split then stops after **one**
  successive-substitution step with an equal-fugacity residual of **2.99e+02**
  (1 MPa) and **3.19e+02** (2 MPa); `flash_tp` raises `ConvergenceError`.
- **The ternary does the same at 3 MPa**, running away to
  `beta = -6.33e+10` and raising the ADR-0016 "vapor fraction outside (0, 1)"
  error. At 5 MPa and above it converges (see (v)).
- This is a gap in the **phi-phi split** for polymer/solvent vapour-liquid
  states, not a property of the polymer support ADR-0022 adds, and it is not
  patched here. Both states are pinned by test so that any future change to
  them is deliberate. `pcsaft-polymer-vle-ethylene` is the named later
  candidate.

- **Tolerance:** asserted 1e-12 absolute on the independent Newton's tie line
  (achieved 1.70e-15); 1e-12 on the mass balance (achieved 2.71e-20); 1e-8 on
  the fugacity residual (achieved 3.41e-13); `dG < 0` and post-split stable;
  1e-8 on FeOs's chemical potentials with matched constants (achieved
  4.55e-13); 1e-8 absolute on the FeOs tie line at 10 MPa (achieved 1.23e-09);
  1e-13 on permutation invariance (achieved 3.47e-15).
- **Independent route:** the two-equation Newton written in the test and again
  in the example; FeOs's chemical potentials; FeOs's own flash where it
  converges.
- **Negative control:** two states that must keep raising (`1 MPa` binary,
  `3 MPa` ternary) - **superseded by Case P-14**, where the same two states are
  the positive result; the 155-state bit-identity fixture, which does not move,
  and the guard instrumentation showing the log-space route is dormant on every
  pre-ADR-0022 state.
- **Test path:** `tests/test_pcsaft_polymer.py` (23 by default, 5 `slow`),
  `tests/validation/test_pcsaft_polymer_vs_feos.py` (5 by default, 1 `slow`),
  `tests/test_component_custom.py`.
- **Script:** `examples/basic/pcsaft_polymer_demo.py` (about 5 s; `--full` for
  the cloud-point bisection and the `Mw = 53000` chain) and
  `examples/validation/20_pcsaft_polymer_vs_feos.py`.

---

## Case P-14: The polyethylene / n-pentane vapour-liquid split, in log mole numbers

- **Source:** a **one-dimensional** equal-fugacity solve written independently
  in the test, in the validation test and in the example (the vapour taken as
  *exactly* pure solvent, one unknown carried as `ln x_solvent`,
  finite-difference Newton started a thousandth away from `flash_tp`'s answer);
  **FeOs 0.10.1** chemical potentials at chemthermo's converged phases and
  densities.
- **Location:** `chemthermo/flash/_log_space.py` (new),
  `chemthermo/flash/_detect.py`; slice `pcsaft-polymer-vle` (ADR-0024). No
  model equation and no stability-solver equation changed.
- **Parameters and provenance:** exactly Case P-13's - `k_ij = -0.006` for
  routes that use one, `k_ij = 0` for every FeOs comparison because feos 0.10.1
  cannot be given one from Python. Nothing here is compared against
  measurement.
- **Assumptions:** monodisperse polymer; 453.0 K throughout; feeds stated as
  polymer **mass** fractions.

### (i) The three vapour-liquid states, and what they are

The feed is `z = (2.3148042237948336e-04, 0.9997685195776205)` (5 wt% polymer).
Below about 2.59 MPa the mixture has two density roots at the feed, the
min-Gibbs one is the **vapour**, and the answer is a solvent vapour over a
solvent-swollen melt. All three converge with `k_seed = "stability-log"` and
`converged_stage = "second-order-log"`; both phases are named from a measured
compressibility.

| P / MPa | melt `x` (polymer, solvent) | melt, wt% solvent | melt phase fraction | `ln y_polymer` | log-space iterations |
| --- | --- | --- | --- | --- | --- |
| 0.5 | `(0.06482279502229682, 0.9351772049777032)` | 5.97 | 0.003570972561424801 | **-466.9816656241** | 7 |
| 1.0 | `(0.02853164733532448, 0.9714683526646756)` | 13.03 | 0.008113111018756114 | **-450.5307940617** | 12 |
| 2.0 | `(0.008772681701263427, 0.9912273182987366)` | 33.20 | 0.026386506459723290 | **-408.9441806100** | 31 |

- Residuals, in the same order: `fugacity_residual` 5.68e-13 / 7.39e-13 /
  2.27e-13 (identical to the stage's own `log_space_residual` here, because no
  mole fraction underflowed); `mass_balance_residual` 2.77e-18 / 1.22e-18 /
  3.52e-19; `delta_g_split_rt` -1.0575e-01 / -1.0172e-01 / -9.1180e-02; every
  phase post-split **stable**.
- The compressibility identity of ADR-0017, computed from the public
  `pressure_Pa` routine rather than from `phase_identity`: vapour `kappa`
  1.0680 / 1.1603 / 1.5282, melt `kappa` 0.0011 / 0.0024 / 0.0068 against the
  0.5 threshold.
- **Before ADR-0024** every one of these raised `ConvergenceError`
  (Case P-13 (vi)).

### (ii) The independent one-dimensional solve

- **Expected outcome:** because the vapour is pure solvent to `1e-196`, the
  melt must satisfy `ln phi_s^V(pure solvent vapour) = ln x_s + ln phi_s^L(x)`
  exactly, in one unknown. Nothing of the flash is used - no stability test,
  no Rachford-Rice, no stage, no phase-count logic.
- **Achieved:** `|x_solvent(flash) - x_solvent(1-D)|` = **2.53e-14** (0.5 MPa),
  **5.55e-15** (1 MPa), **1.22e-15** (2 MPa), against an asserted 1e-10. The
  flash's own vapour-phase solvent fugacity equals the pure-solvent value to
  better than 1e-08 at all three.
- **The polymer's own condition**, which cannot be written in linear mole
  numbers at all: `ln x_PE + ln phi_PE^L` = **-492.6311005272** (0.5 MPa),
  **-504.1563775561** (1 MPa), **-530.5415518401** (2 MPa), against
  `ln y_PE + ln phi_PE^V` to **5.68e-13 / 7.39e-13 / 2.27e-13**.

### (iii) FeOs at chemthermo's phases

- **Achieved** (`k_ij = 0` on both sides, the flash re-run there): the worst
  `|d mu_i / RT|` between the two phases over the three states is
  **2.615e-12** with FeOs's fourteen-figure universal constants substituted
  into chemthermo, and **2.537e-07** as shipped. The difference between those
  two is the constants table, not the model - the same floor Cases P-6 to P-12
  record. Asserted 1e-08 with matched constants.
- FeOs evaluates a state at `y_polymer = 1e-196` without complaint, and
  reproduces `P` at both of chemthermo's densities to 1e-06 relative.

### (iv) An exactly zero mole fraction

- For `Mw = 53000` (`m = 1393.9`) at 2 MPa the vapour's polymer mole fraction
  is `exp(-1314.9852395489)`, which is **not a double**. The result carries
  `y = (0.0, 1.0)`, `diagnostics["log_space_zero_fractions"] = 1` and the
  logarithm in `["log_space_ln_x_min"]` (ADR-0024 decision 3). The melt is
  `(0.0027736965292611865, 0.9972263034707388)` at a phase fraction of
  0.025828116111173305.
- Two consequences, both measured: the material balance is **exact** rather
  than approximate (`1.11e-16`, i.e. round-off in the sum only), because the
  melt holds every mole of polymer the feed had; and `fugacity_residual`
  (2.29e-14) does **not** see the polymer, because it is taken over components
  present in both phases - the stage's own `log_space_residual` (9.09e-13)
  does. `delta_g_split_rt = -9.2125e-02`, post-split stable.
- This state reached the log-space stage through the *third* trigger: the
  K-loop raised `ModelError` because the two phases' `ln phi` differ by more
  than the exponential's range.

### (v) The ternary, and the second entry point

- Polyethylene / n-pentane / n-hexane,
  `z = (2.519896558373394e-04, 0.5441741521284267, 0.45557385821573587)`, 3 MPa:
  before ADR-0024 successive substitution converged on the **trivial**
  solution (`x = y` to twelve figures, `beta = -6.33e+10`) and `flash_tp`
  raised. Handing that to the second-order stage gives two liquids,
  `(1.3272336238799715e-03, 0.5402484470355627, 0.45842431934055733)` at a
  phase fraction of 0.18872173701842598 and
  `(1.863504765793684e-06, 0.5450873602321722, 0.45491077626306187)` at
  0.811278262981574, with `fugacity_residual` 6.82e-13,
  `mass_balance_residual` 5.55e-17, `delta_g_split_rt = -4.5325e-04`,
  post-split stable.
- The stage used is the **linear** one (`converged_stage = "second-order"`, no
  `log_space_*` keys): no composition here is outside machine range. That is
  what makes this state the evidence that the two ADR-0024 entry points are
  separate.

### (vi) Verdict continuity across the whole pressure range

- 25 pressures from 0.3 to 12 MPa at 453 K, 5 wt%, `k_ij = -0.006`, **no
  `ConvergenceError` anywhere**:

  | range | verdict | stage |
  | --- | --- | --- |
  | 0.300 - 2.250 MPa (5 points) | VLE | `second-order-log` |
  | 2.737 - 8.100 MPa (12 points) | LLE | `second-order` |
  | 8.588 - 9.562 MPa (3 points) | LLE | successive substitution |
  | 10.050 - 12.000 MPa (5 points) | single liquid | - |

- The single-phase switch is Case P-13 (i)'s cloud point, **9.7518 MPa**,
  unchanged. Bisecting the VLE / LLE change of character to 1 kPa puts it
  between **2.5868 and 2.5873 MPa**; both sides return `delta_g_split_rt < 0`
  (-4.1536e-03 and -4.1421e-03) and are post-split stable. What changes there
  is the *density root* the incipient phase sits on, and it is a property of
  the model rather than of the route: running the same two states through the
  Wilson-seeded successive-substitution path instead of the log-space one
  reproduces both answers to the printed digits.
- **Only the VLE band is new.** Instrumented pressure by pressure, the LLE band
  takes exactly the path it took before ADR-0024 (the second-order stage from a
  stability seed, or successive substitution alone above 8.5 MPa), which is why
  the 8 MPa numbers of Case P-13 (ii) and the `pcsaft-polymer-lle` benchmark
  hash are unchanged.

### (vii) The two entry points reach the same answer

- Forcing the Wilson K-seed at 1 MPa sends the same state down the other route:
  one successive-substitution step, a linear second-order stage that spends its
  whole 100-iteration budget, and only then the log-space stage seeded from
  that failed iterate. It converges to
  `x_solvent = 0.9714683526646756` and `ln y_PE = -450.5307940616794`, against
  `0.9714683526646756` and `-450.5307940616726` from the stationary-point seed:
  **agreement to the thirteenth figure of a quantity 450 decades below
  machine-representable**, from two seeds and two code paths.
- Seed-independence of the stage itself: the same split comes back from
  `n = 0.5 z`, `0.9 z` and `0.99 z` - three starting points that know nothing
  about the model - to `rel=1e-12` on both quantities.

### Negative controls and what is *not* claimed

- **Dormancy.** No state that converged before ADR-0024 carries a single
  `log_space_*` key: asserted for Peng-Robinson methane/ethane at 240 K /
  3 MPa, PC-SAFT water/n-hexane at 298.15 K / 1 atm and the Case P-13 polymer
  liquid-liquid split at 8 MPa.
- **The fixture did not move.** `refactor_bit_identity_v3.json` (155 states)
  passes unchanged and was **not** regenerated. All nine ADR-0023 benchmark
  cases report **identical result hashes**.
- ~~**A remaining limitation, pinned rather than worked around.** `Mw = 53000`
  at 0.5 and 1 MPa still raises... a **stability trial set** limitation for
  `m = 1393.9`.~~ **RESOLVED by ADR-0025 (Case P-15).** The diagnosis was half
  right. It was upstream of this slice and it was in `stability_tp`, but not in
  the *trial set*: three of its four trials were walking straight at the melt
  and were held 752.2 short of it by a clamp on `ln W` at 700, the melt's
  stationary point being at `ln W_polymer = 1452.21`. With the normalization
  done in logarithms the melt is the minimizer (`tpd_min` -1452.21 / -1346.57)
  and **this stage, unchanged, converges from it**. The reversal is pinned in
  `tests/test_flash_log_space_stage.py::test_the_longest_chain_below_1_mpa_is_now_in_reach`.
- **AMENDED by ADR-0026 (Case P-16).** The stage this case introduced had one
  more weak point than either this case or ADR-0024 recorded: where its Newton
  direction is not a descent direction for the two-phase Gibbs energy it falls
  back to `-r`, whose *length* is `|r|`, so next to the trivial solution - where
  the Gibbs Hessian is indefinite and the direction is therefore always
  rejected - it crawls instead of converging. Three states of the
  `Mw = 53 000` chain measured it (0.3, 2.8 and 2.9 MPa). The stage now retries
  itself once with a curvature safeguard where it finishes above
  `FlashSettings.tol`; **every number in this case is unchanged**, because the
  retry is a second call made only where the first one had already failed.
- **Nothing here is compared against measurement.** The parameters are the
  Case P-12 fixture, the `k_ij` was fitted elsewhere, and the polymer is
  modelled as monodisperse.
- **Tolerance:** asserted 1e-10 absolute on the melt's solvent mole fraction
  against the 1-D solve (achieved 2.53e-14); 1e-12 on the mass balance
  (achieved 2.77e-18); 1e-08 on the fugacity residual and on the log-space
  residual (achieved 7.39e-13 / 9.09e-13); 1e-08 on the polymer's log fugacity
  (achieved 7.39e-13); `dG < 0` and post-split stable at every state; 1e-08 on
  FeOs's chemical potentials with matched constants (achieved 2.615e-12);
  `rel=1e-12` on seed-independence.
- **Independent route:** the one-dimensional equal-fugacity solve, written
  three times over (unit test, validation test, example); FeOs's chemical
  potentials; and the Wilson-seeded path, which is an independent *solver* path
  to the same state.
- **Suite time:** `pytest -q` **224.17 s for 809 tests at `5949a8e` -> 237.18 s
  (3:57) for 830**, uncontended on the one machine; `pytest -q -m slow`
  **694.93 s (11:34) for 34 tests**, against 691.22 s for 28 at the previous
  slice. 27 tests added, none removed; the six marked `slow` are all
  *repetitions* of a map point the default run already covers (another pressure
  on the same three-state map, another model on the same dormancy check).
- **Test path:** `tests/test_flash_log_space_stage.py` (15 by default, 2
  `slow`), `tests/test_pcsaft_polymer.py` (23 by default, 5 `slow`),
  `tests/validation/test_pcsaft_polymer_vle_vs_feos.py` (5 by default, 4
  `slow`).
- **Script:** `python examples/validation/21_pcsaft_polymer_vle.py` (routes 1,
  2, 3 and 5 need no optional dependency; about 5 s; `--full` adds the 25-point
  pressure scan and the bisected VL/LL boundary - and, since ADR-0026, the same
  scan for the `Mw = 53 000` chain, which takes it to about 1 min 15 s).

---

## Case B-1: The benchmark baseline, and three bit-identical optimizations

- **Source:** none, and that is the point. This is not a thermodynamics case:
  it is the *measurement* case. Nothing here is checked against a published
  number; what is checked is that three performance changes moved the clock and
  did not move a single result. The thermodynamic content of every case in the
  workload is already pinned elsewhere in this ledger (F-4, P-8, P-13, L-1,
  R-2, V-1), and the benchmark's `result_hash` is a cross-check on top of
  those, not a substitute for them.
- **Where:** ADR-0023; harness `src/chemthermo/bench/`; records
  `benchmarks/baseline_3ce68df.json` and `benchmarks/after_2ca41bf.json`;
  `benchmarks/README.md`.
- **Assumptions:** wall times are a property of one machine at one moment. All
  numbers below were measured on one Apple M2 Max (12 cores, macOS 26.6.2,
  arm64), CPython 3.11.6, numpy 2.4.2, five timed repeats per case after one
  untimed warm-up, the two records written **47 seconds apart** - the baseline
  from a detached worktree at `3ce68df`, the after record from `2ca41bf`.
- **Components and units:** the nine-case workload of
  `chemthermo/bench/_cases.py`, SI throughout.

### (i) What the profile said, before anything was changed

Measured at `3ce68df` with `cProfile` and with targeted `timeit`:

| observation | number |
| --- | --- |
| `pr-flash-ternary`, share of time in `PengRobinsonEOS.fugacity_coefficients` | 229 calls, ~55 % |
| one `fugacity_coefficients` call | 34.8 us |
| - of which `numpy.roots` | 16.9 us |
| - of which `_mixture_parameters` | 10.5 us |
| - of which the per-component `ln phi` loop and `exp` | 2.9 us |
| - of which input validation | 0.8 us |
| `pcsaft-lle-water-hexane`, share of time in `solve_density_roots` | 291 solves, ~95 % |
| one water / n-hexane density solve at 298.15 K, 1 atm | 5.3 ms |
| - of which the 1599-point scan | 1.28 ms (24 %) |
| - of which ~19 single-point `pressure_and_slope` calls | 3.6 ms (68 %) |
| one `derivatives(1 point, second=True)` | 189 us, of which 100 us association |

And the number that resized the whole slice: over that flash, **291 density
solves were made at 238 distinct `(T, P, x)`**. Broken down by call site: 154
from `_select_surface` (a pinned trial iterating on one root, ADR-0021
decision 3), 88 from the split's `_PhaseRoot`, 24 from
`_select_density_root_surface`, 8 from `_select_min_gibbs`, 15 from
`phase_identity`. The minimum-Gibbs double-solve the slice was scoped around is
the last two rows: **16 avoidable solves of 291, about 5 %**. ADR-0021 had
already removed the bulk of it. The corresponding Peng-Robinson figure is 229
`fugacity_coefficients` calls of which 167 are pinned single-branch and 38 come
from both-branch selectors - 19 avoidable, about 8 %.

**What is therefore *not* claimed:** root reuse is not worth 2x on the PC-SAFT
liquid-liquid flash. It is worth about 5 % there. The larger part of the
measured 1.32x comes from optimizations B and C.

### (ii) The three optimizations, and the evidence each is bit-identical

**A - one root solve serves every branch** (`ln_fugacity_branches`). Identical
by construction: the capability returns `ln phi` *before* the exponential, so
`exp` of it is the same double `fugacity_coefficients` returns, on the same
root, and `eos_branch_terms_all` re-applies the ADR-0022 guard unchanged.
Checked with `==`, not a tolerance, over: the 16-state Case F-4 subset (32
branch comparisons), six water / n-hexane liquid-liquid states, seven
Peng-Robinson grid states at two `k_ij`, the single-root case (both labels must
carry the *same* array, or ADR-0021 decision 4's degeneracy test stops
recognising a one-root region), and the `Mw = 53000` polymer state where
`exp(ln phi)` underflows and `phi is None` must survive the new route. Zero
differences.

**B - the Peng-Robinson companion matrix built directly.** `numpy.roots` on a
monic cubic with a non-zero constant term *is* `eigvals(companion)`; the
bookkeeping around it is skipped. **200 000 random `(A, B)` pairs** over
`A in [1e-6, 1e2]`, `B in [1e-6, 1]`: `np.array_equal` on the eigenvalue array,
**0 mismatches**. 14.7 us -> 8.6 us. The `c3 == 0.0` case, the only one
`numpy.roots` deflates, is delegated back to it.

**C - the density solver stops recomputing what it has.** Three changes, each a
rearrangement of *when*: `a` (which the root solver never reads) is not
computed when only the derivatives are wanted; the bracket's left-edge residual
is carried in from the scan instead of re-evaluated; the refinement's own
`(P, dP/drho)` at the root it returns is handed to the mechanical-stability
filter instead of being recomputed. The second needed a check and got one - a
single-point isotherm evaluation returns the same double the 1599-point grid
evaluation put in that slot, verified over the whole grid for water / n-hexane
at 298.15 K, methane / n-hexane at 300 K and CO2 / n-decane at 240 K
(`np.array_equal`, 0 differences).

### (iii) The measured ratios

`python -m chemthermo.bench --compare benchmarks/baseline_3ce68df.json benchmarks/after_2ca41bf.json`,
exit status 0:

| case | before / s | after / s | ratio |
| --- | ---: | ---: | ---: |
| `pr-flash-ternary` | 0.015261 | 0.013040 | 1.17x |
| `pr-stability-ternary` | 0.004649 | 0.003881 | 1.20x |
| `pr-flash-grid-24` | 0.140889 | 0.117589 | 1.20x |
| `pcsaft-vle-methane-hexane` | 0.240769 | 0.188477 | 1.28x |
| `pcsaft-lle-water-hexane` | 1.626994 | 1.236466 | 1.32x |
| `pcsaft-polymer-lle` | 0.962138 | 0.814358 | 1.18x |
| `nrtl-lle-tessier-p1` *(control)* | 0.043876 | 0.043711 | 1.00x |
| `modified-raoult-vle` *(control)* | 0.020264 | 0.019569 | 1.04x |
| `vlle-364k` *(control)* | 0.099957 | 0.097971 | 1.02x |

**All nine result hashes identical.** The three activity-model cases run none
of the changed code and are the drift control; read the equation-of-state
ratios against their 1.00-1.04x, not against an exact 1.00x. A repeat of the
whole comparison in three alternating rounds at three repeats gave 1.14x, 1.17x,
1.18x, 1.27x, 1.32x, 1.17x on the six equation-of-state cases and 0.98-0.99x on
the three controls, so the ratios above are reproducible to about a point.

Peak allocation (`tracemalloc`, one separate instrumented run) fell where the
density solver stopped allocating an unused `a` and an unused second derivative:
`pcsaft-lle-water-hexane` 0.86 -> 0.78 MiB, `pcsaft-vle-methane-hexane`
0.52 -> 0.49 MiB, `pr-flash-grid-24` 0.04 -> 0.03 MiB.

### (iv) Tried, measured, and rejected

Recorded because each is the obvious next idea and each is wrong for a reason
worth keeping:

- **A Cardano closed-form cubic.** Much faster; does not reproduce these
  doubles. Bit-identity is the gate.
- **Vectorizing the Peng-Robinson per-component `ln phi` loop.** 3.0 us against
  the loop's 1.4 us at n = 3: three numpy operations on a length-3 array cost
  more than the loop.
- **Replacing the association module's `np.einsum` with explicit
  multiply-and-sum.** Not bit-identical for >= 3 sites (checked at 1, 2, 3, 4, 6
  and 8 sites, against both a `.sum(-1)` and a `matmul` form; equal at 1 and 2,
  different from 3 up), *and* slower at the sizes that occur (1.34 us against
  1.00 us).
- **Refining every bracket of one density solve in a single vectorized pass.**
  `solve_site_fractions` tests convergence on `max |residual|` over the whole
  batch, so batching two brackets would give one of them a different number of
  Newton steps and move its last bits.
- **Caching the identity matrix the site-fraction Newton allocates.**
  Implemented, measured, reverted: no gain outside the noise.
- **A last-solve cache on the model instance.** Would catch the 37 further
  duplicate solves the branch capability cannot see (`phase_identity`, the
  post-split re-evaluation). Rejected on design grounds, not measurement: both
  model classes are frozen dataclasses and a stale cache key would be a wrong
  answer rather than a slow one. *(Amended by slice `perf-eos-solve-cache`,
  ADR-0030: those 37 solves are now removed by a memo that lives in the
  **call** rather than on the model, which is not the construction rejected
  here - nothing is written to either dataclass and no key outlives the call
  that made it. See Case B-2.)*

- **Tolerance:** exact. `==` on every comparison above; no tolerance is used
  anywhere in this case, which is what distinguishes a performance change from
  a numerical one.
- **Independent route:** the 155-state pinned fixture
  `tests/fixtures/flash/refactor_bit_identity_v3.json`, **unchanged and not
  regenerated**, is the independent check that the flash answers did not move;
  the branch-by-branch `==` comparisons are the check that the model calls
  underneath them did not either.
- **Negative control:** the three activity-model workload cases, whose code
  none of the three optimizations touches (1.00x / 1.04x / 1.02x); a model
  implementing only the pre-ADR-0023 interface, which must go on being served
  by the per-branch route (asserted by call count);
  `tests/test_bench_harness.py`, which asserts `--compare` **fails** when a
  result hash differs, so the acceptance rule cannot pass vacuously.
- **Suite time:** `pytest -q` **263.80 s for 803 tests at `3ce68df` ->
  224.17 s for 809 at `2ca41bf`** (-15 % while gaining six tests), both runs
  uncontended on the one machine; `pytest -q -m slow` **691.22 s (11:31) for
  28 tests**, against 856.2 s (14:16) for the same 28 at the previous slice (a
  different session, so indicative rather than back to back). 13 tests added
  across the slice (7 harness, 6 branch reuse), none removed, **none marked
  `slow`**.
- **Test path:** `tests/test_eos_branch_reuse.py` (6),
  `tests/test_bench_harness.py` (7), plus the whole existing suite unchanged.
- **Script:** `python -m chemthermo.bench --out record.json` and
  `python -m chemthermo.bench --compare before.json after.json`.

---

## Case B-2: The call-local solve memo, and the structured Hessian that was measured and declined

- **Source:** none, and as with Case B-1 that is the point. This is the
  *measurement* case for slice `perf-eos-solve-cache`. Nothing here is checked
  against a published number; what is checked is that a performance change
  moved the number of solves and did not move a single result. The
  thermodynamic content of every state used below is already pinned elsewhere
  in this ledger (F-4, P-8, P-9, P-10, P-13).
- **Where:** ADR-0030; `src/chemthermo/_eos_memo.py`; the two solve points in
  `src/chemthermo/eos/pcsaft.py` (`_density_roots`) and
  `src/chemthermo/models/peng_robinson.py` (`_compressibility_roots`); records
  `benchmarks/baseline_d8814de.json` and `benchmarks/after_07bf0b3.json`;
  `benchmarks/README.md`.
- **Assumptions:** the machine was **not quiet**. An unrelated process held
  about eight of this machine's twelve cores for the whole session, load
  average 9-11 throughout, and every wall time below was measured under that.
  Machine: Apple M2 Max (12 cores, macOS 26.6.2, arm64), CPython 3.11.6, numpy
  2.4.2. The numbers that carry the claim are therefore the **solve counts**,
  which are integers and a property of the workload; the wall times are
  reported for completeness and are explicitly not the evidence.
- **Components and units:** the nine-case ADR-0023 workload, the Case F-4
  default subset, the ADR-0019 water / n-hexane states and the Case P-10
  water / ethanol / n-hexane tie-triangle states. SI throughout.

### (i) What was duplicated, and by whom

At `d8814de`, counting `PCSAFTEOS._density_roots` calls whose
`(names, T, P, x)` are equal bit for bit over one `invoke` of each case:

| case | solves | distinct | repeats |
| --- | ---: | ---: | ---: |
| `pcsaft-vle-methane-hexane` | 137 | 102 | 35 (25.5 %) |
| `pcsaft-lle-water-hexane` | 275 | 238 | 37 (13.5 %) |
| `pcsaft-polymer-lle` | 601 | 512 | 89 (14.8 %) |

Attributing each of the 37 repeats in the liquid-liquid case to the pair of
call sites that made it (captured from the stack at each solve):

| first call site | second call site | count |
| --- | --- | ---: |
| `_evaluator.ln_fugacity_terms` -> `fugacity_coefficients` | `_reported_terms` -> `identity_label` -> `phase_identity` | 12 |
| `_evaluator.ln_fugacity_terms` -> `fugacity_coefficients` | `_evaluator.branch_terms` -> `ln_fugacity_branches` | 12 |
| `_evaluator.ln_fugacity_terms` -> `fugacity_coefficients` | itself, at the same state | 5 |
| `_split._branch_terms` -> `fugacity_coefficients` | `_detect._name_two_phase_result` -> `phase_identity` | 2 |
| `_split._branch_terms` -> `fugacity_coefficients` | `_evaluator.branch_terms` -> `ln_fugacity_branches` | 2 |
| `_split._branch_terms` -> `fugacity_coefficients` | `stability_tp` -> `identity_label` -> `phase_identity` | 2 |
| `_evaluator.branch_terms` -> `ln_fugacity_branches` | `stability_tp` -> `identity_label` -> `phase_identity` | 1 |
| `_evaluator.branch_terms` -> `ln_fugacity_branches` | itself, at the same state | 1 |

**This is why the memo is where it is.** Six of the 37 are one entry point
called twice; the other 31 are two *different* model methods reaching one root
solve. A memo at the `chemthermo.flash._common` level - the obvious place, and
the one this slice was scoped to first - would have caught the six and none of
the rest, while touching a dozen internal signatures instead of two.

### (ii) Solves removed - the number that is not a wall time

Counting executed `solve_density_roots` calls and executed cubic eigenproblems
over the same prepared payload, with the models' memo lookup returning `None`
and with it active. Integers, and a property of the workload rather than of the
machine:

| case | PC-SAFT density solves off -> on | PR cubic solves off -> on |
| --- | ---: | ---: |
| `pr-flash-ternary` | - | 230 -> 187 (**-18.7 %**) |
| `pr-stability-ternary` | - | 84 -> 73 (**-13.1 %**) |
| `pr-flash-grid-24` | - | 2203 -> 1679 (**-23.8 %**) |
| `nrtl-lle-tessier-p1` (control) | 0 -> 0 | 0 -> 0 |
| `modified-raoult-vle` (control) | 0 -> 0 | 0 -> 0 |
| `vlle-364k` (control) | 0 -> 0 | 0 -> 0 |
| `pcsaft-vle-methane-hexane` | 137 -> 102 (**-25.5 %**) | - |
| `pcsaft-lle-water-hexane` | 275 -> 238 (**-13.5 %**) | - |
| `pcsaft-polymer-lle` | 601 -> 512 (**-14.8 %**) | - |

`pr-flash-grid-24` saves more than the single Peng-Robinson flash does because
it is 24 flashes and each one saves its own repeats; the memo is per call, so
nothing is shared *between* the 24.

### (iii) Bit-identity, and where it comes from

Not from an arithmetic argument. A hit returns **the identical object** the
first solve produced - a `DensityRoots` namedtuple of floats, or a tuple of
compressibility roots - so every double downstream of it is the same double.
What has to be argued is only that no key can name the wrong state:

- the key is the solve routine's own arguments, compared for exact equality;
- `composition_key` gives `-0.0` and `0.0` **different** keys. They compare
  equal and hash equal in Python, so a plain tuple key would merge them. That
  is the one pair of distinct bit patterns that could alias, and it is closed
  by construction rather than by an argument that the models would not care;
- `nan` never compares equal to itself, so a `nan` argument misses;
- `id(model)` is in the PC-SAFT key, and the memo holds a strong reference to
  the model in the same entry, so the id cannot be recycled onto another object
  while the entry lives;
- the 4096-entry FIFO bound can only cause a re-solve, which produces the same
  doubles again. `test_a_memo_bound_of_one_gives_the_same_answer_as_the_default_bound`
  runs a whole associating flash through a **one-entry** memo and compares
  every field with `==`.

Checked, all with `==` and never a tolerance:

| what | states | result |
| --- | ---: | --- |
| Case F-4 default subset, `flash_tp`, memo on vs memo off | 16 | identical, field for field |
| ADR-0019 water / n-hexane, `flash_tp` and `stability_tp` | 4 (x2) | identical |
| ADR-0017 Peng-Robinson grid slice, `flash_tp` and `stability_tp` | 7 (x2) | identical |
| ADR-0022 polymer split (where `exp(ln phi)` underflows) | 1 | identical |
| one-entry memo against the default bound | 1 | identical |
| two `kij` inside one deliberately shared scope | 2 | each equals what it gives alone |
| two models flashed concurrently in two threads | 2 | each equals what it gives alone |
| `refactor_bit_identity_v3.json` (155 states) | 155 | unchanged, not regenerated |
| the nine benchmark `result_hash` values | 9 | identical |
| the 2505-state robustness map, every observable field | 2505 | **0 differences** |

### (iv) The structured multiphase Hessian: measured, and **not shipped**

The assembly is exact as algebra: with the reference phase holding
`N_0 = z - sum_p n_p` and `a_{j,k} = ln x_{j,k} + f_j(x_j)_k`, phase `p+1`
depends only on `n_p`, so `H[(p,k),(q,l)] = delta_pq A^{(p+1)}[k,l] + A^{(0)}[k,l]`
with `A^{(j)}` the per-phase Jacobian. That is `N_phases (NC x NC)` finite
differences instead of `((N_phases - 1) NC)^2`.

Built and compared entry by entry against the Hessian
`_multiphase_second_order` actually forms, at **every** build the stage reaches
over two sets of states: the 20 Case P-9 PC-SAFT water / n-hexane states below
`T3` (`z_water` in `{0.05, 0.3, 0.5, 0.7, 0.95}` x `T3 - {1.0, 0.5, 0.1, 0.01}`
K, 1 atm), which give **8 builds**, all two-phase, **0 of 8 bit-identical**,
largest absolute difference 3.199e-04 and largest relative 7.381e-05; and six
Case P-10 water / ethanol / n-hexane feeds at 1 atm (`(0.1,0.1,0.8)`,
`(0.2,0.2,0.6)`, the centroid and `(0.6,0.2,0.2)` at 333 K, `(0.2,0.6,0.2)` and
`(0.1,0.1,0.8)` at 335 K), which give the **4 builds** below, three of them
three-phase:

| phases | Hessian size | bit-identical | max abs difference | max rel difference | max `|H|` |
| ---: | ---: | --- | ---: | ---: | ---: |
| 2 | 3 | no | 1.557e-07 | 6.208e-08 | 4.633e+01 |
| 3 | 6 | no | 4.079e-04 | 9.675e-05 | 7.709e+05 |
| 3 | 6 | no | 1.018e-04 | 2.522e-05 | 4.416e+04 |
| 3 | 6 | no | 1.047e-04 | 2.536e-04 | 1.051e+04 |

**0 of 12 bit-identical over both sets**, and the reason is structural rather
than incidental.
For a row `(p,k)` and a column `(q,l)` with `p != q`, the shipped code forms
`fl(fl(u - v_-) - fl(u - v_+))` where `u = a_{p+1,k}` is the *same double* in
both gradient evaluations; the structured form computes `fl(v_+ - v_-)`. Those
differ whenever `u` is large next to `v_+ - v_-`, which is the ordinary case
here. Two further differences compound it: the diagonal blocks add their two
contributions in a different order, and the shipped code's reference-phase
perturbation is `fl(z_l - fl(sum_q n_{q,l}))` rather than an exact
`N_{0,l} - h`. A relative difference of 1e-04 in the Hessian moves the Newton
direction, the line search accepts on an Armijo decrease *or* a residual
decrease, so a moved direction moves the accepted iterate. Bit-identity is the
gate; it is not shipped.

**And what it would have bought is smaller than it was**, because the memo
already takes most of it, bit-identically. Counting the model evaluations each
Hessian build makes (distinct `(phase, composition)` pairs is what the memo
reduces the build to):

| phases | Hessian size | evaluations, no memo | with the memo | structured |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 3 | 12 | 12 | 12 |
| 3 | 6 | 36 | 20 | 18 |
| 3 | 6 | 36 | 20 | 18 |
| 3 | 6 | 36 | 21 | 18 |

15 or 16 of the 18 evaluations the structured form would remove are already
removed. The remaining 2 or 3 per build are what a re-audit slice would be
buying, and it would be paying a regenerated bit-identity fixture for them.

### (v) The benchmark pair, and why it is not the evidence

Nine timed repeats each, back to back, `benchmarks/baseline_d8814de.json`
against `benchmarks/after_07bf0b3.json`, every `result_hash` identical:

| case | before / s | after / s | ratio |
| --- | ---: | ---: | ---: |
| `pr-flash-ternary` | 0.016810 | 0.013450 | 1.25x |
| `pr-stability-ternary` | 0.004402 | 0.004119 | 1.07x |
| `pr-flash-grid-24` | 0.158998 | 0.120641 | 1.32x |
| `nrtl-lle-tessier-p1` (control) | 0.057287 | 0.059526 | 0.96x |
| `modified-raoult-vle` (control) | 0.024677 | 0.024191 | 1.02x |
| `vlle-364k` (control) | 0.137144 | 0.133123 | 1.03x |
| `pcsaft-vle-methane-hexane` | 0.216429 | 0.178299 | 1.21x |
| `pcsaft-lle-water-hexane` | 1.595234 | 1.324830 | 1.20x |
| `pcsaft-polymer-lle` | 0.994039 | 0.851818 | 1.17x |

A first attempt at five repeats gave controls of 0.77x / 1.25x / 0.92x and was
discarded as unreadable; nine repeats brought them to 0.96x / 1.02x / 1.03x,
still wider than the 1.00x / 1.04x / 1.02x the same three read for Case B-1 on
a quiet machine.

Because of that, the same workload was also run **in one process with the two
arms interleaved repeat by repeat** - memo active against memo lookup returning
`None`, milliseconds apart, so both arms take the same contention. 11 repeats
each, with the observed per-repeat ratio range:

| case | off / ms | on / ms | ratio | spread |
| --- | ---: | ---: | ---: | --- |
| `pr-flash-ternary` | 13.77 | 13.41 | 1.027x | 1.01x .. 1.03x |
| `pr-stability-ternary` | 4.17 | 4.10 | 1.018x | 1.01x .. 1.02x |
| `pr-flash-grid-24` | 125.67 | 126.70 | 0.992x | 0.94x .. 1.14x |
| `nrtl-lle-tessier-p1` (control) | 52.64 | 54.03 | 0.974x | 0.82x .. 1.33x |
| `modified-raoult-vle` (control) | 22.88 | 24.95 | 0.917x | 0.74x .. 1.13x |
| `vlle-364k` (control) | 133.49 | 132.40 | 1.008x | 0.63x .. 1.12x |
| `pcsaft-vle-methane-hexane` | 190.01 | 146.00 | 1.301x | 1.23x .. 1.31x |
| `pcsaft-lle-water-hexane` | 1363.89 | 1233.82 | 1.105x | 1.01x .. 1.25x |
| `pcsaft-polymer-lle` | 924.09 | 801.81 | 1.153x | 1.05x .. 1.22x |

**What is claimed:** PC-SAFT **1.10x - 1.30x** (the two measurements agree
within their spreads). **What is not claimed:** the Peng-Robinson ratios in the
committed pair. 1.25x / 1.07x / 1.32x there against 1.03x / 1.02x / 0.99x
interleaved, and 1.10x / 1.02x / 1.04x on a second 31-repeat interleaved run.
The arithmetic settles it against the committed pair: Case B-1 measured the
cubic solve at 8.6 us of a 34.8 us `fugacity_coefficients` call, so removing
13-24 % of the solves cannot be worth 30 % of the flash. The Peng-Robinson half
of this change is **1.00x - 1.05x** and is kept for the work it removes, not
for a clock reading. That is a doubt this case records rather than resolves,
and the thing that would resolve it is one uncontended pair.

- **Tolerance:** exact. `==` on every comparison in (iii); no tolerance is used
  anywhere in it. The Hessian comparison in (iv) is the one place a difference
  is *reported* rather than asserted to be zero, and reporting it is the
  finding.
- **Independent route:** the 155-state pinned fixture
  `tests/fixtures/flash/refactor_bit_identity_v3.json`, unchanged and not
  regenerated; the nine benchmark `result_hash` values; and the full 2505-state
  robustness map re-run at `f852726` and diffed field by field against
  `robustness_64831bd.json` with wall time excluded - **0 differences over all
  2505 states and every aggregate field**, and 2500 converged / 5 refused / 0
  invariant violations either way. The sweep read 2017.3 s (33:37) against
  2202.5 s (36:43) at `64831bd`, two loaded machines in different sessions, so
  that pair is a direction and not a ratio. `robustness_64831bd.json` is pruned
  per the `benchmarks/README.md` policy and its `.md` summary kept.
- **Negative control:** the three activity-model workload cases, which make
  zero memoized solves either way because no activity model is memoized - that
  is a deliberate decision (ADR-0030 decision 6) precisely so they stay the
  drift control; and
  `test_two_models_never_read_each_others_entries`, which would fail if the key
  did not carry model identity, together with the assertion inside it that the
  two models actually disagree, so the test cannot pass vacuously.

- **Suite time:** `pytest -q` **901 passed, 98 deselected, 300.74 s
  (5:00) at `d8814de` -> 916 passed, 98 deselected, 286.03 s (4:46)**, the two
  runs back to back on the one machine and both under the same contention.
  The interesting number is the third: the **same 901 tests** (the default run
  with `--ignore=tests/test_eos_memo.py`) read **235.36 s (3:55)**, reproduced
  at 235.46 s on a repeat - **1.28x, -21.7 %**, and *under* the ~240 s ceiling
  brain.md section 10 tracks, on a machine that was not even quiet.
  `tests/test_examples.py` alone reads **101.64 s -> 85.95 s** (1.18x, 43 tests
  either way, the file unchanged).
  `pytest -q tests/test_robustness_map.py -m ""` - the slow full-sweep gate,
  which re-runs all 2505 states in process and checks them against the
  committed record - reads **52 passed in 2028.71 s (33:48)**.
  **The default run as shipped is still over the ceiling, at 286 s**, and the
  reason is this slice's own evidence: `tests/test_eos_memo.py` costs ~50 s
  because it flashes every state **twice**, once with the memo and once
  without, and that is the whole proof. **Nothing was marked `slow`** to hide
  it - the `slow` policy forbids marking the only test of a capability, and a
  bit-identity sweep run at half its states is a weaker proof, not a cheaper
  one. So: the library change carries the pre-existing suite under the ceiling,
  and the slice does not carry the whole default run under it. 15 tests added,
  none removed.
- **Test path:** `tests/test_eos_memo.py` (15), plus the whole existing suite
  unchanged.
- **Script:** `python -m chemthermo.bench --out record.json`,
  `python -m chemthermo.bench --compare before.json after.json`, and
  `python -m chemthermo.bench robustness --out record.json --summary-out record.md`.

---

## Case P-15: The tangent-plane melt a clamp on `ln W` had made unreachable

- **Source:** a **one-dimensional** equal-fugacity solve written independently
  in the example (the vapour taken as *exactly* pure solvent, one unknown
  carried as `ln x_solvent`, finite-difference Newton started a thousandth away
  from `flash_tp`'s answer); **FeOs 0.10.1** chemical potentials at
  chemthermo's converged phases and densities; and, for the stationary point
  itself, Michelsen's own equations (5) and (7) re-derived from `PCSAFTEOS` in
  the test and again in the example.
- **Location:** `chemthermo/stability/tp.py` (`_normalize`, the four removed
  `np.clip(..., -700, 700)`), `chemthermo/stability/results.py`,
  `chemthermo/flash/_detect.py`, `chemthermo/flash/_log_space.py`; slice
  `stability-log-space-sums` (ADR-0025). **No model equation, no split solver,
  no trial set, no initial estimate, no convergence criterion and no tolerance
  changed.**
- **Parameters and provenance:** exactly Case P-13's and Case P-14's -
  polyethylene from `tests/fixtures/pcsaft/martini2009_polymers.json`
  (`m/M = 0.026300 mol/g`, so `m = 1393.9` at `Mw = 53 000`), n-pentane from
  Gross & Sadowski (2001) Table 1, `k_ij = -0.006` except for the FeOs routes,
  which run at `k_ij = 0` on both sides because feos 0.10.1 cannot be given
  one from Python. Nothing here is compared against measurement.
- **Assumptions:** monodisperse polymer; 453.0 K throughout; feed stated as a
  polymer **mass** fraction (5 wt%, so
  `z = (7.163935601491655e-05, 0.9999283606439852)`).

### (i) The defect, measured before it was repaired

`stability_tp`'s deterministic trial set was never the problem. Three of its
four trials - `wilson-liquid` and both `pure-<name>` estimates - were walking
at the melt and could not arrive, because the iteration clamped
`ln W` to `[-700, 700]` and the melt's stationary point is at
`ln W_polymer = 1452.21`.

| P / MPa | `tpd_min` before | minimizing trial before | the three liquid trials before |
| --- | --- | --- | --- |
| 0.5 | -1.0518e-04 | `wilson-vapor`, 8 iterations | `second_order_no_progress` after **51** iterations at `w = (1, 3.67e-303)`, residual **752.2** |
| 1.0 | -3.3766e-04 | `wilson-vapor`, 4 iterations | `second_order_no_progress` after **51** iterations at `w = (1, 6.77e-303)`, residual **646.6** |

`752.2 = 1452.21 - 700` exactly: the residual the clamped trial reports *is*
the distance the clamp was holding it back. `flash_tp` then raised
`ConvergenceError` at both states (Case P-14, "a remaining limitation").

### (ii) The stationary point, after

| P / MPa | `tpd_min` | `ln W` = (polymer, solvent) | iterations | stationarity residual |
| --- | --- | --- | --- | --- |
| 0.5 | **-1452.20680935942** | `(1452.20680935942, 3.616782632957839)` | 3 SSI | 4.26e-14 |
| 1.0 | **-1346.573440835703** | `(1346.5734, 4.2287)` | 3 SSI | 8.53e-14 |

- **Expected outcome:** equation (5), `ln W_i + ln phi_i(w) - d_i = 0` on the
  trial's own surface, and equation (7), `tpd = -ln(sum_i W_i)`.
- **Achieved (re-derived from `PCSAFTEOS` in the example, not read from the
  solver):** equation (5) to **4.3e-14** / 8.5e-14; equation (7) to **exactly
  0.0** at both states.
- `sum_W` itself is `inf` and `tm_at_stationary_point` is `-inf`; `ln_sum_W` is
  1452.2068 / 1346.5734. The normalized `w` is `(1.0, 0.0)` - the solvent's
  share of the melt is `exp(-1448.6)` - which is why `ln W` is reported
  separately and is what the split stage is seeded from.
- The minimizer is the melt (`phase_branch = "liquid"`, `feed_branch = "vapor"`)
  and 3 of 4 trials engage the log-space route at both states.

### (iii) The split it seeds

| P / MPa | melt `x` (polymer, solvent) | melt, wt% solvent | `beta_vapor` | `ln y_polymer` | log-space iterations |
| --- | --- | --- | --- | --- | --- |
| 0.5 | `(0.021187890289570684, 0.9788121097104293)` | 5.916 | 0.996618853739762 | **-1507.9818790382** | 16 |
| 1.0 | `(0.009098397632060337, 0.9909016023679396)` | 12.911 | 0.9921261568342015 | **-1452.6599876925** | 11 |

- `k_seed = "stability-log"`, `converged_stage = "second-order-log"`,
  `phase_label_method = "compressibility"` at both.
- Residuals: `fugacity_residual` 1.30e-13 / 3.27e-14, `log_space_residual`
  1.59e-12 / 4.55e-13, `mass_balance_residual` 1.11e-16 at both,
  `delta_g_split_rt` **-1.0710e-01** / -1.0291e-01, both phases post-split
  **stable**.
- Compressibility identity (ADR-0017), from the public `pressure_Pa` routine:
  vapour `kappa` 1.0680 / 1.1603, melt `kappa` 0.0010 / 0.0023, against the
  0.5 threshold.

### (iv) The independent one-dimensional solve

- **Expected outcome:** the vapour is pure solvent to `exp(-1508)`, so the melt
  must satisfy `ln phi_s^V(pure solvent vapour) = ln x_s + ln phi_s^L(x)`
  exactly, in one unknown. Nothing of the flash or of the stability test is
  used - only `density_roots` and `ln_fugacity_coefficients` at `(T, rho, x)`.
- **Achieved:** `|x_solvent(flash) - x_solvent(1-D)|` = **1.20e-14** (0.5 MPa)
  and **3.55e-15** (1.0 MPa), against an asserted 1e-10.
- The polymer's *own* condition, the one neither the linear split nor the
  clamped stability iteration could write down: `ln f_PE` = **-1591.0434549748**
  (0.5 MPa) and **-1626.3242034592** (1.0 MPa), melt against vapour, closing to
  **1.59e-12** and **1.14e-12**.

### (v) FeOs at chemthermo's phases

- **Expected outcome:** `mu_i^I / RT = mu_i^II / RT` at chemthermo's two
  compositions and two densities, evaluated by an independent implementation.
  The polymer's ideal term is `-inf` on both sides (`y_polymer` is an exact
  `0.0`), so what is compared is the solvent's potential - the equation that
  sets the answer.
- **Achieved:** **1.105e-12** with FeOs's own universal constants substituted
  in, **8.550e-08** as shipped (chemthermo packages the ten figures of Gross &
  Sadowski 2001 as printed; FeOs hard-codes fourteen). Asserted 1e-08 on the
  matched-constants figure.
- FeOs's own `tp_flash` is **not** a reference here, for the reasons Case P-13
  (iii) records; only its potentials are used.

### (vi) What the repair is worth beyond the two pinned states

Over the whole region where the feed is a vapour, 0.4 to 2.0 MPa:

| P / MPa | before | after |
| --- | --- | --- |
| 0.40 | `ConvergenceError` | VLE, melt `x_solvent` 0.9728932929, `tpd_min` -1472.0588 |
| 0.50 | `ConvergenceError` | VLE, 0.9788121097, -1452.2068 |
| 0.75 | `ConvergenceError` | VLE, 0.9868285386, -1400.8038 |
| 1.00 | `ConvergenceError` | VLE, 0.9909016024, -1346.5734 |
| 1.50 | VLE, 0.9950537550 | VLE, 0.9950537550, -1226.8213 |
| 2.00 | VLE, 0.9972263035 | VLE, 0.9972263035, -1081.7506 |

The whole picture, from a 0.3-3.6 MPa sweep at 0.1 MPa steps run against this
commit's parent and against this one, on **both** molar masses (68 states):

| outcome | states |
| --- | --- |
| identical, bit for bit | **48**, including every one of the 34 `Mw = 16 400` states |
| `ConvergenceError` -> a verified split | **13**: 0.4-1.4 MPa (vapour-liquid), 2.6 and 2.7 MPa (liquid-liquid) |
| converged before, moved in the last bits | **7**: 1.5-2.1 MPa |
| converged before, raises now | **0** |

The two liquid-liquid repairs at 2.6 and 2.7 MPa come from the `ln W` hand-over
rather than from the sum: their *verdict* is unchanged, but the minimizing
trial reaches `ln W_polymer = -1187`, so the seed is now built from `ln W`
instead of from a `w` that had rounded to zero.

At 2.2-2.5 MPa the feed's own lowest-Gibbs branch is the liquid root, there is
no melt stationary point distinct from it, and both the verdict and the split
are **identical** to before (2.2 MPa: `tpd_min = -8.4453e-02`, melt
`x_polymer = 0.0021352851683857058`, `beta = 0.9664497477547337`; 2.5 MPa:
`-1.7459e-02`, `0.0012513825804704919`, `0.9427518353436072`). At 0.3 MPa and
2.8-3.2 MPa `flash_tp` raised before and raised after this slice too.
**AMENDED by ADR-0026 (Case P-16):** those six states now converge, and the
"gap in the split stage's budget at a shallow near-critical verdict" recorded
here was two different things - the log-space stage falling back to a step
whose length is `|r|` next to the trivial solution (0.3, 2.8, 2.9 MPa), and a
Rachford-Rice denominator cancelling to an exact zero for a `K` of `1e-18`
(3.0-3.2 MPa). Neither is a budget. Every number in this case is unchanged by
that repair, which engages only where `flash_tp` was about to raise.

### Negative controls and what is *not* claimed

- **Dormancy, measured one evaluation at a time.** The gate runs the
  pre-ADR-0025 expressions wherever every `ln W_i` is inside `[-700, 700]`,
  which is wherever the old `np.clip` returned its argument. Engagement count:
  **0 of 624 trials over the 144-state Peng-Robinson stability grid** (144
  states), **0 of 12** over the three `Mw = 16 400` states of Case P-14 (worst
  `|ln W|` there = **444.2**, comfortably inside the window), and 0 on
  Peng-Robinson methane/ethane at 240 K / 3 MPa, PC-SAFT water/n-hexane at
  298.15 K / 1 atm, and an NRTL n-butanol/water liquid-liquid test.
- **The fixture did not move.** `refactor_bit_identity_v3.json` (155 states)
  passes **unchanged** and was **not** regenerated. All nine ADR-0023 benchmark
  cases report **identical result hashes**.
- **Seven states did move, and they are exactly the ones the gate permits.**
  `Mw = 53 000` between 1.5 and 2.1 MPa had a clamped trial *and* still
  produced an answer; they are now seeded from the melt instead of from a
  shallow vapour-side point and converge on the same answer from a different
  direction. Worst over the seven: `x_polymer` **3.48e-13** relative (2.0 MPa,
  2.7736965292611865e-03 -> 2.7736965292621523e-03; at 1.5 MPa
  0.004946244958517399 -> 0.004946244958516641, 1.53e-13) and `beta_vapor`
  **9.23e-15**. The 2.0 MPa melt composition is pinned at `rel=1e-9` in
  `tests/test_flash_log_space_stage.py` and still passes. **No state that
  converged before raises now.**
- **`ln W` versus the reconstruction it replaces.** The split seed used to be
  built as `ln W = ln w - tpd`. Over the **47 unstable verdicts** of the
  144-state Peng-Robinson grid the two agree to **5.4e-15**. Over *stable*
  verdicts they differ by up to **4.44**, which is not new and not a defect -
  `tpd_min` is measured to the lower envelope of the phase candidates while
  `ln sum_W` is equation (7) on the surface the trial iterated on (ADR-0021) -
  and a stable verdict seeds nothing. This is why `ln W` is handed over only
  where the reconstruction has stopped working.
- **The Newton stage at `|ln W| ~ 1e3` is untested.** Every state this slice
  repairs converges in 3 successive substitutions, so the second-order stage is
  never entered there. Its multiplicative `1e-6` step in `u` is scale-free by
  the same argument ADR-0024 made for its own Hessian, but that is an argument,
  not a measurement.
- **`stable` is still not a global proof.** One family of stationary points
  became reachable; nothing is claimed about stationary points no initial
  estimate approaches.
- **Nothing here is compared against measurement.** The parameters are the Case
  P-12 fixture, the `k_ij` was fitted elsewhere, and the polymer is modelled as
  monodisperse.
- **Tolerance:** asserted 1e-08 on equations (5) and (7) at the stationary
  point (achieved 8.5e-14 and 0.0); 1e-10 absolute on the melt's solvent mole
  fraction against the 1-D solve (achieved 1.20e-14); 1e-12 on the mass balance
  (achieved 1.11e-16); 1e-08 on the fugacity residual, the log-space residual
  and the polymer's log fugacity (achieved 1.30e-13, 1.59e-12, 1.59e-12);
  `dG < 0` and post-split stable at both states; 1e-08 on FeOs's chemical
  potentials with matched constants (achieved 1.105e-12).
- **Independent route:** Michelsen's equations re-derived from the model; the
  one-dimensional equal-fugacity solve; FeOs's chemical potentials; and a
  **synthetic two-branch equation of state** whose stationary point is known in
  closed form (`ln W_i = ln z_i + ln phi_i^V - ln phi_i^L`,
  `tpd = -logsumexp(ln W)`), which exercises the same code path at
  `ln W = 1497.6` with no polymer and no PC-SAFT in sight.
- **Suite time:** `pytest -q` **237.18 s (3:57) for 830 tests at `aa4a666` ->
  250.95 s (4:10) for 852**, uncontended on the one machine; `pytest -q -m slow`
  **694.93 s (11:34) for 34 -> 699.63 s (11:39) for 40**. 22 tests added, none
  removed and one reversed in place; the six marked `slow` are all
  *repetitions* - the 1 MPa half of a two-pressure statement whose 0.5 MPa half
  runs by default, and four further pressures on the same 0.4-1.4 MPa band. The
  new golden path costs about 4.5 s inside `tests/test_examples.py`. 4:10 is a
  little over the ~4 min target and is stated rather than trimmed further, as
  it was at the previous slice.
- **Test path:** `tests/test_stability_log_space.py` (21 by default, 6 `slow`),
  `tests/test_flash_log_space_stage.py` (the former pinned miss, reversed).
- **Script:** `python examples/validation/22_stability_log_space.py` (routes 1,
  2, 3 and 5 need no optional dependency; about 5 s; `--full` adds the
  144-state Peng-Robinson stability grid and the 0.4-2 MPa scan - and, since
  ADR-0026, routes 6 to 8 of Case P-16, which take it to about 1 min 40 s).

---

## Case P-16: The last six refusals of the 0.3-3.6 MPa polymer sweep

- **Source:** for the vapour-liquid state, the **one-dimensional** equal-fugacity
  solve of Case P-14, written independently in the example (the vapour taken as
  *exactly* pure solvent, one unknown carried as `ln x_solvent`,
  finite-difference Newton started a thousandth away from `flash_tp`'s answer);
  for the five liquid-liquid states, a **two-equation Newton** written in the
  example on the two equal-fugacity conditions in logarithms, in which the feed
  never appears, so it knows nothing about Rachford-Rice, the phase count, the
  stability test or the split; **FeOs 0.10.1** chemical potentials at
  chemthermo's converged phases and densities; and, for the Rachford-Rice
  defect, the binary closed form `beta = -(z_1 a_1 + z_2 a_2) / (a_1 a_2)` with
  `a_i = K_i - 1`.
- **Location:** `chemthermo/flash/_log_space.py`
  (`_safeguarded_directions`, `log_space_split(curvature_safeguard=...)`),
  `chemthermo/flash/_split.py` (`_rr_denominators`,
  `_rachford_rice(convex_denominators=...)`,
  `_extended_rachford_rice(convex_denominators=...)`),
  `chemthermo/flash/_detect.py` (`_phi_phi_log_space`'s retry,
  `_flash_tp_tangent_plane`'s last-resort bracket); slice
  `flash-polymer-edge-cases` (ADR-0026). **No model equation, no stability
  solver, no seed, no tolerance and no public signature changed.**
- **Parameters and provenance:** exactly Cases P-13, P-14 and P-15's -
  polyethylene from `tests/fixtures/pcsaft/martini2009_polymers.json`
  (`m/M = 0.026300 mol/g`, so `m = 1393.9` at `Mw = 53 000` and `m = 431.32` at
  `Mw = 16 400`), n-pentane from Gross & Sadowski (2001) Table 1,
  `k_ij = -0.006` except for the FeOs route, which runs at `k_ij = 0` on both
  sides because feos 0.10.1 cannot be given one from Python. Nothing here is
  compared against measurement.
- **Assumptions:** monodisperse polymer; 453.0 K throughout; feed stated as a
  polymer **mass** fraction (5 wt%, so at `Mw = 53 000`
  `z = (7.163935601491655e-05, 0.9999283606439852)`).

### (i) The defect, measured before it was repaired

A 0.3-3.6 MPa sweep at 0.1 MPa steps over both molar masses - 68 states - run
against a worktree at HEAD `584c508`:

| P / MPa | message | measured |
| --- | --- | --- |
| 0.3 | "did not converge the phi-phi split in log mole numbers" | residual `3.901e-05` after 100 iterations |
| 2.8 | same | `1.078e-08` after 100 |
| 2.9 | same | `6.031e+00` after 100 |
| 3.0 | "neither the stability-seeded nor the Wilson K-values bracket a Rachford-Rice root" | `tpd_min = -1.2537e-03` |
| 3.1 | same | `-1.2021e-03` |
| 3.2 | same | `-1.1538e-03` |

All six are on the `Mw = 53 000` chain; all 34 `Mw = 16 400` states converged
before and after. **Two** causes, neither of them about polymers:

**(a) 0.3, 2.8, 2.9 MPa - the log-space stage crawls next to the trivial
solution.** ADR-0024's stage keeps the Newton direction only while it descends
the two-phase Gibbs energy and otherwise falls back to `-r`, whose *length* is
`|r|`. At 0.3 MPa the Newton direction is rejected at every one of the 100
iterations and the stage moves `ln n_solvent` from `-3.18038` to `-3.17648` -
by `4e-03`, where the answer is `3.2` away - with the residual flat at
`3.9e-05`. Both phases stay within a percent of the feed in the solvent: this
is the neighbourhood of the trivial solution, where the Gibbs Hessian in the
mole numbers is **indefinite** (measured eigenvalues `-1.337e-04` and
`+3.967e+06`), which is exactly when the Newton direction is not a descent
direction. Ruled out by measurement, not argument: the finite-difference
Jacobian (step-independent to four figures over `h = 1e-06 .. 1e-02`, and
`dr/dn` symmetric to `1.5e-07` at `h = 1e-04`); damping (the line search
accepts a full step, zero backtracks, at all 100 iterations); and a scaling
accident of the polymer row (the polymer residual converges in 8 iterations; it
is the solvent row that cannot move).

**(b) 3.0, 3.1, 3.2 MPa - the Rachford-Rice denominator cancels to zero.** The
stability seed at 3.0 MPa is `K = (1.11871119e-18, 1.00132622)`, which straddles
one, and `f(0) = +1.2545e-03`, `f(1) = -6.41e+13` bracket a root at
`beta = 0.94591`. `_rachford_rice` forms `t_i = 1 + v (K_i - 1)`; `K_i - 1` is
an **exact** `-1.0` for any `K_i` below the spacing of doubles at one, so at
`v = 1` the sum is an exact `0.0`, the positivity guard fires and `f(1)` is
reported `nan`. The Wilson fallback loses `K` the same way (a non-volatile
component's Wilson estimate is `1e-10`, ADR-0022), so the second attempt fails
for the first attempt's reason and the flash refuses.

### (ii) The six states, after, each against an independent solve

| P / MPa | verdict | route | fugacity residual | mass balance | `dG/RT` | post-split | independent |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.3 | VLE | `stability-log`, safeguarded | 4.865e-13 | 1.11e-16 | -1.0864e-01 | stable | melt `x_C5` **1.4e-14** absolute |
| 2.8 | LLE | `stability-log`, safeguarded | 9.095e-13 | 1.11e-16 | -4.1682e-03 | stable | tie line **2.5e-13** relative |
| 2.9 | LLE | `stability-log`, safeguarded | 9.095e-13 | 1.11e-16 | -3.9940e-03 | stable | **7.1e-13** |
| 3.0 | LLE | `stability`, convex bracket | 7.276e-12 | 1.11e-16 | -3.8290e-03 | stable | **5.7e-13** |
| 3.1 | LLE | `stability`, convex bracket | 3.865e-12 | 1.11e-16 | -3.6725e-03 | stable | **1.0e-12** |
| 3.2 | LLE | `stability`, convex bracket | 2.274e-13 | 1.11e-16 | -3.5238e-03 | stable | **2.8e-13** |

The compositions, which are what the independent solves reproduce:

| P / MPa | polymer-rich `x_polymer` | polymer-lean `ln x_polymer` | polymer-rich phase fraction |
| --- | --- | --- | --- |
| 0.3 (melt) | 0.036808345719672377 | -1528.39130415 (the vapour) | 0.0019462802420003 |
| 2.8 | 9.231250339583601e-04 | -93.81711157873 | 0.07760525755403569 |
| 2.9 | 9.084501098780754e-04 | -90.20777508681 | 0.07885887759376387 |
| 3.0 | 8.940738921079231e-04 | -86.82803360271 | 0.08012688509001775 |
| 3.1 | 8.799780999223999e-04 | -83.65315387491 | 0.08141038512348653 |
| 3.2 | 8.661464766005232e-04 | -80.66231922913 | 0.08271043980469528 |

At 0.3 MPa the melt's solvent mole fraction from the flash is
`0.963191654280328` and from the one-dimensional solve
`0.963191654280342`; the polymer's own equal-fugacity condition in logarithms
agrees to `0.0` (`ln f_polymer = -1577.4489967184` on both sides), which is the
equation no linear parametrization can write down. The five liquid-liquid tie
lines come from a two-equation Newton started **one percent away in `ln x`**
from the flash's answer, so it is not handed its own result.

**FeOs at chemthermo's phases**, all six states at `k_ij = 0` on both sides:
worst `|d mu_i / RT|` between the two phases **6.139e-12** with FeOs's
fourteen-figure dispersion constants substituted in, **3.145e-07** as shipped -
the same ten-versus-fourteen-figure difference reported in Cases P-6 to P-12,
not a disagreement about the equilibrium.

### (iii) The iteration, before and after, at 0.3 MPa

Same seed, same model, same tolerance; only the direction rule differs.

| iteration | ADR-0024 rule | with the curvature safeguard |
| --- | --- | --- |
| 0 (seed) | 1.494998e+03 | 1.494998e+03 |
| 1 | 1.597822e-01 | 7.917254e+01 |
| 2 | 2.772935e-02 | 1.281896e+02 |
| 4 | 1.053965e-03 | 2.597397e+02 |
| 6 | 7.351891e-05 | 7.922556e+01 |
| 8 | 3.898817e-05 | 1.682776e-01 |
| 10 | 3.898861e-05 | 2.094748e-04 |
| 11 | 3.898883e-05 | 9.641553e-09 |
| 12 | 3.898905e-05 | 1.136868e-12 |
| 13 | 3.898926e-05 | 4.865205e-13 |
| 20 | 3.899078e-05 | 4.865205e-13 |
| 100 | 3.900817e-05 | 4.865205e-13 |

The safeguarded sequence gets **worse** before it gets better - iterations 1 to
6 are the modified-Newton direction climbing out of the trivial basin along the
negative-curvature eigenvector - and is then quadratic. The unsafeguarded one
never leaves.

### (iv) The sweep, and bit-identity

The 68-state sweep re-run against a worktree at `584c508`, comparing every
field a caller can observe - phase names, both compositions, phase fractions,
`vapor_fraction` and every diagnostics number - with `==`:

| outcome | states |
| --- | --- |
| identical, bit for bit | **62** |
| `ConvergenceError` -> a verified split | **6** (the table above) |
| converged before, moved in any digit | **0** |
| converged before, raises now | **0** |

Verdict sequences, both chains, 0.3-3.6 MPa: **VLE for 23 states then LLE for
11**, changing character exactly once, between 2.5 and 2.6 MPa for each - and
for `Mw = 16 400` the bisected boundary of Case P-13 is unmoved at
2.5868/2.5873 MPa. Worst equal-fugacity residual over all 68: `9.66e-13`
(`Mw = 16 400`), `1.23e-11` (`Mw = 53 000`). The `-> single-phase` leg of the
sequence is above this band and is checked by the 0.3-12 MPa scans in
`examples/validation/21_pcsaft_polymer_vle.py --full`.

**The 2.5 MPa vapour-liquid state, which Case P-15 left converging through the
log-space stage, is unchanged**: `converged_stage = "second-order-log"`, 33
stage iterations, `log_space_residual = 6.821210263296962e-13`, melt
`x_polymer = 0.0012513825804704919`, `beta = 0.9427518353436072` - the same
doubles as before the slice.

### Negative controls and what is *not* claimed

- **The curvature safeguard is never used on a first attempt.** It is a second
  call to the same function from the same seed, made only where the first call
  finished above `FlashSettings.tol` - which before this slice was a
  `ConvergenceError`. `diagnostics["log_space_curvature_safeguard"]` is `False`
  on every log-space state that converged before (measured on the sweep) and
  `True` on the three it repairs.
- **The convex Rachford-Rice denominator is off by default, and that is
  deliberate.** Switching it on everywhere was implemented and measured: it
  re-routes **15 of the 62** already-converging sweep states - `Mw = 16 400` at
  2.2-2.5 MPa from `stability-log`/`second-order-log` to
  `stability`/`second-order`, and 2.6-3.6 MPa in the last bits - because a seed
  that had had no in-window root suddenly has one. Those answers are as correct
  as before; moving them is a different slice. Asserted instead, with `==`:
  over 400 random K-sets and every K-vector of five real Peng-Robinson split
  iterations, `convex_denominators=True` returns the *same double* wherever the
  naive form had an answer at all.
- **The fixture did not move.** `refactor_bit_identity_v3.json` (155 states)
  passes **unchanged** and was not regenerated; none of its states reaches the
  log-space stage, so none of them gained or lost a diagnostics key. All nine
  ADR-0023 benchmark cases report **identical result hashes**
  (`python -m chemthermo.bench --compare`, the two records measured back to back
  on one machine; per-case ratios 0.99x-1.03x, i.e. noise).
- **A narrow band near 5.2 MPa on the `Mw = 53 000` chain still refuses**, and
  is neither caused nor repaired here. Measured at 5.175, 5.2 and 5.3 MPa
  (5.0, 5.1, 5.4 and 5.5 MPa converge), identically at `584c508`, with the same
  message and residual (`5.8e-04` after 100 successive-substitution and 101
  second-order iterations). It is outside the 0.3-3.6 MPa band this case
  measures and it enters through the *other* log-space entry point
  (`_phi_phi_second_order`, seeded from a failed linear iterate); the symmetric
  retry there was implemented and measured **not** to repair it, so it was not
  shipped (ADR-0002: no unexercised paths).
  `examples/validation/21_pcsaft_polymer_vle.py --full` prints it.
- **Nothing here is a statement about polyethylene.** The polymer parameters
  are one open secondary source citing a paywalled table that was not read, the
  polymer is monodisperse, and no number in this case is compared against
  measurement.
- **Suite time:** `pytest -q` **239.83 s (3:59) for 839 tests at `584c508` ->
  251.77 s (4:11) for 852**, uncontended on the one machine; `pytest -q -m slow`
  **751.87 s (12:31) for 53 -> 818.06 s (13:38) for 56**. 16 tests added, none
  removed; the three marked `slow` are all *repetitions* - three further
  pressures on the same liquid-liquid tie line, of which two (2.8 and 3.0 MPa,
  one per repair) run by default. 4:11 is a little over the ~4 min target and
  is stated rather than trimmed further, as it was at the previous two slices;
  the starting point on this machine was already 3:59.
- **Independent route:** the one-dimensional equal-fugacity solve; a
  two-equation Newton on the tie line in which the feed does not appear; FeOs's
  chemical potentials; the binary closed form of the Rachford-Rice root; and,
  for the direction rule itself, a Jacobian written out by hand in the test
  whose Gibbs Hessian is indefinite.
- **Test path:** `tests/test_pcsaft_polymer.py` (the six states, 3 by default
  and 3 `slow`, plus the K-value the bracket lost),
  `tests/test_flash_log_space_stage.py` (the stall, the repair, the contract of
  `_safeguarded_directions`, the reported flag),
  `tests/test_rachford_rice_extended.py` (the cancellation, the root it hid,
  the bit-identity of the flag being off).
- **Script:** `python examples/validation/22_stability_log_space.py --full`
  (routes 6, 7 and 8: the before/after iteration trace, the six states against
  their independent solves, and the whole 68-state sweep; about 1 min 40 s in
  total) and `python examples/validation/21_pcsaft_polymer_vle.py --full` (the
  two 0.3-12 MPa verdict scans, about 1 min 15 s).

---

## Case R-MAP-1: The robustness map - what every family refuses, and where

- **Source:** none, and like Case B-1 that is the point. This is not a
  thermodynamics case: it is the *coverage* case. Nothing here is checked
  against a published number or an independent implementation; what is measured
  is where `flash_tp` refuses, how it refuses, and whether anything it returns
  breaks its own invariants. The thermodynamic content of the grids is pinned
  elsewhere in this ledger (F-1, F-4, L-1, L-2, P-7..P-9, P-13..P-16, R-1..R-4,
  V-1..V-5); this case is the map drawn over them plus the states none of them
  reaches.
- **Where:** ADR-0027; harness `src/chemthermo/bench/robustness.py`; record
  `benchmarks/robustness_87f0820.json`; summary
  `benchmarks/robustness_87f0820.md`; `benchmarks/README.md`.
- **AMENDED by ADR-0028 (Case P-17).** All **36** refusals below are retired:
  the map at `9adf390` is **2110 of 2110 converging, 0 refusing, 0 violating an
  invariant**. Every diagnosis in section (ii) held - including the lever-rule
  number `4.074164e-04` of class 2, which the repaired solver reproduces as
  `4.074163e-04` - and the two classes recorded there as "not diagnosed
  further" turned out to be one defect with class 2. The `87f0820` JSON record
  was pruned when `9adf390` superseded it (`benchmarks/README.md`); its `.md`
  summary is kept, and every count quoted in this case is that summary's, so
  nothing below has been rewritten.
- **Assumptions:** every state is `flash_tp` at default `FlashSettings`
  (`max_phases=3`, `phase_detection="tangent-plane"`, `post_split_stability=True`).
  Wall times are a property of one machine at one moment - Apple M2 Max
  (12 cores, macOS 26.6.2, arm64), CPython 3.11.6, numpy 2.4.2 - and the tree
  was dirty at measurement because the harness module itself was uncommitted;
  every solver, model and parameter file was at `87f0820`.
- **Components and units:** 2110 states over six families; the grids are fixed
  in the module and each is documented at its definition. SI throughout.
- **Two grids carry values that are not cited, and are labelled so in the
  record:** the two nonzero PC-SAFT `kij` binaries use **illustrative**
  values (0.12 for CO2 / n-decane, 0.045 for methane / n-decane) because
  chemthermo packages no `kij` dataset (ADR-0014), and 7 of the 16 Tessier
  feeds are midpoints this module added for composition coverage. Neither is
  compared against anything. The Peng-Robinson grid is a **reconstruction** of
  the Case F-1 wide scan: Case F-1 names its shape (8 databank mixtures, 11 T,
  13 geometric P) but no file in this repository lists all eight mixtures, so
  the eighth here (Nitrogen / Methane) is this module's choice.

### (i) The map at `87f0820`

| family | states | verdicts | refusal classes | invariant violations | worst mass balance | worst equilibrium residual | worst dG_split/RT | time / s |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `pr-phi-phi` | 1270 | VLE 319, single-liquid 406, single-vapor 545 | - | 0 | 2.12e-13 | 7.79e-08 | -1.20e-05 | 5.6 |
| `pcsaft` | 204 | VLE 132, single-liquid 72 | - | 0 | 1.86e-13 | 3.58e-08 | -3.95e-04 | 72.7 |
| `pcsaft-associating` | 260 | LLE 145, VLE 5, single-liquid 89, single-vapor 21 | - | 0 | 5.07e-13 | 2.83e-09 | -5.34e-07 | 470.9 |
| `modified-raoult` | 108 | LLE 48, VLE 28, VLLE 4, single-liquid 13, single-vapor 15 | - | 0 | 1.11e-16 | 9.23e-13 | -1.07e-06 | 3.8 |
| `gamma-gamma` | 16 | LLE 15, single-liquid 1 | - | 0 | 2.22e-16 | 1.43e-13 | -1.07e-06 | 0.7 |
| `polymer` | 252 | LLE 113, VLE 43, single-liquid 60 | split-non-convergence 34, stability-inconclusive 2 | 0 | 2.22e-16 | 2.09e-07 | -4.12e-07 | 311.9 |

**2074 of 2110 converge, 36 refuse, 0 converge and violate an invariant**, in
865.6 s (14:26) total. The map was run **twice**, about twenty minutes apart on
the one machine: every state count, every verdict, every refusal class and
every worst residual in this table is identical between the two, and only the
clock moved (865.6 s against 918.9 s). That is the reproducibility claim - a
wall time here is not portable, a bucket is. The headline is the shape rather
than the totals:
**every refusal in the whole map is polyethylene / n-pentane.** The 1858 states
of the other five families - the 1270-state Peng-Robinson scan, the 188-state
Case F-4 grid and its two `kij` binaries, 260 associating PC-SAFT states, the
Tessier ternary and quaternary and the water / 1-butanol binary - refuse
nothing.

### (ii) The refusal classes, ranked by count

All 36 are `polymer`; the counts are (class, stage, system):

| # | class | stage | count | system | where |
| --- | --- | --- | ---: | --- | --- |
| 1 | `split-non-convergence` | `phi-phi` | **27** | pe53000 (26) + pe16400 (1) | 1 wt%: 3.6-5.4 and 7.5-8.7 MPa; 15 wt%: 3.9-7.5 and 8.1 MPa; pe16400 5 wt%: 7.5 MPa |
| 2 | `split-non-convergence` | `beta-outside-window` | **4** | pe16400 | 1 wt%, 0.3 / 0.6 / 0.9 / 1.2 MPa |
| 3 | `split-non-convergence` | `log-space` | **3** | pe53000 | 1 wt%, 0.3 / 0.6 / 0.9 MPa |
| 4 | `stability-inconclusive` | `feed` | **2** | pe53000 | 15 wt%, 10.5 and 10.8 MPa |

Six of the eight declared refusal classes are **empty** over 2110 states:
`rr-no-bracket`, `density-root-failure`, `post-split-third-phase`,
`multiphase-solver-failure`, `model-error` and `other-refusal`. That
`other-refusal` is empty is the claim that the classification rules cover what
actually occurs; that `rr-no-bracket` is empty is ADR-0016 and ADR-0026 holding
across a grid neither was measured on.

**Class 1 - the ADR-0026 "What remains" band is 27 states, not three.**
Every message is the `_phi_phi_second_order` one, e.g. `equal-fugacity
residual=3.897e-05 after 100 successive-substitution (max_delta_k=1.539e+89)
and 101 second-order iterations`, which is the entry point roadmap item 1 and
ADR-0026's "What remains" name. What the map adds is its *extent*: measured at
5 wt% the band is the three states near 5.2 MPa that Case P-16 records, and at
**1 wt% and 15 wt% the same defect spans 3.6-8.7 MPa**. This repository's whole
polymer history has swept one weight fraction. **Diagnosis: none attempted
here** - ADR-0026 already recorded that the symmetric curvature retry was
implemented at this entry point and measured *not* to repair it, so the open
question is unchanged and only its size moved.

**Class 2 - `beta = 0` where the lever rule says 4.1e-04.** PE 16400, 1 wt%
polymer (`z_polymer = 4.443385e-05`), 453 K, 0.3-1.2 MPa. `stability_tp` is
emphatic: `tpd_min = -4.550344e+02`, minimizing trial `wilson-liquid` in 3
substitutions, stationary point `w = (1.000000, 1.7297e-197)` - an essentially
pure polyethylene melt. The split then converges to `beta = 0.000000e+00` and
`_detect` refuses it because a single phase contradicts the verdict.
**Evidence that a two-phase answer exists:** the same system at 5 wt%, 453 K,
0.3 MPa converges to melt `x = (0.1090625, 0.8909375)` against vapour
`y = (4.2885e-206, 1.0)` at `beta_melt = 2.122456e-03`. A binary tie line at
fixed `(T, P)` does not depend on the feed, so the lever rule on that tie line
gives `beta_melt = 2.122457e-03` at the 5 wt% feed (reproducing the solver's
own number to 7 digits, which is the check that the lever rule is being applied
correctly) and **`4.074164e-04` at the 1 wt% feed** - strictly inside `(0, 1)`.
**Suspected cause, from that evidence:** the stability seed has the *same
shape* at both feeds (pure melt; `K_solvent` 6.66e-198 at 5 wt%, 1.73e-197 at
1 wt%), so what differs between the converging state and the refusing one is
the size of the incipient phase fraction - 2.1e-03 against 4.1e-04. Stated as
a hypothesis consistent with the measurement, not as a proof: no instrumented
trace of the split iteration was taken.

**Class 3 - three log-space refusals at the dilute end.** PE 53000, 1 wt%,
0.3 / 0.6 / 0.9 MPa, all through `_phi_phi_log_space` from the stability seed,
with residuals `8.009e-02` after 3 iterations, `3.277e-08` after 100 and
`1.657e-06` after 46. The middle one is 3.3e-08 against a `1e-08` acceptance -
close enough that it is worth saying plainly that no tolerance was touched
here. Same 1 wt% feed as class 2 on the other molar mass, which is the reason
to suspect the two share a cause. **Not diagnosed further.**

**Class 4 - a two-state window where no stability trial converges.** PE 53000,
15 wt% (`z_polymer = 2.402e-04`), 10.5 and 10.8 MPa. All four trials end
`second_order_no_progress` after 62-67 iterations (`max_iter = 300`), `tpd` is
`nan` *because* they did not converge, with stationarity residuals **0.685**
and **1.746**. The neighbours settle it: 10.2 MPa converges to `tpd_min =
+2.289939e-04` at a near-trivial minimizer (`K = (1.765e-02, 1.000236)`), and
11.1 MPa to `tpd = 1.4e-15`, the trivial solution. **Two causes are ruled out
by measurement.** It is *not* the ADR-0025 exponential-range problem: the
largest `|ln W|` at the stopping point is **11.09**, nowhere near the 700 that
clamp sat at, and `sum_W` is within 3.9e-04 of exactly 1. It is not a density
root failure either - no root message is raised anywhere in the family.
**Suspected cause:** the second-order stability stage making no progress at a
marginal, near-trivial stationary point, in a window where the verdict is about
to cross from "stable with a small positive tpd" to "trivial". Consistent with
brain.md section 10's standing caveat that a phase count is never better than
the stability test that produced it.

### (iii) The invariants: nothing returned violates one

Checked on every one of the 2074 converged answers: mass balance recomputed
from the returned phase fractions and compositions (`< 1e-10`); each phase's
composition summing to one (`< 1e-10`); phase fractions strictly in `(0, 1)`
and summing to one; the solver's reported equilibrium residual (`< 1e-06`);
`dG_split/RT < 0` on every multiphase answer; the post-split stability verdict
`"stable"` wherever it ran.

- **0 violations.** Worst mass balance anywhere **5.07e-13** (associating
  PC-SAFT), worst equilibrium residual **2.09e-07** (the polymer ternary,
  against 1e-06), worst - i.e. least negative - `dG_split/RT` **-3.95e-04**
  (PC-SAFT).
- **The honest limit:** mass balance and the composition sums are recomputed
  here; the equilibrium residual, `dG_split` and the post-split verdict are
  read from the solver's own `diagnostics`. A defect that made the solver
  mis-report its own residual would pass this check. ADR-0027 records why the
  duplicate verifier was not written.

### (iv) One thing the map surfaced and does not adjudicate

PC-SAFT 2B water / ethanol at `kij = 0` returns **LLE at 19 of its 60 states**
(all of `z_water = 0.8`, and `z_water = 0.5` at 290 K). Example: `z = (0.5,
0.5)`, 290 K, 0.1 MPa gives `beta = (0.039887, 0.960113)` with
`dG_split/RT = -2.442940e-04`, mass balance 7.9e-15, equilibrium residual
6.7e-10 and both phases post-split stable - so the *solver* returned a verified
lowest-Gibbs answer for the model it was given. Water and ethanol are fully
miscible in reality, so this is a statement about a **parameterization** (a
water / alcohol cross pair at `kij = 0`), not about `flash_tp`. Recorded here
because a coverage map is the place a result like that becomes visible;
adjudicating it needs a reference and a `kij`, and is not this case.

### (v) What this case does not establish

- Not correctness: no comparison against a published number or another
  implementation anywhere in it.
- Not a bit-identity baseline: phase compositions are deliberately absent from
  the record (ADR-0027).
- Not portable wall times.
- Not a global stability proof: a `single-*` verdict is only as good as the
  deterministic trial set that produced it.
- Not representative sampling in `--quick`: the 171-state subset is
  cost-bounded so the default suite stays inside its budget, and its offsets
  and strides were chosen from the full run's per-state times.

- **Suite time:** `pytest -q` **257.33 s (4:17) for 852 tests -> 267.75 s
  (4:27) for 891**, measured back to back on the one machine (the pre-slice
  number is this machine *today*; the previous slice recorded 251.77 s for the
  same 852, which is the drift). 39 tests added, none removed or marked `slow`
  except one - the full 2110-state sweep, ~15 min. `--durations` puts the whole
  cost of this slice at **9.12 s**, the module-scoped `--quick` sweep, almost
  all of it PC-SAFT states. 4:27 is over the ~4:20 target and is stated rather
  than trimmed further, as at the previous two slices; a third run of the same
  suite in the same session measured 285.13 s, so read +-7 % into any of these.
  `pytest -q -m slow tests/test_robustness_map.py` passes in **879.98 s
  (14:39)** - the full 2110-state sweep, reproducing the committed record's
  totals, verdict counts and refusal-class counts. The rest of the `slow` set
  was not re-run in this slice; nothing here touches the code it covers.
- **Independent route:** none, by design - see the source note. The one
  independent argument in the case is the lever rule of class 2, which shares
  no code with the split it contradicts.
- **Test path:** `tests/test_robustness_map.py` (the classifier against the
  raised message strings, the quick subset against its committed counts, the
  three pinned refusals, and the full sweep behind `slow`).
- **Script:** `python -m chemthermo.bench robustness --quick` (~9 s) and
  `python -m chemthermo.bench robustness --out record.json --summary-out record.md`
  (~15 min); `--family NAME` runs one family; `--list` prints the grid.

---

## Case P-17: The map's 36 refusals, retired and each one checked

- **Source:** three solvers that share no code with `flash_tp`. For the
  vapour-liquid states, the **one-dimensional** equal-fugacity solve of Case
  P-14 (the vapour taken as exactly pure solvent, one unknown carried as
  `ln x_solvent`, finite-difference Newton). For the liquid-liquid states, the
  **two-equation Newton** of Case P-16 on the two equal-fugacity conditions in
  logarithms, in which the feed never appears, so it knows nothing about
  Rachford-Rice, the phase count, the stability test or the split. For the two
  stability states, a **trial-free scan of equation (2)** over a grid of trial
  compositions, built from `density_roots` and `ln_fugacity_coefficients`
  alone. Plus the **lever rule** on each converged tie line, and **FeOs 0.10.1**
  chemical potentials at chemthermo's converged phases and densities.
- **Location:** `chemthermo/flash/_log_space.py`
  (`lever_rule_phase_fraction`, `stability_seed_ladder`,
  `log_space_seed(phase_fraction=...)`), `chemthermo/flash/_detect.py`
  (`_walk_stability_seed_ladder`, and the two gates that call it),
  `chemthermo/stability/tp.py` (`_substitution_budget_retry` and the
  inconclusive-only retry in `stability_tp`); slice
  `flash-polymer-robustness` (ADR-0028). **No model equation, no tolerance, no
  step rule, no acceptance test, no default setting and no public signature
  changed.**
- **Parameters and provenance:** exactly Cases P-13 to P-16's - polyethylene
  from `tests/fixtures/pcsaft/martini2009_polymers.json`
  (`m/M = 0.026300 mol/g`, so `m = 431.32` at `Mw = 16 400` and `m = 1393.9` at
  `Mw = 53 000`), n-pentane from Gross & Sadowski (2001) Table 1,
  `k_ij = -0.006` except on the FeOs route, which runs at `k_ij = 0` on both
  sides because feos 0.10.1 cannot be given one from Python. Nothing here is
  compared against measurement.
- **Assumptions:** monodisperse polymer; 453.0 K throughout; feeds stated as
  polymer **mass** fractions (1, 5 and 15 wt%, the robustness map's grid).

### (i) The three defects, measured before anything was changed

The 36 refusals of Case R-MAP-1, traced state by state at `87f0820`:

| class | states | what the trace says |
| --- | ---: | --- |
| `split-non-convergence` / `phi-phi` | 27 | successive substitution **diverges** (`max_delta_k` = 1.54e+89 at 3.6 MPa, 5.46e+105 at 5.1 MPa, 2.69e+128 at 8.7 MPa) and both following stages continue from what it left |
| `beta-outside-window` (4) + `log-space` (3) | 7 | the log-space stage converges on, or crawls next to, the **trivial** solution, because its seed's phase fraction is 0.5 and the melt holds 4.07e-04 of the feed |
| `stability-inconclusive` | 2 | every trial ends `second_order_no_progress` at residual 0.685 / 1.746, from a substitution iterate handed over at iteration 50 while it was still travelling |

**Class 1, at 3.6 MPa in full.** The K-loop hands
`_phi_phi_second_order` a pair at `x = (6.45e-90, 1.0)` and
`y = (0.99255, 0.00745)` with `beta = 0.99999` - phase-II mole numbers
`beta * y` exceeding the feed by four orders of magnitude, i.e. not a split.
The linear stage pulls `beta` back to `6.9e-06` (correctly - its box is
`0 < n_i < z_i`) and minimizes from there; the ADR-0024 log-space retry that
follows is seeded from the *same* iterate through `seed_from_iterate`. Both
fail; the state raises at `3.897e-05`. The stationary point at that same state
is `w = (5.61e-09, 1.0)`, `tpd_min = -3.889e-05` - already in the caller's hand
and passed to neither stage.

Measured from the stationary point instead, at the three traced pressures:

| seed / rule | 3.6 MPa | 5.1 MPa | 8.7 MPa |
| --- | --- | --- | --- |
| linear iterate (shipped) | 3.90e-05, 100 it | 1.83e-05, 100 it | 1.03e-06, 45 it |
| stationary point, ADR-0024 rule | 2.05e-03, 100 it | 5.95e-04, 100 it | 2.06e-05, 100 it |
| stationary point + ADR-0026 safeguard | **4.55e-13, 11 it** | **9.09e-13, 10 it** | **5.55e-15, 13 it** |

The 15 wt% half of the band needs less: from the stationary point the
*unsafeguarded* stage converges in 5 to 8 iterations at all fourteen of its
pressures. So the 27 states are two sub-classes of one defect and the seed
alone settles half of them.

**Class 2, at 0.3 MPa on the `Mw = 16 400` chain.** The stage **converges**, to
`8.283e-13`, on `x^I = x^II = z` with `beta` an exact `0.0` - the trivial
solution, whose residual is identically zero and which is therefore an
attractor of the stage's "accept a step that lowers the residual" clause.
`_flash_tp_tangent_plane` refuses it, correctly, against a `tpd_min` of
`-455.03`. The seed's phase fraction is what puts it on that side:

| seed `beta` | outcome at 0.3 MPa |
| --- | --- |
| 0.0001 | trivial solution, 17 it |
| 0.01 | trivial solution, 16 it |
| **0.5 (shipped)** | **trivial solution, 17 it, with or without the safeguard** |
| 0.9 | trivial solution, 15 it |
| 0.99 | **7.21e-13, 8 it**, `beta = 0.9995926` |
| 0.999 | 3.30e-13, 8 it, same answer |
| 0.9999 | 4.55e-13, 6 it, same answer |

`1 - 0.9995926 = 4.074163e-04` is the melt fraction, and Case R-MAP-1 had
already derived `4.074164e-04` for it from the lever rule on the 5 wt% tie
line, before any of this converged. That agreement is the reason the diagnosis
was believed before the fix was written.

**Class 4, the stability stall.** Four things changed one at a time, at
10.5 and 10.8 MPa:

| what was changed | 10.5 MPa | 10.8 MPa |
| --- | --- | --- |
| nothing (shipped) | inconclusive, residual 0.685, 65 it | inconclusive, 1.746, 62 it |
| `second_order_max_iter` 100 -> 1000 | unchanged | unchanged |
| `second_order_max_step` 4.0 -> 0.5 | unchanged (0.685, 63 it) | unchanged (1.746, 64 it) |
| substitution only, `max_iter` 300 | residual **53.4**, `w -> (6.9e-27, 1)` | residual **36.2** |
| substitution only, `max_iter` 3000 | residual **53.4**, identical | residual **36.2** |
| `ssi_iterations` 50 -> 300, same Newton stage | **stable**, 3 of 4 trials converge at iteration 307 | **stable**, same |

So it is not the Jacobian, not the step cap, not the iteration ceiling, and not
"substitution is slow" - on its own substitution walks away. It is the
handover point.

### (ii) The 34 splits, after, each against an independent solve

Every previously refusing split state, at `k_ij = -0.006`, 453 K, with the
independent route run at each one (the 1-D solve for the seven vapour-liquid
states, the two-equation Newton for the twenty-seven liquid-liquid ones):

| quantity, over all 34 | worst |
| --- | --- |
| mass-balance residual | **1.11e-16** (asserted < 1e-12) |
| equal-fugacity residual reported | **3.18e-12** (asserted < 1e-08) |
| relative difference against the independent solve | **2.65e-12** |
| lever rule on the converged tie line, absolute | **2.22e-16** |
| `dG_split/RT`, least negative | **-7.01e-05** (all negative) |
| post-split stability verdict | `stable` at all 34 |

Two cross-checks inside that table are worth naming, because neither can be
produced by a shared code path:

- **The 1 wt% state lands on the tie line a 5 wt% feed already pinned.**
  `Mw = 53 000`, 0.3 MPa: the newly converging 1 wt% feed gives a melt at
  `x_polymer = 0.036808345719631735`, and Case P-16's 5 wt% state - pinned
  since ADR-0026, verified there against the 1-D solve - gives
  `0.036808345719672377`. Relative difference **1.1e-12**, from two different
  ladder entries.
- **1 wt% and 15 wt% meet at every shared pressure.** 3.9 MPa:
  `7.75657320184093e-04` against `7.756573201846113e-04` (6.7e-13); 8.1 MPa:
  `3.647955115718861e-04` against `3.647955115724242e-04` (1.5e-12). The 1 wt%
  states reach those through the safeguarded ladder entry and the 15 wt%
  states through the unsafeguarded one.

Pinned answers, for the record (melt or polymer-rich composition, its phase
fraction, and the lean phase's `ln x_polymer`):

| state | rich `x_polymer` | rich fraction | lean `ln x_polymer` | route |
| --- | --- | --- | --- | --- |
| 16400, 1 wt%, 0.3 MPa | 0.10906250886173996 | 4.0741633069e-04 | -472.8766 | lever-rule seed |
| 16400, 1 wt%, 1.2 MPa | 0.022154085500088714 | 2.0056728217e-03 | -443.2920 | lever-rule seed |
| 53000, 1 wt%, 0.3 MPa | 0.036808345719631735 | 3.7355015626e-04 | -1528.3913 | lever-rule seed |
| 53000, 1 wt%, 3.6 MPa | 8.131909284340129e-04 | 1.6908407133e-02 | -70.2204 | stationary point + safeguard |
| 16400, 5 wt%, 7.5 MPa | 1.129717115748221e-03 | 2.0129111537e-01 | -12.1850 | stationary point + safeguard |
| 53000, 15 wt%, 8.1 MPa | 3.647955115724242e-04 | 6.5834614796e-01 | -20.7155 | stationary point, unsafeguarded |

**FeOs at `k_ij = 0`, matched constants.** Three of the states, flashed at
`k_ij = 0` (where they still come out of the ladder) and handed to FeOs as
composition and density only: worst `|mu_i^I - mu_i^II| / RT` = **4.3e-12**
over the three, against an asserted `1e-07`. `k_ij = 0` is not a choice: feos
0.10.1 raises `RuntimeError: missing field k_ij` for every serialization tried
(Case P-12's note).

### (iii) The two stability states, after

Both return `stable`, on the trivial solution, with
`diagnostics["substitution_budget_retry"] = True` and three of four trials
converging at residuals `1.59e-12` and `1.14e-12`; `flash_tp` returns a single
liquid. Four independent reasons that verdict is the right one:

- **10.2 MPa below** is `stable` at `tpd_min = +2.290e-04` on a non-trivial
  stationary point, and **11.1 MPa above** is `stable` on the trivial one. The
  two states sit exactly where the non-trivial stationary point merges into
  the trivial one, and the merged branch is what the retry finds.
- **A trial-free scan of equation (2)** over 220 trial compositions spanning
  twelve decades of polymer content finds **no negative tangent-plane
  distance**: minimum `+8.03e-09` at 10.5 MPa and `+8.75e-09` at 10.8 MPa,
  both at the feed composition itself. The same scan gives `+6.49e-09`,
  `+7.27e-09` and `+9.45e-09` at 9.9, 10.2 and 11.1 MPa, so the sequence is
  monotone through the window and the two states are not special in it.
- The stationary point the stalled trials were heading for moves smoothly with
  pressure (`w_polymer` = 8.9e-07, 4.2e-06, 1.53e-05, 1.78e-05, 2.40e-04 at
  9.9, 10.2, 10.5, 10.8, 11.1 MPa, the last being the feed), which is the
  merge, not a discontinuity.
- The retry's own arithmetic is the shipped one: `stationarity_met` at
  `1.59e-12`, not a loosened tolerance.

### (iv) The map, re-run, and the dormancy claim

`python -m chemthermo.bench robustness` at `9adf390`, 893.1 s:

| family | states | verdicts | refusals | invariant violations |
| --- | ---: | --- | ---: | ---: |
| `pr-phi-phi` | 1270 | VLE 319, single-liquid 406, single-vapor 545 | 0 | 0 |
| `pcsaft` | 204 | VLE 132, single-liquid 72 | 0 | 0 |
| `pcsaft-associating` | 260 | LLE 145, VLE 5, single-liquid 89, single-vapor 21 | 0 | 0 |
| `modified-raoult` | 108 | LLE 48, VLE 28, VLLE 4, single-liquid 13, single-vapor 15 | 0 | 0 |
| `gamma-gamma` | 16 | LLE 15, single-liquid 1 | 0 | 0 |
| `polymer` | 252 | LLE 140, VLE 50, single-liquid 62 | 0 | 0 |

**2110 of 2110 converge, 0 refuse, 0 converge and violate an invariant.** The
`polymer` family's verdict counts move by exactly the 36: LLE 113 -> 140 (+27),
VLE 43 -> 50 (+7), single-liquid 60 -> 62 (+2), which is class 1, classes 2+3
and class 4 in that order.

**Dormancy, measured rather than argued.** The two records were compared state
by state, on every field the record carries - bucket, verdict, phase names,
phase fractions, `converged_stage`, `phase_set_history`, and each of the five
residuals (`composition_sum`, `mass_balance`, `phase_fraction_sum`,
`equilibrium`, `delta_g_split_rt`) - with `!=`, wall time excluded:
**0 differences across all 2074 previously converged states.** The worst
residuals in every family are the same doubles as at `87f0820`.

### Negative controls and what is *not* claimed

- **Every gate is on a failure path.** `_phi_phi_second_order` walks the ladder
  only above `FlashSettings.tol`, i.e. on the line that raises;
  `_phi_phi_log_space` walks it from entry 3 only above `tol` **or** outside
  `(0, 1)`, both of which ended the flash before; `stability_tp` retries only
  on `inconclusive`. A converged answer cannot reach any of them, which is why
  the 2074-state comparison above is a confirmation rather than the argument.
- **When nothing in the ladder works, nothing is adopted.** The caller keeps
  the exact state it had and raises the message it raised before, so a state
  that still refuses refuses identically. (No state in the map does, but the
  code path is the one that would.)
- **The fixture did not move.** `refactor_bit_identity_v3.json` (155 states)
  passes unchanged and was **not** regenerated. All nine ADR-0023 benchmark
  cases report **identical result hashes** against `after_2ca41bf.json`.
- **`_phi_phi_second_order`'s gate was deliberately left narrower than
  `_phi_phi_log_space`'s.** Adding "or a phase fraction outside (0, 1)" there
  too is symmetric and no measured state exercises it; a synthetic test
  (`tests/test_flash_phi_phi_second_order.py::test_a_converged_vapor_fraction_outside_the_unit_interval_raises`)
  showed what shipping it would do - it turns that injected collapsed split
  into a real two-phase answer and the `beta`-outside-window guard stops being
  reachable from that entry point. Not shipped (ADR-0002).
- **The ladder is three entries because three defects were measured.** It is
  not claimed to be complete, and a fourth geometry may need a fourth entry.
- **Nothing here is a statement about polyethylene.** The polymer parameters
  are one open secondary source citing a paywalled table that was not read, the
  polymer is monodisperse, and no number in this case is compared against
  measurement.
- **Cost.** A refused state used to be cheap. The two stability states are now
  the 4th and 5th slowest in the whole map (8.77 s and 8.44 s, against the
  16.63 s associating-PC-SAFT leader) because they run two full trial sets, and
  the `polymer` family's wall time rises 311.9 s -> 352.1 s (+12.9 %) for
  36 more answers. The other five families move by less than 3 %, which on this
  machine is noise.
- **Suite time: parity, and the ~4:20 target is not reliably met.** Five runs
  on the one machine, the baselines in a detached worktree at `3ab783a` with
  `PYTHONPATH` pointing at its own `src` (the editable install otherwise
  resolves to the working tree, so a worktree would measure this code):

  | tree | `pytest -q` | tests |
  | --- | --- | ---: |
  | `3ab783a` | 262.09 s (4:22) | 891 |
  | this slice | 258.35 s (4:18) | 881 |
  | this slice | 267.06 s (4:27) | 881 |
  | `3ab783a` | **261.67 s (4:21)** | 891 |
  | this slice | **264.22 s (4:24)** | 881 |

  The last two are the pair that was run **immediately back to back**, and they
  are the ones to read: **261.67 s -> 264.22 s, i.e. +1.0 %, which is inside
  this machine's noise** (the same tree read 258.35 s and 267.06 s, a 3.4 %
  spread, and Case R-MAP-1 measured 7 %). So the honest statement is *parity*,
  not a saving, and **the ~4:20 target is met on some runs and missed on
  others - as it was before this slice, the pre-slice suite reading 4:21 and
  4:22 on this machine.** The trim is what bought the parity: it moved 20
  pre-existing parameters (~13 s) into `slow`, which is roughly what this
  slice's own new default coverage costs (~17 s).

  Where the remaining time is, measured: `tests/test_examples.py` is **85.6 s**
  of the ~264 s, 69.7 s of it the 23 validation examples. That is a smoke test
  with one run per example, so none of it is a *repetition* and none of it may
  be marked `slow`; the dev contract's mechanism for it is a `--full` flag on
  the example itself, which is a different change and was not made here.

  23 tests added, none removed; 969 collected against 948. **31 parameters
  moved into `slow`**, 11 of them this slice's own and 20 pre-existing, and
  every one is a *repetition* of something the default run still does:

  | where | newly `slow` | what still runs by default |
  | --- | ---: | --- |
  | `test_flash_phi_phi_second_order` previously-failing states (Case F-4) | 2 | 240 K / 1 MPa, which is also the state the independent Newton comparison uses |
  | `test_modified_raoult_vs_thermo` bubble/dew against `thermo` | 1 | the other feed, all four boundaries |
  | `test_pcsaft_lle_vs_feos` (two tests) | 2 | the high-pressure state of each |
  | `test_pcsaft_flash_vs_teqp` 300 K tie line | 4 | 3 of 7 pressures: both ends and one middle |
  | `test_pcsaft_vs_teqp` coarse saturation remark | 1 | the `pure_VLE_T` comparison at the same state |
  | `test_vlle_water_propanol_butanol` (two tests x 3 T) | 6 | the first feed inside each tie triangle |
  | `test_flash_modified_raoult` bubble / dew scalar equations | 2 | the other composition of each |
  | `test_pcsaft_flash` two-phase verification | 2 | 2 of 4 feeds on the same isotherm |
  | this slice's Case P-17 block | 10 | one state per repair route, plus both dormancy states |
  | this slice's FeOs Case P-17 states | 1 | one state per ladder entry |

  None of them is the only test of a capability. This slice's own additions
  cost **32.6 s** in the default run before the trim and **17.3 s** after: the
  two stability states are 8.7 s each because the retry runs the whole analysis
  twice, so `stability_tp` is asserted by default and `flash_tp`'s
  single-liquid answer - a third analysis - is `slow`.

  `pytest -q -m slow tests/test_robustness_map.py` passes in **925.30 s
  (15:25)** - the full 2110-state sweep against the committed record's totals,
  verdict counts and refusal-class counts (879.98 s at the previous slice, on
  36 fewer answers). The `slow` sets of every module this slice touched pass
  too: `tests/test_pcsaft_polymer.py`, `tests/test_flash_log_space_stage.py`,
  `tests/test_flash_phi_phi_second_order.py`,
  `tests/test_flash_modified_raoult.py` and `tests/test_pcsaft_flash.py`
  together in **99.85 s** for 26 tests, and `tests/validation/` in **166.73 s**
  for 45.
- **Independent route:** the one-dimensional equal-fugacity solve; the
  two-equation tie-line Newton in which the feed does not appear; a trial-free
  scan of equation (2); the lever rule on every converged tie line; FeOs's
  chemical potentials; and the tie line of a 5 wt% feed pinned by an earlier
  slice, met by a 1 wt% feed here.
- **Test path:** `tests/test_pcsaft_polymer.py` (the Case P-17 block: one state
  per repair route by default, the rest `slow`), `tests/test_flash_log_space_stage.py`
  (the lever-rule bound, the seed's box, the ladder's contract),
  `tests/validation/test_pcsaft_polymer_vs_feos.py` (three states against FeOs),
  `tests/test_robustness_map.py` (the empty `PINNED_REFUSALS`, and the full
  sweep behind `slow`).
- **Script:** `python -m chemthermo.bench robustness --quick` (~10 s) and
  `python -m chemthermo.bench robustness --out record.json --summary-out record.md`
  (~15 min).

---

## Case R-MAP-2: Four new families, and the map is hard again

> **Amended by slice `flash-eos-multiphase-robustness` (ADR-0029, Case
> P-18).** Every measurement below stands as taken at `74820b8`; what changed
> afterwards is the answer to the question it asked. All **nine**
> `multiphase-solver-failure` states of (ii) - classes 2, 2 (cont'd) and 3 -
> now converge to verified answers, and the two root causes the diagnoses
> below deliberately stopped short of are in Case P-18 (i) and (ii): an
> irreversible phase removal, and a second-order Hessian that cannot be formed
> on a phase set holding a component at `x ~ 1e-12`. Class 1's five
> `rr-no-bracket` states are untouched and remain by design. The map's totals
> after the repair are in Case P-18 (v); this case's `74820b8` record is kept
> as the *before* measurement and its JSON is pruned per the
> `benchmarks/README.md` policy, its `.md` summary retained.

- **Source:** none, for the same reason as Case R-MAP-1 - this is the
  *coverage* case, not a thermodynamics one. What is measured is where
  `flash_tp` refuses on four grids ADR-0027's own roadmap named as gaps, how it
  refuses, and whether anything it returns breaks an invariant. Every
  diagnosis below is read off evidence already in the record (the message, the
  neighbouring states, a `post_split_stability=False` re-run, a hand-traced
  K-loop), the same standard R-MAP-1 held itself to.
- **Where:** ADR-0027's amendment (slice `robustness-map-coverage`); harness
  `src/chemthermo/bench/robustness.py`; record
  `benchmarks/robustness_74820b8.json` (pruned when ADR-0029's sweep
  superseded it, per the `benchmarks/README.md` policy); summary
  `benchmarks/robustness_74820b8.md` (kept); `benchmarks/README.md`.
- **Assumptions:** every state is `flash_tp` at default `FlashSettings`, as in
  Case R-MAP-1, with one exception the `gamma-phi-legacy` family exists to
  exercise: `flash_mode="gamma-phi"` uses the pre-`flash-auto-phase-detection`
  Wilson heuristic and never runs a stability test (ADR-0007, ADR-0008
  decision 4), so its `post_split_checked` diagnostics key is always `False`
  and the post-split invariant of `_invariant_violations` is silently and
  correctly dormant there - not exempted by a family-specific carve-out, just
  never reached, and stated here because reading the family's table without
  this note could look like an omission. Wall times are Apple M2 Max
  (12 cores), macOS 26.6.2, arm64, CPython 3.11.6, numpy 2.4.2; the tree was
  dirty at measurement (this slice's own harness and test changes were
  uncommitted); every solver, model and parameter file was at `74820b8`.
- **Components and units:** 395 new states over four families, added to the
  2110 of Case R-MAP-1 for **2505 total**; each grid is fixed in the module
  and documented at its definition (`EOS3P_*`, `GAMMA_PHI_LEGACY_*`,
  `NEARCRIT_*`, `ASSOC3_*`). SI throughout.
- **Values not cited:** `NEARCRIT_MEP_P0_PA` (8,076,427.99 Pa) and
  `NEARCRIT_MP_P0_PA` (1,884,640.42 Pa) are **located, not published** - a
  60-iteration bisection of `stability_tp`'s stable/unstable verdict at a
  fixed anchor temperature, run once and cited as a constant the way
  `EOS3P_T3_K` cites Case P-9's Newton solve. Nothing is compared against them.

### (i) The map at `74820b8`

| family | states | verdicts | refusal classes | invariant violations | worst mass balance | worst equilibrium residual | worst dG_split/RT | time / s |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `eos-three-phase` | 117 | LLE 46, LLL 8, VLE 22, VLLE 19, single-liquid 10, single-vapor 3 | multiphase-solver-failure 9 | 0 | 2.11e-12 | 2.80e-08 | -1.25e-04 | 651.6 |
| `gamma-phi-legacy` | 30 | VLE 8, single-liquid 3, single-vapor 14 | rr-no-bracket 5 | 0 | 1.02e-13 | 0.00e+00 | - | 0.0 |
| `pr-near-critical` | 104 | VLE 55, single-liquid 28, single-vapor 21 | - | 0 | 1.47e-13 | 3.33e-08 | -1.94e-07 | 2.1 |
| `pcsaft-associating-ternary` | 144 | LLE 84, VLE 8, VLLE 6, single-liquid 30, single-vapor 16 | - | 0 | 1.74e-13 | 2.98e-09 | -1.30e-04 | 670.6 |

**2491 of 2505 converge, 14 refuse, 0 converge and violate an invariant**, in
2232.8 s (37:13) for the whole map (the original six families cost 893.1 s of
that, matching Case R-MAP-1's own number to within this machine's noise; the
four new families cost the remaining 1339.7 s, almost all of it two families
of PC-SAFT-association three-phase states). `pcsaft-associating-ternary`'s 144
states include **6 verified VLLE** answers - the family was built to sweep
near a cloud point and found one - at no cost in refusals.

**Every state of the original 2110 is unchanged.** Compared field by field
against the committed `9adf390` record (`bucket`, `verdict`, `phases`,
`phase_fractions`, `invariant_violations`, `converged_stage`,
`phase_set_history`, `outcome`, `refusal_class`, `refusal_stage`): **0
differences over all 2110 states**, wall time excepted. This slice touches no
solver, model or parameter file, so that is the expected result, checked
rather than assumed.

### (ii) The refusal classes, ranked by count

All 14 are `multiphase-solver-failure` (9) or `rr-no-bracket` (5); the counts
are (class, stage, system):

| # | class | stage | count | system | where |
| --- | --- | --- | ---: | --- | --- |
| 1 | `rr-no-bracket` | `phi-phi` | **5** | `gammaphi-methane-ethane` | 200 K/3 MPa, 220 K/1 MPa, 240 K/2 MPa, 260 K/4 and 5 MPa |
| 2 | `multiphase-solver-failure` | `collapsed` | **4** | `eos3p-pcsaft-water-n-hexane-t3-scan` | `z_water = 0.05`, `T3 + 0.01` to `T3 + 1.0` K |
| 2 | `multiphase-solver-failure` | `split` | **4** | `eos3p-pr-water-ethanol-n-hexane` | `(0.2, 0.6, 0.2)` at 280 and 300 K; `(0.4, 0.4, 0.2)` and `(0.5, 0.3, 0.2)` at 300 K |
| 4 | `multiphase-solver-failure` | `collapsed` | **1** | `eos3p-pcsaft-water-ethanol-n-hexane` | `(0.1, 0.1, 0.8)`, 333 K |

Six of the eight declared refusal classes are still empty over the full
2505-state map (`stability-inconclusive`, `density-root-failure`,
`post-split-third-phase`, `multiphase-solver-failure`'s wider siblings
(`rachford-rice`, `search`, `feasible-region`), `model-error`,
`other-refusal`). The two classes that do appear were both already declared by
ADR-0027 and already exercised by the original grid on other states -
`multiphase-solver-failure` by the ADR-0011/ADR-0020 search machinery, which
this slice's `eos-three-phase` family is the first grid to stress hard, and
`rr-no-bracket` by ADR-0016's own motivating example, reproduced here on the
one code path ADR-0016 deliberately left unrepaired.

**Class 1 - the legacy path loses its bracket within 0-2 substitutions, every
time.** All five states are Methane/Ethane, NRTL (the packaged synthetic pair)
+ Peng-Robinson, `flash_mode="gamma-phi"`. Traced by hand at 200 K/3 MPa: the
Wilson seed `K = (1.978, 0.0723)` brackets cleanly (`f(0) = 5.8e-12`,
`f(1) = -6.17`), and the very first `K`-update from it -
`K = gamma^L phi^L / phi^V = (1.826, 0.174)` - already has `f(0) = -1.4e-04`
and `f(1) = -2.15`, both negative: the update has compressed the two `K`s
toward 1 (the synthetic NRTL pair and this light-hydrocarbon system are close
to ideal) far enough that the plain, non-extended Rachford-Rice this path is
pinned to (ADR-0016 decision 8: "the legacy `wilson-heuristic` path is
untouched... reproducing pre-ADR-0008 behavior is that path's entire purpose")
has nothing left to bracket. The other four states are the same shape at a
slower clock: 220/240/260 K lose the bracket on the first update as well
(`K` collapsing from Wilson estimates of 9.5/7.0/4.9 down to 8.4/5.3/3.2 in one
step), and 260 K/5 MPa survives two updates before losing it. This is not
ADR-0016's CO2/n-decane oscillation (which grows without bound over many
iterations) - it is a single-step compression, consistent with a near-ideal
system's `K`s having much less room between 1 and their Wilson estimate than a
CO2/n-decane pair does. **Not a defect in the sense the other two classes are**:
this is the deprecated path's documented behaviour (ADR-0007, ADR-0008
decision 4, ADR-0016 decision 8), reproduced rather than discovered, and the
tangent-plane path (`phase_detection="tangent-plane"`, the default since
ADR-0008) converges every one of these five states without incident - checked
directly: `flash_tp` on the same state without `flash_mode="gamma-phi"`
returns a verified two-phase split at all five.

**Class 2 - a hexane-rich water/n-hexane feed above T3 that the
phase-addition/removal search cannot resolve, an extent not previously
measured.** PC-SAFT water/n-hexane, `z_water = 0.05`, 1 atm,
`T = T3 + 0.01` through `T3 + 1.0` K (4 of the T3-scan's 9 offsets at this
feed; `T = T3` and everything at or below it converges to a plain `LLE` pair
with no search entered). With `FlashSettings(post_split_stability=False)` the
direct two-liquid split converges fine (`beta_1 = 0.0280`, `dG_split/RT =
-6.83e-03`) but is doubly unstable - `post_split_status = "unstable"` for
*both* phases - so the search is entered exactly as it is at `z_water = 0.5`
(Case P-9 (i)), where it resolves cleanly as `V -> LV -> LLV -> LL`. Here it
does not resolve within its round budget: the feed's own stability test is
emphatic (`tpd_min = -0.409`, `feed_branch = "liquid"`), so the diagnosis
stops at "the search enters and does not settle" rather than a traced root
cause, in the same spirit R-MAP-1 left its own class 1 - a defect's *extent*
is evidence; its mechanism is a further slice's work, not assumed here.
*(Traced in Case P-18 (i): the `LLV` set recedes, the recession direction
names the hexane-rich liquid - 95 % of the feed - as the phase leaving, and
the water-rich pair left over negative-flashes. The answer is the other
removal, and all four states now return a verified `VLE`.)*

**Class 2 (cont'd) - the one ternary state here is not a new defect.** PC-SAFT
water/ethanol/n-hexane, `(0.1, 0.1, 0.8)`, 333 K, is the *same* feed Case P-10
(i) already recorded: "That feed also raised before this slice... It sits near
the edge of the triangle where the water-rich liquid's amount is tiny." The
wider grid here rediscovers it rather than finding something new - checked:
with `post_split_stability=False` the direct two-liquid split is doubly
unstable (`post_split_status = "unstable"`) exactly as at `z_water = 0.05`
above, the same mechanism at a different system. *(Retired with class 2 by
ADR-0029; it returns a verified `VLE` - two phases, not three, the `LLV` set
at this feed having no Rachford-Rice solution. Case P-18 (i).)*

**Class 3 - the Peng-Robinson three-liquid region moves or narrows between 280
and 300 K, and one feed sits close enough to its edge to refuse at both.**
Water/ethanol/n-hexane, `k_ij = 0`. At 280 K, 11 of 12 swept feeds return the
verified `LLL` triangle of Case P-10 (ii) and one - `(0.2, 0.6, 0.2)`, 60 %
ethanol - refuses with a residual of `1.836e-08`, **1.8 times** the `1e-8`
acceptance, after the full 50-substitution-plus-1-second-order-iteration
budget: caught mid-convergence, not diverging. At 300 K - a temperature Case
P-10 (ii) never swept - the same feed's residual is `1.478e-04`
(**14,780 times** `tol`), and two more feeds refuse alongside it
(`(0.4, 0.4, 0.2)` at `2.04e-05`, `(0.5, 0.3, 0.2)` at `3.31e-06`), all three
converging fine at 280 K. Read together: a feed on the edge of the region
(`(0.2, 0.6, 0.2)`) is marginal at both temperatures and further from
convergence at the higher one, consistent with the three-liquid region having
moved or shrunk between 280 and 300 K rather than with a defect that appears
abruptly. Not traced further than that - which edge, and by how much, needs a
scan across temperature this slice did not run.

*Superseded by the trace in Case P-18 (ii), and the reading above turns out to
have been wrong in its causal half.* The region does not move or narrow: all
three 300 K feeds converge to the **same** tie triangle once the solver can
finish, and the spread of refused residuals (`3.31e-06` to `1.478e-04`) is
simply how far successive substitution had got in its 50 iterations before
handing over. What made the handover useless is that the second-order stage
aborts at its first Hessian column on a phase set whose water-rich liquid
holds n-hexane at `x = 5.3e-14` (280 K) to `1.6e-12` (300 K). The
"caught mid-convergence, not diverging" half of the reading was right, and is
what made the states cheap to repair.

### (iii) The invariants: nothing returned violates one, including the new one

Checked on every one of the 2491 converged answers, by the same rules Case
R-MAP-1 used, plus one this slice adds: a three-phase answer's
`delta_g_vs_two_phase_rt` (present only where the phase addition/removal
search ran and returned three phases) must be negative, mirroring the
existing `delta_g_split_rt < 0` check. **0 violations, including 0 on the new
check** - worst (least negative) `dG_split/RT` **-1.25e-04**
(`eos-three-phase`), which is comfortably negative; every one of the 19 `VLLE`
and 8 `LLL` answers in `eos-three-phase` and the 6 `VLLE` answers in
`pcsaft-associating-ternary` passed it. Worst mass balance anywhere
**2.11e-12** (`eos-three-phase`), worst equilibrium residual **3.33e-08**
(`pr-near-critical`), both comfortably inside their `1e-10` / `1e-6`
tolerances.

### (iv) What this case does not establish

The same five points as Case R-MAP-1 (v), plus: this is not a claim that
`gamma-phi-legacy`'s five refusals are defects to fix - they are the
deprecated path's documented, by-design behaviour, included so the map's
`rr-no-bracket` class is non-vacuous rather than left empty by construction.

- **Suite time:** `pytest -q` **894 passed, 91 deselected, 267.04 s (4:27)**,
  on the same machine the map itself was measured on.
  `pytest -q -m slow tests/test_robustness_map.py` **4 passed, 50 deselected,
  2291.70 s (38:11)** - the full 2505-state sweep against this record's
  totals (2232.8 s of that is the sweep itself, matching (i) above), plus the
  three new tests pinning the three-phase verdicts named in "Test path" below.
- **Independent route:** none, by design, matching Case R-MAP-1 - this case's
  only independent argument is the state-by-state comparison against the
  committed `9adf390` record (0 differences over 2110 states) and the
  `post_split_stability=False` re-runs used to trace classes 2 and 3, both of
  which share code with `flash_tp` itself and are therefore read as diagnostic
  evidence, not as a second solver.
- **Test path:** `tests/test_robustness_map.py` (`EXPECTED_QUICK_FAMILIES`
  gains four rows, `PINNED_REFUSALS` gains the cheap refusals the quick subset
  contains, and two `slow`-marked tests pin the three-phase verdicts Cases P-9
  and P-10 already established independently - `T3 - 0.05 K -> LLE`,
  `T3 + 0.05 K -> VLE`, two of the twelve 333 K tie-triangle feeds `-> VLLE`).
- **Script:** `python -m chemthermo.bench robustness --quick` (~14 s) and
  `python -m chemthermo.bench robustness --out record.json --summary-out record.md`
  (~37 min); `--family NAME` runs one family; `--list` prints the grid.

### Cross-platform (ADR-0032, 2026-09-25): the ladder's *route* is last-bit noise

On the cloud-baseline host (Linux 6.18 x86_64, CPython 3.11.15, numpy 2.4.6)
`test_the_band_the_diverged_k_loop_used_to_end[53000.0-0.15-8100000.0]`
failed on one assertion only: `log_space_seed` was `"linear-iterate"` (the
ADR-0024 retry from the linear stage's iterate) where macOS arm64 and the
GitHub Actions runner give `"stability-w"` (the ADR-0028 ladder). Every other
assertion held: the same tie line to 1.30e-12 relative (asserted 1e-11),
post-split stable, residuals inside their bounds.

Probed on that host over the 17 pressures within +-8 ULP of 8.1 MPa (same
chain, same feed):

| state | `stability-w` | `linear-iterate` | linear stage converged | `ConvergenceError` | worst fraction move among converged |
| --- | --- | --- | --- | --- | --- |
| Mw 53000, 15 wt%, 8.1 MPa | 12 | 3 | 1 | **1** | 4.73e-12 |
| Mw 16400, 5 wt%, 7.5 MPa | 17 (safeguard on) | 0 | 0 | 0 | 9.38e-12 |

Reading: the successive-substitution loop diverges here (`max_delta_k` 1e+89
to 1e+128), so the iterate it leaves - and therefore which rung finishes the
job first - is chaotic in the last bit of the input. The **answer** is not:
every converged outcome is the same tie line. The refusal one ULP-step away
is a **new robustness finding**: the map's "0 refusals in `polymer`" is true
of its grid points on the machines measured, not of their neighbourhoods.
Queued in `.agents/handoffs/continuation-plan.md` (not fixed in the
cross-platform slice, which changes no solver code).

- **Test change:** the route stays pinned to `"stability-w"` on macOS arm64
  and, everywhere, for the 16400 g/mol state (17/17 over its neighbourhood);
  off-platform the 53000 g/mol state admits `"stability-w"` or
  `"linear-iterate"` (`P17_ROUTE_BY_NOISE`), with every numeric assertion
  unchanged.

---

## Case P-18: The map's nine multiphase refusals, retired and each one checked

- **Source:** none for the *answers* - there is no published tie triangle for
  water / ethanol / n-hexane with `k_ij = 0` on either model, and ADR-0027's
  standing warning (the map measures coverage, not correctness) applies to
  every state here. What is established is narrower and is stated that way
  throughout: each returned phase set is an **equilibrium of the model it was
  given**, verified by a solver sharing no code with `flash_tp`, and it is the
  **lowest-Gibbs** one among every admissible alternative that was solved.
- **Where:** slice `flash-eos-multiphase-robustness` (ADR-0029); solver
  `src/chemthermo/flash/_multiphase.py` and the new
  `src/chemthermo/flash/_multiphase_log_space.py`; states from ledger Case
  R-MAP-2 (ii), classes 2 and 3.
- **Assumptions:** default `FlashSettings` throughout (`tol = 1e-8`,
  `ssi_iterations = 50`, `second_order_tol = 1e-12`, `max_phases = 3`,
  `post_split_stability = True`), 1 atm = 101 325 Pa, `k_ij = 0` on both
  models, PC-SAFT with the 2B association scheme for water and ethanol.
  Machine: Apple M2 Max, macOS 26.6.2, CPython 3.11.6, numpy 2.4.2.
- **Components and units:** Water / n-Hexane and Water / Ethanol / n-Hexane,
  mole fractions, SI throughout.
- **Parameters and provenance:** unchanged - PC-SAFT from Gross & Sadowski
  (2001, 2002) as packaged and provenanced in ADR-0014/ADR-0018;
  Peng-Robinson from the packaged critical constants. No parameter file was
  touched by this slice.

### (i) Five `collapsed` states: a removal that could not be taken back

PC-SAFT water / n-hexane, `z_water = 0.05`, 1 atm, at `T3 + 0.01`, `+ 0.1`,
`+ 0.5` and `+ 1.0 K` (`T3 = 334.807826336 K`, Case P-9's Newton), plus the
PC-SAFT ternary `(0.1, 0.1, 0.8)` at 333 K.

**Mechanism, traced.** The route is `L -> LL -> LLV -> LV -> raise`. Above
`T3` a binary at fixed pressure has no three-phase state (Gibbs' phase rule),
so `_multiphase_rachford_rice` correctly reports a receding feasible region
and `_multiphase_ssi` reads `argmin` of the recession rates. (On the ternary
the three-phase region is finite but this feed is outside it, so the same
solve stays inside the Okuno region and returns a converged **negative flash**
- `beta = (-0.761, -4.309, +6.070)` at the first substitution - and
`_multiphase_ssi` reads `argmin` of the fractions instead. Two spellings, one
signal.) At this feed that names the **hexane-rich** liquid - the phase
holding 95 % of the feed. What is left is the water-rich liquid against the
vapour, and that pair negative-flashes: `beta = (-0.20747, +1.20747)` at
`T3 + 0.01 K`. The old code had no way back and raised.

**Why that pair cannot be the answer, computed independently.** The
water-rich vapour-liquid tie line is a perfectly good tie line; it simply does
not contain this feed. Solved here by a 2-equation Newton on
`fugacity_coefficients`, its lever rule gives `beta` of **1.20747, 1.20882,
1.21490, 1.22273** at the four offsets - outside `(0, 1)` at every one.

**The answer, and the ordering.** The other removal leaves the hexane-rich
liquid against the vapour, and that is what `flash_tp` now returns:

| state | x_water (liquid) | x_water (vapor) | beta_vapor | equilibrium residual | dG_split/RT | dG vs two-phase |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `T3 + 0.01 K` | 0.02253935 | 0.21284192 | 0.14429996 | 2.95e-13 | -6.879e-03 | -4.902e-05 |
| `T3 + 0.1 K` | 0.02201259 | 0.21029813 | 0.14864345 | 3.74e-13 | -7.287e-03 | -4.967e-04 |
| `T3 + 0.5 K` | 0.01980545 | 0.19895982 | 0.16853927 | 3.55e-13 | -9.243e-03 | -2.630e-03 |
| `T3 + 1.0 K` | 0.01732053 | 0.18470644 | 0.19523428 | 7.43e-13 | -1.2035e-02 | -5.641e-03 |

Reduced Gibbs energies, every one computed in the test file from
`fugacity_coefficients` and none read out of `flash_tp`:

| state | returned (VL, hexane-rich) | LL (admissible, beta = 0.972) | single liquid | single vapor |
| --- | ---: | ---: | ---: | ---: |
| `T3 + 0.01 K` | **-0.341354043232** | -0.341305025321 (+4.902e-05) | -0.334474632405 | -0.236755238610 |
| `T3 + 0.1 K` | **-0.338912572488** | -0.338415918476 (+4.967e-04) | -0.331625703695 | -0.236726127984 |
| `T3 + 0.5 K` | **-0.328230732139** | -0.325600989922 (+2.630e-03) | -0.318988149938 | -0.236597106096 |
| `T3 + 1.0 K` | **-0.315281870162** | -0.309640821300 (+5.641e-03) | -0.303246938939 | -0.236436648189 |

The margin over the two-liquid pair grows monotonically with temperature,
which is the expected shape: `T3` is where the two coincide. It is also a
**cross-check on the solver's own bookkeeping**: the solver reports
`delta_g_vs_two_phase_rt` of `-4.902e-05`, `-4.967e-04`, `-2.630e-03` and
`-5.641e-03` (the table above), and the independently solved two-liquid pair
sits `+4.902e-05`, `+4.967e-04`, `+2.630e-03` and `+5.641e-03` above the
returned answer - the same four numbers, computed twice by routes that share
no solver code.

The ternary `(0.1, 0.1, 0.8)` at 333 K is the same mechanism on the
tie-triangle's hexane-rich edge, and Case P-10 (i) had already recorded it as
a pre-existing failure. It returns `liquid` `(0.0427106, 0.06676041,
0.89052899)` at `beta = 0.4758653` against `vapor` `(0.15201342, 0.13017844,
0.71780814)`, equilibrium residual **5.05e-13**, `dG_split/RT = -1.397e-02`,
`dG` vs the two-phase candidate **-1.209e-02**; `G/RT = -0.691284237040`
against **-0.677311243024** (single liquid) and **-0.673251317555** (single
vapor). Two phases, not three: the `LLV` set at this feed has no
Rachford-Rice solution.

### (ii) Four `split` states: a Hessian that could not be formed

Peng-Robinson water / ethanol / n-hexane, 1 atm: `(0.2, 0.6, 0.2)` at 280 K
and 300 K, `(0.4, 0.4, 0.2)` and `(0.5, 0.3, 0.2)` at 300 K. All four are
genuine **three-liquid** states.

**Mechanism, traced.** Every one of the four reports `1` second-order
iteration in its refusal message, and that is the whole diagnosis: the
water-rich phase holds n-hexane at `x = 5.32e-14` (280 K) to `1.56e-12`
(300 K), the linear stage perturbs a non-reference mole number by `1e-7`, the
reference phase `z - sum n` goes negative on the first Hessian column, and the
loop breaks having done nothing. Two further measurements decide the fix:

- successive substitution *does* converge these states, in **52, 109, 97 and
  86** iterations against the multiphase budget of
  `min(ssi_iterations, max_iter) = 50`, to residuals of 7.9e-09 to 1.0e-08.
  So the states are reachable and the budget is not the defect - the stage
  that was supposed to take over at 50 is;
- a log-space stage written as a straight generalization of ADR-0024 **does
  not fix it**: keeping the search's own reference phase, it accepts steps on
  an Armijo test whose energy change is `+-9e-16` (floating-point noise in
  `g`) and drives the residual from `1.24e-04` to `1.05e-01` in one step and
  to `1.12e-01` in 100. The reason is structural: `n_hexane` in the reference
  phase is `0.2 - 0.167 - 0.0326 = 2.0e-13`, a cancellation of `O(1)` doubles,
  and it is not a variable in any parametrization of the other phases.

With the reference phase re-chosen to the one whose smallest mole fraction is
largest, the same stage converges from the 50-iteration iterate in **1, 3, 2
and 2** Newton iterations - the `log_space_iterations` the results now carry -
and the phase sets it returns have all-pairs equal-fugacity residuals of
**1.42e-14, 8.88e-15, 1.47e-14 and 1.15e-14**, against an acceptance of
`1e-08`.

**The answers.** At 300 K all three feeds converge to the *same* tie triangle,
which is the evidence that the refused residuals were measuring distance and
not noise:

| state | beta | liquid1 (water-rich) | liquid2 | liquid3 | residual | dG_split/RT | dG vs two-phase | log-space its |
| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| `(0.2, 0.6, 0.2)` 280 K | (0.13802549, 0.54657701, 0.31539750) | (0.998808722, 1.19127789e-03, 5.32400754e-14) | (0.10277099, 0.85987453, 0.03735448) | (0.01891796, 0.41169611, 0.56938593) | 1.42e-14 | -1.5402e-01 | -1.9514e-02 | 1 |
| `(0.2, 0.6, 0.2)` 300 K | (0.12967906, 0.45994477, 0.41037617) | (0.997834456, 2.16554440e-03, 1.55819339e-12) | (0.11615007, 0.81294684, 0.07090309) | (0.04186195, 0.55024766, 0.40789039) | 8.88e-15 | -8.6890e-02 | -3.2215e-03 | 3 |
| `(0.4, 0.4, 0.2)` 300 K | (0.36061156, 0.18042342, 0.45896503) | same as above | same as above | same as above | 1.47e-14 | -4.1244e-01 | -1.1243e-03 | 2 |
| `(0.5, 0.3, 0.2)` 300 K | (0.47607780, 0.04066274, 0.48325946) | same as above | same as above | same as above | 1.15e-14 | -6.3795e-01 | -9.1739e-05 | 2 |

(The phases are reported `liquid1` / `liquid2` / `liquid3` ordered by the
first component's mole fraction, ADR-0019 decision 3; the middle and last
columns above are printed in solve order, which is why `liquid2` is the
ethanol-rich phase at 300 K and the hexane-rich one at 280 K.)

Gibbs ordering, each alternative re-solved in the test file:

| state | three liquids | best admissible pair | single liquid | single vapor |
| --- | ---: | ---: | ---: | ---: |
| `(0.2, 0.6, 0.2)` 280 K | **-3.829037668234** | -3.809524142507 (+1.951e-02) | -3.675022625340 | -0.991524792143 |
| `(0.2, 0.6, 0.2)` 300 K | **-2.683831312913** | -2.680609854845 (+3.221e-03) | -2.596941513508 | -0.984320350053 |
| `(0.4, 0.4, 0.2)` 300 K | **-2.864109965683** | -2.862985629731 (+1.124e-03) | -2.451669289897 | -1.084466353589 |
| `(0.5, 0.3, 0.2)` 300 K | **-2.954249292068** | -2.954157553154 (+9.174e-05) | -2.316301718028 | -1.057078515513 |

At 280 K the third pair (`liquid2 + liquid3`) is also admissible and also
higher (-3.779446163872); at 300 K the feeds `(0.4, 0.4, 0.2)` and
`(0.5, 0.3, 0.2)` have no admissible two-liquid pair other than the one
listed, their lever rules putting the feed outside the others.

### (iii) The independent route

Written in `tests/test_flash_eos_multiphase_robustness.py`, sharing no
*solver* code with `flash_tp` - the model is necessarily shared, through the
public `EquationOfState.fugacity_coefficients` interface, exactly as Cases
P-9 and P-10 share it: the full `NC - 1` equilibrium system - `(N - 1) C`
equal-fugacity equations plus `C - 1` independent mass balances, in
`N (C - 1)` composition degrees of freedom (carried as `ln(x_k / x_last)`, so
that a mole fraction of 5e-14 is representable) plus `N - 1` phase fractions.
That is **3** equations for the binary pair, **5** for the ternary pair and
**8** for the three-liquid set. Started off `flash_tp`'s answer by `1e-3` in
every variable, so convergence back onto it is a result and not an identity.

**Tolerances achieved** over all nine states: Newton residual
`4.89e-15 .. 8.97e-14`; `max |dx|` versus `flash_tp`'s compositions
**2.98e-12** (worst, `T3 + 0.01 K`) and `<= 3.0e-14` on all four
Peng-Robinson states; `max |dbeta|` **5.16e-12** (worst) and `<= 3.3e-14` on
the Peng-Robinson states. Mass balance, recomputed from the returned phases
rather than read from diagnostics: **0.0** on the five PC-SAFT states and
`<= 1.11e-16` on the four Peng-Robinson ones, against the asserted 1e-12.
Post-split stability: `"stable"` on all nine, worst post-split `tpd_min`
**-2.583e-12 (`T3 + 0.1 K`; three of the nine report a small *positive* value, which is what "no negative tangent-plane distance was found" looks like at this scale)**.

### (iv) FeOs, on the PC-SAFT states

FeOs's own chemical potentials evaluated at chemthermo's converged phases and
densities, with FeOs's universal constants substituted for the printed ones
(the standing Case P-6/P-7/P-8/P-9 caveat: as shipped, the constants floor any
such comparison near 2e-06). Worst difference across the phases of a state:
**1.61e-12** and **2.06e-12** on the water-lean binary states at
`T3 + 0.01` and `T3 + 1.0 K`, **3.51e-13** on the ternary `(0.1, 0.1, 0.8)`
at 333 K, against an asserted 1e-8.

### (v) The map after the repair

| family | states | verdicts | refusal classes | invariant violations | worst mass balance | worst equilibrium residual | worst dG_split/RT | time / s |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `eos-three-phase` | 117 | LLE 46, **LLL 12**, **VLE 27**, VLLE 19, single-liquid 10, single-vapor 3 | **-** | 0 | 2.11e-12 | 2.80e-08 | -1.25e-04 | 655.8 |
| `gamma-phi-legacy` | 30 | VLE 8, single-liquid 3, single-vapor 14 | rr-no-bracket 5 | 0 | 1.02e-13 | 0.00e+00 | - | 0.0 |

**2500 of 2505 converge, 5 refuse, 0 converge and violate an invariant**, in
2202.5 s (36:43) at `64831bd` (`benchmarks/robustness_64831bd.json` / `.md`;
tree dirty at measurement - this slice's own documentation was uncommitted -
with every solver, model and parameter file at `64831bd`). The five are
`gamma-phi-legacy`'s, unchanged and by design. `eos-three-phase` refuses
nothing: its `LLL` count goes 8 -> 12 (the four Peng-Robinson three-liquid
states) and its `VLE` count 22 -> 27 (the four water-lean binary states and
the ternary corner).

**Exactly nine states changed, and they are exactly the nine that refused.**
Compared field by field against the committed `74820b8` record on `bucket`,
`verdict`, `phases`, `phase_fractions`, `invariant_violations`,
`converged_stage`, `phase_set_history`, `outcome`, `refusal_class` and
`refusal_stage`: **9 states differ, 2496 are identical**, and every one of the
nine is a `refused -> converged` transition. Dormancy was argued from the code
(both repairs sit on the branch that raised) and this is the measurement of
it.

The two repairs are visible in the record: the five class (i) states carry
`phase_set_history = "L -> LL -> LLV -> LV -> LV"` and
`converged_stage = "second-order"`, the four class (ii) states carry
`converged_stage = "second-order-log"` with `log_space_iterations` of 1, 3, 2
and 2. Wall time 2202.5 s against 2232.8 s at `74820b8` - a 1.4 % difference,
inside this machine's run-to-run noise (`benchmarks/README.md`), and not read
as a speed-up: nine states that used to refuse early now solve.

### (vi) What this case does not establish

- **Not that the model is right at these conditions.** Neither the PC-SAFT nor
  the Peng-Robinson water / ethanol / n-hexane system at `k_ij = 0` is
  validated against a published ternary datum; this repository packages no
  `k_ij` for a water / alcohol or water / alkane cross pair (the same gap
  brain.md's roadmap item 2 records for water / ethanol). The three-liquid
  answers are equilibria **of that model**, nothing more.
- **Not that the removal ranking is now correct**, only that it is no longer
  final. A geometry in which every candidate on the undo stack dead-ends would
  still raise the same message; no state in this repository does.
- **Not a global-minimum proof.** The Gibbs ordering compares the alternatives
  that were *solved* - the single-phase states, every two-phase subset of the
  returned phases, and for the binary the three tie lines its two branches
  admit. A phase set nothing seeded would not appear in it. This is the same
  bound `stability_tp`'s "stable" verdict carries (brain.md section 10).
- **Not a statement about `gamma-phi-legacy`'s five refusals**, which are
  untouched, by design, and pinned as such.

- **Suite time:** `pytest -q` **901 passed, 98 deselected, 267.51 s (4:27)**,
  on the machine named above, against **894 passed, 91 deselected, 267.04 s**
  at `74820b8` (Case R-MAP-2) - +0.5 s, inside this machine's noise. The new
  default-run tests cost 4.1 s of that measured on their own; the slow-marked
  repetitions (three further `T3` offsets, the ternary corner) cost 18.5 s and
  the two new FeOs comparisons 19.2 s, neither of which the default run pays.
  The ~240 s ceiling brain.md section 10 tracks is **still exceeded**, as it
  was before this slice; nothing here narrowed the gap, and roadmap item 1 is
  the only thing on the list that would.
  `pytest -q -m slow tests/test_robustness_map.py` **4 passed, 48 deselected,
  2294.51 s (38:14)** - the full 2505-state sweep reproducing this record's
  totals (2202.5 s of that is the sweep itself), plus the three-phase verdict
  pins.
- **Test path:** `tests/test_flash_eos_multiphase_robustness.py` (nine states,
  the independent Newton, the Gibbs ordering, and the two new rules on their
  own), `tests/validation/test_pcsaft_vlle_water_hexane.py` (the FeOs
  comparisons, `slow`), `tests/test_robustness_map.py` (the quick counts and
  the pins).
- **Script:** `python -m chemthermo.bench robustness --out record.json
  --summary-out record.md`.
