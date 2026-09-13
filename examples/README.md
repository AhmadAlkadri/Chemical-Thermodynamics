# Examples

Run examples from the repo root with `python`.

## Scope Policy

VLLE and PC-SAFT are in scope for Chemical-Thermodynamics. No thermodynamic capability class is categorically out of scope; implementation maturity may vary by module and release.

## Scripts

- `examples/database/01_component_database_demo.py`
  - Inspects canonical packaged runtime DB metadata and sample records.
- `examples/basic/tp_flash_pr_pure.py`
  - Pure-component TP flash with Peng-Robinson EOS.
- `examples/basic/tp_flash_pr_mixture.py`
  - Mixture TP flash with Peng-Robinson EOS.
- `examples/basic/tp_flash_pr_kij_demo.py`
  - Mixture TP flash with a per-pair `kij` mapping (`PengRobinsonEOS(kij={...})`),
    contrasted with the default `kij=0.0` result on the same feed.
- `examples/basic/tp_flash_nrtl_vle.py`
  - Gamma-phi TP flash with NRTL (liquid) + Peng-Robinson (vapor).
- `examples/basic/flash_tp_peng_robinson_demo.py`
  - Existing TP flash demo (Peng-Robinson EOS). Its diagnostics block now also
    prints the tangent-plane phase-detection keys (`phase_detection`,
    `stability_status`, `tpd_min`, `k_seed`, `mass_balance_residual`,
    `fugacity_residual`, `delta_g_split_rt`).
- `examples/basic/flash_tp_auto_phase_demo.py`
  - Automatic 1-vs-2 phase detection in `flash_tp` (ADR-0008): a stable feed, an
    unstable feed whose split is seeded from the tangent-plane minimizer, and
    the Methane / n-Pentane state at 175 K / 1.778 MPa where the legacy Wilson
    heuristic returns a single liquid and the tangent-plane path returns two
    phases. Prints both paths side by side plus the Gibbs-energy evidence.
    Needs no optional dependency; see validation Cases F-1 and F-2.
- `examples/basic/flash_tp_gamma_phi_demo.py`
  - Existing TP flash demo (gamma-phi, NRTL + Peng-Robinson).
  - Uses the packaged NRTL pairs, which are **synthetic illustrative placeholders**
    (not fitted to data, not from any publication); see the `provenance` block in
    `src/chemthermo/parameters/data/activity/nrtl.json`.
- `examples/basic/stability_tp_peng_robinson_demo.py`
  - Michelsen tangent-plane phase stability at two states (unstable and stable).
- `examples/basic/stability_tp_nrtl_lle_demo.py`
  - Liquid-liquid tangent-plane stability with an NRTL activity model
    (`stability_tp(..., activity_model=...)`) on the partially miscible
    n-butanol / water binary: a feed inside the miscibility gap, one outside it,
    and one of the two conjugate liquid phases (marginally stable). Parameters
    are written inline with their citation (Tessier, Brennecke & Stadtherr,
    Chem. Eng. Sci. 55 (2000) 1785, Table 1, pair 2-3) and are **not** the
    packaged synthetic defaults. Needs no optional dependency.
- `examples/basic/flash_tp_nrtl_lle_demo.py`
  - Liquid-liquid TP flash with an NRTL activity model and **no** equation of
    state (`flash_tp(..., activity_model=...)`, mode inferred as
    `gamma-gamma`) on the partially miscible n-butanol / water binary: two feeds
    inside the miscibility gap that return the same tie-line with different
    phase amounts (checked against the lever rule), and one outside it that
    returns a single liquid. Prints the verification residuals and the
    post-split stability block. Parameters are written inline with their
    citation (Tessier, Brennecke & Stadtherr, Chem. Eng. Sci. 55 (2000) 1785,
    Table 1, pair 2-3) and are **not** the packaged synthetic defaults. Needs no
    optional dependency. See validation Cases L-3 and L-4.
- `examples/basic/flash_tp_modified_raoult_demo.py`
  - Low-pressure TP flash with `flash_mode="modified-raoult"`: an NRTL liquid
    with Antoine pure-liquid reference fugacities against an ideal-gas vapor,
    both handed to one tangent-plane test (ADR-0010). Shows 1-propanol / water
    as a subcooled liquid, a vapor-liquid split and a superheated vapor from the
    same call, recomputes `y_i P = x_i gamma_i Psat_i` from the printed numbers,
    and then shows n-butanol / water returning a *liquid-liquid* tie-line from
    the identical call. Prints the Antoine validity window, the candidate labels
    and the post-split block. Parameters are written inline with their citation
    (Tessier, Brennecke & Stadtherr, Chem. Eng. Sci. 55 (2000) 1785, Table 1,
    pairs 1-3 and 2-3) and are **not** the packaged synthetic defaults. Needs no
    optional dependency. See validation Cases R-1 and R-2.
- `examples/validation/00_reference_case.py`
  - Deterministic single-case comparison against `thermo` (optional dependency).
    Unchanged by the tangent-plane phase-detection slice: beta 0.46829044 vs
    thermo 0.46976460, |delta| = 1.474e-03.
- `examples/validation/06_stability_vs_thermo.py`
  - Deterministic stability cross-check against `thermo`'s Michelsen test over 7 states.
- `examples/validation/07_nrtl_tessier_stationary_points.py`
  - Reproduces the published NRTL tangent-plane stationary points and D values of
    Tessier, Brennecke & Stadtherr, Chem. Eng. Sci. 55 (2000) 1785, Table 2
    (n-propanol / n-butanol / water). Needs no optional dependency. Parameters come
    from the cited fixture `tests/fixtures/nrtl/tessier2000_problem1.json`, not from
    the packaged defaults. Two printed D values do not reproduce and are reported as
    `KNOWN-TYPO`; see `.agents/brain/validation-cases.md` Case N-3.
- `examples/validation/08_stability_nrtl_tessier2000.py`
  - Runs `stability_tp` with an NRTL activity model on every published feed of
    Tessier, Brennecke & Stadtherr (2000) Problem 1 (Table 2) and Problem 2
    (Table 5), and compares its verdict, `tpd_min` and minimizing composition
    against an independent damped-Newton refinement of the printed stationary
    points and against the printed five-digit D values. Needs no optional
    dependency. Parameters come from the cited fixtures
    `tests/fixtures/nrtl/tessier2000_problem1.json` and
    `tests/fixtures/nrtl/tessier2000_problem2.json`. One printed D is reported
    as `KNOWN-TYPO`; see `.agents/brain/validation-cases.md` Cases S-6 and S-7.
- `examples/validation/09_lle_tessier2000_tie_lines.py`
  - Runs the liquid-liquid `flash_tp` on every published feed of Tessier,
    Brennecke & Stadtherr (2000) Problem 1 (Table 2) and Problem 2 (Table 5) and
    prints, per feed, the tie-line, the phase fractions, the stage-by-stage
    iteration counts, every verification residual and the post-split stability
    verdict, with PASS/FAIL. The paper publishes stationary points of the
    tangent-plane distance, **not** tie-lines, so each result is compared
    against an independent equal-activity solve written inside the script
    (its own successive substitution plus a damped Newton on the full
    `(x^I, x^II, beta)` system). The stable control feed
    (0.25, 0.25, 0.25, 0.25) must stay one liquid. Needs no optional
    dependency. Parameters come from the cited fixtures. See validation Cases
    L-1 and L-2.
- `examples/validation/10_modified_raoult_water_butanol.py`
  - Water / 1-butanol at 1 atm around its three-phase temperature, with PASS/FAIL
    on every check. Solves the binodal and T3 = 366.2138 K independently inside
    the script (equal activities; then `sum x gamma Psat = P` on one liquid,
    *checked* on the other), prints the shared vapor y = (0.234063, 0.765937),
    and then runs the three negative controls: a stable liquid-liquid split
    2 K below T3, a fully evaporated feed 2 K above it (verified against the
    dew-point equation), and the refusal just below T3, where `flash_tp` raises
    `ConvergenceError` because the converged pair is not a stable phase set.
    Needs no optional dependency. See validation Case R-3.
- `examples/validation/*.py`
  - Optional validation sweeps against `thermo` (requires `pip install -e ".[validation]"`).
  - CSV output is disabled by default; pass `--outdir <dir>` or set `CHEMTHERMO_OUTDIR`.
  - Generated CSV outputs are intentionally not tracked in git.

## CLI quick runs

- Phi-phi TP flash:
  - `chemthermo tp-flash --components Methane,Ethane,Propane --z 0.5,0.3,0.2 --temperature-k 240 --pressure-pa 3000000 --format json`
- Gamma-phi TP flash (NRTL + Peng-Robinson):
  - `chemthermo tp-flash --components Methane,Ethane --z 0.5,0.5 --temperature-k 240 --pressure-pa 3000000 --flash-mode gamma-phi --format json`
- Gamma-phi coverage depends on available NRTL pair data; unsupported pairs return a runtime validation/model error.
- The packaged NRTL pairs are synthetic demo values, not fitted or published parameters. Supply your own via `NRTLParameters.from_pairs(...)` for real work.

## Expected output format

Each script prints:
- Header line describing the model/mode
- Temperature (K) and pressure (Pa)
- Feed composition `z` (mole fractions)
- Phase names and/or vapor fraction (beta)
- Phase compositions (`x` for liquid, `y` for vapor) when present
- Optional K-values (`y/x`) when both liquid and vapor are present
- Diagnostics block (raw keys/values)

All inputs and outputs are SI units. Compositions are mole fractions.
