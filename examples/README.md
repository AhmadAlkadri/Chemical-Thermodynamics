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
