# Examples

Run examples from the repo root with `python`.

## Scope Policy

VLLE and PC-SAFT are in scope for Chemical-Thermodynamics. No thermodynamic capability class is categorically out of scope; implementation maturity may vary by module and release.

## Scripts

- `examples/database/01_component_database_demo.py`
  - Inspects canonical packaged runtime DB metadata and sample records.
- `examples/basic/citation_demo.py`
  - Looks up a property's citation (`chemthermo.cite`) directly and via
    `Component.get_citation`, including the "no citation available" case.
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
- `examples/basic/flash_tp_vlle_demo.py`
  - Three-phase (vapor-liquid-liquid) TP flash, **discovered not assumed**.
    1-propanol / n-butanol / water at 364 K and 1 atm with
    `flash_mode="modified-raoult"`: two feeds inside the tie-triangle return the
    same three phases in different amounts, and a third feed outside it returns
    two, from the identical call. Prints the phase-set history
    (`L -> LV -> LLV`), the stage counts, the verification residuals,
    `delta_g_vs_two_phase_rt` and the post-split block. Parameters are written
    inline with their citation (Tessier, Brennecke & Stadtherr, Chem. Eng. Sci.
    55 (2000) 1785, Table 1) and are **not** the packaged synthetic defaults.
    Needs no optional dependency. See validation Cases V-1 and V-2.
- `examples/basic/pcsaft_properties_demo.py`
  - PC-SAFT (Gross & Sadowski 2001, **non-associating**) residual properties of
    a methane / n-hexane mixture at two densities: `A^res/(R T)`, `Z`, `P` and
    `ln phi_i`, plus the identity `sum_i x_i ln phi_i = A^res/RT + Z - 1 - ln Z`
    recomputed from the printed numbers. Also shows a state *inside* the
    spinodal, where `Z < 0` and the model refuses to produce fugacity
    coefficients rather than returning a `nan`, and the effect of a nonzero
    per-pair `kij`. Parameters come from the packaged
    `src/chemthermo/parameters/data/eos/pcsaft.json` (Gross & Sadowski 2001,
    Table 1; see its `provenance` block). **No density root solving happens
    anywhere**: the state is given as `(T, molar density, x)`. Needs no optional
    dependency. See ADR-0014 and validation Case P-0.
- `examples/basic/flash_tp_pcsaft_demo.py`
  - PC-SAFT in `stability_tp` / `flash_tp` (ADR-0015): methane / n-hexane at
    300 K. Prints the density roots of the feed at three pressures and of pure
    n-hexane at its saturation pressure (two roots, equal fugacity on them, to
    5.6e-12), a two-phase flash at 3 MPa with both phase densities and the
    three verification residuals, the same feed at 8 MPa where the
    tangent-plane test finds one stable phase, and the Peng-Robinson answer
    beside it **for contrast only** - two models, two answers, neither line
    evidence about the other. Needs no optional dependency. See validation
    Cases P-3, P-4 and P-5.
- `examples/basic/pcsaft_association_demo.py`
  - PC-SAFT **with association** (Gross & Sadowski 2002, ADR-0018). Pure water
    term by term - hard chain, dispersion and association side by side, with
    the non-bonded site fraction `X` that the association term is a function of
    (0.9851 in the saturated vapour, 0.0356 in a compressed liquid) - then
    water's saturation state at 373.15 K by equal fugacity on the two density
    roots (100,890.27 Pa, 0.43 % below the 101,325 Pa that *defines* the normal
    boiling point, printed as a remark about the model and asserted nowhere).
    Then a water / n-hexane liquid-liquid split from `flash_tp`, at 1 atm and
    again at 1 MPa - the same tie line, both times named `"liquid1"` /
    `"liquid2"` with `vapor_fraction = None`. The 1 atm state used to raise
    (the split could only pair a vapour-root phase with a liquid-root one) and
    is what ADR-0019 fixed; the script says so where the old caveat used to be.
    `kij = 0`, which is a poor model for water with a hydrocarbon and is said
    so on the page. **Default stops before the 1 MPa section** (about 5 s);
    `--full` runs it. Needs no optional dependency. See ADR-0018, ADR-0019 and
    validation Cases P-6, P-7 and P-8.
- `examples/basic/pcsaft_polymer_demo.py`
  - **A polymer as a PC-SAFT component** (ADR-0022), the golden path for that
    slice. Polyethylene / n-pentane at 453 K: the segments-per-mass record
    (`segments_per_g = 0.0263` mol/g times `Mw = 16400` g/mol gives
    `m = 431.32`, against 2.6896 for the solvent - 160 : 1), a custom
    non-volatile `Component` with **no critical constants at all**, the pure
    melt's single density root over 1-30 MPa, the stability verdict switching
    once with pressure and once with temperature (the LCST-type direction: this
    system demixes on heating), the verified liquid-liquid split at 8 MPa with
    its three residuals, and the effect of `k_ij`. **The polymer parameters are
    a cited test fixture, not packaged data**: as tabulated by Martini et al.
    (2009) citing a paywalled Gross & Sadowski (2002) table that was not read,
    with no second open source found; the polymer is modelled as monodisperse;
    nothing is compared against measurement, and the script says all of this
    before printing anything. About 5 s by default; `--full` bisects the cloud
    point and runs the `Mw = 53000` chain, where `exp(ln phi)` underflows and
    ADR-0022's log-space route carries the flash. Needs no optional dependency.
    See validation Cases P-12 and P-13.
- `examples/basic/flash_tp_pcsaft_lle_demo.py`
  - **Liquid-liquid equilibrium from an equation of state** (ADR-0019), the
    golden path for that slice. Water / n-hexane at 298.15 K: the stability
    verdict and the two density roots that made the state hard, the 1 atm split
    (`liquid1` / `liquid2`, `vapor_fraction = None`, every verification residual
    and a `kappa` recomputed from the public `pressure_Pa` so the liquid
    identities do not rest on `phase_identity` alone), the same tie line from
    three feeds with the lever rule, and the 1 MPa tie line that used to be
    mislabelled `"liquid"` / `"vapor"` by the Wilson-ranking fallback. `kij = 0`
    and the mutual-solubility caveat are printed. **Default is the 1 atm split
    only** (about 5 s); `--full` adds the three feeds and the 1 MPa tie line,
    both of which run in `tests/test_flash_eos_lle.py`. Needs no
    optional dependency. See validation Case P-8.
- `examples/basic/flash_tp_pcsaft_vlle_demo.py`
  - **Three phases from an equation of state** (ADR-0020), the golden path for
    that slice. Three sections: three *liquid* phases from **Peng-Robinson**
    with `kij = 0` (water / ethanol / n-hexane at 280 K, in about a tenth of a
    second - the pure-cubic three-phase state ADR-0019 recorded as not found);
    the water / n-hexane window at 1 atm, where a 4-equation Newton written in
    the script locates the three-phase temperature `T3` and shows that every one
    of the three coexisting compositions has *both* a vapour and a liquid
    density root, followed by `flash_tp` on either side of `T3` - two conjugate
    liquids below it via `V -> LV -> LLV -> LL` (those temperatures used to
    raise), the ordinary vapour-liquid answer above it; and, behind `--full`, a
    PC-SAFT vapour-liquid-liquid tie triangle. Prints the search route, the
    phase amounts and every verification residual. About 8 s by default.
    `kij = 0` and the model caveat are printed, and nothing is compared against
    measurement. Needs no optional dependency. See validation Cases P-9, P-10.
- `examples/validation/00_reference_case.py`
  - Deterministic single-case comparison against `thermo` (optional dependency).
    Unchanged by the tangent-plane phase-detection slice: beta 0.46829044 vs
    thermo 0.46976460, |delta| = 1.474e-03.
- `examples/validation/01_pure_props_grid.py`
  - Peng-Robinson pure-component `Psat(T)` and saturated `Z` against `thermo`
    for methane, 110-180 K; requires `thermo` (exits with a message if absent).
    Optionally writes a CSV to `--outdir` / `CHEMTHERMO_OUTDIR`.
- `examples/validation/02_binary_k_surface.py`
  - Binary methane/ethane TP-flash K-value surface against `thermo` over a
    composition sweep at fixed T, P; requires `thermo`.
- `examples/validation/03_txy_pxy_curves.py`
  - Methane/ethane bubble/dew temperature curve (Txy at fixed P) against
    `thermo`, with RMS temperature differences; requires `thermo`.
- `examples/validation/04_multicomponent_flash_sweep.py`
  - Methane/ethane/propane TP flash over random Dirichlet-sampled
    compositions against `thermo`: phase-count agreement, vapor-fraction and
    composition differences; requires `thermo`.
- `examples/validation/05_robustness_map.py`
  - Convergence and phase-count agreement between chemthermo and `thermo`
    over a 10x10 methane/ethane T-P grid, classified into
    MATCH/PHASE_MISMATCH/CT_FAIL/TH_FAIL/BOTH_FAIL buckets; requires `thermo`.
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
    dew-point equation), and the window just below T3, where the first
    two-phase iterate is the wrong one: `FlashSettings(max_phases=2)` still
    raises there, while the default `max_phases=3` resolves it to the two
    liquids by adding a phase and removing another
    (`phase_set_history = "V -> LV -> LLV -> LL"`).
    Needs no optional dependency. See validation Cases R-3 and V-3.
- `examples/validation/11_vlle_water_propanol_butanol.py`
  - The ternary vapor-liquid-liquid **tie-triangle** of
    1-propanol / n-butanol / water at 1 atm, with PASS/FAIL on every check.
    Solves the triangle at 365, 364 and 363 K from its own six-equation Newton
    iteration (three equal activities, two normalizations, and the bubble
    condition on **one** liquid - that the other liquid also boils, and that
    both share one vapor, are then checked as consequences). Six feeds inside
    the triangle are reproduced to |dx| <= 2.2e-14 and |dbeta| <= 1.6e-13 with
    `G(3 phases) < G(2-phase candidate) < G(feed)` computed in the script; four
    feeds outside it (vapor-liquid region, liquid-liquid region, water-rich
    corner, superheated) are each verified by their own route; and the binary
    refusal window of Case R-3 is shown resolved. Needs no optional dependency.
    See validation Cases V-1, V-2 and V-3. The 363 K feed this script used to
    record as a stability miss is fixed by ADR-0012; it is now the closing
    section of `12_vlle_verdict_map.py`.
- `examples/validation/12_vlle_verdict_map.py`
  - The **verdict map** of the same ternary: how many phases the model has at
    each of 75-76 feeds per temperature, at 363, 364 and 365 K, with PASS/FAIL
    and a printed confusion matrix per temperature. Every feed is classified
    independently inside the script by building *every* state the model admits
    - the tie-triangle, a four-equation vapor-liquid Newton, a seven-equation
    liquid-liquid Newton, and the single phase - and taking the one of lowest
    Gibbs energy, so "inside the tie-triangle means three phases" is checked
    rather than assumed. Zero disagreements at all three temperatures; every
    three-phase answer is the same triangle to |dx| <= 9.5e-12 with phase
    fractions equal to the feed's barycentric weights to 1.6e-11. Ends with the
    feed that validation Case V-2 recorded as a miss, showing the stability
    trials and the surface each ran on. Needs no optional dependency. Runs in
    about 16 s. See validation Case V-5 and ADR-0012.
- `examples/validation/13_pcsaft_vs_teqp.py`
  - PC-SAFT against **teqp** (NIST, MIT), with PASS/FAIL per state. teqp
    implements the same published model but takes every derivative by automatic
    differentiation, where chemthermo writes them analytically - so agreement on
    `Z` and `ln phi` tests the hand-written derivatives, not just the model.
    Compares `A^res/RT`, `Z`, `P` and `ln phi_i` at 14 states (pure, binary,
    ternary, gas-like and liquid-like densities, a nonzero `k_ij`, and one state
    inside the spinodal), runs a negative control (a 1 % change in `sigma` must
    break the agreement), and reproduces teqp's `pure_VLE_T` saturation for
    n-hexane at 300 K and 400 K with a bisection density root finder written
    inside the script. Requires `pip install -e ".[validation]"`; prints a
    message and exits 0 without teqp. See validation Cases P-1 and P-2.
- `examples/validation/14_pcsaft_flash_vs_teqp.py`
  - PC-SAFT **phase equilibrium** against teqp, with PASS/FAIL per check.
    teqp traces its own 300 K methane / n-hexane isotherm and polishes a tie
    line at each of seven pressures; chemthermo reaches the same tie lines
    through the tangent-plane test and the Rachford-Rice split, agreeing to
    `|dx1| <= 2.0e-9` and `|dy1| <= 3.7e-12` with phase densities to 1.3e-9
    relative. The decisive check needs no reference tie line at all: teqp's own
    `get_fugacity_coefficients`, evaluated at chemthermo's converged
    compositions and densities, must make the two phases equal in fugacity
    (achieved 3.8e-9 relative, worst case). Also: pure n-hexane's two
    saturation densities at 300 K and 400 K; bubble pressures found by
    bisecting `stability_tp`'s *verdict* against teqp's `mix_VLE_Tx` (1.4e-8
    relative); and one methane / n-decane case with an **illustrative**
    `kij = 0.03` (not a literature-validated parameter). Requires
    `pip install -e ".[validation]"`; prints a message and exits 0 without
    teqp. Runs in about 25 s. See validation Cases P-3 and P-5.
- `examples/validation/15_flash_split_robustness.py`
  - Phi-phi **split robustness** with PC-SAFT, with PASS/FAIL per check
    (ADR-0016, validation Case F-4). Successive substitution oscillates on some
    feeds - the K-values cross 1 back and forth, the implied vapor fraction
    leaves `[0, 1]` and then no vapor fraction exists at all - and four states
    of the grid below used to end in
    `ConvergenceError("Rachford-Rice failed to bracket a vapor fraction.")`
    despite the tangent-plane test having proved them two-phase. The script
    checks the reference state (CO2 / n-decane, `z = (0.9, 0.1)`, 240 K,
    1.0 MPa) against a damped Newton on the equal-fugacity system **written
    inside the script**, agreeing to `|d beta| = 7.8e-13` and
    `|d composition| = 1.2e-13`; scans the 188-state Case F-4 grid and reports
    the `ConvergenceError` count (0, against 4 before ADR-0016) together with
    the worst invariant on every two-phase answer; and shows that the four
    rescued states still fail with `second_order=False` and on the legacy
    `wilson-heuristic` path, which is deliberate. Needs no optional
    dependency. Runs in about 100 s. See validation Case F-4.
- `examples/validation/16_pcsaft_association_vs_feos.py`
  - PC-SAFT **association** against **FeOs** (feos-org/feos, MIT OR
    Apache-2.0), with PASS/FAIL per check (ADR-0018, validation Cases P-6 and
    P-7). FeOs is used rather than teqp because teqp's `PCSAFT` kind has no
    association term at all; like teqp it gets every derivative by automatic
    differentiation, so no derivative code is shared. Checks eighteen states
    term by term (hard chain, dispersion, association, `A^res/RT`, `Z`,
    `ln phi_i`); settles the `sigma^3`-versus-`d^3` question in the association
    strength numerically, since the 2002 paper is paywalled and was not read;
    re-asserts that non-associating n-hexane is bit-identical to Case P-1's
    pinned values; and runs the equilibrium exam - pure-water saturation
    against FeOs's own `PhaseEquilibrium.pure`, a water/ethanol VLE flash and a
    water/n-hexane liquid-liquid split with FeOs's fugacities evaluated at
    chemthermo's converged phases. Every dispersion-dependent number is
    reported twice, because FeOs hard-codes the 2001 paper's universal
    constants to fourteen figures where the paper prints ten. Requires
    `pip install -e ".[validation]"`; prints a message and exits 0 without
    `feos`. Runs in about 10 s.
- `examples/validation/17_pcsaft_lle_vs_feos.py`
  - The **liquid-liquid tie line** against FeOs, with PASS/FAIL per check
    (ADR-0019, validation Case P-8). Water / n-hexane at 298.15 K, at 1 atm and
    at 1 MPa: the tie line, the phase densities and the phase amounts against
    FeOs's **own** `State.tp_flash` (two solvers, not just two models); FeOs's
    chemical potentials evaluated at chemthermo's phases, which does not depend
    on FeOs's flash converging; the lever rule across three feeds; and a
    negative control (1 % on water's `eps^AB` must move the tie line). Reports
    every FeOs-at-chemthermo's-densities number twice, as shipped and with
    FeOs's fourteen-figure universal constants, and prints the one feed
    (`z = 0.2/0.8`) where FeOs's own flash raises rather than hiding it.
    Requires `pip install -e ".[validation]"`; prints a message and exits 0
    without `feos`. **Default is the 1 atm state only** (about 12 s); `--full`
    adds the 1 MPa tie line, the three feeds and the negative control.
- `examples/validation/18_pcsaft_vlle_water_hexane.py`
  - The **three-phase neighbourhood of water / n-hexane**, with PASS/FAIL per
    check (ADR-0020, validation Cases P-9 and P-10). Four independent routes:
    a 4-equation Newton written in the script locates the three-phase
    temperature `T3` and the three coexisting compositions from
    `fugacity_coefficients` alone; two-equation Newtons give the
    liquid-liquid and vapour-liquid tie lines on either side of `T3` and their
    reduced Gibbs energies say which pair is the equilibrium; FeOs supplies
    chemical potentials at chemthermo's converged phases and densities; and
    FeOs's **own** two-phase flash is run and shown to converge on the
    *metastable* vapour-liquid pair below `T3` - the pair chemthermo's
    post-split stability test refuses. Every FeOs-at-chemthermo's-densities
    number is reported twice, as shipped and with FeOs's fourteen-figure
    universal constants. Routes 1 and 2 run without `feos`; the script says so
    and still exits 0. About 15 s by default; `--full` adds the 41-point scan
    at two feeds and the ternary vapour-liquid-liquid tie triangle (several
    minutes).
- `examples/validation/19_eos_stability_surfaces.py`
  - **Fixed density-root surfaces in the equation-of-state stability trials**
    (ADR-0021, validation Case P-11), with PASS/FAIL per check. A 4-equation
    Newton written in the script locates `T3`; above it, a water-rich
    water / n-hexane feed now returns a vapour-liquid pair where it used to
    return two liquids, and a two-equation Newton plus reduced Gibbs energies
    say which of the two is the equilibrium; the per-trial table shows which
    root surface each trial ran on and which one found the vapour stationary
    point that used to be missed; FeOs supplies chemical potentials at
    chemthermo's converged phases. Runs without `feos` (it says so and still
    exits 0). About 5 s by default; `--full` adds `T3 + 0.05 K` and
    `T3 + 0.5 K`, the matched-universal-constants comparison asserted at 1e-08,
    the two 41-point scans with the verdict boundary bisected to 1e-06 K, and
    the per-surface trial statistics over the 144-state Peng-Robinson grid
    (several minutes).
- `examples/validation/20_pcsaft_polymer_vs_feos.py`
  - **Polymer/solvent PC-SAFT against FeOs** (ADR-0022, validation Cases P-12
    and P-13), with PASS/FAIL per check. Four routes: FeOs's own `A^res/RT`,
    `Z` and `ln phi` at eleven states on chemthermo's own density roots (agreed
    to 5.0e-12 with matched universal constants, 1.4e-06 as shipped); an
    equal-fugacity Newton written in the script that reproduces the 8 MPa tie
    line to 1.7e-15; FeOs's chemical potentials at chemthermo's converged
    phases (4.5e-13 matched, 4.1e-08 as shipped); and FeOs's own `tp_flash` at
    the one pressure where it converges on this system. Two facts about the
    reference are printed rather than hidden - FeOs's flash **raises** here at
    5 and 8 MPa and returns a degenerate pair at 3 MPa, and `k_ij` cannot be
    given to FeOs's PC-SAFT in feos 0.10.1, so every FeOs comparison runs at
    `k_ij = 0` on both sides. The Newton route runs without `feos` (the script
    says so and still exits 0). About 4 s by default; `--full` adds the
    FeOs-flash survey across 3-15 MPa, the `Mw = 53000` chain and the bisected
    cloud point (about 20 s).
- `examples/validation/21_pcsaft_polymer_vle.py`
  - **Polymer/solvent vapour-liquid equilibrium** (ADR-0024, validation Case
    P-14), with PASS/FAIL per check. Below n-pentane's saturation pressure at
    453 K the polyethylene / n-pentane equilibrium is a solvent vapour over a
    solvent-swollen melt, the vapour's polymer mole fraction is `exp(-450)`,
    and `flash_tp` reaches it only by seeding and finishing the split in log
    mole numbers. Five routes: the split itself at 0.5, 1 and 2 MPa with both
    phases named from a measured compressibility and every residual reported;
    a one-dimensional equal-fugacity solve written in the script, with the
    vapour taken as *exactly* pure solvent, which reproduces the melt's solvent
    content to 1.2e-15 - 2.5e-14 absolute; the polymer's own equal-fugacity
    condition checked in logarithms (`ln f` of about -500, closing to 7.4e-13);
    FeOs's chemical potentials at chemthermo's converged phases (2.6e-12 with
    matched universal constants, 2.5e-07 as shipped, at `k_ij = 0` because
    feos 0.10.1 cannot be given one); and the ternary with n-hexane at 3 MPa,
    the second state ADR-0022 pinned as a runaway. Routes 1, 2, 3 and 5 run
    without `feos` (the script says so and still exits 0). About 5 s by
    default; `--full` adds the 25-point pressure scan from 0.3 to 12 MPa - the
    VLE -> LLE -> single-liquid verdict sequence - and the bisected
    vapour-liquid / liquid-liquid boundary (about 30 s).
- `examples/validation/22_stability_log_space.py`
  - **Tangent-plane stability in log mole numbers** (ADR-0025, validation Case
    P-15), with PASS/FAIL per check, and the golden path for that slice. Before
    it, Michelsen's iteration clamped `ln W` to `[-700, 700]`, so a stationary
    point outside that window was not merely inaccurate but *unreachable*: a
    `Mw = 53000` polyethylene melt against a solvent-vapour feed sits at
    `ln W_polymer = 1452`, and the three trials walking at it spent 51
    iterations parked on the clamp with a residual of exactly `1452 - 700`.
    Five routes: the stationary point against Michelsen's own equations (5) and
    (7) re-derived in the script, including that `sum_W` has overflowed to
    `inf` while `ln sum_W` has not; the split it seeds at 0.5 and 1 MPa, with
    both phases named from a measured compressibility and every residual
    reported; a one-dimensional equal-fugacity solve written in the script,
    which reproduces the melt's solvent content to 1.2e-14; FeOs's chemical
    potentials at chemthermo's converged phases (1.1e-12 with matched universal
    constants, 8.6e-08 as shipped, at `k_ij = 0` because feos 0.10.1 cannot be
    given one); and **dormancy**, which is the whole bit-identity argument
    measured rather than argued - the shorter `Mw = 16400` chain, a
    Peng-Robinson state, a PC-SAFT state and an NRTL state all run the
    pre-ADR-0025 arithmetic. Routes 1, 2, 3 and 5 run without `feos` (the
    script says so and still exits 0). About 5 s by default; `--full` adds the
    144-state Peng-Robinson stability grid (0 of 624 trials in log space) and
    the 0.4-2 MPa scan, every point of which either raised or was seeded from
    the wrong stationary point before this slice (about 10 s).
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
