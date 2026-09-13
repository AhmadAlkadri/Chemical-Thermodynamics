# ADR-0017: Phase identity by compressibility, not by tie-break or Wilson ranking

Status: accepted
Date: 2026-09-13

## Context

ADR-0008 decision 3 named the two phi-phi phase labels by convention, because
`EquationOfState` exposed no molar volume at the time:

- **Single-phase results.** `_EOSTangentPlane` holds two candidates, the
  `"vapor"` (largest `Z`) and `"liquid"` (smallest `Z`) compressibility
  branches, and `_select_min_gibbs` keeps whichever minimizes
  `sum_i w_i ln phi_i(w)` - the correct rule when the branches are genuinely
  different roots. When the cubic (or, since ADR-0015, PC-SAFT's density-root
  solver) has a **single** real root, both branches return identical
  fugacity coefficients, the two Gibbs sums tie exactly, and
  `g_res < best_g` (strict) keeps the first candidate tried - `"vapor"`,
  always, by construction order, never by anything physical. A dense
  compressed liquid or a supercritical fluid is exactly this case.
- **Two-phase results.** The successive-substitution loop is symmetric under
  swapping the two phases (`x <-> y`, `beta <-> 1 - beta`, `K <-> 1/K`), so
  something has to break that symmetry before iterating. ADR-0008 decision 3
  chose the Wilson correlation's volatility ranking of the feed components:
  the incipient phase is called vapor-like when it is enriched (relative to
  the feed) in the component with the largest Wilson K relative to the
  smallest. That is a real thermodynamic correlation, but it is built from
  `Tc`, `Pc` and `omega` alone - it knows nothing about the density the split
  actually converged to, and it can be wrong exactly where two components'
  volatilities are close or the split is weakly non-ideal.

Both are documented, in ADR-0008 and in `brain.md`'s convention note, as
*conventions the model cannot avoid* given the `EquationOfState` interface's
one calculation (`fugacity_coefficients`) - not measurement. That was true in
2026-09-13's earlier slices. It stopped being true once PC-SAFT's density-root
solver (ADR-0015) gave the package a genuine state variable to measure from:
`density_roots` at a state is exactly `1/V` on the roots
`fugacity_coefficients` already selects between, and the isotherm module
(`chemthermo.eos._pcsaft_density`) already computes `dP/drho` analytically as
part of finding those roots. Peng-Robinson has the same information latent in
its cubic; it was simply never extracted.

Motivating measured cases at commit `32f693b` (before this slice):

- PC-SAFT, carbon dioxide / n-decane, `z = (0.9, 0.1)`, 230 K, 2.5 MPa: one
  density root (18,676.8 mol/m^3, a dense liquid), reported `"vapor"` with
  `vapor_fraction = 1.0`.
- PC-SAFT, methane / n-hexane, `z = (0.5, 0.5)`, 170 K, 2.0 MPa: one density
  root (12,871.6 mol/m^3, also dense), reported `"vapor"` with
  `vapor_fraction = 1.0`.

Both are liquids by every ordinary sense of the word - a few times the
density of the liquid phase of either binary's own two-phase states on the
same grid - and both were reported `"vapor"` solely because `_select_min_gibbs`
tried the `"vapor"` candidate first.

## Decision

### 1. `kappa = P / (rho * dP/drho)`, evaluated at the root a phase actually converged on

The dimensionless isothermal-compressibility ratio

    kappa = (1 / rho) * (d rho / dP)_T * P  =  -P / (V * (dP/dV)_T)  =  P / (rho * (dP/drho)_T)

is exactly 1 for an ideal gas (`P = rho R T` gives `dP/drho = R T = P/rho`)
and, for a liquid near its own close-packed density, far below 1 - a liquid's
pressure changes enormously for a tiny density change, so `dP/drho` is large
and `kappa` is small. It requires only a first density derivative at fixed
temperature, which both packaged models already have on the branch/root a
phase is evaluated on:

- **Peng-Robinson**: `dP/dV` is the exact derivative of
  `P = RT/(V - b) - a / (V^2 + 2 b V - b^2)`, the same cubic
  `compressibility_factor` already solves for `Z`. No new root-finding.
- **PC-SAFT**: `dP/drho` is the analytic derivative
  `(dP/drho)/RT = 1 + 2 eta a'(eta) + eta^2 a''(eta)` that
  `chemthermo.eos._pcsaft_density.PCSAFTIsotherm.pressure_and_slope` already
  computes for the Newton refinement and the mechanical-stability filter of
  the density-root solver (ADR-0015). No new derivative, no finite
  difference.

Neither model needed a temperature derivative, which is the one thing PC-SAFT
genuinely does not have (see "Alternatives considered").

### 2. The threshold: `KAPPA_LIQUID_THRESHOLD = 0.5`

`kappa < 0.5` names a root `"liquid"`; `kappa >= 0.5` names it `"vapor"`.
Measured at commit `32f693b` plus this slice's code (every number below from
`kappa`, not `Z` or density; PC-SAFT roots via a finite difference of the
public `pressure_Pa(T, rho, x)`, PR roots via the analytic derivative
cross-checked against an independent finite difference built from
`Component.tc_k/pc_pa/omega` alone - see `tests/test_phase_identity.py`):

| grid | class | n | kappa min | kappa max |
|---|---|---:|---:|---:|
| Peng-Robinson, 144-state (`tests/test_flash_phase_detection.py` grid) | two-phase liquid root | 47 | 1.13e-04 | 2.30e-01 |
| Peng-Robinson, 144-state | two-phase vapor root | 47 | 1.01e+00 | 1.48e+00 |
| Peng-Robinson, 144-state | single-phase liquid root | 52 | 1.35e-04 | 1.95e-01 |
| Peng-Robinson, 144-state | single-phase vapor root | 45 | 7.44e-01 | 1.76e+00 |
| PC-SAFT, 188-state (Case F-4 grid, `tests/validation/test_flash_split_robustness_pcsaft.py`) | two-phase liquid root | 123 | 4.88e-04 | 8.50e-03 |
| PC-SAFT, 188-state | two-phase vapor root | 123 | 1.04e+00 | 1.87e+00 |
| PC-SAFT, 188-state | single-phase liquid root | 65 | 1.17e-03 | 3.20e-02 |
| PC-SAFT, 188-state | single-phase vapor root | 0 | - | - |

Over both grids combined: every liquid-root kappa is `<= 0.230`; every
vapor-root kappa is `>= 0.744`. `0.5` sits in the middle of that gap with a
`>= 0.24` margin on each side - not a boundary value chosen to make one state
come out right. (The PC-SAFT grid has no single-phase vapor state at all: its
feed/temperature/pressure ranges - CO2/n-decane and methane/n-hexane, 170-260 K
- never land on a dilute-gas single root; that is a property of the grid, not
of the criterion, and the two-phase vapor roots on the same grid are still
correctly `>= 1.04`.)

**Near-critical states are where the criterion, like any label, stops being
more than a convention**, and this slice does not claim otherwise: a state
whose two converged phases both land on the same side of `0.5` (both
liquid-like or both vapor-like - only ever seen on marginal, weakly unstable
splits) falls back to the historical Wilson-ranking orientation rather than
trusting a `kappa` comparison with no real separation to make. No state in
either validated grid needed that fallback for a *two-phase* result (measured:
`phase_label_method == "compressibility"` on all 47 PR and all 123 PC-SAFT
two-phase states); it is there for the near-critical states this package
cannot yet exhibit on these grids, not a fiction.

### 3. `EquationOfState.phase_identity`, a new method with a documented default

```python
def phase_identity(
    self, *, mixture, temperature_K, pressure_Pa, composition, phase,
) -> str | None:
```

added to `chemthermo.models.base.EquationOfState` as a **concrete, non-abstract**
method whose default implementation returns `None`. `None` is the documented
signal "this model cannot measure an identity here"; every caller must keep
its pre-ADR-0017 behavior when it sees `None` (the fallback is exercised by
`tests/test_phase_identity.py::test_a_model_without_phase_identity_keeps_the_pre_adr_0017_behavior`,
a minimal `EquationOfState` subclass that does not override it). Both
packaged models override it:

- `PengRobinsonEOS.phase_identity`: solves the same cubic
  `compressibility_factor` does, computes `V`, evaluates the analytic
  `dP/dV` above, and returns the threshold verdict. Raises `ModelError` if
  `dP/dV >= 0` (the spinodal branch, which the label selection never reaches
  in practice since only `min`/`max` roots are chosen) rather than reporting
  a meaningless identity.
- `PCSAFTEOS.phase_identity`: reuses `_root_for_phase` (the same root
  `fugacity_coefficients` and `density_roots` use) and
  `PCSAFTIsotherm.pressure_and_slope` for `(P, dP/drho)` at that root's
  `eta`, then the same threshold. Raises `ModelError` under the same
  condition.

`phase` selects which root to identify with the same semantics
`fugacity_coefficients` uses; when the model has one root, `phase="liquid"`
and `phase="vapor"` return the same identity, which is precisely what settles
the case ADR-0008 decision 3 could not.

### 4. Where the label is applied

**Single-phase / `feed_branch` / `phase_branch` (`stability_tp`, both
directly and through `flash_tp`).** `_EOSTangentPlane.identity_label(composition,
label)` calls `eos.phase_identity(..., phase=label)` on the label
`_select_min_gibbs` (or `_select_surface`) already chose, and returns the
measured identity when it is `"liquid"` or `"vapor"`, else the original
label unchanged. It is applied exactly twice per `stability_tp` call, not
inside the successive-substitution loop: once to the feed's own label
(`stability/tp.py`'s `stability_tp`, right after
`evaluator.ln_fugacity_terms(z)`) and once per trial, inside
`_reported_terms`, which already runs only at a trial's convergence (or its
final second-order-stage evaluation), never per iterate. `_ActivityTangentPlane`
and `_ModifiedRaoultTangentPlane` implement `identity_label` as a pass-through
(`return label`): an activity liquid and an ideal gas are two different
models, not two branches of one equation of state competing for the same
root, so there is no tie-break to replace there and this slice does not touch
their labels.

This is a pure relabeling: `ln_fugacity_terms`, `_select_min_gibbs` and the
successive-substitution/Newton iteration are untouched, so no `tpd`, no
`sum_W`, no verdict and no iteration count can move. Only the *string*
`feed_branch` / `phase_branch` can change value.

**Two-phase orientation (`flash_tp`'s phi-phi path only,
`_flash_tp_tangent_plane` in `chemthermo.flash._detect`).** The converged
split already has `x` always evaluated on the model's `"liquid"` branch and
`y` always on its `"vapor"` branch, by construction of
`_split._ln_phi_function` (`phase="liquid"` for `x`, `phase="vapor"` for
`y`), regardless of which orientation the Wilson-ranking seed
(`_incipient_is_vapor_like`, unchanged) used to find that split. After
convergence, `_orient_two_phase_labels` asks `eos.phase_identity(...)` for
`x`'s and `y`'s own root:

- Both identities present and different: the lower-`kappa` one is named
  `"liquid"`, the other `"vapor"` - which may swap the historical
  `x = "liquid"`, `y = "vapor"` assignment, and correspondingly reports
  `vapor_fraction = 1 - beta` instead of `beta` (the two are still the
  *same pair* of numbers; only which one is now called `"vapor"` changed).
- Either identity missing (model does not implement `phase_identity`), or
  both on the same side of the threshold: the historical assignment is kept
  (`x = "liquid"`, `y = "vapor"`, `vapor_fraction = beta`).

`diagnostics["phase_label_method"]` records `"compressibility"` or
`"wilson-ranking"` for a two-phase phi-phi result, and `"compressibility"` or
`"tie-break"` for a single-phase one. **`x`, `y` and `beta` - the converged
compositions and the split fraction the solver actually found - are never
read or written by this function; it only decides which name and which
`vapor_fraction` value go with an already-fixed pair.** This is why the
bit-identity policy below can promise unchanged compositions and fraction
*sets*.

### 5. `chemthermo.__all__` and `cli_schema_version`

`KAPPA_LIQUID_THRESHOLD` is a plain module constant on
`chemthermo.models.base`, not exported from `chemthermo` - it is
implementation detail of `phase_identity`'s two implementations, not a
user-facing tuning knob (matching how `ETA_UNIFORM_STEP` in
`_pcsaft_density.py` is not exported either). `phase_identity` is a new
method on the already-public `EquationOfState`, so no new top-level name is
added. `diagnostics["phase_label_method"]` is a new key inside the existing
free-form `diagnostics` mapping, which ADR-0008 decision 7 already settled
does not need a `cli_schema_version` bump (the mapping's keys are documented
implementation details; nothing existing changed type). `solver.algorithm`
in the CLI output is unchanged.

## Bit-identity policy and what actually changed

`tests/test_flash_refactor_bit_identity.py` pins the whole observable surface
of `flash_tp` (phase names, every phase's composition, phase fractions,
`vapor_fraction`, the full `diagnostics` mapping) for 155 states, bit for bit.
This slice legitimately changes phase *names* - by design, that is the point
- so the fixture was regenerated (`refactor_bit_identity_v2.json`) only after
auditing all 155 states individually against the old fixture
(`refactor_bit_identity_v1.json`, kept in the repository for history). The
audit (script logic now embedded as commentary in the module docstring;
see also validation Case F-5) found:

- **46 of the 144 Peng-Robinson phi-phi grid states changed**, every one a
  single-phase result whose `feed_branch` was the pre-ADR-0017 vapor-first
  tie-break `"vapor"` and is now the measured `"liquid"`, with
  `vapor_fraction` `1.0 -> 0.0` and a new `diagnostics["phase_label_method"]
  = "compressibility"`. Every changed state's kappa is well below the
  threshold (max `0.124`; see Case F-5 for the full 46-row table with every
  mixture, temperature, pressure and kappa).
- **Zero two-phase states changed** on this grid: `phase_label_method` is
  `"compressibility"` on all 47, and `_orient_two_phase_labels` never swapped
  - kappa already agreed with the historical Wilson-ranking orientation
  everywhere on this grid.
- **Zero composition changes, zero fraction-*set* changes, zero other
  diagnostics-number changes**, verified state by state (excluding
  `feed_branch`/`incipient_phase`, which are allowed to carry a different
  *string* and were checked separately above, and the new
  `phase_label_method` key). `diagnostics["phase_state"]` mirrors the phase
  name by construction and changed on the same 46 states, nowhere else.
- The other 11 states (5 gamma-gamma binary feeds, 1 single-component
  gamma-gamma feed, 1 Tessier near-plait feed, 2 gamma-phi cases, 2 legacy
  `wilson-heuristic` cases) are **untouched**: gamma-gamma and the
  Tessier feed never build an `_EOSTangentPlane` at all (no `eos=`);
  gamma-phi stays on the wholly separate legacy path (ADR-0008 decision 6);
  the legacy `wilson-heuristic` phi-phi cases never call `stability_tp`
  (ADR-0008 decision 1) and therefore never reach `identity_label`.

The CLI fixture (`tests/fixtures/cli/tp_flash_v1.json`, methane/ethane/propane
at 240 K / 3 MPa) is a two-phase, well-separated state and is unaffected
(confirmed by `tests/test_cli_tp_flash.py`, unchanged).

## Alternatives considered

- **The Venkatarathnam & Oellrich (2011) `Pi` criterion**
  (`Pi = -(T/Cp) (dP/dT)_v^2 (dV/dP)_T`, *Fluid Phase Equilibria* **301**
  200-203), the criterion actually designed for exactly this
  vapor/liquid-without-a-second-root problem. **Rejected**: it needs
  `(dP/dT)_v` and `Cp`, i.e. temperature derivatives of the equation of
  state. PC-SAFT's implementation in this package (ADR-0014) has none - it is
  built for isothermal density-root solving only - and adding one is a
  materially larger slice (a second analytic derivative chain through every
  PC-SAFT term) for a criterion whose only advantage over `kappa` is
  better-known critical-point behavior, which is not this package's near-term
  need. `kappa` needs only what `_pcsaft_density.py` already computes.
- **A raw compressibility factor or reduced density threshold** (e.g.
  `Z < 0.3`). Rejected: `Z` alone is not dimensionless in the same
  model-independent sense - a `Z` threshold tuned to Peng-Robinson's cubic
  has no reason to transfer to PC-SAFT's very different equation of state,
  and the measured table above shows `Z`-based classes overlapping between
  the two models in ways `kappa` (built from the same physical derivative
  in both) does not. `kappa` is the one quantity with the same value (1) at
  the same physical limit (the ideal gas) in every model.
- **A finite-difference `dP/drho` inside the shipped implementation**, rather
  than the analytic one. Rejected for the same reason ADR-0015 rejected it
  for the density-root solver itself: the quantity being differenced is what
  decides mechanical stability and the identity right where it matters most
  (near the threshold), so cancellation error there is a correctness risk,
  not a rounding nuisance. The analytic derivative is used in the shipped
  code; a finite difference is used only in `tests/test_phase_identity.py`,
  built independently (PR: from `Component.tc_k/pc_pa/omega` directly, never
  the package's private `_mixture_parameters`; PC-SAFT: from the public
  `pressure_Pa(T, rho, x)`, never the private isotherm), specifically so the
  cross-check is not circular.
- **Threshold values other than 0.5.** Any value in `(0.23, 0.74)` separates
  both measured grids identically; `0.5` was chosen because it is the
  natural "closer to the liquid limit than the ideal-gas limit" reading of
  `kappa` and sits roughly in the middle of the measured gap, not because a
  tighter value was needed.
- **Keep the Wilson-ranking orientation for two-phase results and only fix
  the single-phase tie-break.** Considered, since no two-phase state in
  either validated grid actually needed the fix. Rejected: the slice
  declaration and ADR-0008 decision 3 both frame the two-phase orientation as
  the same class of problem (a convention where the model *could* measure and
  chose not to), the mechanism (`phase_identity` on the already-computed
  roots) is the same few lines either way, and leaving it as a pure
  historical-convention fallback - exercised on paper, not on these grids -
  is more honest than asserting it is unreachable.

## Consequences

- Positive: single-root EOS states are now named from a measured property of
  the model, not from candidate-construction order. The two motivating states
  now report `"liquid"`, `vapor_fraction = 0.0`,
  `phase_label_method = "compressibility"`.
- Positive: `stability_tp`'s own `feed_branch` / `phase_branch` - not only
  `flash_tp`'s phase names - are corrected by the same mechanism, so a caller
  reading `StabilityResult` directly, not only `flash_tp`, benefits.
- Positive (validated): 46/144 Peng-Robinson phi-phi grid states and both
  motivating PC-SAFT states relabeled correctly; zero numeric changes anywhere
  in 155 pinned bit-identity states, zero verdict changes, zero composition
  changes. `tests/test_phase_identity.py` cross-checks the analytic PR
  derivative against an independent finite difference to `< 1e-8` relative
  and the PC-SAFT one against the public `pressure_Pa` finite difference
  directly.
- Tradeoff: `phase_identity`'s default returning `None` means a third-party
  `EquationOfState` implementation gets no benefit from this slice until it
  implements the method; that is the documented, intentional fallback, not an
  oversight.
- Tradeoff: two-phase orientation now costs two extra `phase_identity` calls
  per phi-phi split (cheap: PR is a closed-form evaluation on an
  already-solved cubic; PC-SAFT reuses the already-found density root and one
  `pressure_and_slope` call) and single-phase naming costs one to two extra
  calls per `stability_tp` invocation, at `feed_branch` computation and once
  per converged trial - never inside the successive-substitution or
  second-order iteration itself.
- Known limitation, documented rather than hidden: near-critical states where
  the two converged phases fall on the same side of the threshold fall back
  to the Wilson ranking. No state in either validated grid exercises this
  path for a two-phase result; it remains available, untested-on-these-grids,
  for when one does.

## Next slice

None specifically opened by this one. It closes the "Also queued" item in
`brain.md`'s roadmap that named `flash-phase-labels-by-density` (ADR-0015
decision 2, ADR-0016 "Next slice") as the last user-visible phi-phi field
that was a convention rather than a measurement - compressibility, not raw
density, turned out to be the more portable measurement (see "Alternatives
considered" on `Z`/density thresholds), and this slice also carries phase
*densities* implicitly (both models' `phase_identity` computes a root's
density internally already; `PCSAFTEOS.density_roots` /
`PengRobinsonEOS.compressibility_factor` remain the public way to read it).
`pcsaft-association` remains the recommended next increment per ADR-0015 and
`brain.md`'s roadmap; this slice does not change that recommendation.

## Supersedes (optional)

Amends ADR-0008 decision 3 ("Which converged phase is *named* 'vapor' is
decided by volatility ordering"): the Wilson-ranking orientation is now a
documented fallback for the near-critical case, not the rule, and the
"`EquationOfState` exposes no molar volume, so no density-based
identification is available" premise no longer holds for either packaged
model. ADR-0008's decisions 1, 2, 4, 5, 6 and 7 are otherwise unchanged.

## Superseded by (optional)

None.
