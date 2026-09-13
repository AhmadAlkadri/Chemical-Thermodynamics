# ADR-0020: The phase addition/removal search serves an equation of state

Status: accepted
Date: 2026-09-13

## Context

ADR-0011 made the phase *count* an output for the `modified-raoult` path: after
a converged phase set fails its post-split stability test, the incipient phase
found there is added and the set re-solved; a phase whose fraction converges to
zero or below is removed. It was wired to that path only, and said why - "no
state in this repository exercises a third phase" on phi-phi or gamma-gamma,
and ADR-0002 forbids shipping a path nothing exercises.

Two later slices supplied the states. ADR-0018 gave PC-SAFT the association
term, so water and the alkanols became usable; ADR-0019 let each phi-phi phase
sit on its own density root, so a liquid-liquid split became expressible at a
state where a vapour root exists. What was left was exactly the three-phase
case, and ADR-0019 bracketed it rather than guessing: PC-SAFT water / n-hexane
at 1 atm, `z = (0.5, 0.5)`, returned a stable two-liquid result at 322-328 K
and a stable vapour-liquid result at 335-350 K, and **raised
`ConvergenceError` at 330 K and 334 K**, because the post-split stability test
found a converged phase unstable.

A binary at a fixed pressure has no three-phase *region* - Gibbs' phase rule
gives `F = 2 - 3 + 2 = 1`, so three phases meet at a single temperature, and
there the three phase amounts solve an underdetermined system (three unknowns,
two independent balances). The window that raised is therefore not a
three-phase region: it is the ADR-0011 Case R-3 mechanism on an equation of
state. Below `T3` the deepest tangent-plane minimum from the feed is a
**vapour**, the vapour-liquid pair that follows is unstable towards a second
liquid, and the three-phase solve then drives the vapour amount negative. The
answer is two conjugate liquids, reached by adding a phase and removing one.

Genuine three-phase EOS states exist in ternaries, and this slice found two:
PC-SAFT water / ethanol / n-hexane at 1 atm and 333 K has a vapour-liquid-
liquid tie triangle, and **Peng-Robinson with `k_ij = 0`** on the same ternary
at 280 K has a three-*liquid* one - which discharges ADR-0019's "no pure
Peng-Robinson three-phase case found in the databank".

## Decision

### 1. A phase set's model side is one object, and each phase carries its surface

New internal protocol `chemthermo.flash._multiphase._PhaseSetModel`, with two
implementations. The search itself is unchanged; what changed is how it asks
for a phase's tangent-plane fugacity term `t_i(x)`.

| | modified-Raoult (`_ActivityPhaseSet`) | equation of state (`_EosPhaseSet`) |
|---|---|---|
| what a label names | a phase **candidate** (activity liquid / ideal vapour) | a **density root** of one model |
| a phase's surface | that candidate's `ln_fugacity_terms` | its own `_PhaseRoot` (ADR-0019), pinned |
| phase identity | the label itself | `EquationOfState.phase_identity` on the root it converged on (ADR-0017) |
| two liquids ordered by | the order the search created them (roles, ADR-0011) | the first component's mole fraction (ADR-0019 decision 3) |

The pre-slice loop held a `Mapping[label, callable]` and looked each phase's
term function up by label. That cannot express a three-phase EOS set: two of
its phases carry the label `"liquid"` and are *different phases on different
roots*. So the loop now holds a list of surfaces parallel to its list of
labels, one per phase, and `_PhaseSetModel.surface(label)` mints a new one when
a phase is added. The two phases handed over from the two-phase split keep the
very `_PhaseRoot` holders that split used.

**Complexity receipt.** One protocol with four methods replaces a mapping. It
buys: a phase set in which two phases share a label; the ADR-0019 per-phase
root pinning through the multiphase successive substitution, the second-order
stage and the post-split report; and the ADR-0017 measurement of what each
converged phase is. Without it the loop would have to branch on the model
family internally, which ADR-0007 forbids. Cost: two small classes and one
extra `phase_identity` call per phase per round on the EOS path. If it were
omitted, a three-phase EOS answer could not be named or evaluated at all.

### 2. Addition pins the new phase to the branch the stability test found it on

`_PostSplitReport`'s instabilities already carry `branch` - `StabilityResult`'s
`phase_branch`, the ADR-0017-corrected identity of the root the minimizer sits
on. The new phase is pinned to it, which is verbatim what ADR-0019's
`_phi_phi_roots` does for the two phases of a two-phase split, extended to the
third. Removal is unchanged: a converged non-positive phase fraction, or a
recession direction from the multiphase Rachford-Rice, names the departing
phase.

### 3. An incipient phase that duplicates one already present is not retried

**Measured, and the reason this decision exists.** PC-SAFT water / n-hexane,
`z = (0.5, 0.5)`, 1 atm, **330 K**: the post-split test on the `LV` pair
reports *both* phases unstable, and the deepest minimum - the one ADR-0011
adds - is found from the **vapour** and is a water-rich liquid at
`x = (0.99998, 0.00002)`, which is the liquid already in the set at
`(0.99993, 0.00007)`. The multiphase Rachford-Rice removes it on the first
solve, before it can move, and the search runs `LV -> LLV -> LV -> LLV -> ...`
until the round budget.

The rule: when the phase the search *just added* is the phase the solve
removes, the **next** stationary point of the same post-split report is tried
instead, each one at most once. At 330 K that is the minimum found from the
water-rich liquid, `x = (0.01414, 0.98586)` - the hexane-rich liquid the
answer needs - and the search terminates at `LL`.

ADR-0011 rejected "retry from the second-deepest stationary point when the
post-split check fails" as "a heuristic with no stopping rule". This is not
that: it fires only after a solve has *demonstrated* the added phase was a
duplicate, the candidate list is the one finite report, and each entry is used
at most once, so the rule terminates with the phase count. The deepest minimum
is still always tried first.

### 4. Naming an EOS phase set by measurement, not by label

`_ordered_phases` asks `_PhaseSetModel.identity` for every phase, then orders
liquids before the vapour and names them. On the EOS path the identity is
ADR-0017's compressibility measurement on the pinned root; a three-phase set is
`liquid1` / `liquid2` / `vapor` with the liquids ordered by the first
component's mole fraction, `vapor_fraction` = the vapour's fraction and
`phase_regime = "VLLE"`. A set with no vapour reports `vapor_fraction = None`
(the ADR-0011 rule), so `liquid1` / `liquid2` / `liquid3` comes back as
`phase_regime = "LLE"`. One new diagnostics key on a multiphase EOS result,
`phase_label_method` (`"compressibility"`, or `"tie-break"` for a model without
`phase_identity`), matching what a two-phase phi-phi result already reports.

On the modified-Raoult path `identity` returns the label and liquids keep their
creation order, so every ADR-0011 number, name and diagnostics key is
unchanged; the vapour-last sort is stable, which is what makes the two
orderings agree there.

### 5. Where the search is reachable from, and where it still is not

`_flash_tp_tangent_plane` (phi-phi) now calls `_post_split_report` instead of
`_post_split_stability` and hands over to `_flash_tp_phase_addition` when a
phase is unstable, behind `FlashSettings.max_phases` (default 3), exactly as
the modified-Raoult path does. `max_phases = 2` reproduces the pre-slice
refusal, with the ADR-0011 wording ("Raise max_phases ... or pass
`post_split_stability=False`") in place of the old "flash_tp returns at most
two phases in this release". `post_split_stability=False` still returns the
two-phase result and never enters the search.

**`gamma-gamma` is deliberately not wired.** It would be the same two lines,
and no state in this repository exercises a third liquid phase there (Case
L-4), so it would ship a path nothing exercises. The three-liquid Peng-Robinson
state found here is on the *phi-phi* path; an activity-model analogue is still
the trigger that is missing.

### 6. Bit-identity, by non-entry rather than by construction

The search is entered only when the post-split stability test fails. Every
state of the 155-state bit-identity fixture, every state of the 144-state
Peng-Robinson phi-phi grid and every state of the 188-state PC-SAFT grid is
post-split **stable**, so none of them reaches the new branch and every number
and diagnostics key is unchanged. Two independent checks:
`tests/test_flash_refactor_bit_identity.py` passes unchanged against
`refactor_bit_identity_v2.json` (no v3, none needed), and
`tests/test_flash_vlle_eos.py::test_no_bit_identity_fixture_state_carries_a_phase_search_key`
asserts directly on the fixture that no pinned state carries a key the search
produces.

## Alternatives considered

- **Pick the incipient phase that is furthest from every phase already in the
  set**, rather than the deepest minimum with a duplicate-retry. Rejected: it
  would change which phase the modified-Raoult path adds on states whose
  numbers are pinned (Cases V-1, V-3, V-5), for a gain only the 330 K state
  needs, and "furthest" needs a tolerance the package does not have.
- **Detect the duplicate before adding it**, with a looser version of the
  trivial-solution metric. Rejected: at 330 K the two compositions differ by a
  factor of three in the trace component, so `sum_i ln(x_i/x_i')^2 = 1.19`,
  nowhere near `trivial_tol`. Loosening that metric would change the
  stability module's own trivial-solution test.
- **Keep a `Mapping[label, callable]` and key it by `(label, position)`.**
  Rejected: it encodes the phase's identity in a dictionary key instead of in
  the phase, and the second-order stage and the post-split report would each
  have to reconstruct it.
- **Return a three-phase answer for the binary at `T3`.** Rejected as
  arithmetic on an underdetermined system; see "Context". What the tests check
  at `T3` is the geometry - the two-phase answers on either side meet there -
  and the verdict boundary, which locates `T3` to 2.4e-07 K.
- **Wire `gamma-gamma` too, for symmetry.** Rejected; see decision 5.

## Consequences

- Positive: `flash_tp(mixture, temperature_K=..., pressure_Pa=101325.0,
  eos=PCSAFTEOS())` no longer raises anywhere in `[T3 - 1 K, T3 + 1 K]` on
  water / n-hexane (measured on 82 states at two feeds), and the states that
  raised return the two conjugate liquids, matching an independent two-phase
  Newton to 2.3e-13 and the lever rule to 1.2e-13.
- Positive: three-phase equation-of-state answers exist at all. PC-SAFT
  water / ethanol / n-hexane at 333 K returns a `VLLE` tie triangle
  reproduced by an independent 6-equation Newton to 2.6e-11, with amounts from
  the mass balance to 5.5e-11 and `G3 < G2 < G1`; Peng-Robinson on the same
  ternary at 280 K returns three liquids in ~0.1 s.
- Positive: `delta_g_vs_two_phase_rt` on an EOS result is exactly the Gibbs
  energy difference between the answer and the two-phase pair the search
  started from (checked against independent Newton solves to 1e-9).
- Tradeoff: a three-phase PC-SAFT flash costs about 35 s on the development
  machine - 4,725 `fugacity_coefficients` calls, each a full density-root
  scan, half of them inside the second-order stage's finite-difference
  Hessian. No performance work was done here; it is the next-but-one slice
  (`perf-profile-baseline`).
- Tradeoff: the default suite grew and had to be trimmed to pay for it; see
  "Runtime" below and validation Case P-9.
- **Known limitation, measured, unchanged by this slice.** A phase count is
  never better than the stability test that produced it (an invariant since
  ADR-0011). At `z_water = 0.7` and `T > T3`, water / n-hexane returns the two
  liquids, whose reduced Gibbs energy is **higher** than the vapour-liquid
  pair's by 2.5e-03 at 335 K: all four deterministic trials from the
  hexane-rich liquid converge to the trivial solution or to its partner, and
  the vapour stationary point at `tpd = -6.5e-03` is missed. That is
  `_EOSTangentPlane`'s trial set, which ADR-0012 deliberately left with
  per-iterate minimum-Gibbs root selection and no fixed surfaces. Pinned in
  Case P-9 (iv) and in `test_the_window_scan_never_raises[0.7]` so that
  improving the trial set is noticed rather than silently absorbed.
- **Known limitation, measured.** One feed of the 36-point ternary scan at
  333 K, `z = (0.1, 0.1, 0.8)`, raises `ConvergenceError` ("a two-phase set
  converged to a non-positive phase fraction"): the three-phase solve from the
  stability seed does not converge and the removal it signals leaves a
  two-phase set that does not converge either. That feed **also raised before
  this slice** (with `max_phases = 2` it still does), so it is not a
  regression; it is a hard state near the edge of the triangle, reported as a
  failure rather than as a wrong answer.

## Runtime

`pytest -q` was 337 s at the start of this slice (733 tests), against a ~240 s
budget, and
this slice adds a three-phase capability whose flashes are expensive. Two
levers, both already established in the repository:

- Six heavy examples gained a `--full` flag and a cheaper default, as
  `examples/validation/15_flash_split_robustness.py` did in ADR-0017:
  `17_pcsaft_lle_vs_feos` (34.5 s -> 12.1 s), `flash_tp_pcsaft_lle_demo`
  (18.1 -> 4.6), `12_vlle_verdict_map` (16.4 -> 7.0),
  `14_pcsaft_flash_vs_teqp` (16.2 -> 4.6), `16_pcsaft_association_vs_feos`
  (10.2 -> 0.7), `pcsaft_association_demo` (6.7 -> 5.1). Everything a default
  run drops is covered by a test, and `--full` still runs it.
- Six existing tests (and two parameters of a seventh) became `@pytest.mark.slow`,
  each of them a repetition of something the default run still does.

`.agents/dev-contract.md`'s rule for the marker is widened accordingly: a test
may be `slow` when it is a *repetition* of a capability the default run already
covers - a further feed on the same tie line, a further temperature on the same
map - and not only when it is a full grid with a representative subset.

Measured: `pytest -q` **337 s (733 tests) -> 238 s (738 tests)**, on the same
machine, while adding a three-phase capability whose own fast tests cost ~19 s
and whose two golden-path examples cost ~23 s. `pytest -q -m slow` is now **21 tests
in 855 s** (14:14): the 188-state grid (126 s) plus this slice's two 41-point
scans (236 s and 109 s), the PC-SAFT ternary tie triangle (187 s), the
ternary FeOs check (49 s), the verdict-boundary bisection (47 s) and nine
repetitions.

## Supersedes (optional)

Discharges ADR-0011's "What remains: phi-phi and gamma-gamma still stop at two
phases" for **phi-phi**; gamma-gamma is unchanged and the reason is stated in
decision 5. Discharges ADR-0019's "a phi-phi state that needs three phases
still raises" and its "Next slice" note. Amends ADR-0009 decision 4's raise
("a third phase is required ... multiphase flash is the next slice"), which now
applies to `gamma-gamma` only. Records a Peng-Robinson three-phase state, which
ADR-0019 and validation Case P-8 recorded as not found.

## Superseded by (optional)

None.
