# ADR-0011: Phase addition and removal in `flash_tp`

Status: accepted
Date: 2026-09-13

## Context
ADR-0008, ADR-0009 and ADR-0010 made the phase *count* an output of a
tangent-plane stability test rather than an input - up to two. Beyond two, the
package could only refuse:

- ADR-0009 decision 4 re-tests every converged phase and raises
  `ConvergenceError("a third phase is required")` when one of them is unstable.
- ADR-0010 recorded the measured consequence (validation Case R-3): in a window
  about **0.135 K** wide just below the three-phase temperature of water /
  n-butanol at 1 atm, the deepest tangent-plane minimum from the feed is the
  *vapor*, so the split converges a vapor-liquid pair whose liquid lies inside
  the miscibility gap. Refusing it is right - that pair is not the equilibrium -
  but the diagnosis is only half right. Below T3 the answer is not a third
  phase: it is a **different pair of two**, the two conjugate liquids. Reaching
  it needs the second liquid to be *added* and the vapor then *removed*.
- Case L-4 recorded the other half: no state in this repository's
  Peng-Robinson grid genuinely needs a third phase, so the failure path had
  only ever been exercised synthetically.

A ternary does have a genuine three-phase *region*. The independent reference
built for this slice (1-propanol / n-butanol / water, modified Raoult, 1 atm)
has a tie-triangle from about 363 K to the binary T3 = 366.2138 K, and any feed
inside it is a real three-phase state that ADR-0010 could only refuse.

So two things are needed, and they are the same mechanism: a Rachford-Rice that
works for any number of phases, and a loop that reads its answer to decide how
many there are.

## Decision

### 1. Multiphase Rachford-Rice as a constrained convex minimization
New internal module `chemthermo/flash/_multiphase_rr.py`. With a reference
phase `r` and `K_i^j = x_i^j / x_i^r`, material balance gives
`t_i = 1 + sum_{j != r} beta_j (K_i^j - 1)`, `x_i^r = z_i / t_i`,
`x_i^j = K_i^j z_i / t_i`, and the `NP - 1` Rachford-Rice equations
`f_j = sum_i z_i (K_i^j - 1) / t_i = 0`. Their Jacobian is symmetric, so they
are the gradient of

    F(beta) = - sum_i z_i ln( t_i )

(Okuno, Johns & Sepehrnoori, SPE J 15 (2010) 313, eq. 9; Michelsen 1994), whose
Hessian `A^T diag(z / t^2) A` with `A_ij = 1 - K_i^j` is positive semidefinite.
**Solving the multiphase Rachford-Rice is therefore a convex minimization, not
a root find**: one minimum or none, never a spurious root.

The feasible region is Okuno et al.'s `S = { beta : a_i . beta <= b_i }` with
`a_i = (1 - K_i^j)_j` and `b_i = min( 1 - z_i, min_j (1 - K_i^j z_i) )`,
derived from `x_i^j >= 0` alone. Every point of `S` has
`t_i >= max(z_i, max_j K_i^j z_i) > 0`, so **`S` contains no pole**, not even on
its boundary, and a full Newton step to the boundary is safe. The
Leibovici & Neoschil region `t_i >= 0` is deliberately not used: its boundary
*is* the pole set. The initial estimate is the mean of the vertices of
`S` intersected with `beta_j >= 0, sum_j beta_j <= 1`, falling back to the
vertices of `S` alone when that intersection is empty (a negative flash).

Two deviations from the paper, both recorded:
- the stopping test is on `max_j |f_j|` **divided by**
  `sum_i z_i |K_i^j - 1| / t_i` (floored at 1). `f_j` is a cancelling sum of
  terms of that size, so an absolute tolerance is not scale free, and a phase
  set can hold a vapor/liquid `K` of 30 next to a liquid/liquid `K` of 1.1. With
  an absolute 1e-12 the n-butanol / water liquid-liquid set stalls at 2.0e-11
  and the flash fails.
- the line search accepts an Armijo decrease of `F` **or** a decrease of the
  residual. Near the minimum the change in `F` falls below the floating-point
  resolution of `F` and the Armijo test can no longer be met while the residual
  still has digits left. This is the same rule the ADR-0009 second-order stage
  already uses.

### 2. Signs of `beta` are not constrained; that is the removal signal
`S` constrains compositions, not phase amounts, so a converged `beta_j <= 0` is
a meaningful answer - Okuno et al.'s "negative flash" - and it is exactly the
statement *this phase does not exist at this feed*. Phase removal reads that
sign. It is not a failure mode to be guarded against, and no positive threshold
is used: Okuno Example 3 converges to `beta = 2.2e-06`, a real phase.

When the feasible region **recedes** along a descent direction, `F` has no
minimum at all and the phase set has no split of this feed - which is what
Gibbs' phase rule says about three phases in a *binary* at any temperature
other than its single three-phase temperature. The solver reports that
(internal `_NoMultiphaseSolution`) together with the rate of change of every
phase fraction along the recession direction, and the phase whose amount runs
to minus infinity is the one removed.

### 3. Multiphase successive substitution and a multiphase second-order stage
New internal module `chemthermo/flash/_multiphase.py`. Each phase carries a
*phase candidate* (ADR-0010) and its composition. Equal fugacity against the
reference phase is `ln K_i^j = t_i^r(x^r) - t_i^j(x^j)`; the inner solve is
decision 1; the compositions follow from `t_i`. For two phases this is the
existing loop of `_split.py`.

Successive substitution converges linearly and slowly here - 324 iterations to
1e-12 on the validated tie-triangle - so it hands over to a Newton minimization
of the total reduced Gibbs energy in the mole numbers of the non-reference
phases,

    g(n) = sum_j sum_i n_i^j [ ln x_i^j + t_i^j(x^j) ],
    dg / dn_k^j = [ln x_k^j + t_k^j] - [ln x_k^r + t_k^r]

the gradient being exactly the equal-fugacity residual, by the same
Gibbs-Duhem cancellation as ADR-0009. This is that stage generalized from one
non-reference phase to `NP - 1` of them; the derivation is unchanged. From a
20-iteration start it reaches 7e-16 in 3 iterations. The Hessian is central
differences of the gradient, symmetrized and shifted positive definite; the
step is backtracked inside `0 < n_i^j`, `sum_j n_i^j < z_i`.

Minimizing `g` rather than solving the residual system is the ADR-0009 argument
verbatim: a residual system cannot tell the equilibrium from the trivial
solution, and a descent method can.

### 4. The add / remove loop, bounded by `max_phases`
    solve the phase set
      -> a phase fraction <= 0, or two phases collapsed onto each other?
             remove it, re-solve
      -> a converged phase unstable?   add the minimizer found there, re-solve
      -> all fractions positive and every phase stable?   that is the answer

Addition seeds the new phase with the tangent-plane minimizer `w` from the
post-split test on the failing phase. Setting `x^new = w` and computing `K`
from decision 3 reproduces Michelsen's `W`-scaled seed `K_i = W_i / x_i^r`
exactly, because `ln W_i = ln x_i^r + t_i^r - t_i(w)` at a stationary point - so
this is the *same* seeding rule the two-phase paths already use, not a new one.
Collapse of two phases onto one composition is detected with the stability
module's trivial-solution metric, so no new tolerance is introduced.

`ConvergenceError` is raised only when a phase set is unstable at `max_phases`
phases, when a solve fails, or when the round budget (`max_phases + 4`) runs
out.

`post_split_stability=False` still returns the converged two-phase answer with
the failure in diagnostics and never enters the loop. That keyword means "do
not police the phase set", and it is how Case R-3 inspects the pair the search
would otherwise replace.

### 5. Result contract
Phase names: `"vapor"` for the ideal-gas candidate, `"liquid"` for a single
liquid, `"liquid1"` / `"liquid2"` / ... for several, numbered in the order the
search created them. Liquids are reported before the vapor whatever order the
solve held them in. **The liquid numbers are roles, not identities**, exactly as
in ADR-0009; `"vapor"` is the one name with a model-level meaning.
`vapor_fraction` is the fraction of the `"vapor"` phase when the set contains
one and `None` otherwise - the ADR-0009 rule, stated for any phase count.

Diagnostics on a result that went through the search:
`phase_count`, `phase_state`, `phase_regime` (`"VLLE"` joins `"VLE"`, `"LLE"`,
`"single-phase"`), `phase_set_history` (e.g. `"V -> LV -> LLV -> LL"`),
`phases_added`, `phases_removed`, `equilibrium_residual` (max over **all** phase
pairs and components), `mass_balance_residual`, `delta_g_split_rt`,
`delta_g_vs_two_phase_rt`, `ssi_iterations`, `second_order_iterations`,
`rachford_rice_iterations`, `converged_stage`, plus the usual `post_split_*` and
`phase_stability_*` keys.

**Those keys appear only on results that actually entered the search.** Every
single- and two-phase result of the earlier slices is unchanged down to the
last bit, its diagnostics mapping included, which is what
`tests/test_flash_refactor_bit_identity.py` (fixture captured at `e927623`,
before this work) checks and what it would catch.

`FlashSettings` gains `max_phases: int = 3`, validated `>= 1`. `max_phases=2`
reproduces the pre-ADR-0011 refusal exactly. The cap is consulted only for the
*third and further* phases, so `max_phases=1` behaves like `2`: the
one-versus-two decision is thermodynamic, not a setting, and cannot be turned
off by a number.

## Alternatives considered
- **A `flash_mode="vlle"` that assumes three phases.** Rejected by the thesis
  of this whole campaign: the phase count is an output. It would also have no
  answer for the window below T3, where the correct count is two but the naive
  two-phase answer is wrong.
- **Full Gibbs-energy minimization over all phases from scratch, as the
  reference path.** Rejected as the *reference*: the successive-substitution /
  Rachford-Rice first stage is what Michelsen recommends and what the rest of
  this package already is, and the second-order stage here is already a Gibbs
  minimization where it matters. A from-scratch global minimizer would also
  need a starting phase set, which is the part the stability test supplies.
- **Keep raising.** Rejected: the state exists, the model has an answer, and
  ADR-0010 had already recorded the refusal as a half-diagnosis.
- **Retry from the second-deepest stationary point when the post-split check
  fails.** Rejected as a heuristic with no stopping rule; it would also not have
  produced the ternary three-phase answer, only a different pair of two.
- **A positive threshold for "this phase has vanished".** Rejected: Okuno
  Example 3's `beta = 2.2e-06` is a real phase. Existence is decided by the
  post-split stability test, not by a magnitude.
- **Wiring the loop to phi-phi and gamma-gamma as well.** Rejected for this
  slice: see "What remains".

## Consequences
- Positive: `flash_tp(..., flash_mode="modified-raoult")` returns verified
  three-phase states. Validation Case V-1: six feeds inside the tie-triangle at
  363, 364 and 365 K reproduce an independent Newton solve to |dx| <= 2.2e-14
  and |dbeta| <= 1.6e-13, with equilibrium residual <= 1.8e-15, mass balance
  <= 6.9e-18, every phase post-split stable, and G3 < G2 < G1.
- Positive: the ADR-0010 refusal window is resolved, by removal, to the binodal
  tie-line to 3.3e-13 (Case V-3), with the route recorded as
  `"V -> LV -> LLV -> LL"`.
- Positive: the multiphase Rachford-Rice is validated on published data
  independent of this package (Okuno et al. Table 1, Case V-4), including the
  negative flash that phase removal depends on.
- Tradeoff: one new public setting (`max_phases`) and a set of diagnostics keys
  that are present on some results and absent on others. The alternative -
  putting them on every result - would have moved the pinned two-phase numbers.
- Tradeoff: a three-phase solve is roughly 50 successive substitutions plus a
  handful of Newton steps plus ~100-150 Rachford-Rice Newton iterations, on top
  of a two-phase solve and several stability tests. No performance work was
  done.
- **Known limitation, measured.** A thin tie-triangle can hide from the
  deterministic stability trial set. At 363 K a feed weighted towards the two
  close liquid vertices reports `tpd_min = 0.0` and comes back a single liquid,
  which is wrong for this model. The failure is in the *stability* test, not in
  the search - the search is never entered - and it is the honesty note of
  `StabilityResult` made concrete. Pinned in Case V-2 so that a future
  improvement to the trial set is noticed rather than silently absorbed.

## What remains
- **phi-phi and gamma-gamma still stop at two phases**, whatever `max_phases`
  says, and raise as before. The machinery is written generically and both
  could adopt it in a few lines, but no state in this repository exercises a
  third phase on either path (Case L-4), and shipping a path nothing exercises
  is what ADR-0002 forbids. A PR three-phase case in the databank, or a ternary
  with three liquid phases, is the trigger.
- **No PC-SAFT.** Untouched here.
- **Performance.** Unmeasured and unoptimized; see the tradeoff above.
- **The `chemthermo.vlle` plugin boundary.** `chemthermo.vlle` is a public
  plugin boundary (`get_vlle_engine`, `VLLEEngine`, `VLLEResult`) added before
  any of this existed, and `flash_tp(flash_mode="vlle")` still raises a
  `ModelError` telling the caller to install a plugin - for a capability that
  is now in-tree. **This ADR deliberately does not touch it.** Deprecating or
  repurposing a public API is its own decision with its own compatibility
  question (does the boundary become the way to reach a *different* engine, or
  does it go?), and a follow-up ADR, `vlle-plugin-boundary-disposition`, will
  answer it.

## Supersedes (optional)
None. Discharges the "a third phase is required, and flash_tp returns at most
two phases in this release" refusal of ADR-0009 decision 4 and the "Known
limitation, measured" note of ADR-0010, for the `modified-raoult` path only.

## Superseded by (optional)
None.
