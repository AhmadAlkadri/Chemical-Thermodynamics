# ADR-0029: The map's nine multiphase refusals are one irreversible choice and one unusable Hessian

Status: accepted
Date: 2026-09-14

## Context

ADR-0027's amendment (slice `robustness-map-coverage`) grew the robustness map
to 2505 states and made it hard again: **2491 converge, 14 refuse, 0 converge
and violate an invariant** at `74820b8`. Nine of the fourteen are one class,
`multiphase-solver-failure`, in one family, `eos-three-phase`, and brain.md's
roadmap ranked them first because it is the largest class any grid in this
repository has found since ADR-0028. The other five are the deprecated
`gamma-phi` path's `rr-no-bracket` states, which are by design (ADR-0016
decision 8) and are not touched here.

The nine are two defects, not one, and neither is in the thermodynamics. Both
sit in `chemthermo/flash/_multiphase.py`, the ADR-0011 phase addition/removal
search and the phase-set solve underneath it. Three framings were ruled out
before anything was written:

- *"the three-phase sets are near-degenerate and the collapse rule should fire"*
  - no for both classes. In class (i) the two liquids differ by
  `|x_water| = 0.98`; in class (ii) the three liquids are distinct enough that
  every two-phase subset of them is a different tie line.
- *"successive substitution needs a bigger budget"* - half true and not the
  fix. The class (ii) states do converge under successive substitution alone,
  in 52 to 109 iterations against the multiphase budget of
  `min(ssi_iterations, max_iter) = 50`. But the budget is 50 precisely because
  a second-order stage is supposed to take over there, and raising it would
  change the handover point for every multiphase state in the repository. The
  defect is that the stage which was supposed to take over cannot run.
- *"the tolerance is too tight"* - no. Every one of the nine converges below
  `1.5e-14` (class ii) or `7.5e-13` (class i) on an acceptance of `1e-08`.
  Nothing was relaxed.

### Class (i): 5 states where a removal cannot be taken back

`multiphase-solver-failure` / `collapsed`. PC-SAFT (2B water, `k_ij = 0`)
water / n-hexane at 1 atm, `z_water = 0.05`, at `T3 + {0.01, 0.1, 0.5, 1.0} K`
with `T3 = 334.807826336 K`; and the PC-SAFT ternary water / ethanol /
n-hexane feed `(0.1, 0.1, 0.8)` at 333 K, which Case P-10 (i) had already
recorded as a pre-existing failure at the tie-triangle's edge.

Traced at all five. The route is the same every time:

```
L  ->  LL   (the two conjugate liquids; post-split unstable on both phases)
   ->  LLV  (a vapour added, ADR-0020)
   ->  LV   (the three-phase set names a liquid to remove; it is removed)
   ->  raise
```

The third step is where it goes wrong, by two spellings of one signal. On the
binary, a fixed pressure admits three phases at exactly one temperature
(Gibbs' phase rule), so above `T3` the `LLV` set has **no** Rachford-Rice
solution at all: `_multiphase_rachford_rice` raises `_NoMultiphaseSolution`
with the recession direction and `_multiphase_ssi` reads the phase whose
fraction rate is most negative. On the ternary the region is finite but this
feed is outside it, so the solve stays inside the Okuno region and returns a
converged **negative flash** instead - `beta = (-0.761, -4.309, +6.070)` at the
first substitution - and `_multiphase_ssi` reads `argmin` of the fractions. In
both spellings what comes back is the **hexane-rich** liquid, the phase that holds
95 % of the feed. The water-rich pair that is left over then negative-flashes
(`beta = (-0.207, +1.207)` at `T3 + 0.01 K`), and the search, which had no way
to take a removal back, raised "A two-phase set converged to a non-positive
phase fraction".

The answer is the *other* removal. Solved directly, the hexane-rich liquid
against the vapour converges to a residual of `7.4e-13` with
`beta = (0.8557, 0.1443)` and is post-split stable
(`tpd_min = -1.6e-12`). Nothing about the ranking rule is wrong in general -
the recession direction is the right signal, and it is what resolved validation
Case R-3 - but it is a *ranking*, read as if it were a decision.

### Class (ii): 4 states where the second-order stage cannot form a Hessian

`multiphase-solver-failure` / `split`. Peng-Robinson (`k_ij = 0`) water /
ethanol / n-hexane at 1 atm: `(0.2, 0.6, 0.2)` at 280 K and 300 K, and
`(0.4, 0.4, 0.2)` and `(0.5, 0.3, 0.2)` at 300 K. The reported residuals span
`1.836e-08` (1.8x `tol`) to `1.478e-04` (1.5e+4x `tol`), and every one of them
reports *one* second-order iteration.

That "1" is the whole diagnosis. These are three-liquid sets whose water-rich
phase holds n-hexane at `x = 5.3e-14` (280 K) to `1.6e-12` (300 K).
`_multiphase_second_order` minimizes in the **linear** mole numbers of the
non-reference phases with a central difference of `1e-7` and the box
`0 < n_i^j`, `sum_j n_i^j < z_i`. Its first Hessian column moves a
non-reference `n_i` by `1e-7`, which is eleven orders of magnitude more than
the reference phase's entire inventory of n-hexane; the perturbed point leaves
the box, `energy_and_gradient` raises, and the loop `break`s at iteration 1
having done nothing. Successive substitution's unfinished residual is then
reported as a failure.

The obvious repair - carry the stage in `u = ln n`, as ADR-0024 does for the
two-phase split - is **not sufficient on its own**, and measuring that is what
produced decision 2. The quantity that cannot be resolved is
`n_i^r = z_i - sum_j n_i^j`, and the reference phase's mole numbers are not
variables in any parametrization of the others: at 300 K it is
`0.2 - 0.167 - 0.0326 = 2.0e-13`, a difference of `O(1)` doubles whose relative
accuracy is already spent. Run as written, the log-space stage accepts steps on
an Armijo test whose energy change is `+-9e-16` - floating-point noise - and
drives the residual *up*, from `1.2e-04` to `1.3e-01` in one step and nowhere
in 100.

## Decision

### 1. Removal becomes reversible

`_multiphase_ssi` and `_solve_phase_set` now report a *ranking*
(`_MultiphaseSolution.removal_order`) rather than a single index. The first
entry is `argmin` of the fractions, or of the recession rates, filtered to the
non-positive ones - exactly the index the pre-slice code acted on, so the
ordinary path is unchanged by construction. `_flash_tp_phase_addition` pushes
the pre-removal phase set and the untried tail onto an undo stack
(`_RemovalChoice`) each time it removes, and the branch that used to raise "A
two-phase set converged to a non-positive phase fraction" now pops that stack
and takes the next candidate instead. The raise is kept, verbatim, for the case
where the stack is exhausted.

It terminates for the same reason the ADR-0020 duplicate-incipient rule does:
a candidate is consumed when it is taken, the undo path only ever pops (an
exhausted entry is discarded and never pushed again), and the round budget
(`max_phases + 4`) bounds the loop regardless. The repaired route is
one round longer than the old one:

```
L -> LL -> LLV -> LV -> LV      (the second LV is the other pair)
```

Two different two-phase sets print the same compact label in
`phase_set_history`, because that label names the phase *types*. That is
recorded rather than worked around: the history is a route, and the route
genuinely visits two vapour-liquid sets.

### 2. A multiphase log-space stage, with the reference phase re-chosen

New private module `chemthermo/flash/_multiphase_log_space.py`. It is
`chemthermo/flash/_log_space.py` generalized from two phases to any number: the
variables are `u_i^j = ln n_i^j` for each non-reference phase, the
central-difference step is the multiplicative one, the Newton system is written
on the residual `r` rather than on `dg/du = diag(n) r`, the descent test and the
line-search acceptance rule are the same expressions, and the ADR-0026
curvature safeguard (Gill & Murray's modified Newton on
`H = diag(n) J + diag(n r)`) is the same list of directions. The constants are
**imported** from `_log_space` rather than restated, so the two stages cannot
drift apart.

What is new is the reference-phase choice (`reference_phase_index`): the
reference is the phase whose smallest mole fraction over the components present
in the feed is **largest**. That is the phase in which `z - sum n` is a
well-conditioned difference, and it moves every trace composition into the
variables, where a logarithm holds it exactly (`n = 2e-13` is `u = -29`). The
Rachford-Rice reference phase is a free choice - equation (1) of
`_multiphase.py` is symmetric in the phases - so this decides nothing
thermodynamic; the stage undoes the re-ordering before returning, so the
caller's labels, surfaces and compositions still pair positionally.

Measured on the four class (ii) states, from the 50-iteration successive-
substitution iterate: **1 to 3** Newton iterations, and an all-pairs
equal-fugacity residual of `8.9e-15` to `1.5e-14` in the delivered result.
Without the reference re-choice the same code reaches `4.9e-04` to `1.3e-01`
in 100 iterations - it drives the residual *up*.

### 3. The stage runs only where the solve was about to raise

`_solve_phase_set` calls it on exactly one branch: no phase to remove, and the
residual still above `FlashSettings.tol` after successive substitution and the
linear second-order stage - the branch that raised. The result is taken only if
it *lowers* the all-pairs residual `_phase_residual` (the stage's own residual
is measured against its reference phase, which is within a factor two but is
not the number a result reports). A two-entry ladder is walked, unsafeguarded
then safeguarded, in the spirit of ADR-0028's three-entry one.

### 4. Dormancy is by construction, and the record says so

Neither decision can be reached by a state that converges: decision 1 replaces
a `raise`, decision 3 gates on the branch that raises. The new diagnostics key
`log_space_iterations` is conditional on the stage having run, so a result that
converged before this slice carries the same diagnostics mapping key for key.
The claim is checked rather than argued: the full 2505-state sweep is re-run
and every previously converged state is compared field by field against
`benchmarks/robustness_74820b8.json`, `refactor_bit_identity_v3.json` is
unchanged, and all nine ADR-0023 benchmark `result_hash` values are identical.

## Consequences

- **The map's `multiphase-solver-failure` class is empty.** 2500 of 2505
  states converge; the 5 remaining refusals are the deprecated `gamma-phi`
  path's, pinned in `tests/test_robustness_map.py` as by design.
- Four Peng-Robinson three-liquid states and five PC-SAFT vapour-liquid states
  that used to raise now return verified answers (ledger Case P-18).
- One new conditional diagnostics key, `log_space_iterations`, and one new
  value for `converged_stage`, `"second-order-log"` - both names already used
  by the two-phase path (`_detect.py`), not new vocabulary.
- One message text changes, and only on a state that still fails: the
  "multiphase split did not converge" raise gains a clause naming the
  log-space iterations when that stage also ran. The refusal classifier keys on
  "multiphase split", which is before the clause, so the class and stage do not
  move; both forms are pinned in `tests/test_robustness_map.py`.
- **Cost.** A refused state used to be cheap. A state repaired by decision 1
  pays for one extra phase-set solve; one repaired by decision 3 pays for up to
  two log-space stages whose Hessians are `(N-1) C` square (8x8 for a
  three-phase ternary). Neither is on a path a converging state takes.
- **Not claimed.** That the model is right at these conditions - the PC-SAFT
  and Peng-Robinson `k_ij = 0` water / ethanol / n-hexane system is not
  validated against any published ternary datum, and ADR-0027's warning that
  the map measures coverage and not correctness still stands. What is claimed
  is that each returned phase set is an equilibrium *of the model it was given*,
  verified by a Newton sharing no code with `flash_tp`, with the lowest reduced
  Gibbs energy of every admissible alternative that was solved.
- **Not claimed either.** That the removal ranking is now right - only that it
  is no longer final. A geometry where every candidate on the stack dead-ends
  would still raise, and there is no state in this repository that does.

## Alternatives considered

- **Pick a better removal instead of making removal reversible.** The obvious
  candidate is "remove the phase whose composition is furthest from the feed",
  or "iterate the negative flash further before deciding". Both are changes to
  a rule that currently resolves validation Case R-3 and every converging
  three-phase state, i.e. a bit-identity risk taken for a state that can be
  repaired on the failure path alone. Rejected for the same reason ADR-0028
  rejected moving `_SEED_PHASE_FRACTION`.
- **Raise the multiphase successive-substitution budget from
  `ssi_iterations = 50`.** It repairs class (ii) - all four converge in 52 to
  109 iterations - and it moves the handover point for every multiphase state
  in the repository, including the ADR-0011 tie-triangle states whose numbers
  are pinned. Rejected on bit-identity, and it would leave the second-order
  stage still unable to run on a three-liquid set, which is the actual defect.
- **A relative (per-variable) finite-difference step in the linear stage.** It
  makes the first Hessian column formable, and it does not address the
  reference phase, where the cancellation is. Measured to be exactly the
  insufficiency that produced decision 2's second half.
- **Route to the log-space stage on a trace mole fraction, as
  `TRACE_MOLE_FRACTION` does for the two-phase path.** That threshold is
  `1e-30` and these states sit at `1e-12` to `1e-14`; lowering it would re-route
  states that converge today, which is precisely what ADR-0024 and ADR-0028
  both declined to do. Gating on the failure branch costs nothing and risks
  nothing.
- **Choose the reference phase this way in the multiphase Rachford-Rice and the
  successive-substitution loop too.** It is probably better conditioned there
  as well, and it would move every converged multiphase number in the
  repository. Out of scope for a slice whose subject is nine refusals; recorded
  as a narrower open item in brain.md.
- **Accept the class (ii) states at their successive-substitution residual by
  loosening `tol`.** Rejected outright: `1.478e-04` is not an answer, and the
  three feeds at 300 K converge to the *same* tie triangle once solved, which
  is the evidence that the residual was measuring distance and not noise.

## References

- ADR-0011 (the phase addition/removal search and the negative-flash removal
  signal), ADR-0020 (the search on the phi-phi path, and the
  duplicate-incipient rule this one is shaped after), ADR-0019 (per-phase
  density roots), ADR-0024 / ADR-0025 / ADR-0026 (the two-phase log-space
  stage, its `ln W` seed and its curvature safeguard), ADR-0028 (the seed
  ladder and the dormancy argument), ADR-0027 (the robustness map), ADR-0002
  (no unexercised paths).
- Validation Case P-18 in `.agents/brain/validation-cases.md`, with an
  amendment to Case R-MAP-2; Cases P-9 and P-10 for the states themselves.
- Okuno, Johns & Sepehrnoori, *SPE Journal* **15** (2010) 313 (the multiphase
  Rachford-Rice, its feasible region and its recession test); Michelsen,
  *Fluid Phase Equilibria* **9** (1982) 21; Gill, Murray & Wright, *Practical
  Optimization*, sec. 4.4.2 (the modified-Newton direction).
- `src/chemthermo/flash/_multiphase.py`,
  `src/chemthermo/flash/_multiphase_log_space.py`,
  `tests/test_flash_eos_multiphase_robustness.py`,
  `tests/test_robustness_map.py`,
  `tests/validation/test_pcsaft_vlle_water_hexane.py`,
  `benchmarks/robustness_64831bd.json` / `.md`, `benchmarks/README.md`.
