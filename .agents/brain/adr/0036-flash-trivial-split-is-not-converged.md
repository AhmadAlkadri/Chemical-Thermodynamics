# ADR-0036: A trivial split is not a converged split

Status: accepted
Date: 2026-09-25

## Context
ADR-0028's stability-seed ladder rescues phi-phi states whose successive
substitution diverges (the polyethylene / n-pentane band, `max_delta_k` up to
1e+128). Its gates asked for "residual <= tol" (in `_phi_phi_second_order`)
or "residual <= tol and 0 < beta < 1" (elsewhere). The cross-platform slice
(ADR-0032) found the grid point Mw 53000 / 15 wt% / 8.1 MPa converging by a
different ladder rung on each machine, and a scan of the 135 pressures within
+-64 ULP of it (+-1e-8..1e-6 relative included) refused **29** on one Linux
host, in three shapes: "vapor fraction outside (0, 1)", "a two-phase set
converged to a non-positive phase fraction", "the multiphase split did not
converge". Diagnosis: the linear-iterate log-space retry, started from the
diverged K-loop's wreckage, converges on the **trivial** solution - both
phases the feed composition to 2e-15 in `ln x`, `dG_split = 0`, residual
~1e-13. The trivial solution satisfies the equal-fugacity equations exactly
and, since both phases are the feed, the mass balance for *any* `beta`, so
neither half of the existing gate can see it. The ladder was skipped and the
state was refused downstream.

## Decision
1. A split is **physical** only if its residual is at most `tol`, its phase
   fraction is in `(0, 1)`, **and its two phases differ**: some component's
   `ln x` differs by more than `_TRIVIAL_LN_X = 1e-6` (or is present in one
   phase and absent from the other). `_is_trivial_pair` in `flash/_detect.py`.
2. The three ladder gates use it: `_walk_stability_seed_ladder` accepts only a
   physical candidate; `_phi_phi_log_space` and `_phi_phi_second_order` walk
   the ladder whenever the split in hand is not physical.
3. When the ladder finds nothing, the state keeps exactly what it had and
   raises (or proceeds) as before.

## Alternatives considered
- Skip the linear-iterate retry when the K-loop diverged (rejected: 16 of the
  135 neighbours converge correctly through it).
- A tighter solver tolerance or more iterations (rejected: the trivial
  solution is an exact root; no tolerance excludes it).
- Reject trivial splits at the final check and raise (rejected: that is the
  current refusal).

## Consequences
- The 135-pressure neighbourhood of the grid point converges everywhere to
  the pinned tie line (29 refusals -> 0; Linux x86_64).
- Not dormant by construction: a state whose linear stage converged on the
  trivial split *and was then rescued by phase addition* now takes the ladder
  instead. In the scan that is 2 of 135 (answers 1.2e-13 and 6.5e-13 from the
  pin before, 1.1e-12 after, same tie line). The full 2505-state map, run
  before and after on one machine (Linux x86_64) and diffed state by state:
  **0 states differ in any field** (timing excluded); the 155
  `refactor_bit_identity_v3.json` states and the PR stability grid are bit
  for bit unchanged, so the fixture is not regenerated (ledger Case P-17
  "ADR-0036").

- One test changed meaning, and is kept with its original check: the
  synthetic `test_a_converged_vapor_fraction_outside_the_unit_interval_raises`
  (K-loop replaced by a collapsed pair at `beta = -0.2`) is now rescued by the
  ladder to the true methane / ethane VLE (4e-10 from the unpatched flash), so
  it pins the refusal message with the ladder also disabled, and a new test
  pins the rescue.

## Supersedes (optional)
Amends ADR-0028's acceptance rule for ladder candidates.

## Superseded by (optional)
None.
