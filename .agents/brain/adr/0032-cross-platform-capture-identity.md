# ADR-0032: Captured-value guards are exact on their capture platform and bounded elsewhere

Status: accepted
Date: 2026-09-25

## Context
Three kinds of identity guard protect refactors and speed-ups here
(ADR-0017/0021 fixture history, ADR-0023, ADR-0030):

1. **In-process A/B guards**: the same state computed twice, two ways, in one
   run (memo on/off, branch reuse, log-space vs linear). Platform-independent.
2. **Captured-value guards**: a run compared with floats *stored* from one
   machine - `tests/fixtures/flash/refactor_bit_identity_v3.json` (155 flash
   states), the PC-SAFT literals of ledger Cases P-1/P-3, and the per-surface
   count of Case P-11.
3. **Result hashes** of benchmark records (ADR-0023), which the dev contract
   already says to compare only on one machine.

Every capture of kind 2 was made on macOS arm64. The branch's first CI test run
(GitHub Actions `ubuntu-latest`, run 36121300741) failed exactly those three
guards, 3 of 769 tests, and nothing else. Measured on a second Linux x86_64 host
(the 2026-09-25 cloud baseline, CPython 3.11.15, numpy 2.4.6):

- flash fixture: deterministic across two runs; **83 of 155 states** move,
  **688 floats**, **zero discrete fields** (no phase name, key, status,
  verdict, stage or iteration count). Worst move 2.44e-13 absolute on an O(1)
  quantity; compositions and fractions within ~1e-13; the largest *relative*
  moves are on residuals (`max_delta_k`, `mass_balance_residual`, near-zero
  post-split `tpd`), which are differences of O(1) numbers and carry absolute,
  not relative, error.
- PC-SAFT n-hexane literals: `Z` at 7700 mol/m^3 moves 32 ULP (3.5e-15), `ln phi`
  2 ULP; `a_res` is exact here but 3 ULP off on the CI runner. Two Linux x86_64
  machines disagree with each other, not only with macOS: `exp`/`log`/`pow` and
  numpy's CPU-dispatched SIMD kernels are not bit-reproducible across CPUs or C
  libraries.
- Case P-11 surface count: verdicts identical (47 unstable / 97 stable);
  `minimizing_trial_surface` liquid/vapour is 45/32 on macOS, 46/31 on CI,
  47/30 here. On 20 of the 77 states both surfaces reach the **same**
  stationary point (best `tpd` equal to 5e-16, compositions to 6e-13), so the
  "winner" is the last bit.

## Decision
1. **Kind 1 and kind 3 are unchanged**: exact everywhere, and hashes compared on
   one machine.
2. **Kind 2 is exact on the capture platform** (`platform.system() == "Darwin"`
   and `platform.machine() == "arm64"`, `tests/_capture_identity.py`): the
   comparison there is the same `==` as before.
3. **Elsewhere, every discrete field stays exact** (dict keys, sequence lengths,
   strings, ints, bools, `None`, non-finite floats) and a finite float passes
   when `|a - e| <= atol + rtol * max(|a|, |e|)`, with the bound stated at the
   call site and justified against measurement and the solver tolerances:
   - flash fixture: `atol = rtol = 1e-12`. All pinned floats are O(1)-scaled
     (fractions, K-values, `tpd` and `dG` in RT, residuals of those); 4x the
     worst measured move; equal to the tightest solver tolerance
     (`FlashSettings.second_order_tol`) and 1000x tighter than the 1e-09 audit
     tolerance ADR-0021 used when this fixture was regenerated. **Phase
     fractions** (and `vapor_fraction`) get `atol = 1e-12 / delta`, `delta`
     the state's shortest tie line in the fixture (the smallest, over phase
     pairs, of the largest composition difference): the lever rule divides
     composition noise by the tie line's length. Derived after the first CI
     run of this rule (run 36124342571) moved the Tessier near-plait feed's
     fractions by 2.34e-12 (`delta = 0.056`; its bound is 1.8e-11); every
     other multiphase state has `delta >= 0.30`, so its bound stays within
     3.4e-12.
   - PC-SAFT literals: `atol = 5e-14`, `rtol = 0`. 14x the worst measured move;
     the size of the floor at which Case P-1 accepted the same quantities as
     equal to teqp (worst `|dZ|` 2.58e-14, `|d ln phi|` 2.66e-14). That the
     association code adds nothing to a non-associating model stays pinned
     exactly, in-process, on every platform.
4. **A tie is not a result.** A per-surface count is asserted on every platform
   only for states whose winning surface beats the other by more than
   `TIE_MARGIN = 1e-10` (19 vapour, 38 liquid) plus the count of ties (20);
   a decisive winner must equal the reported `minimizing_trial_surface`. The
   margin is 5 orders above the largest tie and 8 below the smallest real gap
   (1.7e-02), and 100x below `tpd_tol`. The exact 45/32 split stays pinned on
   macOS arm64.
5. **A second exact capture per platform is not added.** Two x86_64 Linux hosts
   already disagree, so a Linux capture would make CI's pass depend on which
   runner CPU it is given.
6. Guards added later follow the same rule; a new captured-value guard records
   its capture platform. A guard that passes exact on every measured platform
   is left exact (for example `test_non_associating_density_roots_are_bit_identical`).

## Alternatives considered
- Relax the fixtures to a uniform tolerance everywhere (rejected: gives up the
  bit-identity proof on the machine where refactors are developed).
- Mark the three tests `skipif(not darwin)` or `xfail` on Linux (rejected: CI
  would then check none of the discrete fields or the verdicts).
- Regenerate the fixtures on Linux (rejected: moves the problem to macOS and,
  per the two-host measurement, would still be runner-dependent).
- ULP bounds (rejected for residuals and near-zero `tpd`, where a ULP of the
  value says nothing about the error, which is set by the O(1) quantities the
  residual is a difference of).

## Consequences
- CI on Linux checks the same 155 states, every discrete field exactly, and
  floats to 1e-12; a numerical change the size of a solver tolerance fails on
  every platform. A change smaller than 1e-12 is only caught on macOS arm64,
  which is recorded here rather than hidden.
- A refactor developed on Linux alone gets the bounded guard, not the
  bit-for-bit one; claims of "bit-identical" must name the platform they were
  measured on.
- Evidence: ledger Case P-11, "Cross-platform (ADR-0032)"; negative controls
  for the comparison itself in `tests/test_capture_identity.py`.

## Supersedes (optional)
Amends the bit-identity statements of ADR-0017, ADR-0021, ADR-0023 and ADR-0030
for captured values only; their in-process guards are unchanged.

## Superseded by (optional)
None.
