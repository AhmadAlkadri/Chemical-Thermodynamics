# ADR-0035: On a tpd tie, `stability_tp` reports the better-converged trial

Status: accepted
Date: 2026-09-25

## Context
`stability_tp` reports the trial with the lowest `tpd`. Several trials often
stop at the *same* stationary point; which of them is lowest is then decided
by the last bit, and the tied trials can stop at very different residuals.
Measured (ADR-0032 diagnosis, ledger Case S-9): Tessier (2000) 1-propanol /
water, `z = (0.4, 0.6)`, 361 K - three trials tie (`tpd` equal to 1e-17), one
stopped at stationarity residual 3.3e-11 (inside `tol = 1e-10`), the others
at 1e-16. In one component order the loose trial was reported, so the
reported `trial_composition` differed from the other order's by 8.4e-11, and
by machine.

## Decision
1. After choosing the lowest-`tpd` trial as before, collect the trials whose
   `tpd` lies within `1e-12 * max(1, |tpd|)` of it (`_TIE_TPD`; measured ties
   are <= 5e-16). If the tied trial with the smallest residual has a residual
   at least **1000x** smaller (`_TIE_RESIDUAL_ADVANTAGE`), report it instead.
   The factor makes the rule itself noise-proof: rounding moves residuals by
   O(1) factors, and 1000x is a different stopping point of the iteration.
2. When the rule acts, `diagnostics["minimizing_trial_tie_break"] = "residual"`;
   otherwise the key is absent, so every result the rule does not touch
   carries exactly the mapping it carried before.
3. `tpd_min` is then the reported trial's `tpd`, within the tie window of the
   smallest one; no status can change (the window is 1e-4 of `tpd_tol`).

## Alternatives considered
- Always report the smallest-residual tied trial (rejected: on ties whose
   residuals are all at the noise floor it trades one last-bit choice for
   another and moves captured values on the PR grid).
- Re-polish the reported trial with extra Newton steps (rejected: new
   iterations on every stability call, and it moves every captured number).
- Leave it and keep the test's tie branch (rejected: the reported composition
   is part of the public result and seeds the flash; it should carry the best
   convergence the trial set reached).

## Consequences
- Dormant everywhere measured, bit for bit on one machine (Linux x86_64):
  all 155 `refactor_bit_identity_v3.json` flash states, the 144-state PR
  stability grid (49 states with tied trials, 0 switched) and every state of
  `robustness --quick`. The fixture is therefore not regenerated; the macOS
  exact check (ADR-0032) re-confirms this on the capture platform.
- The Tessier case now reports a 1e-16-converged trial in both orders; the
  reported compositions agree to 4.4e-16 (Case S-9).

## Supersedes (optional)
None. Discharges continuation-plan item A3.

## Superseded by (optional)
None.
