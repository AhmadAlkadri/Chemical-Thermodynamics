# ADR-0016: Negative-flash Rachford-Rice and a second-order stage for the phi-phi split

Status: accepted
Date: 2026-09-13

## Context
The phi-phi split of `flash_tp` is successive substitution on `K`, with the
vapor fraction re-solved from Rachford-Rice after every `K` update
(`chemthermo.flash._split._solve_k_loop`). Until this ADR that inner solve
searched `[0, 1]` only, and reported failure when

    f(beta) = sum_i z_i (K_i - 1) / (1 + beta (K_i - 1))

did not change sign between `beta = 0` and `beta = 1`. The loop then raised
`ConvergenceError("Rachford-Rice failed to bracket a vapor fraction.")`.

That is a failure of the *iterate*, not of the *state*. Measured case, PC-SAFT
(`PCSAFTEOS()`), carbon dioxide / n-decane, `z = (0.9, 0.1)` at 240 K and
1.0 MPa:

- the tangent-plane test is unambiguous: `tpd_min = -5.668e-02`, feed branch
  liquid, incipient composition `w = (1 - 4.2e-08, 4.2e-08)`;
- the seeded `K` brackets properly (`f(0) = +0.058`, `f(1) = -2.2e+05`), so the
  split starts;
- successive substitution then *oscillates*: `K_CO2` swings above and below 1
  with growing amplitude (1.176, 1.095, 1.196, 1.073, 1.225, 1.042, 1.258,
  1.009, 0.991, ...). At iterations 1, 3, 5 and 7 the root of `f` is negative;
  at iteration 9 every `K` is below 1 and there is no root at all;
- the equilibrium exists and is unremarkable: an independent damped Newton on
  the equal-fugacity system in vapor mole numbers (residual `6.7e-13`) gives
  `beta = 0.1835824343`, `x = (0.87751368, 0.12248632)`,
  `y = (1 - 7.1e-08, 7.07e-08)`, `dG_split/RT = -5.5909e-03`.

So a state the package can *prove* is two-phase, and whose answer is neither
near-critical nor delicate, was returned as an error because one intermediate
`K` set had no in-window root. Four states of the Case F-4 grid failed this
way, all on the CO2 / n-decane binary, on both the tangent-plane and the legacy
path.

Two separate things are wrong. The `beta` update is *over-constrained* (the
physical range of the answer is being imposed on every iterate), and the
first-order iteration has *no fallback* when it oscillates - unlike the
liquid-liquid and modified-Raoult splits, which have had the ADR-0009
second-order stage since that ADR and which deliberately left phi-phi alone so
no phi-phi number would move.

## Decision

1. **The two-phase Rachford-Rice update may leave `[0, 1]` during iteration:
   solve on the Leibovici-Neoschil window instead** (the "negative flash" of
   Whitson & Michelsen, *Fluid Phase Equilibria* **53** (1989) 51-71).
   Non-negative phase compositions require every denominator
   `t_i = 1 + beta (K_i - 1)` to be positive, which is

       1 / (1 - K_max)  <  beta  <  1 / (1 - K_min)

   whenever `K_max > 1 > K_min` - the two-phase case of the `t_i > 0` region
   already discussed in `chemthermo.flash._multiphase_rr`. `f` is strictly
   decreasing there (`f' = -sum_i z_i (K_i - 1)^2 / t_i^2 < 0`), so the window
   contains **at most one** root, and exactly one when both endpoints are poles.
   `chemthermo.flash._split._extended_rachford_rice` finds it by safeguarded
   Newton inside a maintained bracket. When all `K` fall on one side of 1 there
   is no root and no window; that is the existing single-phase signal and is
   reported as such.

2. **Bit-identity is by construction, not by measurement.**
   `_extended_rachford_rice` *calls* the unchanged `_rachford_rice` first and
   returns its answer verbatim whenever it has one. The extended window is
   reached only where the old code raised. Every `(z, K)` on which the
   pre-ADR-0016 loop made progress therefore produces the same double, to the
   last bit, and no state that converged before this slice can move. This was
   chosen over reusing `_multiphase_rachford_rice` at `NP = 2`, which is a
   different algorithm (a constrained Newton minimization from a vertex-mean
   start) and would have re-derived every in-window root by a different
   arithmetic path.

3. **The window bounds are taken over all components, including those absent
   from the feed.** A component with `z_i = 0` contributes nothing to `f`, so
   its bound can only shrink the window and can never introduce a spurious
   root. What it buys is that every composition the loop forms inside the
   window is non-negative componentwise, with no `z_i = 0` special case in
   `_solve_k_loop`.

4. **The ADR-0009 second-order stage now finishes an unconverged phi-phi
   split.** It is the *same* function
   (`chemthermo.flash._second_order._second_order_split`), a damped Newton
   minimization of the two-phase Gibbs energy in the phase-II mole numbers; it
   takes the phases' tangent-plane terms as callables, and the phi-phi path
   supplies `ln phi` on the liquid root branch for phase I and on the vapor
   root branch for phase II - exactly the branches successive substitution
   evaluated. The derivation in that module is unchanged: its equation (2) is
   Gibbs-Duhem at fixed `T, P` applied phase by phase, which `ln phi` satisfies
   for the same reason `ln gamma` does. The stage is what makes the minimum,
   rather than any stationary point, the thing being converged to - which is
   why the oscillating iterate above lands on the equilibrium and not on the
   trivial solution.

5. **The stage runs only where the old path would have failed.** Successive
   substitution is given the **full `settings.max_iter` budget it had before
   this slice** and the stage runs only if it did not converge in it, or if an
   updated `K` had no admissible vapor fraction at all (all `K` on one side of
   1), in which case there is no next iterate to form and the loop stops there.
   `FlashSettings.ssi_iterations` is deliberately **not** consulted on the
   phi-phi path, unlike the liquid-liquid and modified-Raoult paths which hand
   over at it: a handover at 50 iterations would change every currently
   converging state that needs more, and decision 2's guarantee is worth more
   than the symmetry. `FlashSettings(second_order=False)` reproduces the
   pre-ADR-0016 first-order split.

6. **The extra diagnostics keys appear only when the stage actually ran.**
   A phi-phi result that converged in the first stage carries exactly the
   mapping it carried before this slice; one that went through the stage
   additionally carries `ssi_iterations`, `second_order_iterations`,
   `converged_stage` and `negative_flash_steps`. This asymmetry is deliberate
   and is the same trade ADR-0011 made: `tests/test_flash_refactor_bit_identity.py`
   compares whole `diagnostics` mappings, so an unconditional key is a pinned
   number changing. A caller must therefore use `.get()` for these four keys.

7. **A converged split whose vapor fraction is outside `(0, 1)` is an error,
   not a single-phase answer.** A negative flash is legitimate *during*
   iteration; as a result it says the phase set collapses to one phase, which
   contradicts the tangent-plane verdict that started the split. The phi-phi
   path raises `ConvergenceError` naming `beta` and `tpd_min`. (Contrast
   `chemthermo.flash._multiphase`, where a converged `beta_j <= 0` *is* the
   answer - it is the phase-removal signal for phases three and above. The
   difference is that there the phase count is being searched; here it was
   decided by a stability test.)

8. **The legacy `phase_detection="wilson-heuristic"` path is untouched.** It
   still raises "Rachford-Rice failed to bracket a vapor fraction." on the four
   Case F-4 states, because reproducing pre-ADR-0008 behavior is that path's
   entire purpose (ADR-0008, ADR-0009). The gamma-gamma and modified-Raoult
   splits also keep the in-window solver: they already hand over to the
   second-order stage at `ssi_iterations`, no state in this repository fails
   their Rachford-Rice, and the change would be untested breadth.

## Alternatives considered
- **Damp the `K` update** (`FlashSettings(damping=...)`). It does suppress this
  oscillation, but the damping factor is a tuning knob with no thermodynamic
  content, it slows every well-behaved state that shares the setting, and it
  only postpones the failure: a sufficiently stiff feed oscillates at any fixed
  damping. Rejected as a fix; still available as a user setting.
- **Retry the split from the Wilson `K` when the stability seed's iteration
  fails.** Already the fallback for a seed that cannot be bracketed *at the
  start* (`_flash_tp_tangent_plane`), and it does not help here: the seed is
  fine, and the Wilson start walks into the same oscillation.
- **Raise a better error.** Honest, and it was the status quo's only virtue.
  Rejected because the state is two-phase, is provably so, and its answer is
  reachable by a method the package already owns.
- **Reuse `_multiphase_rachford_rice` at `NP = 2`.** Rejected under decision 2:
  it would re-derive in-window roots by different arithmetic and put the
  bit-identity fixture at risk for no gain.
- **Hand over to the second order stage at `ssi_iterations`, as the
  liquid-liquid path does.** Rejected under decision 5.

## Consequences
- Positive: the four Case F-4 failures become verified two-phase results. At
  240 K / 1.0 MPa the returned split matches the independent Newton to
  `3.6e-13` in `beta` and `5.5e-14` in every mole fraction, with
  `fugacity_residual = 5.1e-13`, `mass_balance_residual = 0`,
  `delta_g_split_rt = -5.5909e-03` and both phases stable.
- Positive: over the 188-state Case F-4 PC-SAFT grid the `ConvergenceError`
  count on the tangent-plane path goes 4 -> 0, and all 184 states that
  converged before are bit-identical.
- Positive: over a 1248-state Peng-Robinson grid, 1247 states are bit-identical
  and the one state that previously exhausted the iteration limit
  (methane/ethane/propane at 290 K, 8 MPa - the weakly unstable near-critical
  state recorded as open in `brain.md`) now converges through the stage to a
  `fugacity_residual` of `1.8e-15` with `dG/RT = -1.5e-05`. That open item is
  discharged.
- Tradeoff: the four keys of decision 6 are conditional, so `diagnostics` is no
  longer a fixed key set on the phi-phi path. Documented in `flash_tp`.
- Tradeoff: a starved iteration budget is no longer a reliable way to *force* a
  phi-phi non-convergence, because the stage rescues it. Three tests that used
  `max_iter=1` as a non-convergence fixture now pass `second_order=False` (or
  use the legacy path) and one new test pins the rescue.
- Tradeoff: `FlashSettings(second_order=False)` disables the *stage* but not
  the extended window, which is unconditional on the phi-phi tangent-plane
  path. `phase_detection="wilson-heuristic"` remains the full pre-slice
  reproduction.
- Unchanged: the legacy path, gamma-phi, gamma-gamma, modified-Raoult, the
  multiphase path, every stability solver, and every pinned number in the
  suite.

## Next slice
This is also the mechanism a future **`flash-phase-addition-eos`** slice needs.
Wiring the ADR-0011 addition/removal search to the phi-phi path means solving
phase sets that do *not* exist at the feed, and reading the sign of a converged
phase fraction as the removal signal; both require an inner solve that treats a
fraction outside `[0, 1]` as information rather than as a failure. Decision 1
supplies that for `NP = 2` and `_multiphase_rr` already supplies it for
`NP > 2`, so what remains for that slice is the search itself - and, per
ADR-0011 "What remains" and ADR-0015, a real three-phase EOS state to exercise
it. The recommended intermediate step is `flash-phase-labels-by-density`
(ADR-0015 decision 2), which would let a phi-phi result name its phases from
the density root rather than from a Wilson volatility ranking.

## Supersedes (optional)
Amends ADR-0009's "the phi-phi and gamma-phi splits never call it" (the
second-order stage) to "the gamma-phi split never calls it"; the phi-phi split
now does, under decision 5. ADR-0009 is otherwise unchanged, and gamma-phi -
which has no stability test (ADR-0007) and stays on the legacy path - is
unchanged.

## Superseded by (optional)
None.
