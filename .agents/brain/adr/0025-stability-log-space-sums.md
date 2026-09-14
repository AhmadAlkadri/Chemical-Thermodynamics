# ADR-0025: Michelsen's trial mole numbers are summed in logarithms, not clamped

Status: accepted
Date: 2026-09-13

## Context

ADR-0024 gave the two-phase split a log-space stage and retired the
`Mw = 16 400` polymer/solvent VLE defects. It left one thing pinned as still
broken, and said in so many words that it was **not** a split problem:
`Mw = 53 000` polyethylene (`m = 1393.9`) in n-pentane at 453 K and 0.5 or
1 MPa. There, `stability_tp` reported the feed unstable with
`tpd_min = -1.05e-04` / `-3.38e-04` - a shallow *vapour-side* stationary point -
and never found the melt, so ADR-0024's stage was seeded 1300 orders of
magnitude from the answer and spent its budget. Case P-14 recorded it as "a
**stability trial set** limitation for `m = 1393.9`".

That diagnosis was half right. The trial set was fine. Three of its four trials
- `wilson-liquid` and both `pure-<name>` estimates - were walking straight at
the melt and could not get there.

**What was actually in the way.** Michelsen's iteration carries the
*unnormalized* mole numbers `W`, and equation (8) produces them as logarithms:
`ln W_i <- d_i - ln phi_i(w)`. The only place `W` was ever needed as a double
is the normalization `w = W / sum_j W_j`. `chemthermo/stability/tp.py` guarded
that normalization with a clamp,

```python
ln_w_new[active] = np.clip(d[active] - ln_f[active], -700.0, 700.0)
```

and three more `np.clip(u, -700, 700)` in the Newton stage - `-700` and `+700`
being, to two figures, where `exp` stops existing as a double.

A melt of a 1393.9-segment chain against a solvent-vapour feed has
`ln phi_polymer = -1549.71` and `d_polymer = -97.50`, so its stationary point
sits at `ln W_polymer = 1452.21` (`ln W_solvent = 3.62`). The clamp put it 752
out of reach. What was measured before this slice: all three liquid-surface
trials ran out their 50 successive substitutions, entered the Newton stage,
and ended `second_order_no_progress` after 51 iterations parked at
`w = (1, 3.7e-303)` with a stationarity residual of **752.2** - which is
exactly `1452.2 - 700`, the distance the clamp was holding them back.

So the clamp was not protecting a result. It was keeping `exp` in range, and by
doing so it made a whole family of stationary points *unreachable* rather than
merely inaccurate. `logsumexp` removes the need for it.

Three facts decided the shape of the change:

1. **The verdict never needed `sum_W` as a double.** Equation (7) is
   `tpd = -ln(sum_i W_i)` and the sign test is `sum_W > 1 <=> tpd < 0`. Both
   are statements about the *logarithm*. Only `tm* = 1 - sum_W`, a reported
   diagnostic, wants the sum itself, and at `ln sum_W = 1452` it is `-inf`.
2. **The normalized `w` cannot carry such a point either.** At the melt,
   `w_solvent = exp(3.62 - 1452.21) = exp(-1448.6)`, which is an exact `0.0`.
   That is the honest answer for a mole fraction (ADR-0024 decision 3 already
   settled the convention), but it means a consumer that rebuilds `ln W` from
   `w` gets `log(0)`. The magnitude has to travel separately.
3. **Everything else in the repository stays inside `[-700, 700]`.** Measured:
   0 of 624 trials over the 144-state Peng-Robinson stability grid, 0 of 12
   over the three `Mw = 16 400` states of Case P-14 (worst `|ln W| = 444.2`),
   and none on the Peng-Robinson, PC-SAFT water/n-hexane or NRTL controls.

## Decision

### 1. The normalization is done in logs where, and only where, the clamp would have engaged

`chemthermo.stability.tp._normalize` is the single place `w`, `sum_W` and
`ln sum_W` are formed, and it has two branches:

- every `ln W_i` inside `[-700, 700]`: form `W = exp(ln W)`, sum it, divide.
  This is the pre-ADR-0025 expression character for character, and it runs
  wherever the old `np.clip` returned its argument untouched - which is to say
  wherever the old clamp was dormant. **Every result that existed before this
  decision record is bit-identical by construction**, not by measurement: the
  same operations in the same order on the same doubles, and `ln sum_W` is
  `math.log(sum_W)` so that `tpd_from_sum_W` is the same double it was.
- otherwise: `ln sum_W = logsumexp(ln W)` and `w_i = exp(ln W_i - ln sum_W)`.
  `w_i` may underflow to an exact `0.0` and `sum_W` may overflow to `inf` or
  underflow to `0.0`; `ln sum_W` is the one that is always meaningful.

The four clamps are gone. The Newton stage's variables were already `u = ln W`,
so removing them there is removing a box from a variable that never needed one;
its line search keeps a **finiteness** test (reject an `inf`/`nan` candidate
and halve the scale) in place of the box, which is what the clamp was doing for
such a candidate anyway.

The gate is per-*evaluation*, not per-trial: a trial that leaves the window and
comes back runs the old expressions again on the way back. That is the most
conservative reading of "only where a clip would have engaged", and it is what
makes the dormancy claim checkable one evaluation at a time.

### 2. `ln W` is reported, and it is what a log-space consumer is seeded from

`StabilityTrial` gains `ln_W` (unnormalized, at the final point), `ln_sum_W`
and `log_space`; `StabilityResult` gains `trial_ln_W` (the minimizer's) and the
diagnostics keys `ln_sum_W`, `minimizing_trial_log_space`,
`log_space_trial_count` and `log_space_trials`. The last three are written
**only when the route engaged**, so their absence is the statement that
nothing left the old arithmetic.

`tm_at_stationary_point` is reported as `-inf` where `sum_W` has overflowed,
rather than omitted: the sign is the verdict and the magnitude is in
`tpd_from_sum_W`.

ADR-0024's seed (`flash/_log_space.py::log_space_seed`, and the K-seed in
`flash/_detect.py`) rebuilt Michelsen's `W` from the normalized `w` as
`ln W = ln w - tpd`. That is exact where `w` is representable and meaningless
where it is not, so both now take an optional `ln_capital_w`, and `_detect`
passes it **only** when the minimizing trial engaged log space. Where it is not
passed the pre-ADR-0025 expression runs unchanged, so no converged split is
re-seeded.

Measured, because the equality is load-bearing: over the 47 unstable verdicts
of the 144-state Peng-Robinson grid the two roads agree to **5.4e-15**. Over
*stable* verdicts they can differ by whole units - up to 4.44 on that grid -
and that is not new and not a defect: `tpd_min` is measured to the lower
envelope of the phase candidates while `ln sum_W` is equation (7) on the
surface the trial iterated on, and they part company exactly where a pinned
trial stops above the other candidate (ADR-0021, and the note in `_summarize`).
A stable verdict seeds nothing, so nothing reads either number there. The gap
is recorded rather than hidden, so that widening the gate would have to
confront it.

### 3. The physics is untouched

No model equation, no split solver, no trial set, no initial estimate, no
convergence criterion and no tolerance changed. The trial count is still
`n + 2`; the surfaces are still ADR-0021's; the stationarity residual is still
measured on the trial's own surface and the tangent-plane distance still on the
lowest-Gibbs candidate. This decision record is about one normalization.

## Consequences

**What it buys.** `Mw = 53 000` polyethylene / n-pentane at 453 K, 5 wt%:

| P / MPa | `tpd_min` before | `tpd_min` after | melt `ln W` after | trials |
| --- | --- | --- | --- | --- |
| 0.5 | -1.0518e-04 | **-1452.2068** | `(1452.2068, 3.6168)` | 3 SSI |
| 1.0 | -3.3766e-04 | **-1346.5734** | `(1346.5734, 4.2287)` | 3 SSI |

The three liquid-surface trials converge in **3 successive substitutions**
where they used to spend 51 iterations failing, and `flash_tp` returns the
vapour-liquid split ADR-0024's stage was written for: melt `x_solvent`
0.9788121097 / 0.9909016024 (5.92 / 12.91 wt% solvent), `beta_vapor`
0.9966188537 / 0.9921261568, `ln y_polymer` -1507.98 / -1452.66,
`fugacity_residual` 1.3e-13 / 3.3e-14, `mass_balance_residual` 1.1e-16,
`delta_g_split_rt` -1.071e-01 / -1.029e-01, both phases post-split **stable**.
The melt agrees with an independently written one-dimensional equal-fugacity
solve to **1.2e-14** / 3.6e-15 and with FeOs's chemical potentials to
**1.1e-12** with matched constants (8.6e-08 as shipped). Case P-14's pinned
miss is retired; Case P-15 records the reversal.

It is wider than the two pinned states. A 0.3-3.6 MPa sweep at 0.1 MPa steps,
run against this commit's parent and against this one, on both molar masses -
68 states - says:

| | states |
| --- | --- |
| identical, bit for bit | **48** (every one of the 34 `Mw = 16 400` states) |
| `ConvergenceError` -> a verified split | **13** (0.4-1.4 MPa vapour-liquid; 2.6 and 2.7 MPa liquid-liquid) |
| converged before, moved in the last bits | **7** (1.5-2.1 MPa) |
| converged before, raises now | **0** |

The two liquid-liquid repairs at 2.6 and 2.7 MPa are decision 2 rather than
decision 1: their *verdict* is unchanged, but the minimizing trial reaches
`ln W_polymer = -1187`, so the seed is now built from `ln W` instead of from a
`w` that had rounded to zero.

**What moved that had an answer before.** Seven states, all the same system
between 1.5 and 2.1 MPa, and all of them states where the old clamp engaged
*and* a result still came out. They used to be seeded from the shallow
vapour-side stationary point and are now seeded from the melt, and they
converge on the same answer from a different direction: worst **3.48e-13
relative** on the melt composition (2.0 MPa: `x_polymer`
2.7736965292611865e-03 -> 2.7736965292621523e-03) and **9.23e-15** on the phase
fraction. The 2.0 MPa melt composition is pinned at `rel=1e-9` in
`tests/test_flash_log_space_stage.py` and still passes. Nothing else in the
repository moved: `refactor_bit_identity_v3.json` (155 states) passes
**unchanged** and was not regenerated, all nine ADR-0023 benchmark result
hashes are **identical**, the whole suite passes with one test reversed rather
than adjusted, and the 144-state Peng-Robinson grid engages the new route 0
times in 624 trials.

**The cost.** One branch in one function, three new fields on `StabilityTrial`,
one on `StabilityResult`, four conditional diagnostics keys, and an optional
argument on two seed builders. `logsumexp` costs one extra pass over an
`n`-vector, and only on the branch that needs it.

**What remains.**

1. **The multiphase stability path has no log-space form.** `flash/_multiphase.py`
   runs its own post-split stability tests through the same `stability_tp`, so
   it inherits this; what it does *not* have is ADR-0024's log-space split, so
   a three-phase state whose compositions leave machine range would still fail.
   No state in this repository needs it, so it is not written (ADR-0002).
2. **The Newton stage's conditioning at `|ln W| ~ 1e3` is untested.** Every
   state this slice repairs converges during successive substitution, in three
   steps, so the second-order stage is never entered there. Its central
   difference step is `1e-6` in `u`, which is a *multiplicative* step on `W`
   and so is scale-free - the same argument ADR-0024 made for its own Hessian -
   but that is an argument, not a measurement. A state that needs the Newton
   stage *and* sits outside the window has not been found.
3. **Six states of this system still raise, and this slice did not touch
   them.** 0.3 MPa and 2.8, 2.9, 3.0, 3.1, 3.2 MPa raised before and raise
   now on the 0.1 MPa sweep, with the same
   message: the tangent-plane verdict is a shallow near-critical one and the
   split does not converge from it. That is the ADR-0024 stage's own budget,
   not a stability miss - `ln W` there is inside the window - so it is recorded
   rather than reached for.
4. **The far end of the window is exercised but not stressed.** `ln W` below
   `-700` no longer clamps either, which is the right thing (`sum_W = exp(-799)`
   underflows to `0.0` and the composition is still exact in logs). The states
   that reach it here - the same polymer at 2.2-2.7 MPa, whose minimizing trial
   reaches `ln W_polymer = -1282` to `-1187` - return an unchanged `tpd_min`,
   because a component at `exp(-1282)` and one at `exp(-700)` are both an exact
   zero as a mole fraction. Only the *seed* notices, which is what repairs 2.6
   and 2.7 MPa.
5. **The verdict is still not a global proof.** This makes one family of
   stationary points reachable; it says nothing about stationary points no
   initial estimate approaches. The honesty note on `StabilityResult` stands.

## Alternatives considered

- **Widen the clamp to `[-5000, 5000]`.** Rejected: `exp(5000)` is still `inf`,
  so the normalization would return `nan` rather than a wrong answer, and every
  number in the window would have had to be re-audited. The clamp's value was
  never its width.
- **Always normalize in logs.** Rejected: `logsumexp` and `sum(exp(.))/x` differ
  in the last bits, so every pinned number, the 155-state fixture and the nine
  benchmark hashes would have had to be regenerated for a change that alters no
  physics. The gate buys bit-identity by construction, which is the same
  argument ADR-0022's `log_fugacity_coefficients` guard and ADR-0023's
  `ln_fugacity_branches` made.
- **Detect the polymer and branch on it.** Rejected on principle: nothing in
  the solver knows what a polymer is, and a magnitude test is the honest
  trigger for a magnitude problem.
- **Seed the split from `ln W` everywhere.** Rejected, with evidence: on stable
  verdicts `ln W` and `ln w - tpd` differ by up to 4.44 over the Peng-Robinson
  grid (decision 2), and while a stable verdict seeds nothing today, making the
  seed depend on which surface a trial stopped on would be a real change with
  no state asking for it.
- **Report `tm_at_stationary_point` as `nan`, or omit it, when `sum_W`
  overflows.** Rejected: `-inf` is the correct limit and carries the sign,
  which is the part of `tm*` that means anything.

## References

- M. L. Michelsen, "The isothermal flash problem. Part I. Stability",
  *Fluid Phase Equilibria* **9** (1982) 1-19 - equations (3) to (7) of
  `chemthermo/stability/tp.py`'s module docstring.
- ADR-0022 (`log_fugacity_coefficients`, the same "guard, not a rewrite"
  pattern one level down), ADR-0024 (the log-space split stage this feeds),
  ADR-0021 (pinned trial surfaces; the source of the `tpd_min` versus
  `ln sum_W` distinction in decision 2).
- Validation Case P-15 (`.agents/brain/validation-cases.md`), and the
  amendment to Case P-14.
- Golden path: `python examples/validation/22_stability_log_space.py`
  (`--full` adds the 144-state grid and the 0.4-2 MPa scan).
