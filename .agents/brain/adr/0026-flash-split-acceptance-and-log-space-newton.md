# ADR-0026: The split's last two refusals are a curvature safeguard and a denominator

Status: accepted
Date: 2026-09-14

## Context

ADR-0024 gave the two-phase split a log-space stage; ADR-0025 removed the
`ln W` clamp that had been hiding the polymer melt from the stability test.
Together they turned most of the polyethylene / n-pentane band at 453 K into
verified two-phase answers. They did not turn all of it: a sweep from 0.3 to
3.6 MPa at 0.1 MPa steps over both molar masses - 68 states - had **six**
raising `ConvergenceError` at HEAD `584c508`, all of them on the
`Mw = 53 000` chain (`m = 1393.9`, 5 wt%, `k_ij = -0.006`, feed
`z_polymer = 7.164e-05`). The roadmap carried them as item 1.

The six are two different defects, and neither is about polymers.

### (a) 0.3, 2.8 and 2.9 MPa: the stage crawls next to the trivial solution

ADR-0024's stage takes its step on `r(u) = 0` and keeps it only while it
descends the two-phase Gibbs energy `g`; the documented fallback where it does
not is `-r`, the log-space successive-substitution step, which is a descent
direction for `g` always. The fallback is safe and it is *slow*: its length is
`|r|`, so once the residual is small the step is small, whatever distance is
left to travel.

What was measured at 0.3 MPa, from the shipped seed: the Newton direction is
rejected at **every one** of the 100 iterations, and the stage spends its whole
budget moving `ln n` of the solvent from `-3.1804` to `-3.1765` - by `4e-03`,
where the answer is `3.2` away - with the residual flat at `3.9e-05` and
*rising* in the tenth digit. Both phases sit within a percent of the feed: this
is the neighbourhood of the **trivial** solution `x^I = x^II = z`, and the
reason the descent test rejects Newton there is that the Gibbs Hessian in the
mole numbers is **indefinite** - measured eigenvalues `-1.337e-04` and
`+3.967e+06`. The rejection is correct. The replacement is the problem: `-r`
cannot leave a saddle at any useful speed.

Three things were ruled out first, by measurement rather than argument:

- **The finite-difference Jacobian is not the problem.** At the stalled iterate
  every entry agrees to four figures across `h = 1e-06`, `1e-04` and `1e-02`
  (the `h = 1e-08` column shows the expected round-off, 2% in the smallest
  entry), and `dr/dn = J diag(n)^-1` is symmetric - as the Hessian of `g` must
  be - to `1.5e-04` relative at the shipped step and `1.5e-07` at `h = 1e-04`.
- **Damping is not the problem.** The line search accepts a full step
  (`scale = 1.0`, zero backtracks) at every one of the 100 iterations.
- **It is not a scaling accident of the polymer row.** The polymer's residual
  converges in eight iterations and then stops mattering; it is the *solvent*
  row, weighted by `n_solvent = 0.042`, that has nowhere to go.

2.8 and 2.9 MPa are the same defect on a liquid-liquid geometry. At 2.8 MPa the
crawl consumes 97 of the 100 iterations; the Newton direction becomes
admissible at the 98th and the residual then falls `8.58e-02 -> 4.37e-04 ->
1.078e-08` in three steps, so the state runs out of budget roughly one
iteration short of converging. At 2.9 MPa the crawl consumes all of them and
the residual ends at `6.03e+00`.

The near-miss at 2.8 MPa is worth naming because it invites the wrong fix. It
is not a tolerance that is too strict: the acceptance rule already *is*
`FlashSettings.tol` (`_phi_phi_log_space` raises on
`refined.residual > settings.tol`, not on `second_order_tol`), and
`1.078e-08` is above it. Loosening it would have retired one of the six and
shipped an answer whose iteration was still in motion, while 0.3 MPa
(`3.9e-05`) and 2.9 MPa (`6.0e+00`) went on refusing. **No tolerance was
changed by this ADR**; with the repair below, 2.8 MPa converges to `9.09e-13`
in twelve iterations.

### (b) 3.0, 3.1 and 3.2 MPa: a denominator that cancels to zero

These raised

    Feed is unstable (tpd_min=-1.2537e-03) but neither the stability-seeded nor
    the Wilson K-values bracket a Rachford-Rice root.

and the stability seed at 3.0 MPa is `K = (1.11871119e-18, 1.00132622)`, which
straddles one. By hand, `f(0) = +1.254e-03` and `f(1) = -6.4e+13`: a root sits
at `beta = 0.9459`, and both `_rachford_rice` and the extended solver's first
call reported none.

The cause is one floating-point expression. `_rachford_rice` forms its
denominator as `t_i = 1 + v (K_i - 1)`. For any `K_i` below the spacing of
doubles at one, `K_i - 1` rounds to **exactly** `-1.0`, so at `v = 1` the sum
`1 + (-1)` is an exact `0.0`, the positivity guard fires, `f(1)` is reported
`nan`, and the sign test cannot bracket. Nothing about the physics is
degenerate; the equation lost its `K` to a cancellation. The Wilson fallback
loses it the same way - a non-volatile component's Wilson K is `1e-10`
(ADR-0022) - so the second attempt fails for the same reason as the first and
the flash refuses.

`t_i = (1 - v) + v K_i` is the same quantity written as the convex combination
it is: strictly positive for every `K_i > 0` and `0 <= v <= 1`, exact at both
endpoints, no cancellation.

## Decision

### 1. The log-space stage gains a curvature safeguard, used only on a retry

`chemthermo.flash._log_space.log_space_split` takes
`curvature_safeguard: bool = False`. `False` is the ADR-0024 iteration,
character for character - the same Jacobian, the same single direction, the
same `slope >= 0` swap to `-r`, the same line search and the same acceptance
test.

`True` replaces the single `-r` fallback with an **ordered list** of
directions, built by `_safeguarded_directions`, of which the line search tries
each in turn and keeps the first it accepts:

1. the Newton direction on `r = 0`, while it descends `g` - which near the
   solution is the ordinary quadratic step and is unchanged;
2. the **modified-Newton** direction of the Gibbs Hessian itself,
   `H = diag(n) J + diag(n r)` (`J = dr/du` already carries the second
   `diag(n)`), symmetrized, with every eigenvalue replaced by its magnitude and
   floored at `1e-12` of the largest. This is Gill & Murray's modification
   (*Practical Optimization*, sec. 4.4.2): where `H` is positive definite it is
   Newton's own step, and where it is not, a negative-curvature eigenvector is
   followed **downhill** instead of being discarded. Its length is set by the
   curvature, not by `|r|`, which is the whole point;
3. `-r`, exactly as before.

A direction whose slope is not negative is dropped rather than repaired, so
every direction offered is a descent direction and the list always ends with
one.

**The flag is never set on a first attempt.** `_phi_phi_log_space` runs the
stage as it always has; only if the result is still above `FlashSettings.tol`
does it call the stage a *second* time, from the *same seed*, with the
safeguard on, and keep the better of the two. A second call cannot move a
number the first call returned, so bit-identity is by construction rather than
by audit. `diagnostics["log_space_curvature_safeguard"]` records which rule
produced the answer.

### 2. Rachford-Rice gains a cancellation-free denominator, used only where it would otherwise refuse

`_rachford_rice` and `_extended_rachford_rice` take
`convex_denominators: bool = False`, threaded through one helper,
`_rr_denominators(v, K, convex=...)`. It returns the naive `1 + v (K - 1)`
whenever that is admissible - so every value that was ever finite is the same
double - and decides only what happens where it is not: the pre-ADR-0026
`None` (and with it a `nan` and a "no root" verdict), or one recomputation as
`(1 - v) + v K`.

The convex form is **not** the default. It is not the same double as the naive
form for an ordinary `K` below `0.5`, and this module's answers are pinned bit
for bit; switching it on everywhere was tried and measured to re-route 15 of
the 62 states of the sweep that already converged - four of them from the
log-space stage to successive substitution - because a seed that had had no
in-window root suddenly had one. Those answers are as correct as before, but
moving them is not this slice's business.

The one caller that switches it on is
`_flash_tp_tangent_plane`, at the point where it is about to raise: the
stationary point's K-values have been reported rootless, the Wilson fallback
has been reported rootless, and the stationary point's K-values are asked once
more with the convex denominator and on the wider (negative-flash) window.
The seed that comes back is still the stationary point's, so
`diagnostics["k_seed"]` stays `"stability"` and the split keeps the per-phase
branch pinning that goes with it (ADR-0019); what changed is only how the
bracket was found, and that is its own key,
`diagnostics["rachford_rice_convex_denominators"]`.

## Consequences

**All 68 states of the 0.3-3.6 MPa sweep converge, and the other 62 are
bit-identical.** Every field a caller can observe - phase names, both
compositions, phase fractions, `vapor_fraction` and every diagnostics number -
compared with `==` against the same sweep run from a worktree at `584c508`.

The six, each verified by a solve that shares no code with `flash_tp`:

| P / MPa | verdict | route | fugacity residual | independent check |
| --- | --- | --- | --- | --- |
| 0.3 | VLE | `stability-log`, safeguarded | 4.87e-13 | melt `x_C5` to 1.4e-14 by the Case P-14 1-D solve |
| 2.8 | LLE | `stability-log`, safeguarded | 9.09e-13 | tie line to 2.5e-13 relative by a 2-equation Newton |
| 2.9 | LLE | `stability-log`, safeguarded | 9.09e-13 | 7.1e-13 |
| 3.0 | LLE | `stability`, convex bracket | 7.28e-12 | 5.7e-13 |
| 3.1 | LLE | `stability`, convex bracket | 3.87e-12 | 1.0e-12 |
| 3.2 | LLE | `stability`, convex bracket | 2.27e-13 | 2.8e-13 |

Mass balance is below `1.2e-16` at all six, `dG_split/RT < 0` at all six, both
phases post-split stable at all six, and FeOs's chemical potentials at
chemthermo's compositions and densities agree to `6.1e-12` with FeOs's own
dispersion constants substituted in (`3.1e-07` as shipped, the known
ten-versus-fourteen-figure difference of Cases P-6 to P-12), at `k_ij = 0` on
both sides.

The 0.3 MPa iteration, which is the clearest statement of what changed:

| iteration | ADR-0024 rule | with the safeguard |
| --- | --- | --- |
| 1 | 1.598e-01 | 7.917e+01 |
| 4 | 1.054e-03 | 2.597e+02 |
| 8 | 3.899e-05 | 1.683e-01 |
| 10 | 3.899e-05 | 2.095e-04 |
| 11 | 3.899e-05 | 9.642e-09 |
| 12 | 3.899e-05 | 1.137e-12 |
| 100 | 3.901e-05 | 4.865e-13 |

The safeguarded sequence gets *worse* before it gets better - that is what
climbing out of the trivial basin looks like - and is then quadratic.

**Nothing else moved.** `python -m chemthermo.bench --compare` reports every
one of the nine cases' `result_hash` identical; the bit-identity fixture v3
(155 states) is untouched and no state gained or lost a diagnostics key there,
because none of them reaches the log-space stage; `pytest -q` and
`pytest -q -m slow` pass.

**Two new diagnostics keys**, both only where the routes they describe ran:
`log_space_curvature_safeguard` (with the other `log_space_*` keys) and
`rachford_rice_convex_denominators` (only when the last-resort bracket was
used). No public signature changed.

### What remains

- **A narrow band near 5.2 MPa on the `Mw = 53 000` chain still refuses.**
  Measured at 5.175, 5.2 and 5.3 MPa; 5.0, 5.1, 5.4 and 5.5 MPa converge. It
  refuses identically at HEAD `584c508`, with the same message and the same
  residual (`5.8e-04` after 100 successive-substitution and 101 second-order
  iterations), so it is neither caused nor repaired here. It is outside the
  0.3-3.6 MPa band this ADR measures, it enters through the *other* log-space
  entry point (`_phi_phi_second_order`, seeded from a failed linear iterate),
  and the symmetric retry there was implemented and measured **not** to repair
  it - so it was not shipped, because it would have been an unexercised path
  (ADR-0002). `examples/validation/21_pcsaft_polymer_vle.py --full` prints the
  band rather than hiding it.
- **The safeguard is wired at one entry point only.** See above: the second
  one has no state in this repository that it would help.
- **`convex_denominators` is a keyword, not the shape of the module.** The
  honest end state is one denominator everywhere; getting there means
  re-auditing the 15 sweep states it re-routes, plus whatever else it touches
  across the validation grids, which is a slice of its own.

## Alternatives considered

- **Accept the split above `settings.tol` when the budget is exhausted, instead
  of raising.** This is what 2.8 MPa's `1.078e-08` seems to ask for, and it
  would have retired exactly one of the six while leaving 0.3 MPa (`3.9e-05`)
  and 2.9 MPa (`6.0e+00`) refusing. Rejected for three reasons: the acceptance
  rule already *is* `tol` and `1.078e-08` is above it, so this is a request to
  loosen a documented tolerance for one state; the number is not a converged
  residual but an iteration caught in motion (it had fallen four orders of
  magnitude in the previous step); and the stall it papers over is the same one
  that makes the other two states hopeless. The tolerances in
  `FlashSettings` are unchanged.
- **Take the Newton step on `r = 0` regardless of descent, with a line search
  on `|r|` alone.** Simpler, and wrong: measured from the shipped seed it
  converges to the **trivial** solution at 2.8 and 2.9 MPa (`beta` driven to
  `0.0`, `u` to `-8e+04`) and stalls at `8.5e+01` at 0.4 MPa, a state that
  converges today. `r = 0` has a whole line of trivial roots; the descent test
  on `g` is what keeps the stage away from them, and it has to stay.
- **Change `_SEED_PHASE_FRACTION` (the seed's `beta = 0.5`).** The stall is
  partly a seed placement - at 0.3 MPa the seed puts 4.2% of the system in
  phase II where the answer is 99.8% - and a different constant does repair
  that state. It also moves every state ADR-0024 and ADR-0025 converged, since
  they all start from the same expression. Rejected on bit-identity, and it
  would have been a fix for one geometry rather than for the iteration.
- **A trust region instead of a modified Hessian.** Equivalent in effect here
  and more machinery: the line search already present provides the
  globalization, and all that was missing was a direction with the right
  scaling in a non-convex region.
- **Repair `_rachford_rice` unconditionally.** See decision 2: measured to
  re-route 15 converged states.

## References

- ADR-0016 (the extended/negative-flash window), ADR-0019 (per-phase density
  roots), ADR-0022 (non-volatile components and the `1e-10` Wilson estimate),
  ADR-0024 (the log-space split stage), ADR-0025 (`ln W` summed in logs).
- Validation Case P-16 in `.agents/brain/validation-cases.md`, with amendments
  to Cases P-14 and P-15.
- Gill, Murray & Wright, *Practical Optimization*, Academic Press 1981,
  sec. 4.4.2 (modified Newton methods for indefinite Hessians).
- Leibovici & Neoschil, *Fluid Phase Equilibria* **74** (1992) 303 (the
  Rachford-Rice window); Whitson & Michelsen, *Fluid Phase Equilibria* **53**
  (1989) 51 (the negative flash).
- `examples/validation/22_stability_log_space.py --full` (routes 6 to 8),
  `examples/validation/21_pcsaft_polymer_vle.py --full`,
  `tests/test_flash_log_space_stage.py`, `tests/test_rachford_rice_extended.py`,
  `tests/test_pcsaft_polymer.py`.
