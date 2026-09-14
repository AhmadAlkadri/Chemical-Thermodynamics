# ADR-0024: A phase split whose compositions leave machine range is solved in log mole numbers

Status: accepted
Date: 2026-09-13

## Context
ADR-0022 made a polymer a PC-SAFT component, and validation Case P-13 pinned
two states where the resulting flash **fails**: polyethylene (Mw = 16 400,
m = 431.32) in n-pentane at 453 K and 1 or 2 MPa, and the ternary with
n-hexane at 3 MPa. Those are not exotic corners - n-pentane is subcritical at
453 K, so *everything* below about 2.6 MPa is in that regime, and
polymer/gas systems live there entirely.

Three measurements say what the defect is, and none of them is about the
equation of state:

1. **The equilibrium is outside machine range.** The vapour over the melt is
   essentially pure solvent: at 1 MPa the melt holds `x_pentane = 0.97147` and
   the vapour holds `ln y_polymer = -450.53`, i.e. `y_polymer = 2.2e-196`. At
   Mw = 53 000 the same number is `exp(-1315)`, which is not a double at all.
   `_second_order_split` parametrizes the split by the linear mole numbers of
   one phase with a floor of `1e-300`, a box `0 < n_i < z_i` and a
   finite-difference step of `1e-7`: it can neither represent the answer nor
   step towards it, and it exhausted its 100 iterations on both pinned states.

2. **Successive substitution has nowhere to start.** `stability_tp` is right -
   the feed is unstable with `tpd_min = -412` at 1 MPa - and its minimizer is a
   genuine stationary point: an essentially pure polymer melt,
   `w = (1, 2.35e-178)`, converged to a stationarity residual of `1.4e-14`.
   The K-values that stationary point implies are `(2.5e-183, 0.0467)`. **Every
   one of them is below one**, so the Rachford-Rice function has no root in
   `[0, 1]` *and* the Leibovici-Neoschil window `1/(1 - K_max) < beta <
   1/(1 - K_min)` of ADR-0016 is empty. The documented Wilson fallback then
   produces one substitution step with an equal-fugacity residual of `3e+02`
   and hands a hopeless iterate to the linear stage.

3. **The ternary is a different failure wearing the same clothes.** At 3 MPa
   successive substitution *converges* - on the trivial solution, `x = y` to
   twelve figures with `beta = -6.3e+10`. Nothing was wrong with the iteration;
   the answer was refused downstream, by the "vapor fraction outside (0, 1)"
   check, and `flash_tp` raised without ever offering the state to the
   second-order stage.

The constraint on any fix is the one every slice here works under: the bit
identity fixture v3, the pinned Peng-Robinson and PC-SAFT numbers, the Case
P-13 liquid-liquid answers and the ADR-0023 benchmark result hashes must not
move.

## Decision

### 1. A second-order stage written in `u = ln n`
`chemthermo.flash._log_space.log_space_split` minimizes the same two-phase
reduced Gibbs energy `_second_order_split` minimizes, over the same variable,
under the substitution `n_i = exp(u_i)`:

- mole numbers are positive for every finite `u`, so the floor and the
  lower half of the box are gone and `u ~ -450` is an ordinary double;
- the second phase's mole fractions are formed as
  `ln y_i = u_i - logsumexp(u)`, so a mole fraction below the exponential's
  range still has a finite, accurate logarithm, and the ADR-0022 log-space
  `ln phi` route means no model evaluation needs `exp` either;
- the gradient is unchanged. The residual `r_k = ln(y_k phi_k^II) -
  ln(x_k phi_k^I)` is `dg/dn_k` (ADR-0009 equation 3), so `dg/du_k = n_k r_k`
  by the chain rule and the two parametrizations have the same stationary
  points.

The **Newton system is written on `r(u) = 0`, not on `dg/du = 0`**. With
`J = dr/dn` the true Hessian in `u` is `diag(n) J diag(n) + diag(n r)`, whose
condition number is of order `n_max / n_min` - `1e+190` on this state, i.e.
numerically singular - while `dr/du = J diag(n)` has entries of order one
there. The two give the same direction in exact arithmetic; one of them
survives in doubles. `dr/du` is built by central differences in `u` (a
*multiplicative* step on `n`, which is what makes one step size usable across
190 decades). Descent is still measured on the energy, through
`(n r) . du`; a direction that is not a descent direction is replaced by `-r`,
which is successive substitution in log space. The line search and the Armijo
constant are `_second_order_split`'s.

The stage runs **only** where the existing path has already failed or cannot
start (decision 2), so no state that converged before this ADR reaches it and
every pinned number is unchanged by construction, not by measurement.

### 2. Three triggers, all of them on states that used to raise
On the phi-phi tangent-plane path:

- **The seed trigger.** When the stability minimizer has a component at or
  below `1e-30` (`TRACE_MOLE_FRACTION`) **and** its K-values bracket no
  Rachford-Rice root, the K-loop is skipped entirely and the split is solved in
  log space from a seed built out of the same stationary point. Both conditions
  are required: a trace component with a workable bracket is ordinary, and a
  collapsed bracket at an ordinary stationary point is what the Wilson fallback
  is for. `diagnostics["k_seed"] = "stability-log"`.

  The seed is `n_i = beta K_i z_i / ((1 - beta) + beta K_i)` at a **fixed**
  `beta = 0.5`, formed in logarithms from `ln K_i = ln W_i - ln z_i` (or its
  reciprocal). Two properties make it usable where Rachford-Rice is not: the
  denominator is a convex combination of `1` and `K_i` and so is strictly
  positive - no pole to step over, no window to bracket - and `n_i / z_i` lies
  strictly in `(0, 1)` for every positive `K`, so the seed satisfies the
  two-phase box componentwise whatever the K-values are. It is only a starting
  point: the stage is a descent method, and the same answer comes back from
  `n = 0.5 z`, `0.9 z` and `0.99 z` to twelve figures.

- **The hand-over trigger.** When the linear second-order stage finishes above
  `settings.tol`, the same minimization is retried in log space from the best
  iterate it produced (`log_space_seed = "linear-iterate"`). Reached only on a
  state that would otherwise raise.

- **The unrepresentable-K trigger.** When successive substitution raises
  `ModelError` because the two phases' `ln phi` differ by more than the
  exponential's range, the log-space stage is started from the stationary
  point instead. This is the Mw = 53 000 chain below the solvent's saturation
  pressure.

And one change that is not about log space at all: **a first stage that
converged on the trivial solution is offered to the second-order stage**
instead of being refused outright. That is the ternary of measurement 3; the
"vapor fraction outside (0, 1)" refusal still stands after the stage, so the
only states affected are ones that raised before.

### 3. An exactly zero mole fraction is a legitimate phase composition
When a converged phase's mole fraction is outside the exponential's range,
`0.0` is the nearest double there is and the result carries it.
`Composition` already accepts zeros (`validate_fractions` requires
non-negative, not positive), so nothing in the public contract changes.

The number is not lost. `diagnostics["log_space_ln_x_min"]` is the smallest log
mole fraction of the split's second phase and
`["log_space_ln_x_min_component"]` names the component;
`["log_space_zero_fractions"]` counts how many underflowed. Two consequences
are recorded rather than hidden: the **material balance becomes exact** (the
other phase then holds every mole of that component the feed had), and
`fugacity_residual` - taken over components present in both phases - cannot see
that component, so the stage's own log-space residual is reported as
`["log_space_residual"]`.

### 4. Everything downstream is unchanged
The log-space stage returns the same `_SplitSolution` the K-loop returns, so
verification, compressibility-based naming (ADR-0017/0019), post-split
stability, phase addition and assembly are the code that was already there.
The polymer VLE states come back named `"liquid"` / `"vapor"` from a measured
`kappa` (1.16 and 0.0024 at 1 MPa), post-split stable, with
`delta_g_split_rt < 0`.

## Consequences

**What now works.** `flash_tp` returns a verified vapour-liquid split for
polyethylene(16 400) / n-pentane at 5 wt%, 453 K and every pressure from 0.3 to
2.5 MPa; the melt's solvent content agrees with an independently written
one-dimensional equal-fugacity solve to 1.2e-15 - 2.5e-14 absolute, and FeOs's
chemical potentials at those phases agree to 2.6e-12 (matched universal
constants) / 2.5e-07 (as shipped). The 453 K pressure scan from 0.3 to 12 MPa
is continuous in verdict - VLE to 2.25 MPa, LLE from 2.74 to 9.56 MPa, one
liquid above the Case P-13 cloud point of 9.7518 MPa - with no
`ConvergenceError` anywhere. The ternary at 3 MPa converges. Both Case P-13 (vi)
defects are retired and pinned the other way round as Case P-14.

**What did not move.** 807 tests pass unchanged, the bit-identity fixture v3 is
untouched, and all nine ADR-0023 benchmark cases report identical result
hashes at 0.96x-1.01x - including `pcsaft-polymer-lle`, whose whole
liquid-liquid regime goes through the unchanged path.

**The cost.** One more module (about 300 lines) and one more branch in
`_flash_tp_tangent_plane`. The stage is *more* expensive per iteration than the
linear one for the same problem size - the finite-difference Jacobian is the
same `n` gradient evaluations, but the line search may take a very long first
step - which is the second reason it is gated rather than made the default; the
first is bit identity. On the states that use it a whole flash is 0.3 - 1.5 s,
against the `ConvergenceError` it used to be.

**Why not always log space.** Two reasons, in this order. (i) Bit identity:
`ln n` and `n` do not produce the same doubles, so switching the default would
move every pinned two-phase number in the repository, and this slice would have
had to re-derive rather than preserve them. (ii) Conditioning in the other
direction: `dr/du` is well conditioned when the mole numbers span many decades
and *worse* conditioned than the linear system when they do not, because the
substitution multiplies each column by `n_j`. The linear stage remains the
right tool for an ordinary split.

**Alternatives considered.**

- *A better Rachford-Rice bracket.* There is no bracket to find: the window is
  provably empty when every `K_i` is on one side of one, and that is a fact
  about the stationary point, not about the solver.
- *Seeding from Michelsen's `W` directly, unscaled.* `sum W = 9.1e+178` at this
  state, so `W` is not a phase of one mole of feed and clipping it to `z`
  gives the melt the whole feed, i.e. the trivial split. The fixed-`beta`
  Rachford-Rice denominator is what turns `ln K` into an admissible split.
- *Seeding from the shallow vapour-side stationary point the `wilson-vapor`
  trial also finds* (`tpd = -1.8e-04` at 1 MPa). It does have a Rachford-Rice
  bracket, but it is a near-trivial stationary point and the substitution from
  it is the one that was already measured to stall.
- *A `sqrt(n)` parametrization*, which is symmetric, well conditioned and keeps
  a true Hessian. It represents `1e-196` but not `exp(-1315)`, and the 53 000
  chain is the state that decided it.

**What remains.**

1. **Multiphase.** `_multiphase.py` is untouched: a three-phase state whose
   compositions leave machine range would fail the same way the two-phase one
   did. No state in this repository needs it yet, so it is not written
   (ADR-0002).
2. **The longest chain below 1 MPa.** Mw = 53 000 at 0.5 and 1 MPa still does
   not converge, and the reason is upstream of this ADR: `stability_tp`'s
   deepest stationary point there is a shallow vapour-side one
   (`tpd ~ -1e-04`), so the melt is never found and the stage is seeded 1300
   orders of magnitude away from the answer. That is a **stability trial set**
   problem for `m = 1393.9`, pinned in
   `tests/test_flash_log_space_stage.py`.
3. **Performance.** The stage has never been profiled and is not in the
   ADR-0023 workload. Adding a polymer-VLE case to the benchmark would change
   the committed records, so it is left to whoever next touches them.
