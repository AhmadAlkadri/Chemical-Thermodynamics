# ADR-0009: Liquid-liquid `flash_tp` from an activity model, and post-split stability

Status: accepted
Date: 2026-09-13

## Context
After ADR-0008, phi-phi `flash_tp` decides one phase versus two from Michelsen's
tangent-plane test and seeds the split from the minimizer. Two gaps remained,
both recorded there as "next slices":

1. **No liquid-liquid flash.** `stability_tp(..., activity_model=...)` could
   already *detect* a liquid-liquid split (ADR-0007), and the Tessier et al.
   (2000) global minima were reproduced, but there was no way to obtain the
   resulting tie-line. `flash_tp` required an `eos`.
2. **No post-split verification.** ADR-0008 decision 4 verified a converged
   two-phase result for material balance, equal fugacities and a negative Gibbs
   energy change, but explicitly did *not* re-test the converged phases for
   stability: "a three-phase state will still come back as two phases".

## Decision

### 1. `eos` becomes optional and a third mode exists
`flash_tp(mixture, *, temperature_K, pressure_Pa, eos=None, activity_model=None,
flash_mode=None, settings=None)`.

- `flash_mode` defaults to `None`, meaning *infer from the models supplied*:
  an activity model with no EOS is `"gamma-gamma"`, anything else is
  `"phi-phi"` (the historical default). Naming a mode explicitly always wins,
  so `flash_mode="phi-phi"` without an `eos` is still a `ModelError` rather
  than a silent reinterpretation, and `flash_mode="gamma-gamma"` with an `eos`
  is a `ModelError` for the same reason `stability_tp` refuses both models
  (ADR-0007): combined gamma-phi equilibrium is `flash_mode="gamma-phi"`, which
  is unchanged.
- `eos` gains a default of `None` so the declared call
  `flash_tp(mixture, temperature_K=..., pressure_Pa=..., activity_model=...)`
  is legal. Passing no model at all is still a `ModelError`.

### 2. The liquid-liquid path is the tangent-plane path, with `gamma` for `phi`
Both phases are liquids at the same pure-liquid reference state, so `mu_i^0`
cancels and the equilibrium condition is equality of activities,
`x_i^I gamma_i^I = x_i^II gamma_i^II`. The flow is identical to phi-phi:

    flash_tp -> stability_tp(feed, activity_model) -> single phase | seeded split

- stable -> one phase named `"liquid"`, `phase_fractions = {"liquid": 1.0}`,
  `vapor_fraction = None`, `termination_reason = "feed_stable_tangent_plane"`;
- unstable -> the split is seeded exactly as phi-phi seeds it, from Michelsen's
  unnormalized mole numbers `W = w exp(-tpd)` with `K_i = W_i / z_i` (ADR-0008
  decision 2), and the same Rachford-Rice / successive-substitution loop runs
  with the K update `K_i = gamma_i^I / gamma_i^II`;
- inconclusive -> `ConvergenceError`.

`vapor_fraction` is `None` for every gamma-gamma result: neither phase is a
vapor, and reporting a number there would be fiction.

**Phase names are roles, not identities.** `"liquid1"` is the phase the split
started from as feed-like, `"liquid2"` the one started from the tangent-plane
minimizer. Nothing distinguishes two liquids the way volatility distinguishes a
vapor from a liquid (ADR-0008 decision 3 had at least a Wilson ranking to lean
on; there is no analogue here), so no attempt is made to name them by
composition. Two feeds on the same tie-line can return the same pair of
compositions with the labels swapped - measured, for n-butanol/water at
z1 = 0.10 versus z1 = 0.20. Callers must compare the phase *set*.

Diagnostics mirror the phi-phi path: `phase_detection`, `stability_status`,
`tpd_min`, `stability_trials`, `k_seed`, `mass_balance_residual`,
`delta_g_split_rt`, plus `equilibrium_residual`
(`max_i |ln(x_i^I gamma_i^I) - ln(x_i^II gamma_i^II)|`, the gamma-gamma name for
what phi-phi calls `fugacity_residual`), `ssi_iterations`,
`second_order_iterations` and `converged_stage`.

### 3. The second-order stage minimizes the two-phase Gibbs energy
Successive substitution on the equal-activity condition converges linearly with
a ratio close to one near a plait point: on the four Tessier Problem 1 feeds it
needs 536, 794, 1315 and 3922 iterations from the stability seed, against a
default budget of 100. A second-order stage is therefore not optional.

The stage is a damped Newton **minimization** of the reduced two-phase Gibbs
energy in the phase-II mole numbers `n_i` (`0 < n_i < z_i`):

    g(n) = sum_i (z_i - n_i) ln(x_i^I gamma_i^I) + sum_i n_i ln(x_i^II gamma_i^II)

whose gradient collapses, by Gibbs-Duhem plus `sum_i N_i d ln x_i = 0`, to

    dg / dn_k = ln(x_k^II gamma_k^II) - ln(x_k^I gamma_k^I)

i.e. **the gradient is exactly the equal-activity residual the result reports**.
The Hessian is central-differenced from that gradient (no model derivatives, no
`scipy`), symmetrized, and shifted by a multiple of the identity when it is not
positive definite; the line search backtracks inside the box on an Armijo
decrease of `g` or a decrease of the residual.

**The alternative was implemented, measured and rejected.** Newton on the
residual system in `(ln K, beta)` - the first option the design named - solves
three of the four Problem 1 feeds but walks into the trivial branch on
`z = (0.12, 0.05, 0.83)` (`beta -> -8`, the two phases merging) and stalls at a
residual of 1.4e-07; it only converges there if given 200 or more substitutions
first. A residual system cannot distinguish the equilibrium from the trivial
solution - both are roots. Minimizing `g` can: the trivial solution
(`n_i = beta z_i`) is a stationary ridge at `g = g(feed)`, and an unstable feed
has `g < g(feed)` at the true split, so monotone descent cannot reach it. The
descent method converges that feed in 7 iterations to a residual of 4.4e-16.

The stage targets `FlashSettings.second_order_tol` (default 1e-12), tighter than
`tol`, because it converges quadratically (the extra step is cheap) and because
the equal-activity residual is what a caller verifies. A split is *accepted* as
soon as it meets `tol`.

**It applies to the liquid-liquid split only.** The phi-phi and gamma-phi splits
are untouched by this release, which is what keeps every validated phi-phi
number identical: the worst phi-phi split on the in-repo grid uses 32
substitutions, so even a shared 50-iteration budget would have left the stage
dormant, but not running it at all makes the claim exact rather than measured.
Accelerating the phi-phi split (one state in the 1144-state scan of ADR-0008
still exhausts its budget) stays open.

### 4. Every two-phase result is re-tested for stability
Each converged phase is fed back into `stability_tp` with the same model, and
the outcome is always recorded: `post_split_checked`, `post_split_stable`,
`post_split_status`, `post_split_tpd_min`, `phase_stability_<name>` and
`phase_stability_tpd_min_<name>`.

- **Converging onto the partner is not an instability.** Two coexisting phases
  share one tangent plane, so a stability test on either finds the other with
  `tpd = 0` (validation Case S-3). Up to the split's own tolerance that zero is
  a small negative number - worst measured -7.0e-09 over the in-repo
  Peng-Robinson grid, inside the default `tpd_tol = 1e-8` by a factor of only
  1.4. A minimizer that *is* the partner phase is therefore classified
  `"marginal"`. "Is the partner" uses the stability module's own
  trivial-solution measure, `sum_i ln(w_i / x_i^partner)^2 < trivial_tol`,
  rather than a new tolerance.
- **Anything else negative raises.** `ConvergenceError` states that the
  two-phase solution is not a stable phase set and that a third phase is
  required. It raises **by default** because returning a two-phase answer that
  the package itself has proved wrong is the failure mode this whole campaign
  exists to remove: ADR-0008 shipped "the converged phases are not re-tested",
  and a silent two-phase answer for a three-phase state is exactly the kind of
  plausible-looking wrong number ADR-0007 refused to produce for gamma-phi
  stability. `FlashSettings(post_split_stability=False)` returns the result
  anyway with the failure in `diagnostics`; the flag gates the *raise*, not the
  computation, so the diagnostics are there either way.
- **Two paths cannot run it and say so.** `gamma-phi` has no stability test at
  all (ADR-0007) and the legacy `phase_detection="wilson-heuristic"` path exists
  to reproduce pre-ADR-0008 behavior unchanged. Both report
  `post_split_checked = False` with a `post_split_skipped_reason`.

### 5. `cli_schema_version` stays 1
Same reasoning as ADR-0008: the new keys live inside the free-form `diagnostics`
mapping, which removes nothing and changes no type. Asserted in
`tests/test_cli_tp_flash.py::test_cli_tp_flash_json_diagnostics_carry_the_post_split_keys`.

## Alternatives considered
- **A separate `flash_lle(...)` function.** Rejected: the phase *count* is an
  output of this package, not an input. A user who calls `flash_lle` has already
  decided there are two liquids; `flash_tp` with an activity model returns one
  phase when the feed is stable, which is the answer a separate entry point
  would make awkward to express. It would also duplicate the stability call, the
  split loop, the verification and the post-split check.
- **An "assume two liquids" mode** that skips the stability test and solves the
  split directly. Rejected for the same reason, and because it contradicts the
  campaign thesis that a split must be *discovered and verified*, never assumed:
  it would happily return a converged fixed point for a single-phase feed.
- **Naming the two liquids by composition** (for example "aqueous" / "organic",
  or by which is richer in component 1). Rejected: any such rule is a heuristic
  about the chemistry, not a thermodynamic statement, and it would break under
  component reordering or a different system. Roles plus an explicit warning are
  honest; a name that looks physical and is not would be worse.
- **Post-split check off by default.** Rejected: see decision 4.
- **A composition-distance tolerance for the partner test.** Rejected in favour
  of reusing `trivial_tol`, so the flash does not invent a second notion of "the
  same phase" alongside the stability module's.
- **GDEM (dominant-eigenvalue) acceleration** instead of a second-order stage.
  Not chosen: it accelerates the same fixed point and inherits the same
  trivial-solution attraction, and it would not have produced the
  round-off-level `equilibrium_residual` that makes the verification meaningful.

## Consequences
- Positive: a verified liquid-liquid tie-line is available from a single call
  with no EOS, and the split is *discovered* by the stability test rather than
  assumed. Validation Cases L-1, L-2, L-3.
- Positive: every two-phase answer on the reference paths now carries a
  statement about its own phase set, and a state needing a third phase is
  reported instead of silently returned. Case L-4.
- Positive: the equal-activity residual of a liquid-liquid split is at round-off
  (<= 1.8e-14 over all nine validated feeds), so `delta_g_split_rt` and the
  material balance are meaningful to the digits printed.
- Tradeoff: every two-phase flash on the tangent-plane paths now pays for two
  extra stability analyses. Measured: 0.9 s for the 47 two-phase states of the
  in-repo grid, roughly doubling their cost.
- Tradeoff: five more `FlashSettings` fields, and a `flash_mode` default that is
  now `None` rather than `"phi-phi"` (observationally identical for every
  existing call).
- Tradeoff: the post-split check is only as sharp as the split it checks. With a
  deliberately loosened `tol` (1e-3) a genuinely two-phase state can be reported
  as needing a third phase, because the minimizer is then too far from the
  partner to be recognised as it. At default tolerances the margin is eight
  orders of magnitude.
- Known limitation: `flash_tp` still returns at most two phases. A three-phase
  state is now *detected* and refused, not solved.
- Known limitation: `"stable"` remains bounded by `stability_tp`'s deterministic
  trial set, and that bound now applies to the post-split verdict too - a third
  phase that no trial reaches would not be found.

## Next slices
- `flash-phase-addition`: multiphase Rachford-Rice / Michelsen multiphase flash,
  so a `FlashResult` can carry more than two phases and the post-split failure
  becomes an answer instead of an error.
- `vlle`: once more than two phases exist, wire the multiphase split through the
  `chemthermo.vlle` boundary. This also needs the pure-liquid reference fugacity
  that ADR-0007 records as missing, since a VLLE state mixes an activity liquid
  with a vapor.

## Supersedes (optional)
None. Discharges ADR-0008's "test both converged phases for stability now" and
"phase addition/removal and LLE are the next slice" notes, except for the
multiphase result contract, which moves to `flash-phase-addition`.

## Superseded by (optional)
None.
