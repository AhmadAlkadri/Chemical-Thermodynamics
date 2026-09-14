# ADR-0030: The repeated equation-of-state solve is memoized in the call, not on the model

Status: accepted
Date: 2026-09-14

## Context
ADR-0023 built the measurement instrument and took three bit-identical
optimizations against it. It also recorded, in its own "Alternatives
considered", the one it would not take:

> **Cache the last root solve on the model instance.** It would catch more than
> the branch capability does - `phase_identity` and the post-split
> re-evaluation ask for the same state again - and it is bit-identical too.
> Rejected: both model classes are frozen dataclasses, a hidden mutable cache
> makes them effectively stateful and thread-unsafe, and a stale-key bug in it
> would be a *wrong thermodynamic answer* rather than a slow one.

Every clause of that rejection is still true of a cache *on the model*. None of
it is true of a memo that lives in the **call**, and that is what this ADR
adopts. The duplication it catches is real and was measured again at `d8814de`,
counting `PCSAFTEOS._density_roots` calls whose `(names, T, P, x)` are equal
bit for bit:

| case | solves | distinct states | repeat |
| --- | ---: | ---: | ---: |
| `pcsaft-vle-methane-hexane` | 137 | 102 | 25.5 % |
| `pcsaft-lle-water-hexane` | 275 | 238 | 13.5 % |
| `pcsaft-polymer-lle` | 601 | 512 | 14.8 % |

and, attributing each repeat to the pair of call sites that made it, it is four
pairs, three of which ADR-0023's roadmap named by name:

| pair, `pcsaft-lle-water-hexane` (37 repeats in all) | count |
| --- | ---: |
| `fugacity_coefficients` then `phase_identity` at the same state | 16 |
| `fugacity_coefficients` then `ln_fugacity_branches` at the same state | 14 |
| the same entry point twice (the post-split test re-evaluating a split composition) | 6 |
| `ln_fugacity_branches` then `phase_identity` at the same state | 1 |

**31 of the 37 are two different model methods reaching one root solve**, and
no single caller can see any of those, which is exactly why the ADR-0023
capability could not remove them: `ln_fugacity_branches` deduplicates *within*
one caller's need for two branches, and says nothing about a second caller
arriving at the same state a moment later. The remaining 6 are one entry point
called twice, and they are the only ones a memo at the caller seam would have
caught.

## Decision

### 1. The memo is a call scope, not model state
`chemthermo._eos_memo` (internal, ADR-0001) holds a bounded FIFO
`EosSolveMemo` and a `ContextVar` naming the one that is active.
`chemthermo.flash_tp` and `chemthermo.stability_tp` each wear the module's
`scoped` decorator, which installs a memo for the duration of the call and
resets the context variable in a `finally`. A `stability_tp` call made *inside*
a `flash_tp` call finds the flash's memo already active and keeps it, so one
flash is one memo rather than one per nested stability test.

Nothing is written to a model. `PengRobinsonEOS` and `PCSAFTEOS` stay the
frozen, stateless dataclasses they were; they *read* the active memo and take
the pre-ADR-0030 path when there is none, which is what happens for every
direct model call outside a flash or a stability test. A `ContextVar` is
per-thread and per-async-task, so two concurrent calls sharing one model
instance never see each other's memo. The memo cannot outlive its call, so
there is no lifetime in which a key could go stale relative to anything.

### 2. Two solve points, chosen because they are where the time is
- `PCSAFTEOS._density_roots`, keyed on
  `("pcsaft-density-roots", id(model), names, T, P, composition_key(x))`. Every
  `(T, P, x)` entry point of that class funnels through it -
  `fugacity_coefficients`, `log_fugacity_coefficients`,
  `ln_fugacity_branches`, `density_roots`, `molar_volume` and `phase_identity`
  - and it is 1.461 s of the 1.588 s a `pcsaft-lle-water-hexane` flash takes
  (92 %), essentially all of it inside `solve_density_roots`.
- `PengRobinsonEOS._compressibility_roots`, keyed on
  `("pr-compressibility-roots", A, B)`. The cubic's coefficients are functions
  of `A` and `B` alone, so the key needs neither the mixture nor the model
  instance; `fugacity_coefficients`, `ln_fugacity_branches`,
  `compressibility_factor` and `phase_identity` all reach it with the same
  `(A, B)` at one state.

`id(model)` appears in the PC-SAFT key because two instances may differ in
parameters or `k_ij` while agreeing on every argument the method receives. The
memo stores the instance alongside the value purely to hold a strong reference,
so the id cannot be recycled onto a different object while the entry lives.

### 3. A hit returns the identical object, which is what makes it bit-identical
Not an equal value - the same object. The first solve's `DensityRoots` (a
`NamedTuple` of floats) is handed back, so every double downstream of it is the
same double, and the bit-identity argument is a statement about object
identity rather than about floating-point arithmetic. The Peng-Robinson entry
stores a tuple and returns a fresh `list` of it, so no caller can mutate what a
later hit will read; the doubles are unchanged either way.

The key is the solver's own arguments compared for **exact equality**, with one
piece of care taken deliberately: `composition_key` gives `-0.0` and `0.0`
different keys. They compare equal and hash equal in Python while being
different bit patterns, and that is precisely the class of stale key ADR-0023
was worried about. Nothing else in a key can compare equal without being the
same double; `nan` never compares equal to itself, so a `nan` argument simply
misses and is recomputed.

### 4. The bound is a memory bound, never a numerical one
4096 entries, FIFO. Evicting an entry can only cause the solver to redo a solve
that produces the same doubles again, so no result anywhere depends on
`max_entries` or on the eviction order -
`tests/test_eos_memo.py::test_a_memo_bound_of_one_gives_the_same_answer_as_the_default_bound`
runs a whole associating flash through a **one-entry** memo and compares every
field with `==`. The largest memo any state in this repository fills is 512
entries (the polymer case).

### 5. The counters are not in `diagnostics`
The memo counts `hits` and `misses`, and they stay out of every result. A
`FlashResult.diagnostics` mapping is compared **whole** by
`tests/test_flash_refactor_bit_identity.py`, so a new key there is a changed
answer by this repository's own rule (that fixture's own history records "no
state gained or lost a diagnostics key" as the thing being proved). No
`FlashSettings` / `StabilitySettings` field was added either: those fields are
recorded verbatim in every benchmark record as the convergence criteria in
force, and a new one would change what a record says about a run whose numerics
did not move.

A caller that wants the counters opens the scope itself, which works because
`scoped` reuses an already-active memo rather than nesting a second one:

```python
with chemthermo._eos_memo.activated() as memo:
    chemthermo.flash_tp(...)
memo.hits, memo.misses
```

That is the context object, and it is how `tests/test_eos_memo.py` reads them.
Widening the scope that way also widens the memo - it then spans every call
inside the `with`, still exact because the key is exact, but no longer the
per-call lifetime the library itself uses. That is documented in the module and
is the only reason the wider scope is reachable at all.

### 6. The activity path is deliberately untouched
No activity-coefficient evaluation is memoized. The reason is not that it would
be wrong - it would be the same construction - but that `nrtl-lle-tessier-p1`,
`modified-raoult-vle` and `vlle-364k` are the benchmark's **drift control**
(ADR-0023 decision 3). They must run the same code before and after or there is
nothing to read the equation-of-state ratios against. Their measured solve
counts are 0 with the memo and 0 without it, which is the control working.

## What was measured

### Solves removed - machine-independent, and the honest headline
Counting executed `solve_density_roots` calls and executed cubic eigenproblems,
with the models' memo lookup returning `None` and with it active. These are
integer counts of the same fixed workload and do not depend on the machine:

| case | PC-SAFT density solves off -> on | PR cubic solves off -> on |
| --- | ---: | ---: |
| `pr-flash-ternary` | - | 230 -> 187 (**-18.7 %**) |
| `pr-stability-ternary` | - | 84 -> 73 (**-13.1 %**) |
| `pr-flash-grid-24` | - | 2203 -> 1679 (**-23.8 %**) |
| `nrtl-lle-tessier-p1` (control) | 0 -> 0 | 0 -> 0 |
| `modified-raoult-vle` (control) | 0 -> 0 | 0 -> 0 |
| `vlle-364k` (control) | 0 -> 0 | 0 -> 0 |
| `pcsaft-vle-methane-hexane` | 137 -> 102 (**-25.5 %**) | - |
| `pcsaft-lle-water-hexane` | 275 -> 238 (**-13.5 %**) | - |
| `pcsaft-polymer-lle` | 601 -> 512 (**-14.8 %**) | - |

### Wall times - and a caveat that has to be read first
The benchmark pair `benchmarks/baseline_d8814de.json` /
`benchmarks/after_07bf0b3.json` was taken back to back on one machine at nine
timed repeats, as ADR-0023 decision 3 requires, and **the machine was not
quiet**: an unrelated process held about eight of its twelve cores throughout
the session. All nine `result_hash` values are identical, which is the part of
the comparison a loaded machine cannot corrupt.

| case | before / s | after / s | ratio |
| --- | ---: | ---: | ---: |
| `pr-flash-ternary` | 0.016810 | 0.013450 | 1.25x |
| `pr-stability-ternary` | 0.004402 | 0.004119 | 1.07x |
| `pr-flash-grid-24` | 0.158998 | 0.120641 | 1.32x |
| `nrtl-lle-tessier-p1` (control) | 0.057287 | 0.059526 | 0.96x |
| `modified-raoult-vle` (control) | 0.024677 | 0.024191 | 1.02x |
| `vlle-364k` (control) | 0.137144 | 0.133123 | 1.03x |
| `pcsaft-vle-methane-hexane` | 0.216429 | 0.178299 | 1.21x |
| `pcsaft-lle-water-hexane` | 1.595234 | 1.324830 | 1.20x |
| `pcsaft-polymer-lle` | 0.994039 | 0.851818 | 1.17x |

The three controls at 0.96x / 1.02x / 1.03x are tighter than they were at five
repeats (0.77x / 1.25x / 0.92x on the first attempt), which is what the extra
repeats bought; they are still wider than ADR-0023's 1.00x / 1.04x / 1.02x on a
quiet machine.

A second measurement was therefore made, and it is the better estimate: the
same workload in **one process with the two arms interleaved repeat by
repeat**, the only difference between them being whether the models' memo
lookup returns the call's memo or `None`. Both arms take the same contention,
milliseconds apart; 11 repeats each:

| case | off / ms | on / ms | ratio | per-repeat spread |
| --- | ---: | ---: | ---: | --- |
| `pr-flash-ternary` | 13.77 | 13.41 | 1.027x | 1.01x .. 1.03x |
| `pr-stability-ternary` | 4.17 | 4.10 | 1.018x | 1.01x .. 1.02x |
| `pr-flash-grid-24` | 125.67 | 126.70 | 0.992x | 0.94x .. 1.14x |
| `nrtl-lle-tessier-p1` (control) | 52.64 | 54.03 | 0.974x | 0.82x .. 1.33x |
| `modified-raoult-vle` (control) | 22.88 | 24.95 | 0.917x | 0.74x .. 1.13x |
| `vlle-364k` (control) | 133.49 | 132.40 | 1.008x | 0.63x .. 1.12x |
| `pcsaft-vle-methane-hexane` | 190.01 | 146.00 | 1.301x | 1.23x .. 1.31x |
| `pcsaft-lle-water-hexane` | 1363.89 | 1233.82 | 1.105x | 1.01x .. 1.25x |
| `pcsaft-polymer-lle` | 924.09 | 801.81 | 1.153x | 1.05x .. 1.22x |

**What is claimed, and what is not.** The PC-SAFT cases agree between the two
measurements to within their spreads: **1.10x - 1.30x**, and that is the claim.
The Peng-Robinson half is close to break-even and the committed pair overstates
it: 1.25x / 1.07x / 1.32x there against the interleaved 1.03x / 1.02x / 0.99x,
and 1.10x / 1.02x / 1.04x on a second 31-repeat interleaved run. The arithmetic
says the interleaved figure is the right one - ADR-0023 measured the cubic
solve at 8.6 us of a 34.8 us `fugacity_coefficients` call, so removing a
quarter of the solves cannot be worth 30 % - so the Peng-Robinson claim is
**1.00x - 1.05x, with 13 % to 24 % fewer cubic solves**, and it is kept for the
work it removes rather than for a clock reading. A single uncontended pair is
the one piece of evidence this ADR is short of.

### The multiphase second-order stage
The stage this slice was also asked to look at builds its Hessian by perturbing
the whole non-reference mole-number vector and recomputing the *whole*
gradient, which evaluates every phase's fugacity terms for every column. Most
of those evaluations are at compositions bit-identical to the base point - only
the perturbed phase and the reference phase actually move - so the memo serves
them. Measured over the four builds the PC-SAFT water / ethanol / n-hexane
tie-triangle states reach:

| phases | Hessian size | model evaluations per build | with the memo | a structured Hessian would need |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 3 | 12 | 12 | 12 |
| 3 | 6 | 36 | 20 | 18 |
| 3 | 6 | 36 | 20 | 18 |
| 3 | 6 | 36 | 21 | 18 |

so on the three-phase builds the memo already removes **15 or 16 of the 18**
evaluations a structured Hessian would remove, and it removes them
bit-identically.

## The structured Hessian: measured, and **not shipped**
The assembly is exact as algebra. With `n_{p}` the mole numbers of non-reference
phase `p`, the reference phase holding `N_0 = z - sum_p n_p`, and
`a_{j,k} = ln x_{j,k} + f_j(x_j)_k`, the gradient is
`g_{p,k} = a_{p+1,k} - a_{0,k}` and phase `p+1` depends only on `n_p`, so

    H[(p,k),(q,l)] = delta_{pq} A^{(p+1)}[k,l] + A^{(0)}[k,l]

with `A^{(j)}` the Jacobian of phase `j`'s log-activity with respect to that
phase's own mole numbers. That is `N_phases x (NC x NC)` finite differences
instead of `((N_phases - 1) NC)^2`, which is the 36 -> 18 in the table above.

It is **not bit-identical**, and the reason is structural rather than
incidental. For a column `(q,l)` and a row `(p,k)` with `p != q`, the shipped
code forms `fl(fl(u - v_-) - fl(u - v_+))` where `u = a_{p+1,k}` is the same
double in both gradient evaluations and `v_±` are the reference phase's two
perturbed values; the structured form computes `fl(v_+ - v_-)`. Those differ
whenever `u` is large next to `v_+ - v_-`, which is the ordinary case here.
Two further differences compound it: on the diagonal blocks the two
contributions are divided by `2h` and added in a different order, and the
shipped code's reference-phase perturbation is `fl(z_l - fl(sum_q n_{q,l}))`
rather than an exact `N_{0,l} - h`.

Measured at **every** build the stage reaches, over two sets of states - the
20 Case P-9 PC-SAFT water / n-hexane states below `T3` and the six Case P-10
water / ethanol / n-hexane feeds:

| | 20 Case P-9 binary states | 6 Case P-10 ternary feeds |
| --- | --- | --- |
| builds compared | 8 (all two-phase) | 4 (one two-phase, three three-phase) |
| bit-identical | **0 of 8** | **0 of 4** |
| largest absolute difference in a Hessian entry | 3.20e-04 | 4.08e-04 |
| largest relative difference in a Hessian entry | 7.38e-05 | 2.54e-04 |
| largest `|H|` entry at those builds | 5.65e+04 | 7.71e+05 |

A relative difference of 1e-4 in the Hessian moves the Newton direction, and
this stage's line search accepts on an Armijo decrease *or* a residual
decrease, so a moved direction moves the accepted iterate and with it the
converged compositions in their last bits. Bit-identity is the gate, so the
structured Hessian is **not shipped**. It is left for a re-audit slice that is
willing to regenerate the bit-identity fixture and audit every moved state, in
the shape ADR-0017 and ADR-0021 used - and it is worth noting that after this
slice, what such a re-audit would buy on a three-phase build is 20 evaluations
down to 18, not 36 down to 18.

## Alternatives considered
- **A cache on the model instance** (ADR-0023's own rejected option). Rejected
  again, for ADR-0023's reasons, and this ADR is the way to get the same
  duplicate elimination without them.
- **A memo threaded through every internal signature** (evaluator, `_PhaseRoot`,
  `_verify`, `_multiphase`, `_detect`). Considered first and rejected on
  measurement rather than on taste: the duplication is *inside the model* -
  `phase_identity` and `fugacity_coefficients` at one state are two different
  model methods reaching one root solve - so a memo at the
  `chemthermo.flash._common` level catches only 6 of the 37 repeats in
  `pcsaft-lle-water-hexane` and none of the `phase_identity` ones, while
  touching a dozen signatures instead of two.
- **Serving `fugacity_coefficients(phase=P)` from a cached
  `ln_fugacity_branches` result** as `exp(branches[P])`. It would remove 14 more
  repeats. Rejected: it is a *different code path*, not a repeat of the same
  one - the ADR-0022 guard behaves differently when `exp` underflows - and the
  equality it relies on is ADR-0023's construction claim rather than object
  identity. The memo's whole safety argument is "the same object", and this
  would be the one entry that is not.
- **An `EquationOfState` method that returns a call-scoped copy of the model**
  (`with_solve_cache(memo)`, built with `copy.copy` plus
  `object.__setattr__`). Rejected: it adds a method to the public ABC (ADR-0001)
  and puts a scoped clone of the user's model into solver code, to buy exactly
  what the context variable buys with neither.
- **Memoizing `PCSAFTEOS._isotherm` as well.** Measured and dropped:
  `_density_roots` is 1.461 s of the 1.588 s `pcsaft-lle-water-hexane` flash and
  the isotherm construction is not visible next to it.
- **Memoizing the activity models.** Rejected on purpose; see decision 6.
- **Adding a three-phase equation-of-state case to the benchmark workload.**
  Wanted, and rejected on cost: one such state costs 13-40 s, so at the
  workload's five timed repeats plus a warm-up it would turn a half-minute
  instrument into a five-minute one. The three-phase evidence is the Hessian
  table above instead, measured directly and recorded in ledger Case B-2.
- **Putting `eos_evaluations` / `eos_memo_hits` in `diagnostics`.** Rejected;
  see decision 5.

## Consequences
- One new internal module, `src/chemthermo/_eos_memo.py`. Nothing is
  re-exported from `chemthermo`; the public API surface is unchanged, and so is
  the `EquationOfState` interface - this slice adds no method to it, unlike
  ADR-0022 and ADR-0023.
- `flash_tp` and `stability_tp` are wrapped by `functools.wraps`, so their
  names, docstrings and signatures are unchanged to `inspect` and to callers.
- A model called **directly**, outside `flash_tp` / `stability_tp`, sees no
  memo and runs exactly as before. Every test that calls a model directly is
  therefore untouched by construction.
- No fixture was regenerated and no pinned number moved:
  `refactor_bit_identity_v3.json` is unchanged, all nine benchmark
  `result_hash` values are identical, and the 2505-state robustness map re-run
  at `f852726` is identical to `robustness_64831bd.json` **field for field over
  all 2505 states and every aggregate field**, wall time excluded - 0
  differences, 2500 converged / 5 refused / 0 invariant violations either way.
  `benchmarks/robustness_64831bd.json` is pruned and its `.md` summary kept,
  per the `benchmarks/README.md` policy;
  `tests/test_robustness_map.py::COMMITTED_RECORD` now names the new file.
- `tests/test_eos_memo.py` adds 15 tests. Two of them are the Case F-4 and
  associating-state bit-identity sweeps and they are the two most expensive
  tests in the default suite, because each state is flashed **twice** - once
  with the memo and once without. Neither is marked `slow`: they are the
  slice's own evidence, the `slow` policy forbids marking the only test of a
  capability, and a bit-identity sweep run at half its states is a weaker
  proof rather than a cheaper one.
- **Suite time, and the part of the declaration this slice does not meet.**
  `pytest -q` reads **901 passed, 300.74 s (5:00) at `d8814de` -> 916 passed,
  286.03 s (4:46)**, the two runs back to back under the same contention. The
  same **901** tests after the change - the default run with
  `--ignore=tests/test_eos_memo.py` - read **235.36 s (3:55)**, reproduced at
  235.46 s: **1.28x, -21.7 %, and under the ~240 s ceiling** brain.md section
  10 tracks, on a machine that was not quiet.
  `pytest -q tests/test_robustness_map.py -m ""` - the slow full-sweep gate -
  reads **52 passed in 2028.71 s (33:48)**. The default run **as shipped is
  still over the ceiling at 286 s**, because this slice's own evidence costs
  ~50 s. So the library change carries the pre-existing suite under the
  ceiling and the slice as a whole does not carry the default run under it;
  that is stated rather than engineered away.
- `tests/test_examples.py` is unchanged, and measured: **101.64 s at `d8814de`
  -> 85.95 s** (1.18x), 43 tests either way, the two runs back to back on the
  one machine. Every heavy example already defaults to its cheap subset and
  takes `--full` for the rest (ADR-0017, ADR-0020), so there is no `--quick`
  to pass and nothing to weaken; the examples got faster because the library
  did, which is the only way this slice was allowed to move that number.

## What remains
- **One uncontended benchmark pair.** The committed pair was taken on a loaded
  machine; the hashes are conclusive and the wall times are not.
- **The 14 `fugacity_coefficients`-then-`ln_fugacity_branches` repeats** in the
  liquid-liquid case are still paid, because serving one from the other is a
  different code path (see Alternatives). Closing them means deciding whether
  ADR-0023's construction claim is strong enough to key a memo on, which is its
  own argument.
- **The structured multiphase Hessian**, with the numbers above, as a re-audit
  slice rather than a performance slice.
- **The memo is per call, so a sweep pays for every state once per call.** The
  robustness map and the benchmark grid case flash the same feed at neighbouring
  states and share nothing between them. That is the correct default - a
  cross-call memo is a cache with a lifetime, and a lifetime is where staleness
  lives - but a caller that knows its states repeat can open `activated()`
  itself, and nothing in this repository does.

## Supersedes (optional)
None. It **reverses one line** of ADR-0023's "Alternatives considered" - the
rejection of a last-solve cache - for a construction that is not the one
rejected: the memo is in the call, not on the model. ADR-0023's measurement
contract, its workload and its bit-identity gate are unchanged and are what
this slice was held to.
