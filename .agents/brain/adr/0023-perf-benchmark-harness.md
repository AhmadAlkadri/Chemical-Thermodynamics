# ADR-0023: A benchmark record is the evidence a performance change is accepted on

Status: accepted
Date: 2026-09-13

## Context
Every slice so far has been measured against *answers*: a pinned fixture, a
published table, an independent implementation. Runtime has been measured too,
but informally - a suite time in a commit message, a "76.6 s -> 44.3 s" in
ADR-0021 - with no artefact anyone can re-run, no record of the machine it was
taken on, and no way to tell a real gain from an afternoon's thermal drift.
This slice needed exactly that, because two of the three things it changes are
worth between 1 % and 3 % each and the fourth is worth 25 %, and there is no
way to know which is which by reading the diff.

The other half of the problem is the one this repository already takes
seriously everywhere else: **a faster run that moved a number is not an
optimization.** A performance slice has to prove the answer held, and "the
tests still pass" is not that proof when the tests carry tolerances.

The profile that motivated the work (Apple M2 Max, Python 3.11.6, numpy 2.4.2):

- the Peng-Robinson ternary flash at 240 K and 3 MPa spends about 55 % of its
  time in `PengRobinsonEOS.fugacity_coefficients`, and half of *that* is
  `numpy.roots` - 16.9 us of a 34.8 us call;
- the PC-SAFT water / n-hexane liquid-liquid flash at 298.15 K and 1 atm spends
  95 % of its time in `solve_density_roots`, and 68 % of a single 5.3 ms solve
  is the ~19 single-point `pressure_and_slope` evaluations the safeguarded
  Newton makes;
- across that whole flash, 291 density-root solves were made at 238 distinct
  states. The double-solving the minimum-Gibbs rule causes is real but it is
  **53 solves of 291**, not half of them: ADR-0021 already removed the bulk of
  it by pinning each trial iteration to one root.

That last number is why this ADR records a measurement rule rather than a
prediction. The slice was scoped expecting root reuse to be worth roughly 2x on
PC-SAFT; it is worth about 5 %, and the harness is what said so before the code
was written.

## Decision

### 1. `chemthermo.bench` produces the record
An internal package (ADR-0001: nothing re-exported from `chemthermo`), run as
`python -m chemthermo.bench --out record.json` or `python tools/bench.py`. It
runs a **fixed** nine-case workload - Peng-Robinson flash, stability and a
24-state grid; NRTL gamma-gamma; modified-Raoult VLE and the 364 K VLLE feed;
PC-SAFT vapour-liquid, associating liquid-liquid and the ADR-0022 polymer - and
records per case:

| field | why it is in the record |
| --- | --- |
| model, components, `(T, P)`, overall composition | what was run |
| phase count, phase names, compositions, fractions | what came out |
| convergence criteria (every `FlashSettings` / `StabilitySettings` field that can move an answer) | a wall time at a different tolerance is a different measurement |
| initialization (`k_seed`, `phase_detection`, `feed_branch`, `incipient_phase`) | the same state reached two ways costs two different amounts |
| iteration counts, from the solver's own diagnostics | distinguishes "each step got cheaper" from "there were fewer steps" |
| derivative mode (`"analytic"`) | constant here, and worth saying: no model in this library differentiates by finite difference |
| median / min / max wall time over `--repeats` timed runs after one untimed warm-up | see decision 4 |
| peak `tracemalloc` allocation of one separate instrumented run | the tracer roughly triples the time it observes, so it never runs inside a timed repeat |
| refusal status and message for any state the solver declined | a "speedup" from a state that stopped converging has to be visible |
| platform, CPU string, core count, Python and numpy versions, git commit and dirty flag | a wall time means nothing without them |

Workload *preparation* - databank reads, parameter resolution, mixture
construction - happens outside the timed region. What is timed is `flash_tp`
and `stability_tp`.

### 2. Bit-identity is the gate, and the result hash is the check
An optimization in this repository may not move a thermodynamic result. Not
"within tolerance": not at all. The three changes this slice ships are
bit-identical **by construction** - each is a rearrangement of *when* an
expression is evaluated, never of the expression - and the evidence is:

- the 155-state pinned fixture `refactor_bit_identity_v3.json`, floats compared
  with `==`, unchanged and not regenerated;
- `tests/test_eos_branch_reuse.py`, which compares the one-solve branch route
  against the per-branch route with `==` over the Case F-4 subset, the
  water / n-hexane liquid-liquid states, the Peng-Robinson grid, the
  single-root case and the polymer state where `exp(ln phi)` underflows;
- every other pinned test in the suite.

The record's `result_hash` - phase names, compositions and fractions at twelve
significant digits, plus the state and the refusal status - is a *check* on top
of that, not a substitute for it: it is what
`python -m chemthermo.bench --compare before.json after.json` asserts, and that
command **exits 1** when any case's hash differs.

### 3. A performance change is accepted on a baseline, an after, and a ratio
The rule, in full:

> A committed baseline record, a committed after record, identical
> `result_hash` on every shared case, and a measured improvement on the cases
> the change was about - with no case materially worse.

Both records live in `benchmarks/`, named for the commit they were measured at,
and both must be measured **back to back on one machine**. That is not
pedantry: the first baseline taken in this slice and a run twenty minutes later
differed by 5-10 % on cases whose code had not changed, which is larger than
two of the three optimizations here. The workload's three activity-model cases
touch none of the code this slice changes and are the drift control: they
measured 1.00x, 1.04x and 1.02x in the accepted pair.

### 4. Median of timed repeats, after a warm-up
The first call through any of these paths pays costs an optimization does not
move. The median rather than the minimum because this machine produces
occasional slow repeats and no fast ones, so the minimum flatters; min and max
are recorded so the spread stays visible.

### 5. The workload is fixed
Adding a case is a normal change - older records simply lack it and `--compare`
says "after only". Editing or removing one invalidates every committed baseline
for that case and has to be said in `benchmarks/README.md`.

## The three optimizations this ADR is adopted alongside

**A. One root solve serves every branch.**
`EquationOfState.ln_fugacity_branches` (optional, `None` by default) returns
`ln phi` on every admissible branch from one root solve. The minimum-Gibbs rule
(ADR-0005), the root-count measurement (ADR-0021 decision 4) and the split's
lowest-Gibbs fallback (ADR-0019) each have to compare the branches and so asked
the model twice at one `(T, P, x)`. Bit-identity is by construction because the
capability returns `ln phi` *before* the exponential: `exp` of it is the same
double `fugacity_coefficients` returns, so `flash._common.eos_branch_terms_all`
reconstructs exactly what the per-branch route produced, ADR-0022 guard
included. A pinned trial's iteration still asks for one branch (ADR-0021
decision 3) and is untouched - asking for both there would undo that saving.

**B. Peng-Robinson reaches LAPACK without `numpy.roots`' bookkeeping.**
For a monic polynomial with a non-zero constant term, `numpy.roots` *is*
"build the companion matrix, take its eigenvalues"; the companion matrix is now
built directly and `numpy.linalg.eigvals` called on it, skipping the
`atleast_1d`, the non-zero scan, the trimming, the dtype check, the division by
a leading coefficient that is exactly `1.0` and the `hstack`. 14.7 us -> 8.6 us,
and the eigenvalues are `==` over 200 000 random `(A, B)` pairs. A constant term
of exactly zero, the one case `numpy.roots` deflates instead, is handed back to
it.

**C. The PC-SAFT density solver stops recomputing what it already knows.**
`A^res/RT` itself is never read by the root solver - `pressure` needs `a'`,
`pressure_and_slope` needs `a'` and `a''` - and it is the only part of
`PCSAFTIsotherm.derivatives` that takes a logarithm; `value=False` skips it. The
refinement's bracket already carries the residual at its left edge, which the
refinement used to recompute, and the refinement already has `(P, dP/drho)` at
the root it returns, which the mechanical-stability filter used to recompute.
The scan grid, the bracketing and the refinement arithmetic are untouched.

Measured, back to back at 5 repeats, Apple M2 Max / Python 3.11.6 / numpy
2.4.2, `baseline_3ce68df.json` against `after_2ca41bf.json`:

| case | before | after | ratio |
| --- | ---: | ---: | ---: |
| `pr-flash-ternary` | 15.26 ms | 13.04 ms | 1.17x |
| `pr-stability-ternary` | 4.65 ms | 3.88 ms | 1.20x |
| `pr-flash-grid-24` | 140.9 ms | 117.6 ms | 1.20x |
| `pcsaft-vle-methane-hexane` | 240.8 ms | 188.5 ms | 1.28x |
| `pcsaft-lle-water-hexane` | 1627.0 ms | 1236.5 ms | 1.32x |
| `pcsaft-polymer-lle` | 962.1 ms | 814.4 ms | 1.18x |
| `nrtl-lle-tessier-p1` (control) | 43.88 ms | 43.71 ms | 1.00x |
| `modified-raoult-vle` (control) | 20.26 ms | 19.57 ms | 1.04x |
| `vlle-364k` (control) | 99.96 ms | 97.97 ms | 1.02x |

All nine result hashes identical.

## Alternatives considered
- **Cache the last root solve on the model instance.** It would catch more than
  the branch capability does - `phase_identity` and the post-split re-evaluation
  ask for the same state again - and it is bit-identical too. Rejected: both
  model classes are frozen dataclasses, a hidden mutable cache makes them
  effectively stateful and thread-unsafe, and a stale-key bug in it would be a
  *wrong thermodynamic answer* rather than a slow one. The capability is
  stateless and cannot be stale.
- **A closed-form (Cardano) cubic solver.** Far faster than `eigvals` and it
  does not reproduce these doubles. Bit-identity is the gate; rejected, and the
  code says so where someone will next be tempted.
- **Vectorize the Peng-Robinson per-component `ln phi` loop.** Measured at
  3.0 us against the loop's 1.4 us at the two- and three-component sizes these
  reference paths use: three numpy operations on a length-3 array cost more than
  the loop they replace. Rejected on measurement, and the loop carries the
  number.
- **Replace the association module's `np.einsum` calls with explicit
  multiply-and-sum.** *Not* bit-identical for three or more sites (checked at 1,
  2, 3, 4, 6 and 8 sites against both a `.sum(-1)` and a `matmul` form), and
  slower at the sizes that occur (1.34 us against 1.00 us). Rejected twice over.
- **Refine every bracket of a density solve in one vectorized pass.** Rejected:
  `solve_site_fractions` tests convergence on `max |residual|` over the whole
  batch, so batching two brackets together would give one of them a different
  number of Newton steps and move its last bits.
- **Cache the identity matrix the site-fraction Newton step allocates.**
  Implemented, measured, reverted: no gain outside the noise, and it puts a
  shared array where a later edit could write to it.
- **Make the harness a public API.** Rejected: it is a maintainer's instrument.
  It is also the only place in the package that shells out to `git` and reads
  `sysctl`, neither of which belongs behind a public import.

## Consequences
- `benchmarks/baseline_3ce68df.json` and `benchmarks/after_2ca41bf.json` are
  committed, with `benchmarks/README.md` explaining how to compare them.
- `EquationOfState` gains one optional method. Implementations that do not
  provide it are unaffected: the default returns `None` and every caller keeps
  its per-branch behaviour, which `tests/test_eos_branch_reuse.py` pins with a
  model that implements only the pre-ADR-0023 interface.
- `PCSAFTIsotherm.derivatives` and `AssociationIsotherm.derivatives` gain a
  `value` keyword. Both are internal (ADR-0001); with `value=False` the `a` slot
  is `nan` rather than a plausible-looking zero.
- `_refine` returns a `_RefinedRoot(eta, pressure, slope)` instead of a float,
  and takes the bracket's left residual as an argument.
- No fixture was regenerated and no pinned number moved.
- The default suite got **faster while gaining six tests**: `pytest -q`
  **263.80 s for 803 tests at `3ce68df` -> 224.17 s for 809 at `2ca41bf`**
  (-15 %), both runs uncontended on the one machine. `pytest -q -m slow` is
  **691.22 s (11:31) for 28 tests**, against the 856.2 s (14:16) recorded for
  the same 28 at the previous slice - a different session, so read that pair as
  indicative rather than back to back. **No test was marked `slow` in this
  slice** and none was deleted: the 803 -> 809 is the six of
  `tests/test_eos_branch_reuse.py`, and the 796 -> 803 before it was the seven
  of `tests/test_bench_harness.py`.

## Supersedes (optional)
None. Adds a decision rule that did not exist. Extends the `EquationOfState`
interface of ADR-0001 with one optional capability, in the same shape ADR-0022
added `log_fugacity_coefficients`, and leaves ADR-0005, ADR-0012, ADR-0019,
ADR-0021 and ADR-0022 semantically unchanged - which is the point: this slice
is allowed to move the clock and nothing else.
