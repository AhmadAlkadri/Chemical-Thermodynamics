# ADR-0027: A robustness map is the evidence a *solver* slice is picked on

Status: accepted
Date: 2026-09-14

## Context
ADR-0023 gave this repository an instrument for the clock. There was none for
coverage. "Where does chemthermo still refuse?" was answered by recall - the
last thing a slice happened to hit - and the answers were scattered: the
188-state Case F-4 grid lives in one test, the 68-state polymer sweep in an
example, the 1144-state Peng-Robinson scan only as prose in Case F-1 (no file
in this repository lists its eight mixtures). Three roadmap items in a row were
chosen because a polymer state was in front of someone, not because a count
said the polymer path was where the refusals are.

Two consequences of that were visible before this slice started. The roadmap's
item 1 named a **three-state** band near 5.2 MPa on one molar mass of one
polymer, which is a real defect and is also, by count, unlikely to be the
largest one. And nothing anywhere counted the outcome that is strictly worse
than a refusal: a flash that *returns* and violates an invariant. A refusal is
loud; a wrong answer with a positive `dG_split` is not.

## Decision

### 1. `chemthermo.bench.robustness` produces the map
An internal module inside the existing internal harness package (ADR-0001,
ADR-0023: nothing re-exported from `chemthermo`), run as

```bash
python -m chemthermo.bench robustness --out benchmarks/robustness_<sha>.json
```

It sweeps six families over **fixed** state and composition grids - 2110 states
- and records, per state: family, system, components, feed, `(T, P)`, wall
time, and exactly one classification bucket.

`--family NAME` restricts the sweep to one family, so it partitions and
resumes. `--quick` runs a fixed ~150-state subset across all six, which is what
the default test suite runs. `--list` prints the grid.

The subcommand is dispatched on the first argument rather than through
`add_subparsers`, because the timing CLI's flags were its entire contract
before this slice existed (`benchmarks/README.md`, `tools/bench.py`,
`tests/test_bench_harness.py`) and `python -m chemthermo.bench --out x.json`
has to keep meaning what it meant.

### 2. Every state lands in exactly one bucket
A **verdict** (`single-liquid`, `single-vapor`, `VLE`, `LLE`, `VLLE`, `LLL`,
`other`) when `flash_tp` returned and every invariant held;
`converged-invariant-violated` when it returned and one did not; or one of
eight **refusal classes** when it raised:

| class | what it means |
| --- | --- |
| `stability-inconclusive` | no tangent-plane trial converged (feed or post-split), or an unstable feed produced no minimizer |
| `rr-no-bracket` | the K-values do not bracket a Rachford-Rice root, on either seed |
| `density-root-failure` | no admissible density/compressibility root, or only mechanically unstable ones |
| `post-split-third-phase` | the converged set is not stable and `max_phases` forbids another |
| `multiphase-solver-failure` | the multiphase Rachford-Rice, split or add/remove search failed |
| `split-non-convergence` | a split ran out of iterations; the *stage* is recorded (`phi-phi`, `log-space`, `gamma-gamma`, `modified-raoult`, `beta-outside-window`, ...) |
| `model-error` | `ModelError` / `PropertyNotFoundError` / range or composition errors |
| `other-refusal` | anything the rules do not name; an empty bucket is the claim that the rules cover what occurs |

The class is read off the exception **type and message**, by an ordered rule
list in one function. That is a deliberate tradeoff: it is coupled to message
text and will need editing when a message is reworded, and in exchange it needs
no change to any solver, which this slice is forbidden from touching. The exact
message and state are recorded next to every classification, so a class is
where a diagnosis starts and never what it concludes.

### 3. Invariants are checked on every converged answer, and a violation outranks a refusal
Mass balance (`< 1e-10`) and the phase composition sums are **recomputed** here
from the returned phase fractions and compositions. The equilibrium residual
(`< 1e-6`), `dG_split/RT < 0` and the post-split verdict are read from the
solver's own diagnostics - recomputing equal fugacities would need per-model
code this module deliberately does not carry, and that limit is stated rather
than hidden. Phase fractions must lie strictly in `(0, 1)` and sum to one.

A state that converges and violates one of these is counted in its own bucket,
not among the verdicts, because it is worse than a refusal: a refusal tells the
caller something went wrong.

### 4. Every refusal class found is pinned by a test that expects today's failure
`tests/test_robustness_map.py` runs the `--quick` subset and asserts the state
count, the per-family bucket counts, and that each pinned state still refuses
with the same class. The tests are written so that a **fix** breaks them. That
is the point: this repository has repeatedly rediscovered a refusal it had
already seen, and a pinned expectation converts "it still fails" into a fact
with a date on it.

## What this measures, and what it does not
- **Not a correctness check.** Nothing here is compared against a published
  number or an independent implementation; the ledger (F-1..F-5, L-1..L-4,
  P-0..P-16, R-1..R-4, V-1..V-5) is where that lives. What this checks is
  internal consistency, which catches a wrong answer only when it is wrong in
  one of the five ways listed above.
- **Not a bit-identity fixture.** Phase *compositions* are deliberately absent
  from the record. `benchmarks/*.json` (ADR-0023) and
  `refactor_bit_identity_v3.json` are the baselines; a second one that had to
  be regenerated whenever a last bit moved would make a coverage map expensive
  to keep.
- **Not portable wall times.** Per-state times are recorded so the *shape* of
  the cost is visible (which family, which state), not so two machines can be
  compared - `benchmarks/README.md` says why.
- **Never better than the stability test underneath it.** A state counted
  `single-*` is one where the deterministic trial set found no negative
  tangent-plane distance (brain.md section 10). The map inherits that bound and
  cannot detect a miss.
- **Two grids carry values that are not cited.** The two nonzero PC-SAFT `kij`
  binaries use *illustrative* values (chemthermo packages no `kij` dataset,
  ADR-0014), and seven of the sixteen Tessier feeds are midpoints this module
  added for composition coverage. Both are labelled at their definitions and in
  the record. Neither is compared against anything.

## How it is meant to be used to pick a slice
1. Run it at HEAD, commit `benchmarks/robustness_<sha>.json` and the `.md`
   summary next to it.
2. Rank the refusal classes by count, then by how many *families* each spans -
   a class that appears in one system of one family is a state, a class that
   appears in three families is a mechanism.
3. Diagnose the top class from evidence already in the record (the message, the
   stage, the neighbouring states that converged), then from a stability trial
   or a Gibbs comparison at one of its states. Only then write the slice.
4. When the slice lands, re-run and diff: the class it targeted should shrink
   and nothing else should grow.

The map does **not** rank by importance. A class with 30 states in one polymer
system may still be the right next slice if those 30 states are the only ones
anybody asks for. It ranks by *evidence*, and the judgment stays with the
maintainer.

## Alternatives considered
- **Extend `chemthermo.bench`'s existing case list instead.** Rejected: the
  ADR-0023 workload is *fixed* by decision 5 of that ADR, because a benchmark
  whose workload drifts measures nothing, and a coverage sweep wants to grow
  whenever a new system is interesting. Different lifetime, different artefact,
  same package.
- **Classify from exception types alone.** Too coarse to be useful:
  `ConvergenceError` covers the stability verdict, both Rachford-Rice failures,
  four split stages, the multiphase search and the post-split test - which is
  every interesting distinction collapsed into one bucket.
- **Add structured fields to the exceptions and classify on those.** The right
  long-term answer, and out of scope here: it is a public-API change to a
  package whose error types are documented (ADR-0001), and this slice may not
  touch solver code. Recorded as the follow-up if message-matching proves
  brittle.
- **Recompute the equilibrium residual here from the models.** Would make the
  invariant check independent of the solver's own reporting. Rejected for this
  slice: it needs a per-model fugacity/activity evaluation path that duplicates
  `flash/_verify.py`, and duplicating a verifier to check the verifier is how
  two of them drift apart.
- **Run the full sweep in CI.** It is minutes, not seconds. The `--quick`
  subset runs by default and the full sweep is a committed record, which is the
  same shape as the `slow` marker policy in `.agents/dev-contract.md`.
- **Make the sweep parallel.** Rejected here: wall time per state is one of the
  things recorded, and a process pool would make those numbers a property of
  the scheduler. The sweep partitions by `--family` instead.

## Consequences
- New internal module `src/chemthermo/bench/robustness.py` and one subcommand.
  No public API, no solver, no model and no existing test changed.
- `benchmarks/robustness_<sha>.json` and `benchmarks/robustness_<sha>.md` are
  committed, and `benchmarks/README.md` gains a section.
- `tests/test_robustness_map.py` runs the quick subset by default (a few
  seconds) and the full sweep behind `slow`.
- The brain roadmap's item 1 is replaced by the evidence-ranked list, and
  validation Case R-MAP-1 carries the table.
- The classification rules are coupled to exception message text. When a
  message is reworded, `other-refusal` will grow and the pinned tests will
  notice - which is the failure mode this design accepts, chosen over touching
  the solvers.

## Supersedes (optional)
None. Extends the ADR-0023 harness package with a second artefact; leaves the
benchmark workload, its records and its acceptance rule untouched.
