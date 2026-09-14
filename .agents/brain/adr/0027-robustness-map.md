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

## Amendment (slice `robustness-map-coverage`)

Status: accepted
Date: 2026-09-14

### Context

ADR-0028 retired every refusal the original 2110-state map found and, in doing
so, discharged its own roadmap item 1 with an honest reading rather than a
victory: "a map that refuses nothing is not a solved solver" (brain.md section
10). The map's own structure said where it had stopped being hard - six of
eight refusal classes empty, no three-phase equation-of-state window,
`gamma-phi` not swept at all, the polymer family a single temperature at a
single `k_ij` - and this amendment is that list, worked.

### Decision

Four families join the six of decision 1, all in the same internal module
(`src/chemthermo/bench/robustness.py`), none touching a solver, model,
parameter file or public API:

| family | states | what it closes |
| --- | ---: | --- |
| `eos-three-phase` | 117 | "no three-phase EOS window is in the grid" - PC-SAFT water/n-hexane around `T3` (Case P-9), PC-SAFT and Peng-Robinson water/ethanol/n-hexane tie-triangles (Case P-10) |
| `gamma-phi-legacy` | 30 | "`gamma-phi` is not swept at all" - the deprecated `flash_mode="gamma-phi"` path (ADR-0010), NRTL-packaged Methane/Ethane, over the CLI's own contract state and a small T/P grid |
| `pr-near-critical` | 104 | near-critical Peng-Robinson states (a boundary located by bisecting `stability_tp`, not assumed) plus the CO2/n-decane and Methane/n-pentane windows that produced ADR-0016's and Case F-2's defects |
| `pcsaft-associating-ternary` | 144 | the polymer family was "one temperature at a single `k_ij`" for association coverage too - a second associating ternary (water/1-propanol/n-hexane) and a finer feed grid on both, aimed at cloud points |

Total grid: **2505 states** (2110 + 395), `--quick` **224 states** (171 + 53,
~14.2 s against the 15 s budget). Each grid is documented at its definition in
the module, the same way as the original six; nothing here changes how a state
is classified (`classify_refusal`, `_REFUSAL_RULES`) - the four new families
exercise rules that were already declared and already dormant
(`multiphase-solver-failure`, `rr-no-bracket`), which is itself evidence for
decision 2's claim that the rule list covers what occurs.

One check is added to `_invariant_violations`: a three-phase answer's
`delta_g_vs_two_phase_rt` (ADR-0020, present only on a result that entered the
phase addition/removal search and returned three phases) must be negative,
mirroring the existing `delta_g_split_rt < 0` check one level up. It is
dormant on every family but `eos-three-phase` and the rare state elsewhere
that lands on three phases by chance.

### The `gamma-phi-legacy` family is read differently from the other nine

`flash_mode="gamma-phi"` has had no stability test since ADR-0007 and stays on
the pre-`flash-auto-phase-detection` Wilson heuristic by design (ADR-0008
decision 4); `_invariant_violations`'s post-split check is conditioned on
`diagnostics["post_split_checked"]`, which this path never sets
(`_legacy.py`), so it is silently and correctly excluded rather than exempted
by a family-specific carve-out. Its `rr-no-bracket` refusals are not a defect
this slice diagnoses further: they are ADR-0016's own motivating example,
reproduced on the one code path ADR-0016 deliberately left unrepaired ("the
legacy `phase_detection="wilson-heuristic"` path is untouched... because
reproducing pre-ADR-0008 behavior is that path's entire purpose").

### What the new grids found

**2491 of 2505 converge, 14 refuse, 0 converge and violate an invariant**, in
2232.8 s (37:13); every one of the original 2110 states is unchanged, checked
field by field against the committed `9adf390` record. All 14 refusals are in
two of the four new families: `eos-three-phase` (9, `multiphase-solver-failure`
- the ADR-0011/ADR-0020 phase addition/removal search either collapsing a
converged set to a non-positive phase fraction or leaving the multiphase split
short of `tol`) and `gamma-phi-legacy` (5, `rr-no-bracket` - the legacy Wilson
heuristic's own documented failure mode, ADR-0016 decision 8, reproduced by
hand-trace as a single-step `K` compression rather than ADR-0016's diverging
oscillation). `pr-near-critical` and `pcsaft-associating-ternary` refuse
**nothing** - 248 states between them, the latter finding 6 verified `VLLE`
answers near a cloud point at no cost. One of the nine `eos-three-phase`
refusals is not new: `(0.1, 0.1, 0.8)` at 333 K is the same feed Case P-10 (i)
already recorded as a pre-existing, non-regression failure. See ledger Case
R-MAP-2 for the ranked table, every example state and message, and the
diagnosis of each class from evidence already in the record.

### Consequences

- `benchmarks/robustness_74820b8.json` / `.md` regenerated at the new grid
  size; the superseded `9adf390` JSON pruned per the existing policy, its
  `.md` kept.
- `tests/test_robustness_map.py`'s `EXPECTED_QUICK_FAMILIES` gains four rows,
  `PINNED_REFUSALS` gains the eight cheap refusals the quick subset contains
  (one binary-scan collapse, two Peng-Robinson-ternary splits, five legacy
  `rr-no-bracket` states), and two new tests pin the three-phase verdicts
  Cases P-9 and P-10 already established independently (`T3 - 0.05 K -> LLE`,
  `T3 + 0.05 K -> VLE`, two of the twelve 333 K tie-triangle feeds `-> VLLE`),
  marked `slow` as repetitions of capability `test_flash_vlle_eos.py` already
  covers by default or behind `slow`.
- The brain roadmap is re-ranked again from this result; see brain.md section
  9 and validation Case R-MAP-2.
- **Cost.** A PC-SAFT ternary VLLE state costs 13-40 s (ADR-0020's own "~35 s"
  measurement, reproduced here: the family's 117 states cost 651.6 s, and the
  associating-ternary family's 144 states cost 670.6 s). The full sweep grew
  from 893.1 s to 2232.8 s accordingly. The `--quick` subset avoids every such
  state but one cheap single-phase pin per PC-SAFT-ternary system, which is
  why the family's quick bucket counts are not representative of its
  full-sweep verdict mix.

## Superseded by (optional)
None.
