# Benchmarks

Reproducible benchmark records for the reference flash paths. The harness is
`chemthermo.bench` (internal, not public API); the decision it implements is
ADR-0023.

## Running it

```bash
python -m chemthermo.bench --out benchmarks/mine.json
# or, from a source checkout without an install:
python tools/bench.py --out benchmarks/mine.json
```

Useful flags:

| flag | meaning |
| --- | --- |
| `--repeats N` | timed repeats per case after one untimed warm-up (default 5) |
| `--case ID` | run one case (repeatable); `--list` prints the ids |
| `--out PATH` | write the JSON record (the summary table always prints) |
| `--compare BEFORE AFTER` | compare two records instead of running |

The whole workload takes about half a minute at the default repeat count.

## Comparing two records

```bash
python -m chemthermo.bench --compare benchmarks/baseline_<sha>.json benchmarks/after_<sha>.json
```

prints one line per case - median wall time before, after, and the ratio - and
**exits 1** if any case's `result_hash` differs. That is the acceptance rule
for a performance change in this repository:

> Same `result_hash` on every case, a better median on the cases the change
> was about, and no case materially worse.

The hash is taken over the *accepted thermodynamic answer* - phase names, every
phase composition and every phase fraction, at twelve significant digits - plus
the state and the refusal status. It deliberately does **not** cover wall time,
memory or iteration counts, so a record whose only change is the clock hashes
the same.

A hash match is a check, not the bit-identity proof. Bit-identity is proved by
`tests/test_flash_refactor_bit_identity.py` (the 155-state pinned fixture,
floats compared with `==`) and, for the equation-of-state branch capability, by
`tests/test_eos_branch_reuse.py`.

## What a record contains

Per case: model, components, state `(T, P)`, overall composition, phase count
obtained, the phase names / compositions / fractions, the convergence criteria
in force (every `FlashSettings` and `StabilitySettings` field that can move an
answer), the initialization the solver reported (`k_seed`, `phase_detection`,
`feed_branch`, `incipient_phase`), the iteration counts it reported, the
derivative mode (`"analytic"` - chemthermo differentiates no model by finite
difference), the median / min / max wall time over the timed repeats, the peak
`tracemalloc` allocation of one separate instrumented run, the `result_hash`,
and the refusal status of any state the solver declined.

Per record: schema tag, UTC timestamp, the git commit and whether the tree was
dirty, and the environment (platform, CPU string, core count, Python and numpy
versions).

Workload preparation - databank reads, parameter resolution, mixture
construction - happens **outside** the timed region. What is timed is
`flash_tp` / `stability_tp`.

## The committed records

| file | what it is |
| --- | --- |
| `baseline_3ce68df.json` | the reference paths before the `perf-baseline-and-root-reuse` optimizations |
| `after_2ca41bf.json` | the same workload after them |

Both were measured on the machine named in their own `environment` block
(Apple M2 Max, Python 3.11.6, numpy 2.4.2), five timed repeats per case, **47
seconds apart** - the baseline from a detached worktree at `3ce68df`, the after
record from `2ca41bf`. That matters: a baseline taken twenty minutes earlier in
the same session read 5-10 % slower on cases whose code had not changed, which
is larger than two of the three optimizations being measured. Compare records
taken back to back on one machine, or compare nothing.

```
case                           before / s    after / s   speedup  result
------------------------------------------------------------------------------
modified-raoult-vle              0.020264     0.019569     1.04x  identical
nrtl-lle-tessier-p1              0.043876     0.043711     1.00x  identical
pcsaft-lle-water-hexane          1.626994     1.236466     1.32x  identical
pcsaft-polymer-lle               0.962138     0.814358     1.18x  identical
pcsaft-vle-methane-hexane        0.240769     0.188477     1.28x  identical
pr-flash-grid-24                 0.140889     0.117589     1.20x  identical
pr-flash-ternary                 0.015261     0.013040     1.17x  identical
pr-stability-ternary             0.004649     0.003881     1.20x  identical
vlle-364k                        0.099957     0.097971     1.02x  identical
```

The three activity-model cases are the **drift control**: none of the three
optimizations touches the code they run, so their 1.00x / 1.04x / 1.02x is what
this machine's noise looks like over one pair of runs. Read the
equation-of-state ratios against that, not against 1.00x exactly.

**Wall times are not portable** between machines. The `result_hash` *is*, and
is what a comparison on a different machine can still assert.

## Adding a case

Adding a case to `chemthermo/bench/_cases.py` is a normal change: existing
records simply do not carry it, and `--compare` reports it as "after only".
Editing or removing a case invalidates every committed baseline for that case;
say so here and regenerate.

## The robustness map (ADR-0027)

A second artefact lives here, and it measures coverage rather than speed:

```bash
python -m chemthermo.bench robustness --out benchmarks/robustness_<sha>.json \
                                      --summary-out benchmarks/robustness_<sha>.md
```

| flag | meaning |
| --- | --- |
| `--family NAME` | run one of `pr-phi-phi`, `pcsaft`, `pcsaft-associating`, `modified-raoult`, `gamma-gamma`, `polymer` - the sweep partitions and resumes |
| `--quick` | the 171-state cost-bounded subset across all six families (~9 s), which is what `tests/test_robustness_map.py` runs |
| `--list` | the grid: every system, its state count and its quick count |
| `--out` / `--summary-out` | the JSON record / the Markdown summary table |
| `--quiet` | no per-system progress line |

It flashes 2110 fixed states and puts each one in exactly **one** bucket: a
phase verdict (`single-liquid`, `single-vapor`, `VLE`, `LLE`, `VLLE`, `LLL`),
`converged-invariant-violated`, or one of eight refusal classes read off the
exception type and message. The exact state, the exact message and the
invariant residuals are recorded per state.

**It is not a correctness check** - nothing here is compared against a
published number or another implementation; `.agents/brain/validation-cases.md`
is where that lives. It is not a bit-identity fixture either: phase
*compositions* are deliberately not recorded, because `baseline_*.json` /
`after_*.json` and `refactor_bit_identity_v3.json` are the baselines and a
second one would have to be regenerated whenever a last bit moved. Read
ADR-0027 before using it to justify anything.

| file | what it is |
| --- | --- |
| `robustness_87f0820.json` | the full 2110-state sweep at `87f0820` |
| `robustness_87f0820.md` | its summary table |

At `87f0820`: **2074 of 2110 states converge, 36 refuse, 0 converge and
violate an invariant.** Every refusal is in the `polymer` family; the whole
Peng-Robinson, PC-SAFT, associating-PC-SAFT, modified-Raoult and
`gamma-gamma` sweep (1858 states) refuses nothing. See ledger Case R-MAP-1 for
the ranked classes and the diagnoses.

Adding a system or a state to `chemthermo/bench/robustness.py` is a normal
change - unlike the timing workload above, this grid is *meant* to grow -
but it moves the committed counts, so regenerate the record and the
`tests/test_robustness_map.py` expectation in the same commit. The JSON is
about 1.8 MB (one entry per state, phase compositions excluded); regenerating
it at a new commit means a new file, so prune the superseded one rather than
accumulating a sweep per slice.

## The polymer case

`pcsaft-polymer-lle` reads its parameters from
`tests/fixtures/pcsaft/martini2009_polymers.json`, because chemthermo packages
no polymer PC-SAFT parameters (ADR-0022). Outside a source checkout that file
is absent and the case records `status: "skipped"` instead of failing the run.
