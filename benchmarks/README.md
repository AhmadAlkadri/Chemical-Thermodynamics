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
| `baseline_<sha>.json` | the state of the reference paths before the `perf-baseline-and-root-reuse` optimizations |
| `after_<sha>.json` | the same workload after them |

Both were measured on the machine named in their own `environment` block.
**Wall times are not portable**: compare records measured on the same machine,
back to back, or compare nothing. The `result_hash` *is* portable and is what
a comparison on a different machine can still assert.

## Adding a case

Adding a case to `chemthermo/bench/_cases.py` is a normal change: existing
records simply do not carry it, and `--compare` reports it as "after only".
Editing or removing a case invalidates every committed baseline for that case;
say so here and regenerate.

## The polymer case

`pcsaft-polymer-lle` reads its parameters from
`tests/fixtures/pcsaft/martini2009_polymers.json`, because chemthermo packages
no polymer PC-SAFT parameters (ADR-0022). Outside a source checkout that file
is absent and the case records `status: "skipped"` instead of failing the run.
