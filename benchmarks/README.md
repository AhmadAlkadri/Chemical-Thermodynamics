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
| `baseline_d8814de.json` | the same workload before the ADR-0030 call-local solve memo |
| `after_07bf0b3.json` | the same workload after it |

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

## The ADR-0030 pair, and a measurement that had to be made differently

`baseline_d8814de.json` / `after_07bf0b3.json` measure the call-local
equation-of-state solve memo of ADR-0030. They were taken back to back on one
machine, as the rule requires, and **the machine was not quiet**: an unrelated
process held about eight of its twelve cores for the whole session. Read the
three activity-model controls first, because they are what that looks like -
they touch none of the code this slice changes and they still scatter across
1.02x, 0.96x and 1.03x at nine timed repeats, against the
1.00x / 1.04x / 1.02x the same three read for ADR-0023 on a quiet machine. A
first attempt at five repeats gave controls of 0.77x / 1.25x / 0.92x and was
discarded as unreadable.

```
case                           before / s    after / s   speedup  result
------------------------------------------------------------------------------
modified-raoult-vle              0.024677     0.024191     1.02x  identical
nrtl-lle-tessier-p1              0.057287     0.059526     0.96x  identical
pcsaft-lle-water-hexane          1.595234     1.324830     1.20x  identical
pcsaft-polymer-lle               0.994039     0.851818     1.17x  identical
pcsaft-vle-methane-hexane        0.216429     0.178299     1.21x  identical
pr-flash-grid-24                 0.158998     0.120641     1.32x  identical
pr-flash-ternary                 0.016810     0.013450     1.25x  identical
pr-stability-ternary             0.004402     0.004119     1.07x  identical
vlle-364k                        0.137144     0.133123     1.03x  identical

result hashes identical across every shared case
```

A second measurement was made because of that, and it is the better estimate of
what the change is worth: **the same workload run in one process with the two
arms interleaved repeat by repeat**, the only difference between them being
whether the models' memo lookup returns the call's memo or `None` (the
pre-ADR-0030 path). Both arms then take the same contention, milliseconds
apart, 11 repeats each:

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

Where the two disagree, believe the interleaved one. The committed pair's
Peng-Robinson ratios (1.25x and 1.32x) are larger than the interleaved
measurement supports (1.03x and 0.99x, and 1.10x / 1.04x on a second
31-repeat run): the cubic solve is about a sixth of a Peng-Robinson
`fugacity_coefficients` call - ADR-0023 measured 8.6 us of 34.8 us - so
removing a quarter of the solves cannot be worth 30 %. The PC-SAFT cases agree
between the two measurements to within their spreads: **1.10x - 1.30x**, and
that is the honest range.

**All nine result hashes are identical**, which is the part of the comparison a
loaded machine cannot corrupt, and it is the part the acceptance rule is
actually about.

Because the wall times are not trustworthy here, this slice's ratio is stated
from something that is: the **number of solves** the same fixed workload makes,
counted with the models' memo lookup returning `None` and with it active. These
are integers, they are a property of the workload rather than of the machine,
and they reproduce exactly on any hardware.

| case | PC-SAFT density solves off -> on | PR cubic solves off -> on |
| --- | ---: | ---: |
| `pr-flash-ternary` | - | 230 -> 187 (-18.7 %) |
| `pr-stability-ternary` | - | 84 -> 73 (-13.1 %) |
| `pr-flash-grid-24` | - | 2203 -> 1679 (-23.8 %) |
| `nrtl-lle-tessier-p1` (control) | 0 -> 0 | 0 -> 0 |
| `modified-raoult-vle` (control) | 0 -> 0 | 0 -> 0 |
| `vlle-364k` (control) | 0 -> 0 | 0 -> 0 |
| `pcsaft-vle-methane-hexane` | 137 -> 102 (-25.5 %) | - |
| `pcsaft-lle-water-hexane` | 275 -> 238 (-13.5 %) | - |
| `pcsaft-polymer-lle` | 601 -> 512 (-14.8 %) | - |

The three controls making zero solves either way is the control working: no
activity-coefficient evaluation is memoized, deliberately, so that they go on
running the same code before and after (ADR-0030 decision 6).

A **three-phase equation-of-state case was considered for the workload and not
added**: one such state costs 13-40 s, so at five timed repeats plus a warm-up
it would turn this half-minute instrument into a five-minute one. ADR-0030
measures the three-phase effect directly instead, and ledger Case B-2 records
it.

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
| `--family NAME` | run one of `pr-phi-phi`, `pcsaft`, `pcsaft-associating`, `modified-raoult`, `gamma-gamma`, `polymer`, `eos-three-phase`, `gamma-phi-legacy`, `pr-near-critical`, `pcsaft-associating-ternary` - the sweep partitions and resumes |
| `--quick` | the 224-state cost-bounded subset across all ten families (~14 s), which is what `tests/test_robustness_map.py` runs |
| `--list` | the grid: every system, its state count and its quick count |
| `--out` / `--summary-out` | the JSON record / the Markdown summary table |
| `--quiet` | no per-system progress line |

It flashes 2505 fixed states and puts each one in exactly **one** bucket: a
phase verdict (`single-liquid`, `single-vapor`, `VLE`, `LLE`, `VLLE`, `LLL`),
`converged-invariant-violated`, or one of eight refusal classes read off the
exception type and message. The exact state, the exact message and the
invariant residuals are recorded per state. Four families -
`eos-three-phase`, `gamma-phi-legacy`, `pr-near-critical` and
`pcsaft-associating-ternary` - were added by slice `robustness-map-coverage`
(ADR-0027 amendment) to close the gaps ADR-0027's own roadmap named: no
three-phase equation-of-state window, `gamma-phi` unswept, no near-critical
Peng-Robinson coverage, and one associating ternary out of the packaged four.

**It is not a correctness check** - nothing here is compared against a
published number or another implementation; `.agents/brain/validation-cases.md`
is where that lives. It is not a bit-identity fixture either: phase
*compositions* are deliberately not recorded, because `baseline_*.json` /
`after_*.json` and `refactor_bit_identity_v3.json` are the baselines and a
second one would have to be regenerated whenever a last bit moved. Read
ADR-0027 before using it to justify anything.

| file | what it is |
| --- | --- |
| `robustness_761fd57.json` | the full 2505-state sweep at `761fd57` (after ADR-0035/0036; Linux x86_64, CPython 3.11.15, numpy 2.4.6) |
| `robustness_761fd57.md` | its summary table |
| `robustness_f852726.md` | the summary of the superseded `f852726` sweep (macOS arm64, ADR-0030) |
| `robustness_64831bd.md` | the summary of the superseded `64831bd` sweep, the *before* measurement for ADR-0030 |
| `robustness_74820b8.md` | the summary of the superseded `74820b8` sweep, the *before* measurement for ADR-0029 |
| `robustness_9adf390.md` | the summary of the superseded `9adf390` (2110-state) sweep |
| `robustness_87f0820.md` | the summary of the superseded `87f0820` sweep |

At `761fd57` (the first record taken on Linux): **2500 of 2505 states
converge, 5 refuse (the deprecated gamma-phi path's `rr-no-bracket`, by
design), 0 converge and violate an invariant** - the same totals and, matched
by `(system, state_index)`, the same bucket for every state as the macOS
`f852726` record; a before/after pair on one machine showed ADR-0036 changes
no state (ledger Case P-17 "ADR-0036"). Records from different machines must
be matched by index: generated grid pressures can differ by 1 ULP.

At `f852726`: **2500 of 2505 states converge, 5 refuse, 0 converge and violate an invariant**, unchanged from `64831bd` - the ADR-0030 call-local solve
memo is a performance change and the whole sweep was re-run to say so. The
comparison is **0 differences over all 2505 states and every aggregate field**,
wall time excluded, which is the strongest of this repository's three
bit-identity gates because it is the widest. The 9 `multiphase-solver-failure`
states `eos-three-phase` refused at `74820b8` were repaired by ADR-0029 (ledger
Case P-18), so what is left is the deprecated `gamma-phi-legacy` path's 5
`rr-no-bracket` states, which are by design (ADR-0016 decision 8) and are
pinned as such rather than as a defect queue.

The sweep itself read **2017.3 s (33:37)** against 2202.5 s (36:43) at
`64831bd`. Both were taken on a loaded machine in different sessions, so read
that pair as indicative of the direction and not as a ratio.

`robustness_64831bd.json`, `robustness_74820b8.json` and
`robustness_9adf390.json` were pruned when their successors superseded them,
per the pruning policy at the end of this section; their `.md` summaries are
kept, because the counts in them are what ADR-0027, ADR-0028, ADR-0029,
ADR-0030 and Cases R-MAP-1 / R-MAP-2 / P-17 / P-18 quote.

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
