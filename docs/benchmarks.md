# Benchmarks and the robustness map

How performance and robustness are measured (ADR-0023, ADR-0027). Records live in `benchmarks/`.

`python -m chemthermo.bench --out record.json` (or `python tools/bench.py`)
runs a fixed nine-case workload over the reference flash paths - Peng-Robinson
flash, stability and a 24-state grid, NRTL `gamma-gamma`, modified-Raoult VLE
and the 364 K VLLE feed, PC-SAFT vapour-liquid, associating liquid-liquid and
the polymer split - and writes a record carrying, per case, the model, state,
composition, phase count, convergence criteria, initialization, iteration
counts, derivative mode, median wall time over timed repeats, peak allocation,
the machine and Python/numpy versions it was measured on, and a hash of the
accepted thermodynamic answer.
`python -m chemthermo.bench --compare before.json after.json` prints the
per-case speedup and **exits 1 if any result hash moved**, which is the whole
acceptance rule for a performance change here: same answer, better clock. Two
records are committed under `benchmarks/`, with `benchmarks/README.md`
explaining how to read them; see ADR-0023 for the policy and for the three
optimizations measured against them (Peng-Robinson flash 1.17x, PC-SAFT
liquid-liquid 1.32x, every result hash identical). The harness is internal: it
is not importable from `chemthermo` and is not part of the public API.

Since ADR-0030 each `flash_tp` / `stability_tp` call also carries a bounded
**call-local memo** of the density/compressibility root solves it makes, so a
state solved once in a call is not solved again in it. Nothing is cached on a
model - both model classes stay frozen and stateless, the memo lives in the
call and is discarded when it returns - and a hit returns the identical object
the first solve produced, which is what makes it bit-identical rather than
merely close. Measured as integer counts of the same fixed workload, and so
independent of the machine: **25.5 %, 13.5 % and 14.8 % fewer PC-SAFT density
solves** on the three PC-SAFT cases and **18.7 %, 13.1 % and 23.8 % fewer
cubic solves** on the three Peng-Robinson ones, with the three activity-model
cases unchanged at zero either way (they are the drift control). Every result
hash is identical and no fixture was regenerated.

`python -m chemthermo.bench robustness --out record.json` is the coverage
instrument next to that speed one: it sweeps all ten model families over fixed
state and composition grids (2505 states) and writes a classified record -
phase verdict or refusal class per state, with the exact state, message and
invariant residuals - so the next solver slice is chosen from counts rather
than recall. It found 36 refusals at `87f0820`, every one of them
polyethylene / n-pentane, and **0 at `9adf390`** after ADR-0028 retired them;
growing the grid to four more families (three-phase equation-of-state windows,
the deprecated `gamma-phi` path, near-critical Peng-Robinson states,
associating ternaries) found **14 more at `74820b8`**, none of them in the
original 2110 states. Nine of those fourteen were one class,
`multiphase-solver-failure`, and ADR-0029 repaired all nine - a removal the
phase-addition search can now take back, and a multiphase log-space stage for
a three-liquid set holding a component at `x ~ 1e-12` - leaving **5 at
`64831bd`**, all of them the deprecated `gamma-phi` path's own by-design
behaviour - and **5 again at `f852726`**, where the whole record is identical
field for field to the `64831bd` one because ADR-0030 is a performance change
and re-running the sweep is how that is proved over 2505 states. See
`benchmarks/README.md`, ADR-0027 (and its `robustness-map-coverage`
amendment), ADR-0028, ADR-0029 and ADR-0030.
