# Command line (`chemthermo`)

The `chemthermo` command (also `python -m chemthermo`). Contract: ADR-0003, ADR-0004, ADR-0033 in `.agents/brain/adr/`.

Run TP flash without writing Python:

```bash
chemthermo tp-flash --components Methane,Ethane,Propane --z 0.5,0.3,0.2 --temperature-k 240 --pressure-pa 3000000 --format json
```

Run gamma-phi mode (NRTL liquid activity + Peng-Robinson EOS):

```bash
chemthermo tp-flash --components Methane,Ethane --z 0.5,0.5 --temperature-k 240 --pressure-pa 3000000 --flash-mode gamma-phi --format json
```

Module execution is also supported:

```bash
python -m chemthermo tp-flash --components Methane,Ethane --z 0.5,0.5 --temperature-k 240 --pressure-pa 3000000
```

Notes:
- `--flash-mode` defaults to `phi-phi`; valid choices are `phi-phi` and `gamma-phi`.
- Gamma-phi mode currently uses `NRTL()` for the liquid activity model with the
  packaged **synthetic** pair parameters (see "NRTL activity coefficients").
- NRTL pair coverage is data-dependent; missing pair data returns a runtime validation/model error.

Test a feed for phase stability (tangent-plane analysis, ADR-0033):

```bash
chemthermo stability-tp --components Methane,n-Hexane --z 0.5,0.5 --temperature-k 300 --pressure-pa 2000000
chemthermo stability-tp --components Methane,n-Hexane --z 0.5,0.5 --temperature-k 300 --pressure-pa 2000000 --eos pc-saft --format json
```

`status` is `stable`, `unstable` or `inconclusive`. **`stable` means no
negative tangent-plane distance was found from a deterministic trial set - not
a global proof** - and every JSON payload says so in
`result.stability_scope = "bounded-trial-set"`. An unstable feed reports its
incipient phase (`trial_composition`, `k_values = w / z`, `phase_branch`).

Choose the equation of state and the phase budget of `tp-flash`:

```bash
# PC-SAFT (packaged parameters, kij = 0, phi-phi only)
chemthermo tp-flash --components Methane,n-Hexane --z 0.5,0.5 --temperature-k 300 --pressure-pa 2000000 --eos pc-saft --format json
# three liquids (Peng-Robinson water / ethanol / n-hexane); --max-phases defaults to 3
chemthermo tp-flash --components Water,Ethanol,n-Hexane --z 0.2,0.4,0.4 --temperature-k 280 --pressure-pa 101325 --max-phases 3
```

More in `examples/cli/`. Contract (ADR-0003, 0004, 0033):
- `cli_schema_version` is `1` for every command. New commands and flags are
  additive: an invocation that worked before prints exactly the same bytes, and
  keys a new flag adds (`solver.settings.max_phases`) appear only when it is given.
- Any number of phases fits the v1 layout: `result.phases` and
  `result.phase_fractions` are keyed by phase name (`liquid1`, `liquid2`, ...),
  and `vapor_fraction` is `null` when there is no vapour.
- Exit codes: `0` success; `1` validation or model error (unknown component,
  missing parameters); `2` usage error (including `--eos pc-saft` with
  `--flash-mode gamma-phi`, or `--max-phases 0`); `3` solver non-convergence -
  for `stability-tp`, an `inconclusive` verdict, with the JSON still printed.
- Not exposed yet: `modified-raoult` / `gamma-gamma` modes, `kij` or parameter
  files, polymers.
