# Examples

Run examples from the repo root with `python`.

## Scope Policy

VLLE and PC-SAFT are in scope for Chemical-Thermodynamics. No thermodynamic capability class is categorically out of scope; implementation maturity may vary by module and release.

## Scripts

- `examples/database/01_component_database_demo.py`
  - Inspects canonical packaged runtime DB metadata and sample records.
- `examples/basic/tp_flash_pr_pure.py`
  - Pure-component TP flash with Peng-Robinson EOS.
- `examples/basic/tp_flash_pr_mixture.py`
  - Mixture TP flash with Peng-Robinson EOS.
- `examples/basic/tp_flash_nrtl_vle.py`
  - Gamma-phi TP flash with NRTL (liquid) + Peng-Robinson (vapor).
- `examples/basic/flash_tp_peng_robinson_demo.py`
  - Existing TP flash demo (Peng-Robinson EOS).
- `examples/basic/flash_tp_gamma_phi_demo.py`
  - Existing TP flash demo (gamma-phi, NRTL + Peng-Robinson).
- `examples/basic/stability_tp_peng_robinson_demo.py`
  - Michelsen tangent-plane phase stability at two states (unstable and stable).
- `examples/validation/00_reference_case.py`
  - Deterministic single-case comparison against `thermo` (optional dependency).
- `examples/validation/06_stability_vs_thermo.py`
  - Deterministic stability cross-check against `thermo`'s Michelsen test over 7 states.
- `examples/validation/*.py`
  - Optional validation sweeps against `thermo` (requires `pip install -e ".[validation]"`).
  - CSV output is disabled by default; pass `--outdir <dir>` or set `CHEMTHERMO_OUTDIR`.
  - Generated CSV outputs are intentionally not tracked in git.

## CLI quick runs

- Phi-phi TP flash:
  - `chemthermo tp-flash --components Methane,Ethane,Propane --z 0.5,0.3,0.2 --temperature-k 240 --pressure-pa 3000000 --format json`
- Gamma-phi TP flash (NRTL + Peng-Robinson):
  - `chemthermo tp-flash --components Methane,Ethane --z 0.5,0.5 --temperature-k 240 --pressure-pa 3000000 --flash-mode gamma-phi --format json`
- Gamma-phi coverage depends on available NRTL pair data; unsupported pairs return a runtime validation/model error.

## Expected output format

Each script prints:
- Header line describing the model/mode
- Temperature (K) and pressure (Pa)
- Feed composition `z` (mole fractions)
- Phase names and/or vapor fraction (beta)
- Phase compositions (`x` for liquid, `y` for vapor) when present
- Optional K-values (`y/x`) when both liquid and vapor are present
- Diagnostics block (raw keys/values)

All inputs and outputs are SI units. Compositions are mole fractions.
