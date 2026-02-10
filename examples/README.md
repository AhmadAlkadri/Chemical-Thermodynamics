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
- `examples/validation/00_reference_case.py`
  - Deterministic single-case comparison against `thermo` (optional dependency).
- `examples/validation/*.py`
  - Optional validation sweeps against `thermo` (requires `pip install -e ".[validation]"`).
  - CSV output is disabled by default; pass `--outdir <dir>` or set `CHEMTHERMO_OUTDIR`.
  - Generated CSV outputs are intentionally not tracked in git.

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
