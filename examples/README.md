# Examples

Run examples from the repo root with `python`.

## Scripts

- `examples/tp_flash_pr_pure.py`
  - Pure-component TP flash with Peng-Robinson EOS.
- `examples/tp_flash_pr_mixture.py`
  - Mixture TP flash with Peng-Robinson EOS.
- `examples/tp_flash_nrtl_vle.py`
  - Gamma-phi TP flash with NRTL (liquid) + Peng-Robinson (vapor).
- `examples/flash_tp_peng_robinson_demo.py`
  - Existing TP flash demo (Peng-Robinson EOS).
- `examples/flash_tp_gamma_phi_demo.py`
  - Existing TP flash demo (gamma-phi, NRTL + Peng-Robinson).
- `examples/flash_tp_vlle_demo.py`
  - Existing VLLE scaffold demo.

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
