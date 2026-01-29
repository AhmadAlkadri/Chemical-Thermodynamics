# Validation notebooks

These notebooks mirror the automated validation tests in `tests/validation/` so you can compare
chemthermo against reference libraries interactively.

## Setup

From the repo root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
# Validation dependencies (thermo)
pip install -e ".[validation]"
# Notebook environment
pip install notebook
```

If you plan to run any CoolProp-based checks in the future, install it separately:

```bash
pip install CoolProp
```

## What the notebooks compare

- **EOS flash validation:** chemthermo TP flash (Peng-Robinson) versus `thermo` PR flash.
- **Activity model validation:** chemthermo NRTL activity coefficients versus `thermo`'s NRTL helper.

Each notebook prints compact tables for inputs, results, and differences. The final cell contains
assertions that use the same tolerances as `tests/validation/test_flash_vs_thermo.py`. A failure
means the comparison exceeded those tolerances; a skip message indicates a single-phase result
where a two-phase comparison was expected.

## Contents

- `notebooks/validation/01_flash_pr_vs_thermo.ipynb`
- `notebooks/validation/02_flash_nrtl_vs_thermo.ipynb`

Run with:

```bash
jupyter notebook
```
