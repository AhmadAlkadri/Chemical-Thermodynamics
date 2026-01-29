# Notebooks

These notebooks demonstrate TP flash calculations using the current public API.

## Setup

From the repo root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
# Optional validation dependencies
pip install -e ".[validation]"
# Notebook environment
pip install notebook
```

## Notebook output stripping

This repo uses `nbstripout` as a git filter to strip notebook outputs on commit
while keeping outputs visible locally.

Install the tooling (once):

```bash
pip install -e ".[dev]"
python -m nbstripout --install --attributes .gitattributes
```

Notes:
- Running a notebook will add outputs and `git status` will show the file as modified.
- Committing will automatically strip outputs, keeping the repository clean.

## Contents

- `notebooks/01_tp_flash_peng_robinson.ipynb`
- `notebooks/02_tp_flash_nrtl.ipynb`
- `notebooks/validation/README.md`
- `notebooks/validation/01_flash_pr_vs_thermo.ipynb`
- `notebooks/validation/02_flash_nrtl_vs_thermo.ipynb`

Run with:

```bash
jupyter notebook
```
