---
name: chemthermo-change-loop
description: Use when implementing or reviewing changes in this Chemical-Thermodynamics repository. Apply the repo contract by reading .agents/brain/brain.md first, preserving SI-unit and public-API invariants, shipping thin vertical slices with a runnable golden path, running CI-equivalent quality gates (ruff format, ruff format --check, ruff check, pyright, pytest), and updating ADR/brain docs when architecture or public API changes.
---

# Chemthermo Change Loop

## Inputs
- User request and acceptance criteria
- Files expected to change
- Whether public API or architecture boundaries are affected
- Local environment status (Python 3.11+, editable install)

## Workflow

1. Establish scope and constraints.
- Read `README.md` and `.agents/brain/brain.md` before structural edits.
- If the change is architectural or API-facing, read relevant ADRs in `.agents/brain/adr/`.
- Confirm non-goals remain intact: VLLE and PC-SAFT are out of scope.

2. Classify the change.
- Treat edits touching `src/chemthermo/__init__.py`, `src/chemthermo/eos/__init__.py`, or `src/chemthermo/vlle/__init__.py` as public API-sensitive.
- Treat changes to solver/model internals as invariant-sensitive (SI units, composition validation, deterministic flash behavior).
- Treat notebook edits as hygiene-sensitive (strip outputs, keep runnable narrative).

3. Implement a thin vertical slice.
- Prefer the smallest end-to-end usable increment over scaffold-only work.
- Keep interfaces and behavior coherent across code, tests, and examples.
- Avoid opportunistic refactors unrelated to the requested change.

4. Validate in escalating order.
- Run focused tests for touched modules first.
- Run quality gates (must pass) before finalizing:
  - `ruff format src tests`
  - `ruff format --check src tests`
  - `ruff check src tests`
  - `pyright`
  - `pytest -q`
- Run a golden path script for behavior sanity:
  - `python examples/flash_tp_peng_robinson_demo.py`
- If validation extras are installed, run optional reference checks under `tests/validation/`.

5. Update docs and decision records.
- Update README/examples/notebooks docs when user-visible behavior changes.
- Update `.agents/brain/steering-brief.md` when architecture context changes.
- Add or update an ADR for public API or architectural decisions.

6. Report outcomes clearly.
- Summarize what changed, what was validated, and any known risks.
- If any checks were skipped, state exactly which checks and why.

## Notebook-specific checklist
- Keep notebooks in `notebooks/` runnable top-to-bottom.
- Preserve concise markdown context for each code block.
- Ensure outputs are stripped before commit (nbstripout workflow).
- Keep supporting utilities in `notebooks/_nb_utils.py` or regular modules when reuse is needed.

## Guardrails
- Do not invent thermodynamic constants, component data, or external-reference results.
- Do not loosen CI or invariant checks unless explicitly requested.
- Keep commits small and logically grouped.
- Keep the working tree clean at handoff.
