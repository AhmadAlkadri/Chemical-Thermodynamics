# ADR-0002: Thin vertical slices and golden-path proof

Status: accepted
Date: 2026-02-10

## Context
The repository needs an execution model that prevents partial scaffolding from
looking like delivered capability. Contributors and agents must ship changes
that are immediately usable end-to-end and verifiable from a clean checkout.

## Decision
We adopt a thin-vertical-slice philosophy for feature development.

A thin vertical slice is a change that:
1) Delivers one complete, usable user capability now (for example, "run a
   documented binary TP-flash flow") rather than partial scaffolding (for
   example, "added model class stubs only").
2) Touches all required layers for that capability: interface (API/CLI), core
   logic/model, data/parameters, tests, and documentation.
3) Includes runnable end-to-end validation via at least one golden-path
   command or test proving the intended user flow.

We explicitly avoid horizontal-only layering where broad infrastructure lands
first and usable behavior lands later.

## Non-negotiable criteria (hard gate)
A change is incomplete and must not be handed off until all are present:
1) Slice declaration: "After this change, user can X by running Y."
2) User-visible capability is runnable end-to-end.
3) Golden path command/test and outcome are reported.
4) Focused tests for touched behavior and outcomes are reported.
5) CI-equivalent checks and installability smoke checks are reported.
6) User-facing docs are updated when behavior changes.

## Consequences
- Positive: each task delivers tangible value and keeps API behavior honest.
- Positive: review quality improves because evidence is explicit and runnable.
- Tradeoff: initial slices may be narrowly scoped, which is preferred over
  broad but non-usable scaffolding.

## Supersedes (optional)
None.

## Superseded by (optional)
None.
