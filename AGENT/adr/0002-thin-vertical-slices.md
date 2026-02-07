# 2. Adopt Thin Vertical Slices for Feature Development

Date: 2026-02-06
Status: Accepted

## Context
The `chemthermo` repository previously contained several "scaffolded" or "partial" features (e.g., PC-SAFT, VLLE boundary) that existed in code but offered no end-to-end capability to the user. This created a misleading public API surface and complicated maintenance.

## Decision
We will adopt a "Thin Vertical Slice" philosophy for all future development. 

A **Thin Vertical Slice** is defined as a change that:
1.  **Delivers a complete, usable capability** to the end-user (e.g., "I can flash Methane with PC-SAFT", not "I added the PC-SAFT class").
2.  **Touches all necessary layers**: User Interface (API/CLI) -> Logic/Model -> Data -> Tests -> Documentation.
3.  **Is validated end-to-end**: Must include a runnable "Golden Path" example or test case demonstrating the user flow.

We will **avoid** horizontal layering (e.g., "Implement the parameter database for all 500 components" before the model works). Instead, we will implement the model for *one* component, verify it works, and *then* expand the database.

## Consequences
*   **Positive**: Every PR / Task delivers tangible value. "Work in Progress" is minimized. The API is always "honest" (if it's there, it works).
*   **Negative**: Initial implementations might look "hardcoded" or "limited" (e.g., PC-SAFT only working for Methane). This is acceptable and preferred over "generic but broken".
