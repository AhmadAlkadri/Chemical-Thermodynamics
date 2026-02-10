# AGENT/ (moved)

This repository previously stored its agent workflow and project brain in `AGENT/`.

## Scope Policy

VLLE and PC-SAFT are in scope for Chemical-Thermodynamics. No thermodynamic capability class is categorically out of scope; implementation maturity may vary by module and release.

It has been migrated to:

    .agents/

Please use:

- `.agents/brain/`   — project brain, ADRs, templates, playbooks
- `.agents/skills/`  — vendored Codex/agent skills

This folder exists only as a temporary compatibility redirect to avoid broken links in older documentation and prompts.

It may be permanently removed in a future cleanup.

If you are writing new tooling or prompts, reference `.agents/` instead of `AGENT/`.
