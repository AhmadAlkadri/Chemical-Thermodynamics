# AGENTS

This file is the authoritative delivery-policy surface for this repository.
All other docs may reference this file, but should not redefine these rules.

## Delivery Contract (Hard Gates)

### 1) Clean handoff (universal)
At handoff, the repository must have no tracked changes.

- Required check: `git status --porcelain`
- Pass condition: command prints no output.

### 2) Commit cadence (universal)
When work is delivered in slices, each slice must have at least one commit
containing a Git trailer in the commit message body:

- `Slice: <slug>`

Slug format:

- Regex: `^[a-z0-9][a-z0-9_-]*$`
- Examples: `db`, `cli`, `validation`, `policy`, `docs`, `ci`

### 3) Evidence commands
Run and record these commands for compliance evidence:

- `git status --porcelain`
- `git log --oneline -n 20`
- `git log -n 50 --format=%B | rg "^Slice:"`

## Override Policy
No overrides. These hard gates apply to all contributors (humans and agents).
