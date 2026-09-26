# ADR-0039: One owner-directed rewrite of commit attribution; tags follow their content

Status: accepted
Date: 2026-09-26

## Context
Commits made by automated coding sessions carried the automation's author and
committer identity, co-author trailers and session-link trailers. Before the
first PyPI release the owner decided that the published history should record
the owner as the author of this work. ADR-0031 item 5 says a published tag is
never moved, and the release packets said never to force-push.

## Decision
1. **One rewrite, metadata only.** On 2026-09-26 the owner directed a rewrite of
   `dev/sprint` and of the tags `v0.2.0b1` and `v0.3.0b1`: the automation
   identity becomes the owner's, the automation's trailers are removed, and
   the signatures made invalid by that are dropped. Trees, parent order, dates
   and every other part of each message are unchanged. `main` and `v0.1.0`
   contained nothing to rewrite and were not touched.
2. **A migrated tag names a commit with the same tree as its original
   target.** That is the only sense in which a tag may move, and only under a
   decision like this one; ADR-0031 item 5 otherwise stands, and the new
   `v0.4.0` tag is never moved.
3. **Historical records keep the ids they were produced on.** CI runs,
   benchmark records and validation tables are not edited; the map in
   `.agents/reports/history-rewrite-2026-09-26-map.tsv` connects them to the
   new ids. Active instructions (release packets, checkout targets) use the
   new ids or a tag.
4. Commits are authored by the owner's identity with no automation trailers
   from now on.
5. Published GitHub release assets and any PyPI files are not rebuilt or
   replaced because of this rewrite: their trees are unchanged.

## Alternatives considered
- `.mailmap` only (rejected by the owner: it changes display, not history).
- Leave the history and change only future commits (rejected by the owner).
- Re-cut the old versions as new tags (rejected: the content did not change,
  so a new version would be misleading).

## Consequences
- Every commit id after `v0.1.0` from the stability work onward changed (139
  commits); clones made before 2026-09-26 must re-fetch (report section
  "Other clones").
- GitHub's read-only pull-request refs and old commit pages still reach the old
  objects; the rewrite does not claim otherwise.
- Details and the equivalence proof: `.agents/reports/history-rewrite-2026-09-26.md`.

## Supersedes (optional)
Amends ADR-0031 item 5 for the two migrated tags only.

## Superseded by (optional)
None.
