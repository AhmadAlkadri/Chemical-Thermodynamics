# History rewrite, 2026-09-26: attribution metadata only

Owner-directed (ADR-0039). Commit identities and attribution trailers written by
automated sessions were replaced or removed across `dev/sprint` and the tags
`v0.2.0b1` and `v0.3.0b1`. **No file content changed**: every rewritten commit
has the same tree as the original it replaces, the same (mapped) parents in the
same order, and the same author and committer dates and timezone offsets.

The full map, old SHA to new SHA with what changed in each, is
`.agents/reports/history-rewrite-2026-09-26-map.tsv` (139 rows).

## What changed

| change | commits |
| --- | --- |
| automation author and committer identity replaced by `Ahmad Alkadri <ahmadalkadri@berkeley.edu>` (dates kept); co-author and session trailers removed; the SSH signature on the object dropped | 22 |
| session-link trailer removed | 112 |
| nothing but the parent id (a descendant of a rewritten commit) | 5 |
| unchanged (all of `main`, `v0.1.0` and the first 9 commits after it) | 107 |

A signed object cannot keep its signature once its content changes, so the 22
signatures (all made by the automation's signing key, on the automation's own
commits) were removed, not re-created. The two GitHub-signed commits on `main`
were not touched and keep their signatures. No other message text changed: the
only difference in any message is the removed trailer lines (and the blank line
that separated them).

## Refs

| ref | before | after |
| --- | --- | --- |
| `refs/heads/main` | `a7a8ca742cadf41ad228ffa8c23d95d3252f5ce9` | unchanged |
| `refs/tags/v0.1.0` (tag object / commit) | `eb8c663d` / `a7a8ca74` | unchanged |
| `refs/tags/v0.2.0b1` (tag object) | `965acdf4967c8a2b1418dc625c7fed521ad5916e` | `0ad9d3c4d0e30b52bda45203c4b5b83466b9ba43` |
| `v0.2.0b1^{}` | `bd370681056baf487f0419a6dd05213b2d83932c` | `b5956d8ef8bc45576618a2932f949a4bafc16750` |
| `refs/tags/v0.3.0b1` (tag object) | `1edb10da6081dc56617409ada66245ede4ffe797` | `3a6df574d9913c225a7a25c0b11788a301affaef` |
| `v0.3.0b1^{}` | `64130a2d5373e88cc65c28fdd047807f04a5daf4` | `40630c312171820e8b9bdc9d4d378eb111c06eb6` |
| `refs/heads/dev/sprint` (before the bookkeeping commits) | `e95ae69a0af714e70d713ada802ee91b141e9156` | `ad28c142689bbbc03acb9bf31af34369053686c5` |

The new tag objects keep the original tagger, date and message; only the
`object` line differs. Each moved tag names a commit whose tree is identical to
its original target, so **each version's files are exactly what was released**,
and the GitHub release assets (built from the original commits) are unchanged.

## Commits named in the handoffs

Records written before this date name the original ids. They are historical
evidence and were not edited; read them through this table (or the map).

| original | rewritten | what it is |
| --- | --- | --- |
| `5041dd7` | `1b07b8f` | campaign checkpoint (first push) |
| `bd37068` | `b5956d8` | `v0.2.0b1` |
| `5864ab0` | `ff7cd6c` | cloud baseline |
| `64130a2` | `40630c3` | `v0.3.0b1` |
| `f852726` | `f030711` | macOS robustness map record |
| `761fd57` | `41cd3fd` | last solver change; Linux robustness map record |
| `bc92c5b` | `f9c5f3d` | packaging, library code as released |
| `4aa4753` | `f2165cd` | 0.4.0 release commit tested by the cloud session |
| `e95ae69` | `ad28c14` | 0.4.0 release packet |

CI runs, benchmark file names (`benchmarks/robustness_761fd57.*`) and
validation tables that name an original id ran on that original commit, whose
tree is identical to the rewritten one; they did not run on the new id.

## Equivalence proof (run before any push)

A script compared every one of the 246 commits reachable from the original
`main`, `dev/sprint` and tags with its image: tree equal; parent list equal to
the mapped original list; mapping one-to-one; `rev-list --topo-order` of each
branch equal element by element under the map; non-automation author and
committer lines byte-identical; automation lines differing only in name and
email; messages equal after removing exactly the two trailer patterns; no
other header changed except the dropped signatures; each moved tag's tree
equal to its original target's. Result: **0 differences**. After the rewrite,
no author, committer, tagger or commit message in the published history names
the automation.

## Limits

GitHub keeps the old head of pull request #1 in its read-only `refs/pull/1/*`
refs and in the PR timeline, and old commit pages stay reachable by id; CI
records name the old ids. Those were not, and cannot be, rewritten here. The
owner keeps a private rollback bundle of the original objects off the
repository.

## Other clones

A clone made before this date has the old `dev/sprint` and tags. Do not merge,
rebase or cherry-pick the old commits into the new history. To update:

```bash
git fetch origin --prune
git fetch origin --force 'refs/tags/v0.2.0b1:refs/tags/v0.2.0b1' 'refs/tags/v0.3.0b1:refs/tags/v0.3.0b1'
git checkout dev/sprint && git reset --keep origin/dev/sprint   # only with no local-only commits
```

Local-only commits made on the old history can be moved with
`git rebase --onto <new-sha> <old-sha> <branch>`, mapping `<old-sha>` through
the map above.
