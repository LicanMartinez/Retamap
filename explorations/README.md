# explorations/ — side explorations (not part of the canonical pipeline)

Self-contained exploratory branches of the pipeline that were **run, measured
against the independent validation, and then set aside**. Nothing here produces
the maps used in the paper; the canonical pipeline is stages `01`–`08` in
[`../`](../).

**Why this folder is under version control** (unlike
[`../experimental/`](../experimental/), which is gitignored): these explorations
were carried through to a validated result, so the scripts *and* their measured
outcome are worth keeping backed up and runnable in the GEE Code Editor in case
the question is revisited. `experimental/` is for scratch forks that were never
evaluated.

## Convention

Each exploration lives in its own subfolder with a README stating **what was
tried, what the validation said, and why it was not adopted** — so that a future
reader can tell in one minute whether it is worth re-opening, without re-running
anything.

Scripts here keep their original numeric prefixes (e.g. `05y`, `06y`) so their
relationship to the canonical stage is obvious, and write assets under their own
suffix so they can never collide with the canonical ones.

| Subfolder | Question it asked | Verdict |
|---|---|---|
| [`yearlyRF2/`](yearlyRF2/) | Does training one RF2 **per year** beat the single pooled RF2? | **Not adopted** (2026-08-09) — worse retama omission in both validated years (McNemar p = 0.016 in 2025), and it breaks the map's zero-commission property in the far background. |

## Related

- [`../experimental/`](../experimental/) — scratch forks, never evaluated, gitignored.
- [`../../legacy/`](../../legacy/) — retired material that is not expected to run again.
