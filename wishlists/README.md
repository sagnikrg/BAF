# Simulation wishlists

One file per chat, named after its nickname (`Vision.md`, `EDITH.md`, `FRIDAY.md`, `JARVIS.md`, `Ultron.md`). Each chat keeps only its own file and commits it whenever an entry changes. Tony sets the priorities across the files; Vision merges them into one table when asked.

## Status of an entry

| Status | Meaning |
|---|---|
| planned | on the wishlist, not approved yet |
| approved | Tony said yes, not scripted or submitted yet |
| queued | scripted and submitted (cluster) or in the laptop job queue; give the job name / cluster id |
| final | done, data in place (give the path) |

Entries move down the list; keep an entry when it is final, so the file is the record of what was run.

## One entry

```markdown
### <short id>: <what, one line>
- **Status:** planned | approved | queued (<job / cluster id>) | final (<data path>), <date>
- **Where:** cluster (BAF) | laptop (job queue)
- **Parameters:** L, theta, g (or epsilon) grid, h, realisations per point
- **Cost:** core-hours in total, longest single job, number of jobs
- **For:** Fig. X.Y `fig:label` (figure number from the latest full draft, plus the \label, since numbers move), and what it settles there
- **Notes:** dependencies, quick or long, anything Tony should know for the priority
```

Short entries are fine; the figure reference and the cost are the parts Tony needs to set the priority.
