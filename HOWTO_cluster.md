# How to use BAF well (for every chat)

Summary of the public BAF slides (PI IT team, Bonn, 29 Nov 2023; user tutorial https://unibonn.github.io/htcondor-bonn/),
and what follows from them for our jobs. Tony, 2026-11-13: use this for every cluster job.

## How BAF hands out cores

- **Fair share.** Every user has an *effective priority* that follows the recent usage (cores x priority factor 1000)
  exponentially, with a **half-life of 24 h**. Higher usage = worse priority; it recovers within a few days.
- **Only queued jobs get a share.** The free slots are split, proportionally to priority ("pie slices"), among the
  accounts that **have jobs waiting**. An empty queue gets nothing. So keep a continuous queue: submit the next batch
  before the current one drains (e.g. V1.3 right behind V2), rather than waiting for one batch to finish.
- **Tetris.** Jobs with different resource requests leave the cluster fragmented in many dimensions (cores, memory, disk,
  time). Small jobs fit into the gaps that big ones cannot. Hence our design: **one core per job, many short-ish jobs,
  memory requested close to the real peak.**
- `/etc/profile` in the job sets the `*_NUM_THREADS` from `Request_cpus`; we still export them to 1 in run.sh.
- About 7000 logical cores (2023). Nodes in the HRZ room are `CephFS_IO` "medium" (1 Gbit/s); our jobs say
  `+CephFS_IO = "low"`, so they can go everywhere.
- A node with iowait > 15 % is taken out as unhealthy: keep the I/O of a job small (our 133 MB tarball copy per job is
  the biggest I/O we do; fine so far).

## Rules for our jobs

1. **Batch submit, always.** One HTCondor submit file per batch with `Queue <vars> from items.txt`, never one jdl and
   one run file per job (thousands of files, endless `condor_q`, slow `mopthefloor`). Template:
   `BAF/mbl_dtc_polfed/cat_batch.sh` (or `BAF/mbl_dtc/cat_batch.sh`): one batch per (set, L, theta) in
   `~/local/batches/<B>/{run.sh, job.jdl, items.txt, log/}`, the julia script snapshot in `$BUDDY/batches/<B>/`,
   `JobBatchName` naming the batch, one `condor_submit` per batch. Values containing commas (epsilon lists) are passed
   with `:` and converted back in run.sh. `run.sh` starts with `cd /jwd || exit 1` (it ends with `rm -rf *`).
   `mopthefloor_BAF.sh` then removes every batch with no job left in `condor_q`, by itself, at the next submission.
2. **Memory = measured peak x ~1.3**, not a round guess. After the first jobs:
   `condor_history $USER -limit 20 -af JobBatchName MemoryUsage RequestMemory` (MemoryUsage in MB), and lower the
   request of the remaining batches. Too low gets the job held; too high strands cores.
3. **Job length: ~1-10 h.** Split long runs into parts (as V2 splits L=16 into 4 epsilon parts); `+MaxRuntimeHours`
   with margin.
4. **Keep the queue fed.** Submit the next batch before the previous one drains.
5. **Test in an interactive job**, not in a batch test that waits in the queue: same container and environment as the
   jobs, dedicated nodes, minutes instead of hours:
   ```
   condor_submit -interactive -append '+ContainerOS="Debian12"' -append 'request_memory=4GB'
   # then in the job shell: unpack the tarball as run.sh does, run one short job by hand, check the output
   ```
6. **Dry run first:** `bash <folder>/BAF.sh --dry` writes the batches and lets `condor_submit -dry-run` expand them.

## Useful commands

| what | command |
|---|---|
| my priority vs others | `condor_userprio -allusers` |
| my jobs, compact | `condor_q -batch` (one line per batch) |
| any attribute | `condor_q -af JobBatchName JobStatus RequestMemory` |
| finished jobs | `condor_history $USER -limit 50 -af JobBatchName ExitCode MemoryUsage RemoteWallClockTime` |
| why idle | `condor_q -better-analyze <id>` |
| follow a job | `condor_tail -f <id>` |
