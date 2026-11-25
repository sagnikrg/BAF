# condor_ql (BAF/tools copy)

Marvin Lenk's `condor_ql` (https://github.com/marvinlenk/condor_ql, MIT, see LICENSE) with a batch view and `batch_time_left` added.

## Setup on BAF (copy and paste, once)

```bash
curl -fsSL https://raw.githubusercontent.com/sagnikrg/BAF/main/tools/condor_ql/setup.sh | bash && source ~/.bashrc
```

The script:
- puts `condor_ql` and `batch_time_left` in `~/.local/condor_ql`;
- adds that folder to the PATH in `~/.bashrc`;
- makes sure `~/.bash_profile` (when it exists) loads `~/.bashrc`, so the commands still work in a new ssh session.

Run the same line again to update.

## condor_ql

- `condor_ql`: with no flags it is `condor_ql -b -u $USER`. It goes to `~`, clears the screen, and shows one line per batch (`JobBatchName`) of your jobs across the whole cluster (all schedds): schedd, submitted, done, running, idle, held, total, a progress bar, the longest running wall time, and the peak RAM against the request. `condor_userprio` follows below.
- `condor_ql -b -u <user>`: the same for another user.
- `condor_ql -j`, or any of the original flags (`-s 2`, `-hi`, `-u`, `-hn`): the original per job view.

DONE = the cluster's `TotalSubmitProcs` minus the jobs still in the queue. A batch whose jobs have all left the queue disappears from the view.

## batch_time_left

`batch_time_left [pattern]` (asks for the pattern when not given; part of the batch name, case insensitive) prints, per matching batch in the queue (all schedds), done / running / idle / held and the time until the last job ends:
- some jobs finished: the job time is the mean wall time of the finished jobs (condor_history of the batch's schedd and clusters);
- none finished yet: the job time is the requested MaxRuntimeHours ("raw", an upper bound).

Running jobs need the job time minus what they already ran; idle jobs start as the currently running slots free up. Held jobs are not counted.
