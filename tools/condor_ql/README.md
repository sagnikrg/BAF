# condor_ql (BAF/tools copy)

Marvin Lenk's `condor_ql` (https://github.com/marvinlenk/condor_ql, MIT, see LICENSE) with a batch view added.

- `condor_ql`: with no flags it is `condor_ql -b -u sghosh`. It goes to `~`, clears the screen, and shows one line per batch (`JobBatchName`) of sghosh across the whole cluster: submitted, done, running, idle, held, total, a progress bar, the longest running wall time, and the peak RAM against the request. `condor_userprio` follows below.
- `condor_ql -b -u <user>`: the same for another user.
- `condor_ql -j`, or any of the original flags (`-s 2`, `-hi`, `-u`, `-hn`): the original per job view.

DONE = the cluster's `TotalSubmitProcs` minus the jobs still in the queue. A batch whose jobs have all left the queue disappears from the view.

## Setup on BAF (once)

```bash
cd ~/local/BAF && git pull
chmod +x ~/local/BAF/tools/condor_ql/condor_ql
cat >> ~/.bashrc <<'EOF'
export PATH=~/local/BAF/tools/condor_ql:$PATH
condor_ql() { cd ~ && command condor_ql "$@"; }   # cd ~ in the shell itself, from wherever it is run
EOF
source ~/.bashrc
condor_ql
```

Remove the older line `export PATH=/home/sghosh/physik/local/condor_ql:$PATH` from `~/.bashrc`, so that this copy is the one that runs.

## batch_time_left

`batch_time_left [pattern]` (asks for the pattern when not given; part of the batch name, case insensitive) prints, per matching batch in the queue (all schedds), done / running / idle / held and the time until the last job ends:
- some jobs finished: the job time is the mean wall time of the finished jobs (condor_history of the batch's schedd and clusters);
- none finished yet: the job time is the requested MaxRuntimeHours ("raw", an upper bound).
Running jobs need the job time minus what they already ran; idle jobs start as the currently running slots free up. Held jobs are not counted. It is on the PATH with condor_ql (same folder); `chmod +x` it once after the pull.
