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
