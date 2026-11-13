#!/usr/bin/env bash
#########
# BAF.sh of symmetry_breaking (Chapter 7, J.A.R.V.I.S.; one per folder, run from anywhere):
#   bash symmetry_breaking/BAF.sh                 all jobs of jobs.txt
#   bash symmetry_breaking/BAF.sh jobs_test.txt   the ~2 minute test (3 short jobs)
#   add --dry to write the batches into ~/local/batches without submitting (and let condor_submit -dry-run expand them)
# Pulls both BAF copies, mops the batches whose jobs are all gone from condor_q, writes one batch per (kind, L, theta)
# with cat_batch.sh (BAF/HOWTO_cluster.md), submits.
#########
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
jobs=jobs.txt; dry=0
for a in "$@"; do [[ "$a" == "--dry" ]] && dry=1 || jobs="$a"; done
[[ -f "$here/$jobs" ]] || { echo "ERROR: no $here/$jobs" >&2; exit 1; }
export JOBS="$here/$jobs"

cd "$here"/../              # the BAF root of this copy
if [[ $dry == 0 ]]; then
  git pull
  ( cd "$BUDDY/BAF/" && git pull )
fi
# mop OFF (Tony 2026-11-13: it removed the files of queued V2 batches); run mopthefloor_BAF.sh by hand only
# bash "$here/../mopthefloor_BAF.sh" --yes || echo "(mop skipped)"

bash "$here/cat_batch.sh"
if [[ $dry == 0 ]]; then
  bash "$here/condor_submit.sh"
else
  if command -v condor_submit >/dev/null; then
    for B in $(cat ~/local/batches/.last); do
      ( cd ~/local/batches/$B && condor_submit -dry-run dry.txt job.jdl >/dev/null && echo "dry-run $B: $(grep -c '^Args\|^Arguments' dry.txt) jobs expanded, first: $(grep -m1 '^Arguments\|^Args' dry.txt)" )
    done
  fi
  echo "dry run: batches in ~/local/batches, nothing submitted"
fi
echo "Done."
