#!/usr/bin/env bash
#########
# BAF.sh of mbl_dtc (one per folder, run from anywhere):
#   bash mbl_dtc/BAF.sh                   all jobs of jobs.txt
#   bash mbl_dtc/BAF.sh jobs_v1.3.txt     V1.3
#   bash mbl_dtc/BAF.sh jobs_test.txt     the small test (sets test/testgrid/testscr, their own file names)
#   add --dry to write the batches into ~/local/batches without submitting (and let condor_submit -dry-run expand them)
# Pulls both BAF copies, mops the batches whose jobs are all gone from condor_q, writes one batch per (set, L, theta)
# with cat_batch.sh, submits.
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
bash "$here/../mopthefloor_BAF.sh" --yes || echo "(mop skipped)"

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
