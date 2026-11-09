#!/usr/bin/env bash
#########
# BAF.sh of vertical_phases (one per folder, run from anywhere):
#   bash vertical_phases/BAF.sh                 all jobs of jobs.txt
#   bash vertical_phases/BAF.sh jobs_test.txt   the small test (2 short jobs)
#   add --dry to write the run and job files into ~/local without submitting
# Pulls both BAF copies, writes ~/local/run_file_*.sh and job_file_*.jdl, submits.
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
mkdir -p ~/local/log

bash "$here/cat_run.sh"
bash "$here/cat_env.sh"
bash "$here/cat_jl.sh"
if [[ $dry == 0 ]]; then
  bash "$here/condor_submit.sh"
else
  echo "dry run: files in ~/local, nothing submitted"
fi
echo "Done."
