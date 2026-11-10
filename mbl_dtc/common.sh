#!/bin/bash
# Shared by cat_run.sh, cat_env.sh, cat_jl.sh, condor_submit.sh: reads jobs.txt once, so the four never disagree.
# each_job <function>  calls  <function> L theta itr job eps mem tag   for every job of every line of jobs.txt
# optional 6th column 'from': first job number of the line (default 1), to resubmit the tail of a line

MBL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOBS="${JOBS:-$MBL_DIR/jobs.txt}"

each_job() {
  local fn=$1 L th itr nj mem from eps i
  while read -r L th itr nj mem from _; do
    [[ -z "${L:-}" || "$L" == \#* ]] && continue
    if [[ ! -f "$MBL_DIR/eps/theta_${th}.txt" ]]; then
      echo "ERROR: no epsilon list eps/theta_${th}.txt" >&2; exit 1
    fi
    eps=$(tr -d ' \r\n' < "$MBL_DIR/eps/theta_${th}.txt")
    [[ "${from:-}" =~ ^[0-9]+$ ]] || from=1
    for (( i=from; i<=nj; i++ )); do
      "$fn" "$L" "$th" "$itr" "$i" "$eps" "$mem" "mbldtc_win_L${L}_theta${th}_${i}"
    done
  done < "$JOBS"
}

count_jobs() { local n=0; _c() { n=$((n+1)); }; each_job _c; echo $n; }
