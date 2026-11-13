#!/bin/bash
# Shared by cat_batch.sh, condor_submit.sh and check.sh: reads jobs.txt once, so they never disagree.
# each_job <function>  calls  <function> L theta itr job eps mem tag   for every job of every line of jobs.txt
# optional 6th column 'from': first job number of the line (default 1), to resubmit the tail of a line
# optional 7th column 'set': win (default, eps/theta_<theta>.txt), grid (the old 0.025 grid, eps/grid.txt)
#   or scr (the scrambling reference point, eps/scr.txt, V1.3);
#   test / testgrid: the same lists for jobs_test.txt, with their own file names, so a test never overwrites data;
#   it names the files mbldtc_<set>_L<L>_theta_<theta>_<job>.hdf5, so the two sets never collide

MBL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOBS="${JOBS:-$MBL_DIR/jobs.txt}"

each_job() {
  local fn=$1 L th itr nj mem from set eps i ef
  while read -r L th itr nj mem from set _; do
    [[ -z "${L:-}" || "$L" == \#* ]] && continue
    set=${set:-win}
    case $set in
      win|test)          ef="$MBL_DIR/eps/theta_${th}.txt" ;;
      grid|testgrid)     ef="$MBL_DIR/eps/grid.txt" ;;
      scr|testscr)       ef="$MBL_DIR/eps/scr.txt" ;;
      *) echo "ERROR: set must be win, grid, scr, test, testgrid or testscr, not $set" >&2; exit 1 ;;
    esac
    if [[ ! -f "$ef" ]]; then
      echo "ERROR: no epsilon list $ef" >&2; exit 1
    fi
    eps=$(tr -d ' \r\n' < "$ef")
    [[ "${from:-}" =~ ^[0-9]+$ ]] || from=1
    for (( i=from; i<=nj; i++ )); do
      "$fn" "$L" "$th" "$itr" "$i" "$eps" "$mem" "mbldtc_${set}_L${L}_theta${th}_${i}" "$set"
    done
  done < "$JOBS"
}

count_jobs() { local n=0; _c() { n=$((n+1)); }; each_job _c; echo $n; }
