#!/bin/bash
# Local runner of eigen_laz.jl (Tony 2026-11-15): the jobs of jobs_laz_local.txt on NLANES single threaded lanes,
# raw data written directly to E:/mbl_dtc/rawdata/lazarides/L<L>/ (logs in .../logs). A job whose log ends in
# "maxrss_MB" is complete and skipped on a restart; an unfinished job (eigen_laz.jl opens its file with "cw")
# is removed and rerun from scratch.
#   bash run_laz_local.sh [NLANES] [LIST]   (default 6, jobs_laz_local.txt; lines "#R ..." are the later runs, skipped)
# A line "WAIT <file stem>" holds its lane until that job (started elsewhere) has finished. Calls jobq done laz_L8_L10 at the end.
NL=${1:-6}
LIST="${2:-$(cd "$(dirname "$0")" && pwd)/jobs_laz_local.txt}"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=/e/mbl_dtc/rawdata/lazarides
RAW=/e/mbl_dtc/rawdata/eigen   # output in $RAW/L<L>/theta_<th> (2026-11-19 layout), logs stay in $OUT/logs
mkdir -p "$OUT/logs"
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 HERE OUT RAW
one() {
  if [[ $1 == WAIT* ]]; then until grep -q maxrss_MB "$OUT/logs/${1#WAIT }.log" 2>/dev/null; do sleep 30; done; return 0; fi
  read -r L th itr job eps <<< "$1"
  f="mbldtc_laz_L${L}_theta_${th}_${job}"; log="$OUT/logs/$f.log"
  grep -q maxrss_MB "$log" 2>/dev/null && return 0
  mkdir -p "$RAW/L$L/theta_$th"; rm -f "$RAW/L$L/theta_$th/$f.hdf5"
  cd "$RAW/L$L/theta_$th" && julia +1.9.4 --startup-file=no "$HERE/eigen_laz.jl" $L $th $itr $job $eps laz >| "$log" 2>&1
  echo "$(date '+%F %T') $f exit $?" >> "$OUT/logs/progress.txt"
}
export -f one
grep -v '^\s*#' "$LIST" | grep -v '^\s*$' | xargs -d '\n' -P "$NL" -I{} bash -c 'one "$@"' _ {}
echo "$(date '+%F %T') all done" >> "$OUT/logs/progress.txt"
powershell -NoProfile -ExecutionPolicy Bypass -File 'C:\Users\Ritam\.jobq\jobq.ps1' done laz_L8_L10
