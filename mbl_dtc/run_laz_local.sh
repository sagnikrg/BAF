#!/bin/bash
# Local runner of eigen_laz.jl (Tony 2026-11-15): the jobs of jobs_laz_local.txt on NLANES single threaded lanes,
# raw data written directly to E:/mbl_dtc/rawdata/lazarides/L<L>/ (logs in .../logs). A job whose log ends in
# "maxrss_MB" is complete and skipped on a restart; an unfinished job (eigen_laz.jl opens its file with "cw")
# is removed and rerun from scratch.
#   bash run_laz_local.sh [NLANES]   (default 6); calls jobq done laz_L8_L10 at the end
NL=${1:-6}
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=/e/mbl_dtc/rawdata/lazarides
mkdir -p "$OUT/L8" "$OUT/L10" "$OUT/logs"
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1 HERE OUT
one() {
  read -r L th itr job eps <<< "$1"
  f="mbldtc_laz_L${L}_theta_${th}_${job}"; log="$OUT/logs/$f.log"
  grep -q maxrss_MB "$log" 2>/dev/null && return 0
  rm -f "$OUT/L$L/$f.hdf5"
  cd "$OUT/L$L" && julia +1.9.4 --startup-file=no "$HERE/eigen_laz.jl" $L $th $itr $job $eps laz > "$log" 2>&1
  echo "$(date '+%F %T') $f exit $?" >> "$OUT/logs/progress.txt"
}
export -f one
grep -v '^\s*#' "$HERE/jobs_laz_local.txt" | grep -v '^\s*$' | xargs -d '\n' -P "$NL" -I{} bash -c 'one "$@"' _ {}
echo "$(date '+%F %T') all done" >> "$OUT/logs/progress.txt"
powershell -NoProfile -ExecutionPolicy Bypass -File 'C:\Users\Ritam\.jobq\jobq.ps1' done laz_L8_L10
