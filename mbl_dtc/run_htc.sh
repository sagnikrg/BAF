#!/bin/bash
# HTC dispatcher of eigen_win.jl jobs on the laptop (SCRIPT, default eigen_win.jl: histograms + the full diag/pi_diag lists, Tony 2026-11-19; eigen_laz.jl only on request) (Tony 2026-11-19, T1.2): like run_laz_local.sh, but a job starts only when
#   - fewer jobs run than the lane limit, read each round from ~/.jobq/lanes_<BATCH> (default 6; Vision sets 7 after Tony's "bye"), and
#   - the RAM in use after the start (in use now + PEAK_GB) stays <= the cap of ~/.jobq/ramcap_<BATCH> (percent, default 75; Tony 2026-11-19: 90).
# Waits first until the jobq entry WAITFOR (default T1_L8) is gone, so the batches run one after the other.
#   bash run_htc.sh LIST JOBQNAME [PEAK_GB=1.2] [WAITFOR=T1_L8]
LIST="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; JQ=$2; PEAK=${3:-1.2}; WAITFOR=${4:-T1_L8}; SCRIPT=${SCRIPT:-eigen_win.jl}
HERE="$(cd "$(dirname "$0")" && pwd)"; OUT=/e/mbl_dtc/rawdata/lazarides; RAW=/e/mbl_dtc/rawdata/eigen
LANES=/c/Users/Ritam/.jobq/lanes_$JQ; RAMCAP=/c/Users/Ritam/.jobq/ramcap_$JQ; JULIA=/c/Users/Ritam/.julia/juliaup/julia-1.9.4+0.x64.w64.mingw32/bin/julia.exe
JQPS='C:\Users\Ritam\.jobq\jobq.ps1'
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$OUT/logs"; P="$OUT/logs/progress.txt"
while [ -n "$WAITFOR" ] && [ -f "/c/Users/Ritam/.jobq/queue/$WAITFOR.json" ]; do sleep 60; done
echo "$(date '+%F %T') $JQ dispatcher starts" >> "$P"
powershell -NoProfile -ExecutionPolicy Bypass -File "$JQPS" request "$JQ" -Cores 6 -Ram 8 -Disk E -Minutes 2400 -Chat Vision > /dev/null
one() {
  read -r L th itr job eps <<< "$1"
  f="mbldtc_laz_L${L}_theta_${th}_${job}"; log="$OUT/logs/$f.log"
  mkdir -p "$RAW/L$L/theta_$th"; rm -f "$RAW/L$L/theta_$th/$f.hdf5"
  cd "$RAW/L$L/theta_$th" && "$JULIA" --startup-file=no "$HERE/$SCRIPT" $L $th $itr $job $eps laz >| "$log" 2>&1
  echo "$(date '+%F %T') $f exit $?" >> "$P"
}
ramok() {   # in use + PEAK <= cap% of total (MB from Windows)
  local c; c=$(cat "$RAMCAP" 2>/dev/null || echo 75)
  read -r tot av <<< "$(powershell -NoProfile -Command "[int]((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1MB), [int](Get-Counter '\Memory\Available MBytes').CounterSamples.CookedValue" | tr -d '\r' | tr '\n' ' ')"
  awk -v t="$tot" -v a="$av" -v p="$PEAK" -v c="$c" 'BEGIN{exit !((t-a)+p*1024 <= c/100*t)}'
}
declare -a pids=()
grep -v '^\s*#' "$LIST" | grep -v '^\s*$' > "$OUT/logs/$JQ.todo"
while read -r line; do
  read -r L th itr job eps <<< "$line"
  grep -q maxrss_MB "$OUT/logs/mbldtc_laz_L${L}_theta_${th}_${job}.log" 2>/dev/null && continue
  while :; do
    alive=(); for p in "${pids[@]}"; do kill -0 "$p" 2>/dev/null && alive+=("$p"); done; pids=("${alive[@]}")
    n=$(cat "$LANES" 2>/dev/null || echo 6)
    if [ "${#pids[@]}" -lt "$n" ] && ramok; then break; fi
    sleep 30
  done
  one "$line" & pids+=($!)
  sleep 20   # let the new job reach its memory before the next check
done < "$OUT/logs/$JQ.todo"
wait
echo "$(date '+%F %T') $JQ all done" >> "$P"
powershell -NoProfile -ExecutionPolicy Bypass -File "$JQPS" done "$JQ" > /dev/null
