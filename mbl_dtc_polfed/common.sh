#!/bin/bash
# Shared by cat_batch.sh, condor_submit.sh and check.sh: reads the job list once, so they never disagree.
#
# Job list line:   L  theta  realisations  parts  memory_GB  hours  [from]  [set]
#   realisations : realisations r = from..realisations (from defaults to 1)
#   parts        : jobs per realisation; the epsilon list eps/L<L>_theta_<theta>.txt is cut into this many contiguous parts
#   hours        : MaxRuntimeHours of the jobs
#   set          : polfed (default), test, or v21 (V2.1, its own eps lists eps/v21/L<L>_theta_<theta>.txt): names the files
#                  mbldtc_<set>_..., so a test or V2.1 never overwrites V2 data
#   part         : only this part of each realisation (default all; check.sh writes one line per missing job with it)
#
# each_job <function>  calls  <function> L theta r part eps mem hours tag set parts   for every job

POL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOBS="${JOBS:-$POL_DIR/jobs.txt}"

# the julia tarball of the jobs (built by Tony on BAF, julia 1.11.1 (no 1.11.7 module on BAF) with ITensors, ITensorMPS, Arpack, ToeplitzMatrices, HDF5)
TARBALL="${TARBALL:-julia-1.11.1-polfed.tar.gz}"
JULIA_MODULE="${JULIA_MODULE:-julia/1.11.1}"

each_job() {
  local fn=$1 L th nr parts mem hours from set only ef ne per j a b eps r
  while read -r L th nr parts mem hours from set only _; do
    [[ -z "${L:-}" || "$L" == \#* ]] && continue
    from=${from:-1}; [[ "$from" =~ ^[0-9]+$ ]] || from=1
    set=${set:-polfed}
    [[ $set == polfed || $set == test || $set == v21 ]] || { echo "ERROR: set must be polfed, test or v21, not $set" >&2; exit 1; }
    ef="$POL_DIR/eps/L${L}_theta_${th}.txt"; [[ $set == v21 ]] && ef="$POL_DIR/eps/v21/L${L}_theta_${th}.txt"
    [[ -f "$ef" ]] || { echo "ERROR: no epsilon list $ef" >&2; exit 1; }
    IFS=',' read -r -a E <<< "$(tr -d ' \r\n' < "$ef")"
    ne=${#E[@]}; per=$(( (ne + parts - 1) / parts ))
    for (( r=from; r<=nr; r++ )); do
      for (( j=1; j<=parts; j++ )); do
        [[ -n "${only:-}" && "$only" != "$j" ]] && continue
        a=$(( (j-1)*per )); (( a >= ne )) && break
        b=$(( a+per < ne ? a+per : ne ))
        eps=$(IFS=','; echo "${E[*]:a:b-a}")
        "$fn" "$L" "$th" "$r" "$j" "$eps" "$mem" "$hours" "mbldtc_${set}_L${L}_theta${th}_r${r}_${j}" "$set" "$parts"
      done
    done
  done < "$JOBS"
}

count_jobs() { local n=0; _c() { n=$((n+1)); }; each_job _c; echo $n; }
