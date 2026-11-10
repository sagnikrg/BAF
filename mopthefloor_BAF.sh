#!/usr/bin/env bash
#########
# mopthefloor_BAF.sh: clean up the files the BAF scripts generated, once the queue is empty.
# Run by hand on the cluster:   bash mopthefloor_BAF.sh          (asks for a key phrase, lists, asks, then deletes)
#                               bash mopthefloor_BAF.sh --logs   (also empties ~/local/log)
# Deletes only generated files, never the scripts or anything tracked in git:
#   ~/local/run_file_*.sh, ~/local/job_file_*.jdl
#   $BUDDY/*.jl written by the cat_jl/create_file scripts (top level of $BUDDY only):
#     mbldtc_L*.jl, eigen_win.jl, anderson_metastability_L*.jl, Ilattice_L_*.jl, opo_fp_*.jl
# Nothing inside the git copies ~/local/BAF or $BUDDY/BAF is touched.
#########
set -uo pipefail
logs=0; [[ "${1:-}" == "--logs" ]] && logs=1
: "${BUDDY:?BUDDY is not set}"

shopt -s nullglob
files=( ~/local/run_file_*.sh ~/local/job_file_*.jdl
        "$BUDDY"/mbldtc_L*.jl "$BUDDY"/eigen_win.jl "$BUDDY"/anderson_metastability_L*.jl
        "$BUDDY"/Ilattice_L_*.jl "$BUDDY"/opo_fp_*.jl )
(( logs )) && files+=( ~/local/log/err.* ~/local/log/out.* ~/local/log/log.* )
ex=(); for f in "${files[@]}"; do [[ -f "$f" ]] && ex+=( "$f" ); done   # eigen_win.jl is no glob: keep only what exists
files=( ${ex[@]+"${ex[@]}"} )

# only the files of one run: their names carry its key phrase (e.g. sb, vertical, mbldtc); empty = all
read -r -p "Key phrase of the run to mop (e.g. sb; empty = everything): " key
if [[ -n "$key" ]]; then
  keep=(); for f in "${files[@]}"; do [[ "$(basename "$f")" == *"$key"* ]] && keep+=( "$f" ); done
  files=( ${keep[@]+"${keep[@]}"} )
fi

if (( ${#files[@]} == 0 )); then echo "Floor is already clean${key:+ for \"$key\"}."; exit 0; fi

if command -v condor_q >/dev/null && [[ -n "$(condor_q -nobatch -format '%d\n' ClusterId 2>/dev/null)" ]]; then
  echo "WARNING: you still have jobs in the queue (condor_q); running jobs need their run_file and \$BUDDY/*.jl."
fi
echo "Will delete ${#files[@]} files:"
printf '  %s\n' "${files[@]:0:5}"; (( ${#files[@]} > 5 )) && echo "  ... and $(( ${#files[@]} - 5 )) more"
for d in ~/local "$BUDDY" ~/local/log; do
  n=0; for f in "${files[@]}"; do [[ "$(dirname "$f")" == "$d" ]] && ((n++)); done
  (( n )) && echo "  $n in $d"
done
read -r -p "Delete them? [y/N] " ok
[[ "$ok" == "y" || "$ok" == "Y" ]] || { echo "Nothing deleted."; exit 0; }
rm -f -- "${files[@]}"
echo "Mopped: ${#files[@]} files deleted."
