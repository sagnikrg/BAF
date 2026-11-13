#!/usr/bin/env bash
#########
# restore_batches.sh: rewrites the files of batches that still have jobs in condor_q but whose files are gone
# (2026-11-13: BAF.sh ran mopthefloor_BAF.sh, which removed the files of the queued V2 batches; the jobs went on hold).
#   bash ~/local/BAF/restore_batches.sh            lists what is missing and rewrites it
#   bash ~/local/BAF/restore_batches.sh --release  the same, then releases the held jobs of these batches
# Per batch directory (taken from the Cmd of the jobs, exactly as condor knows it): log/, run.sh, job.jdl, and
# $BUDDY/batches/<B>/<julia script>. Existing files are never overwritten. items.txt is not needed by queued jobs.
# Folder by JobBatchName: "mbldtc_polfed: ..." -> mbl_dtc_polfed (polfed_bands.jl), "mbldtc_<set>" -> mbl_dtc (eigen_win.jl).
#########
set -uo pipefail
release=0; [[ "${1:-}" == "--release" ]] && release=1
: "${BUDDY:?BUDDY is not set}"
here=$(cd "$(dirname "$0")" && pwd)
RB="$BUDDY/batches"

# one line per batch directory: dir memMB hours JobBatchName
if ! q=$(condor_q "$USER" -af Cmd RequestMemory MaxRuntimeHours JobBatchName 2>/dev/null); then echo "condor_q does not answer"; exit 1; fi
lines=$(grep '/batches/[^/]*/run.sh ' <<< "$q" | awk '{ d=$1; sub(/\/run\.sh$/,"",d); $1=d; print }' | sort -u -k1,1)
[[ -z $lines ]] && { echo "no batch jobs in condor_q"; exit 0; }

fixed=(); nmiss=0
while read -r dir mem hours name; do
  B=$(basename "$dir")
  if [[ $name == mbldtc_polfed:* ]]; then F=mbl_dtc_polfed; jl=polfed_bands.jl; else F=mbl_dtc; jl=eigen_win.jl; fi
  L=$(sed -E 's/.*_L([0-9]+)_theta.*/\1/' <<< "$B"); th=$(sed -E 's/.*_theta([0-9.]+)_.*/\1/' <<< "$B")
  set=$(sed -E 's/^mbldtc_(.*)_L[0-9]+_theta.*/\1/' <<< "$B")
  memgb=$(( (${mem%.*} + 1023) / 1024 )); hours=${hours%.*}
  miss=()
  [[ -d $dir/log ]]        || miss+=(log/)
  [[ -f $dir/run.sh ]]     || miss+=(run.sh)
  [[ -f $dir/job.jdl ]]    || miss+=(job.jdl)
  [[ -f $RB/$B/$jl ]]      || miss+=("\$BUDDY/batches/$B/$jl")
  (( ${#miss[@]} )) || { echo "ok       $B"; continue; }
  echo "missing  $B: ${miss[*]}"; nmiss=$((nmiss+1))
  # the generators of this folder's cat_batch.sh (make_run, make_jdl), with its common.sh settings
  ( source "$here/$F/common.sh"
    eval "$(grep -E '^(TARBALL|JULIA_MODULE)=' "$here/$F/cat_batch.sh")"
    eval "$(sed -n '/^make_run()/,/^}/p; /^make_jdl()/,/^}/p' "$here/$F/cat_batch.sh")"
    set -e
    mkdir -p "$dir/log" "$RB/$B"
    [[ -f $RB/$B/$jl ]]   || cp "$here/$F/$jl" "$RB/$B/$jl"
    [[ -f $dir/run.sh ]]  || { make_run "$B" > "$dir/run.sh.new"; chmod +x "$dir/run.sh.new"; mv "$dir/run.sh.new" "$dir/run.sh"; }
    if [[ ! -f $dir/job.jdl ]]; then
      if [[ $F == mbl_dtc_polfed ]]; then make_jdl "$B" "$L" "$th" "$memgb" "$hours" > "$dir/job.jdl.new"
      else make_jdl "$B" "$L" "$th" "$memgb" "$hours" "$set" > "$dir/job.jdl.new"; fi
      mv "$dir/job.jdl.new" "$dir/job.jdl"
    fi ) || { echo "  FAILED to rewrite $B"; continue; }
  fixed+=("$dir")
done <<< "$lines"

echo "restore: $nmiss of $(grep -c . <<< "$lines") batches rewritten"
if (( release )); then
  # the held jobs of every batch that has all its files now (only these; other held jobs are left alone)
  while read -r dir _; do
    [[ -f $dir/run.sh && -d $dir/log ]] && condor_release -constraint "Owner == \"$USER\" && JobStatus == 5 && Cmd == \"$dir/run.sh\"" 2>&1 | tail -1
  done <<< "$lines"
else
  echo "check: ls ${fixed[0]:-<batch dir>}; then: bash $0 --release   (releases only the held jobs of these batches)"
fi
