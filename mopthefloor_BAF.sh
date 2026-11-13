#!/usr/bin/env bash
#########
# mopthefloor_BAF.sh: removes what the BAF scripts generated for jobs that are no longer in condor_q.
#   bash mopthefloor_BAF.sh          shows what it would remove, asks once
#   bash mopthefloor_BAF.sh --yes    removes without asking (BAF.sh calls it this way before every submission)
#   add --keep-logs to keep the log/ directories of finished batches (they are removed with the batch by default)
#
# New layout (cat_batch.sh): one directory per batch, ~/local/batches/<B>/ and $BUDDY/batches/<B>/.
#   A batch is in use while any job of condor_q has its executable in ~/local/batches/<B>/; all others are removed whole.
# Old layout (one file per job): ~/local/run_file_*.sh, ~/local/job_file_*.jdl, ~/local/log/{err,out,log}.*,
#   and the julia scripts in the top level of $BUDDY; every file except those of the old style jobs still in condor_q
#   (the shared julia scripts only once no old style job is left).
# If condor_q does not answer, nothing is removed. Nothing tracked in git (~/local/BAF, $BUDDY/BAF) is ever touched.
#########
set -uo pipefail
yes=0; keeplogs=0
for a in "$@"; do case $a in --yes) yes=1 ;; --keep-logs) keeplogs=1 ;; esac; done
: "${BUDDY:?BUDDY is not set}"
LB=~/local/batches; RB="$BUDDY/batches"

# executables of all of my jobs still in the queue (idle, running, held)
if ! cmds=$(condor_q -af Cmd 2>/dev/null); then echo "mop: condor_q does not answer, nothing removed"; exit 1; fi

# new layout: batches not in use
gone=()
if [[ -d $LB ]]; then
  for d in "$LB"/*/; do
    [[ -d $d ]] || continue
    B=$(basename "$d")
    grep -qF "$LB/$B/" <<< "$cmds" || gone+=( "$B" )
  done
fi

# old layout: only if no old style job is left
# old layout, file by file: everything except the files of old style jobs still in the queue
oldlist=$(mktemp); keep=$(mktemp); trap 'rm -f "$oldlist" "$keep"' EXIT
grep -o 'run_file_[^/]*\.sh$' <<< "$cmds" | sed 's/^run_file_//; s/\.sh$//' | while read -r t; do
  printf '%s\n' "run_file_$t.sh" "job_file_$t.jdl" "err.$t" "out.$t" "log.$t"
done > "$keep"
{ find ~/local -maxdepth 1 \( -name 'run_file_*.sh' -o -name 'job_file_*.jdl' \) 2>/dev/null
  (( keeplogs )) || find ~/local/log -maxdepth 1 \( -name 'err.*' -o -name 'out.*' -o -name 'log.*' \) 2>/dev/null
  # the julia scripts of the old layout are shared by all jobs of a run: only once no old style job is left
  [[ -s $keep ]] || find "$BUDDY" -maxdepth 1 \( -name 'mbldtc_L*.jl' -o -name 'eigen_win.jl' -o -name 'polfed_bands.jl' -o -name 'anderson_metastability_L*.jl' -o -name 'Ilattice_L_*.jl' -o -name 'opo_fp_*.jl' \) 2>/dev/null
} | awk -F/ 'NR==FNR { k[$0]=1; next } !($NF in k)' "$keep" - > "$oldlist"
old=$(grep -c . "$oldlist")
(( $(grep -c . "$keep") )) && echo "mop: keeping the files of $(( $(grep -c . "$keep") / 5 )) old style jobs still in the queue"

if (( ${#gone[@]} == 0 && old == 0 )); then echo "mop: floor is clean ($(grep -c . <<< "$cmds") jobs in the queue)"; exit 0; fi
echo "mop: ${#gone[@]} finished batches, $old old style files"
for B in "${gone[@]:0:10}"; do echo "  $B"; done; (( ${#gone[@]} > 10 )) && echo "  ... and $(( ${#gone[@]} - 10 )) more"
if (( ! yes )); then
  read -r -p "Remove them? [y/N] " ok
  [[ "$ok" == "y" || "$ok" == "Y" ]] || { echo "Nothing removed."; exit 0; }
fi

for B in "${gone[@]}"; do
  if (( keeplogs )); then find "$LB/$B" -mindepth 1 -maxdepth 1 ! -name log -exec rm -rf {} +; else rm -rf "${LB:?}/$B"; fi
  rm -rf "${RB:?}/$B"
done
(( old )) && tr '\n' '\0' < "$oldlist" | xargs -0 rm -f --
echo "mop: removed ${#gone[@]} batches and $old old style files"
