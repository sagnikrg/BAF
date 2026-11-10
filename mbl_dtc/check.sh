#!/bin/bash
# Which jobs of a job list are done, still in the queue, or lost (failed / never submitted)?
#
#   bash mbl_dtc/check.sh [joblist]          (default jobs.txt; run on the submit node)
#
#   done   : julia_<tag>.log in data/mbl_dtc/logs ends in "julia exit code 0"
#   queue  : run_file_<tag>.sh is the executable of a job in condor_q (idle, running or held)
#   failed : a log with another exit code, and not in the queue
#   lost   : no log and not in the queue (never submitted, or died before the log was copied)
#
# Writes jobs_missing.txt (one line per failed or lost job, same format as jobs.txt with 'from'), so
#   bash mbl_dtc/BAF.sh jobs_missing.txt
# resubmits exactly those. A failed job may have left a partial hdf5 in data/mbl_dtc; the rerun writes its
# file in a fresh /jwd and its cp overwrites the partial one, so there is nothing to clean.
# Jobs not submitted yet count as lost: check a list only after its submit has finished.

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export JOBS="$here/${1:-jobs.txt}"
source "$here/common.sh"

LOGS=/cephfs/user/sghosh/data/mbl_dtc/logs
out="$here/jobs_missing.txt"

# tags of the jobs in the queue (from the executable run_file_<tag>.sh), and the held ones
queued=$(condor_q "$USER" -af Cmd JobStatus 2>/dev/null)
inq=$(echo "$queued" | awk '{n=split($1,p,"/"); f=p[n]; sub(/^run_file_/,"",f); sub(/\.sh$/,"",f); print f}')
held=$(echo "$queued" | awk '$2==5' | wc -l)

declare -A Q; for t in $inq; do Q[$t]=1; done
nd=0; nq=0; nf=0; nl=0
{
  echo "# failed or lost jobs of $(basename "$JOBS"), from check.sh on $(date '+%F %T')"
  echo "# L   theta  itr_per_job  njobs   memory_GB  from"
} > "$out"

check() {
  local L=$1 th=$2 itr=$3 i=$4 mem=$6 tag=$7 log="$LOGS/julia_$7.log"
  if [[ -f "$log" ]] && tail -n 1 "$log" | grep -q "julia exit code 0"; then nd=$((nd+1)); return; fi
  if [[ -n "${Q[$tag]:-}" ]]; then nq=$((nq+1)); return; fi
  if [[ -f "$log" ]]; then nf=$((nf+1)); echo "failed: $tag ($(tail -n 1 "$log"))"
  else nl=$((nl+1)); echo "lost:   $tag"; fi
  printf "%3s   %-6s %-12s %-7s %-10s %s\n" "$L" "$th" "$itr" "$i" "$mem" "$i" >> "$out"
}
each_job check

echo "$(basename "$JOBS"): done $nd, in queue $nq (held $held), failed $nf, lost $nl"
if (( nf + nl > 0 )); then echo "resubmit with: bash mbl_dtc/BAF.sh jobs_missing.txt"
else rm -f "$out"; echo "nothing to resubmit"; fi
