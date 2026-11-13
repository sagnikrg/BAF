#!/bin/bash
# Which jobs of a job list are done, still in the queue, or lost (failed / never submitted)?
# Also the way to finish a submit that was cut off: everything not done and not in the queue is written out again.
#
#   bash mbl_dtc/check.sh [joblist] [--since "YYYY-MM-DD HH:MM"]     (default jobs.txt; run on the submit node)
#   bash mbl_dtc/BAF.sh jobs_missing.txt                             (submits exactly the missing ones)
#
#   done   : julia_<tag>.log in data/mbl_dtc/logs ends in "julia exit code 0" (and, with --since, is newer than that)
#   queue  : a job in condor_q (idle, running or held) with this tag: the last of its Arguments (batches of cat_batch.sh),
#            or its executable run_file_<tag>.sh (the old layout of one file per job)
#   failed : a log with another exit code, and not in the queue
#   lost   : no (new enough) log and not in the queue: never submitted, cut off, or died before the log was copied
#
# --since: logs older than this do not count as done, e.g. for jobs rerun under an old name (V1.2 reruns
#          job 1 of L=8, 10 theta 0.3, whose old test logs say exit code 0).
# jobs_missing.txt keeps the 'set' column (win / grid), one line per job. A failed job may have left a partial hdf5;
# the rerun overwrites it. If condor_q does not answer, the script stops, so queued jobs are never submitted twice.

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
list=jobs.txt; since=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --since) since=$2; shift 2 ;;
    *) list=$1; shift ;;
  esac
done
export JOBS="$here/$list"
source "$here/common.sh"

LOGS=/cephfs/user/sghosh/data/mbl_dtc/logs
out="$here/jobs_missing.txt"

# tags of the jobs in the queue (last argument of a batch job, or the old executable run_file_<tag>.sh), and the held ones
if ! queued=$(condor_q "$USER" -af Cmd JobStatus Arguments); then
  echo "ERROR: condor_q did not answer, nothing written (try again in a few minutes)" >&2; exit 1
fi
inq=$(echo "$queued" | awk 'NF{n=split($1,p,"/"); f=p[n]; if (f ~ /^run_file_/) { sub(/^run_file_/,"",f); sub(/\.sh$/,"",f); print f } else print $NF}')
held=$(echo "$queued" | awk '$2==5' | wc -l)

# logs newer than --since
declare -A NEW
if [[ -n "$since" ]]; then
  while read -r f; do NEW[$(basename "$f")]=1; done < <(find "$LOGS" -name 'julia_*.log' -newermt "$since")
fi

declare -A Q; for t in $inq; do Q[$t]=1; done
nd=0; nq=0; nf=0; nl=0
{
  echo "# failed or lost jobs of $list, from check.sh on $(date '+%F %T')${since:+ (logs since $since)}"
  echo "# L   theta  itr_per_job  njobs   memory_GB  from  set"
} > "$out"

check() {
  local L=$1 th=$2 itr=$3 i=$4 mem=$6 tag=$7 set=$8 log="$LOGS/julia_$7.log" fresh=1
  [[ -n "$since" && -z "${NEW[julia_$7.log]:-}" ]] && fresh=0
  if (( fresh )) && [[ -f "$log" ]] && tail -n 1 "$log" | grep -q "julia exit code 0"; then nd=$((nd+1)); return; fi
  if [[ -n "${Q[$tag]:-}" ]]; then nq=$((nq+1)); return; fi
  if (( fresh )) && [[ -f "$log" ]]; then nf=$((nf+1)); echo "failed: $tag ($(tail -n 1 "$log"))"
  else nl=$((nl+1)); fi
  printf "%3s   %-6s %-12s %-7s %-10s %-5s %s\n" "$L" "$th" "$itr" "$i" "$mem" "$i" "$set" >> "$out"
}
each_job check

echo "$list: done $nd, in queue $nq (held $held), failed $nf, lost / not submitted $nl"
if (( nf + nl > 0 )); then echo "submit the rest with: bash mbl_dtc/BAF.sh jobs_missing.txt"
else rm -f "$out"; echo "nothing to submit"; fi
