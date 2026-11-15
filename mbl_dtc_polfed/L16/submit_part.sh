#!/bin/bash
# submits one L=16 batch, line $1 of ~/local/batches/.L16 (written by prepare.sh), with the guard of condor_submit.sh:
# dry-run expansion, 7 arguments per job and one job per line of items.txt, then condor_submit with 3 retries 60 s apart
set -uo pipefail
LB="${LB:-$HOME/local/batches}"
B=$(sed -n "${1}p" "$LB/.L16" 2>/dev/null); [[ -n "$B" && -d "$LB/$B" ]] || { echo "STOP: no batch on line $1 of $LB/.L16, run prepare.sh first" >&2; exit 1; }
cd "$LB/$B" || exit 1
condor_submit -dry-run log/dry.txt job.jdl >/dev/null || { echo "STOP: dry-run failed for $B" >&2; exit 1; }
bad=$(grep -i '^Arguments\|^Args' log/dry.txt | sed -E 's/^[^=]*= *"?//; s/"$//' | awk 'NF!=7' | head -1)
n=$(grep -ci '^Arguments\|^Args' log/dry.txt); m=$(wc -l < items.txt)
[[ -z "$bad" && "$n" == "$m" ]] || { echo "STOP: $B expands to $n jobs for $m items, or a job without 7 arguments: '$bad'" >&2; exit 1; }
echo "guard ok: $B, $n jobs"
for t in 1 2 3 4; do
  condor_submit job.jdl && { echo "submitted $B"; exit 0; }
  [[ $t -lt 4 ]] && { echo "retry $t: $B" >&2; sleep 60; }
done
echo "FAILED: $B" >&2; exit 1
