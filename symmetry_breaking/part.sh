#!/bin/bash
# Fine run of Chapter 7 (J.A.R.V.I.S., Tony 2026-11-17), one part at a time (every part < 10000 jobs, < 1 GB output except p4 ~9 GB):
#   bash ~/local/BAF/symmetry_breaking/part.sh <1|2|3|4|rerun>   (rerun = jobs_fine_rerun.txt, the lost output of 2026-11-23)
#   p1 L=14 scan | p2 L=12 scan, fine g | p3 L=12 full + kb | p4 L=8, L=10
# Download and clear the output of a part from /cephfs/user/sghosh/data/symmetry_breaking before the one after next goes in.
# 1. stops if this part is already in condor_q (JobBatchName "... (part <n>)")
# 2. pulls both copies, writes the batches of jobs_fine_p<n>.txt with cat_batch.sh (PART=<n> in the batch names)
# 3. guard: dry-run of every batch, 7 arguments per job, one job per item; then condor_submit.sh (3 retries per batch)
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd); P="${1:?part 1..4}"
J="$here/jobs_fine_p$P.txt"; [[ -f "$J" ]] || J="$here/jobs_fine_$P.txt"; [[ -f "$J" ]] || { echo "STOP: no jobs_fine_p$P.txt or jobs_fine_$P.txt" >&2; exit 1; }

if condor_q -batch 2>/dev/null | grep -q "^.*sb .*(part $P)"; then echo "STOP: part $P is already in condor_q"; exit 0; fi

( cd "$here/.." && git pull ); ( cd "$BUDDY/BAF/" && git pull )
JOBS="$J" PART="$P" bash "$here/cat_batch.sh" || { echo "STOP: cat_batch failed" >&2; exit 1; }

tot=0
for B in $(cat ~/local/batches/.last); do
  cd ~/local/batches/$B || exit 1
  condor_submit -dry-run log/dry.txt job.jdl >/dev/null || { echo "STOP: dry-run failed for $B" >&2; exit 1; }
  bad=$(grep -i '^Arguments\|^Args' log/dry.txt | sed -E 's/^[^=]*= *"?//; s/"$//' | awk 'NF!=7' | head -1)
  n=$(grep -ci '^Arguments\|^Args' log/dry.txt); m=$(wc -l < items.txt)
  [[ -z "$bad" && "$n" == "$m" ]] || { echo "STOP: $B expands to $n jobs for $m items, or a job without 7 arguments: '$bad'" >&2; exit 1; }
  echo "guard ok: $B, $n jobs"; tot=$((tot+n))
done
echo "part $P: $tot jobs pass the guard, submitting"
bash "$here/condor_submit.sh"
