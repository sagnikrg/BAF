#!/bin/bash
# V2.1 L=16, theta=0.0, part 1/3 alone (10000 jobs: r=201..533 at all 30 eps, r=534 at the first 10), the rest queued later.
#   bash ~/local/BAF/mbl_dtc_polfed/L16/p1_theta0.0.sh
# 1. stops if part 1 of L=16 theta=0.0 is already in condor_q
# 2. reuses a complete p1 batch which an earlier cat already wrote (~/local/batches/mbldtc_v21_L16_theta0.0_*_p1 with
#    run.sh, job.jdl, 10000 items and its polfed_bands.jl snapshot in $BUDDY/batches); otherwise cats part 1 only (ONLY_PART=1)
# 3. guard (dry-run, 7 arguments per job, one job per item) and condor_submit with 3 retries
set -uo pipefail
here=$(cd "$(dirname "$0")/.." && pwd); LB="${LB:-$HOME/local/batches}"; TH=0.0; P=1

if condor_q -batch 2>/dev/null | grep -q "L=16, theta=$TH (part $P/3)"; then echo "STOP: part $P of L=16 theta=$TH is already in condor_q"; exit 0; fi

B=""
for d in $(ls -dt "$LB"/mbldtc_v21_L16_theta${TH}_*_p$P 2>/dev/null); do
  b=$(basename "$d")
  if [[ -f "$d/run.sh" && -f "$d/job.jdl" && -f "$d/items.txt" && $(wc -l < "$d/items.txt") == 10000 && -f "$BUDDY/batches/$b/polfed_bands.jl" ]]; then B=$b; break; fi
done
if [[ -n "$B" ]]; then
  echo "reusing the complete batch $B"
else
  ( cd "$here/.." && git pull ); ( cd "$BUDDY/BAF/" && git pull )
  JOBS="$here/jobs_v2.1_L16_theta$TH.txt" ONLY_PART=$P LB="$LB" bash "$here/cat_batch.sh" || exit 1
  B=$(head -1 "$LB/.last")
  [[ "$B" == *_p$P && $(wc -l < "$LB/$B/items.txt") == 10000 ]] || { echo "STOP: the cat did not give a 10000 job part $P: $B" >&2; exit 1; }
fi

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
