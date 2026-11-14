#!/bin/bash
# submits the batches cat_batch.sh wrote last (~/local/batches/.last), one condor_submit per batch
# (a submit file is all or nothing: if one fails after 3 retries 60 s apart we stop there; the batches before it are in,
#  resubmit the rest by giving their job list lines alone)

LB=~/local/batches
mapfile -t BS < "$LB/.last"
# guard (Tony 2026-11-14, V1.3 died because an empty queue variable shifted the julia arguments): expand every batch
# with condor_submit -dry-run first (into log/dry.txt, overwritten each time); every job must get exactly 7 arguments
# and the batch must expand to one job per line of items.txt, else nothing is submitted
for B in "${BS[@]}"; do
  ( cd "$LB/$B" && condor_submit -dry-run log/dry.txt job.jdl >/dev/null ) || { echo "STOP: dry-run failed for $B, nothing submitted" >&2; exit 1; }
  bad=$(grep -i '^Arguments\|^Args' "$LB/$B/log/dry.txt" | sed -E 's/^[^=]*= *"?//; s/"$//' | awk 'NF!=7' | head -1)
  n=$(grep -ci '^Arguments\|^Args' "$LB/$B/log/dry.txt"); m=$(wc -l < "$LB/$B/items.txt")
  if [[ -n "$bad" || "$n" != "$m" ]]; then echo "STOP: $B expands to $n jobs for $m items, or a job without 7 arguments: '$bad'; nothing submitted" >&2; exit 1; fi
done
echo "guard: all batches expand to 7 arguments per job"
echo "submitting ${#BS[@]} batches ..."
for B in "${BS[@]}"; do
  ok=0
  for t in 1 2 3 4; do
    if ( cd "$LB/$B" && condor_submit job.jdl ); then ok=1; break; fi
    [[ $t -lt 4 ]] && { echo "retry $t: $B" >&2; sleep 60; }
  done
  (( ok )) || { echo "FAILED, stopping at batch $B; the batches before it are submitted" >&2; exit 1; }
  sleep 2
done
echo "condor_submit: done"
