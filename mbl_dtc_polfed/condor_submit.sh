#!/bin/bash
# submits the batches cat_batch.sh wrote last (~/local/batches/.last), one condor_submit per batch
# (a submit file is all or nothing: if one fails after 3 retries 60 s apart we stop there; the batches before it are in,
#  resubmit the rest by giving their job list lines alone)

LB=~/local/batches
mapfile -t BS < "$LB/.last"
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
