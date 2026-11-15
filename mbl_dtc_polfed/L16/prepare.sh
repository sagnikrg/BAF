#!/bin/bash
# V2.1 L=16, step 1 (run once): pulls both BAF copies and writes the 6 batches (theta 0.0, 0.1 x parts p1-p3, 10000 jobs each)
# into ~/local/batches, nothing is submitted. The batch names go to ~/local/batches/.L16 (lines 1-3 theta 0.0 p1-p3, 4-6 theta 0.1 p1-p3).
#   bash ~/local/BAF/mbl_dtc_polfed/L16/prepare.sh
# then submit each part on its own:  bash ~/local/BAF/mbl_dtc_polfed/L16/submit_theta0.0_p1.sh   (etc.)
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
if [[ "${NOPULL:-0}" != 1 ]]; then ( cd "$here/.." && git pull ); ( cd "$BUDDY/BAF/" && git pull ); fi
LB="${LB:-$HOME/local/batches}"
JOBS="$here/jobs_v2.1_L16.txt" LB="$LB" bash "$here/cat_batch.sh"
cp "$LB/.last" "$LB/.L16"
echo "L16 batches (in $LB/.L16):"; cat -n "$LB/.L16"
