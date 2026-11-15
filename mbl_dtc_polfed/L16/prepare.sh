#!/bin/bash
# V2.1 L=16, step 1: pulls both BAF copies and writes the 3 batches of ONE theta (parts p1-p3, 10000 jobs each)
# into ~/local/batches, nothing is submitted (one cat per theta, the cat of both took too long).
# The batch names go to ~/local/batches/.L16_theta<theta> (lines 1-3 = p1-p3).
#   bash ~/local/BAF/mbl_dtc_polfed/L16/prepare.sh 0.0      then submit_theta0.0_p1.sh .. p3.sh
#   bash ~/local/BAF/mbl_dtc_polfed/L16/prepare.sh 0.1      then submit_theta0.1_p1.sh .. p3.sh
set -euo pipefail
th="${1:?give theta: 0.0 or 0.1}"
here=$(cd "$(dirname "$0")/.." && pwd)
if [[ "${NOPULL:-0}" != 1 ]]; then ( cd "$here/.." && git pull ); ( cd "$BUDDY/BAF/" && git pull ); fi
LB="${LB:-$HOME/local/batches}"
JOBS="$here/jobs_v2.1_L16_theta$th.txt" LB="$LB" bash "$here/cat_batch.sh"
cp "$LB/.last" "$LB/.L16_theta$th"
echo "L16 theta=$th batches (in $LB/.L16_theta$th):"; cat -n "$LB/.L16_theta$th"
