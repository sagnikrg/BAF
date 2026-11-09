#!/bin/bash
# puts the one julia script of all jobs into $BUDDY (it takes L theta itr job eps as arguments,
# so no script per job is needed any more; the old per-job generators are in git history)

MBL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cp "$MBL_DIR/eigen_win.jl" "$BUDDY/eigen_win.jl"
echo "cat_jl: eigen_win.jl -> \$BUDDY"
