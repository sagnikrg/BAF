#!/usr/bin/env bash
# missing.sh [E dir] : lines of jobs.txt whose output is not in <E dir>/L<L>/ (default /e/symmetry_breaking)
# output is in jobs.txt format (one job per block), ready to be used as a rerun list
E=${1:-/e/symmetry_breaking}
grep -v '^#' jobs.txt | awk 'NF' | while read -r kind L th g r0 r1 block cpus mem hrs; do
  for ((a=r0; a<=r1; a+=block)); do b=$((a+block-1)); ((b>r1)) && b=$r1
    [[ -f "$E/L$L/sb_L${L}_${kind}_theta${th}_g${g}_r${a}-${b}.hdf5" ]] || echo "$kind $L $th $g $a $b $block $cpus $mem $hrs"
  done
done
