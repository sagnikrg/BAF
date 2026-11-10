#!/bin/bash
here=$(cd "$(dirname "$0")" && pwd)
cd ~/local/
grep -v '^\s*#' "${JOBS:-$here/jobs.txt}" | grep -v '^\s*$' | while read -r kind L th g r0 r1 block cpus mem hours; do
for (( a=r0; a<=r1; a+=block )); do
b=$(( a+block-1 )); (( b>r1 )) && b=$r1
condor_submit job_file_sb_${kind}_L${L}_th${th}_g${g}_r${a}-${b}.jdl
done
done
