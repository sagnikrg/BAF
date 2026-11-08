#!/bin/bash
here=$(cd "$(dirname "$0")" && pwd)
cd ~/local/
grep -v '^\s*#' "$here/jobs.txt" | grep -v '^\s*$' | while read -r L g r0 r1 block cpus mem hours thetas; do
for (( a=r0; a<=r1; a+=block )); do
b=$(( a+block-1 )); (( b>r1 )) && b=$r1
condor_submit job_file_vertical_L${L}_g${g}_r${a}-${b}.jdl
done
done
