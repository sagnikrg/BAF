#!/bin/bash
# Writes one HTCondor job_file per job of jobs.txt into ~/local
here=$(cd "$(dirname "$0")" && pwd)
cd ~/local/
grep -v '^\s*#' "${JOBS:-$here/jobs.txt}" | grep -v '^\s*$' | while read -r L g r0 r1 block cpus mem hours thetas; do
for (( a=r0; a<=r1; a+=block )); do
b=$(( a+block-1 )); (( b>r1 )) && b=$r1
name=vertical_L${L}_g${g}_r${a}-${b}
cat <<EOS >job_file_${name}.jdl

#Job Script to be submitted using HTCondor (same layout as mbl_dtc)

Executable              = run_file_${name}.sh
JobBatchName            = vertical: L=${L}, g=${g}, r=${a}-${b}
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=(SUBHOST);

Arguments               = 8

Universe                = vanilla

Transfer_executable     =
Transfer_input_files    =
Transfer_output_files   =


Error                   = log/err.${name}
Output                  = log/out.${name}
Log                     = log/log.${name}

Request_memory          = ${mem} GB
Request_cpus            = ${cpus}

Request_disk            = 2 GB

+CephFS_IO              = "low"

+MaxRuntimeHours        = ${hours}

+ContainerOS            = "Debian11"

# Submit job
Queue 1
EOS
done
done
