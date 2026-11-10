#!/bin/bash
# creates one HTCondor job_file_<tag>.jdl per job in ~/local (jobs from the job list, see common.sh)

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
cd ~/local/
mkdir -p log

make_jdl() {
local L=$1 th=$2 r=$3 j=$4 eps=$5 mem=$6 hours=$7 tag=$8
cat <<EOF >job_file_${tag}.jdl

#Job Script to be submitted using HTCondor

Executable              = run_file_${tag}.sh
JobBatchName            = mbldtc_polfed: L=${L}, theta=${th}
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=(SUBHOST);

Arguments		= 8

Universe                = vanilla

Transfer_executable     =
Transfer_input_files    =
Transfer_output_files   =


Error                   = log/err.${tag}
Output                  = log/out.${tag}
Log                     = log/log.${tag}

Request_memory          = ${mem}  GB
Request_cpus            = 1

Request_disk            = 2 GB

+CephFS_IO  = "low"

+MaxRuntimeHours	= ${hours}

+ContainerOS            = "Debian12"

# Submit job
Queue 1

EOF
}

each_job make_jdl
echo "cat_env: $(count_jobs) jdl files in ~/local"
