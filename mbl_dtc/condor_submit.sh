#!/bin/bash
# submits every job of jobs.txt (see common.sh)

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
cd ~/local/

submit() { condor_submit "job_file_$7.jdl" > /dev/null || echo "FAILED: $7" >&2; }

echo "submitting $(count_jobs) jobs ..."
each_job submit
echo "condor_submit: done"
