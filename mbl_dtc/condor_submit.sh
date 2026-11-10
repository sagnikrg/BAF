#!/bin/bash
# submits every job of jobs.txt (see common.sh)

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
cd ~/local/

# gentle on the schedd: a short pause after every job, 3 retries 60 s apart, stop at the first job that
# still fails (submission is in job order, so everything before it is in; rerun from there with 'from')
submit() {
  local t
  for t in 1 2 3 4; do
    condor_submit "job_file_$7.jdl" > /dev/null && { sleep 0.3; return 0; }
    [[ $t -lt 4 ]] && { echo "retry $t: $7" >&2; sleep 60; }
  done
  echo "FAILED, stopping: $7 (L=$1 theta=$2 job $4); resubmit this line with from=$4" >&2; exit 1
}

echo "submitting $(count_jobs) jobs ..."
each_job submit
echo "condor_submit: done"
