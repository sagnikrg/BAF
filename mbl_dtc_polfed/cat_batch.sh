#!/bin/bash
# Writes one HTCondor batch per (set, L, theta) of the job list, instead of one run file and one jdl per job:
#   ~/local/batches/<B>/run.sh      the one run script of the batch (takes the job's numbers as arguments)
#   ~/local/batches/<B>/job.jdl     one submit file, "queue ... from items.txt"
#   ~/local/batches/<B>/items.txt   one line per job: L theta r part eps(: separated) set tag
#   ~/local/batches/<B>/log/        err/out of the jobs and one cluster log
#   $BUDDY/batches/<B>/polfed_bands.jl   the julia script of this batch (a snapshot, so a later git pull never changes running jobs)
# <B> = mbldtc_<set>_L<L>_theta<theta>_<stamp> (set polfed or test); mopthefloor_BAF.sh removes a batch once none of its jobs is in condor_q.
# The list of the batches written is in ~/local/batches/.last (read by condor_submit.sh).

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
STAMP="${STAMP:-$(date +%Y%m%d-%H%M%S)}"
LB=~/local/batches; RB="$BUDDY/batches"
mkdir -p "$LB" "$RB"
: > "$LB/.last"
declare -A MEM

make_run() {   # $1 = batch
cat <<EOF
#!/bin/bash
# one job of batch $1, arguments: L theta r part eps(: separated) set tag
L=\$1; th=\$2; r=\$3; j=\$4; eps=\${5//:/,}; set=\$6; tag=\$7

# Goes to the job work directory (never go on anywhere else, the end of this script empties it)
cd /jwd || exit 1
source /etc/profile

# Copies, extracts and removes the Julia tarball
cp \$BUDDY/BAF/mbl_dtc_polfed/.julia/${TARBALL} ./
tar -xf ${TARBALL}
rm -f ${TARBALL}

lscpu --json | grep "Model name" | awk -F '"' '{print \$8}' > model_name.txt

# Loads Julia
module load ${JULIA_MODULE}

# the script of this batch and the header (itensor and polfed only)
cp \$BUDDY/batches/$1/polfed_bands.jl ./
cp -rf \$BUDDY/BAF/mbl_dtc_polfed/.header ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1
export JULIA_PKG_OFFLINE=true

# Do the real thing here
julia polfed_bands.jl \$L \$th \$r \$j \$eps 64 \$set > julia_\$tag.log 2>&1
ec=\$?; echo "julia exit code \$ec" >> julia_\$tag.log
rm polfed_bands.jl

# Clean up the working directory
mkdir -p /cephfs/user/sghosh/data/mbl_dtc_polfed/logs
cp julia_\$tag.log /cephfs/user/sghosh/data/mbl_dtc_polfed/logs/.
cp *hdf5 /cephfs/user/sghosh/data/mbl_dtc_polfed/.
rm -rf *
exit \$ec   # the job exit code is julia's (not that of the last cp/rm)
EOF
}

make_jdl() {   # $1 = batch, L theta mem hours
cat <<EOF

#Job Script to be submitted using HTCondor (one batch: all jobs of items.txt)

Executable              = run.sh
JobBatchName            = mbldtc_polfed: L=$2, theta=$3
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=(SUBHOST);

Arguments               = "\$(L) \$(th) \$(r) \$(j) \$(eps) \$(set) \$(tag)"

Universe                = vanilla

Transfer_executable     =
Transfer_input_files    =
Transfer_output_files   =


Error                   = log/err.\$(tag)
Output                  = log/out.\$(tag)
Log                     = log/cluster.log

Request_memory          = $4  GB
Request_cpus            = 1

Request_disk            = 2 GB

+CephFS_IO  = "low"

+MaxRuntimeHours	= $5

+ContainerOS            = "Debian12"

# Submit the jobs, one per line of items.txt
Queue L,th,r,j,eps,set,tag from items.txt

EOF
}

add_job() {
  local L=$1 th=$2 r=$3 j=$4 eps=$5 mem=$6 hours=$7 tag=$8 set=$9
  local B="mbldtc_${set}_L${L}_theta${th}_${STAMP}" d
  d="$LB/$B"
  if [[ -z "${MEM[$B]:-}" ]]; then
    MEM[$B]="$mem/$hours"
    mkdir -p "$d/log" "$RB/$B"
    cp "$POL_DIR/polfed_bands.jl" "$RB/$B/polfed_bands.jl"
    make_run "$B" > "$d/run.sh"; chmod +x "$d/run.sh"
    make_jdl "$B" "$L" "$th" "$mem" "$hours" > "$d/job.jdl"
    : > "$d/items.txt"
    echo "$B" >> "$LB/.last"
  elif [[ "${MEM[$B]}" != "$mem/$hours" ]]; then
    echo "ERROR: two lines for $set L=$L theta=$th with different memory/hours; give them one" >&2; exit 1
  fi
  echo "$L $th $r $j ${eps//,/:} $set $tag" >> "$d/items.txt"
}

each_job add_job
for B in $(cat "$LB/.last"); do echo "cat_batch: $B  $(wc -l < "$LB/$B/items.txt") jobs"; done
echo "cat_batch: $(count_jobs) jobs in $(wc -l < "$LB/.last") batches in ~/local/batches"
