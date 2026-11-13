#!/bin/bash
# Writes one HTCondor batch per (set, L, theta) of the job list, instead of one run file and one jdl per job:
#   ~/local/batches/<B>/run.sh      the one run script of the batch (takes the job's numbers as arguments)
#   ~/local/batches/<B>/job.jdl     one submit file, "queue ... from items.txt"
#   ~/local/batches/<B>/items.txt   one line per job: L theta itr job eps(: separated) set tag
#   ~/local/batches/<B>/log/        err/out of the jobs and one cluster log
#   $BUDDY/batches/<B>/eigen_win.jl   the julia script of this batch (a snapshot, so a later git pull never changes running jobs)
# <B> = mbldtc_<set>_L<L>_theta<theta>_<stamp>; mopthefloor_BAF.sh removes a batch once none of its jobs is in condor_q.
# The list of the batches written is in ~/local/batches/.last (read by condor_submit.sh).

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
STAMP="${STAMP:-$(date +%Y%m%d-%H%M%S)}"
LB=~/local/batches; RB="$BUDDY/batches"
TARBALL=julia-1.9.4-08-08-24.tar.gz
mkdir -p "$LB" "$RB"
: > "$LB/.last"
declare -A MEM

make_run() {   # $1 = batch
cat <<EOF
#!/bin/bash
# one job of batch $1, arguments: L theta itr job eps(: separated) set tag
L=\$1; th=\$2; itr=\$3; i=\$4; eps=\${5//:/,}; set=\$6; tag=\$7

# Goes to the job work directory (never go on anywhere else, the end of this script empties it)
cd /jwd || exit 1
source /etc/profile

# Copies, extracts and removes the Julia tarball
cp \$BUDDY/BAF/mbl_dtc/.julia/${TARBALL} ./
tar -xf ${TARBALL}
rm -f ${TARBALL}

lscpu --json | grep "Model name" | awk -F '"' '{print \$8}' > model_name.txt

# Loads Julia
module load julia/1.9.4

# the script of this batch and the header
cp \$BUDDY/batches/$1/eigen_win.jl ./
cp -rf \$BUDDY/BAF/mbl_dtc/.header ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Do the real thing here
julia eigen_win.jl \$L \$th \$itr \$i \$eps \$set > julia_\$tag.log 2>&1
ec=\$?; echo "julia exit code \$ec" >> julia_\$tag.log
rm eigen_win.jl

# Clean up the working directory
mkdir -p /cephfs/user/sghosh/data/mbl_dtc/logs
cp julia_\$tag.log /cephfs/user/sghosh/data/mbl_dtc/logs/.
cp *hdf5 /cephfs/user/sghosh/data/mbl_dtc/.
rm -rf *
exit \$ec   # the job exit code is julia's (not that of the last cp/rm)
EOF
}

make_jdl() {   # $1 = batch, L theta mem hours set
cat <<EOF

#Job Script to be submitted using HTCondor (one batch: all jobs of items.txt)

Executable              = run.sh
JobBatchName            = mbldtc_$6: L=$2, theta=$3
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=(SUBHOST);

Arguments               = "\$(L) \$(th) \$(itr) \$(i) \$(eps) \$(set) \$(tag)"

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
Queue L,th,itr,i,eps,set,tag from items.txt

EOF
}

add_job() {
  local L=$1 th=$2 itr=$3 i=$4 eps=$5 mem=$6 tag=$7 set=$8 hours=24
  local B="mbldtc_${set}_L${L}_theta${th}_${STAMP}" d
  d="$LB/$B"
  if [[ -z "${MEM[$B]:-}" ]]; then
    MEM[$B]="$mem/$itr"
    mkdir -p "$d/log" "$RB/$B"
    cp "$MBL_DIR/eigen_win.jl" "$RB/$B/eigen_win.jl"
    make_run "$B" > "$d/run.sh"; chmod +x "$d/run.sh"
    make_jdl "$B" "$L" "$th" "$mem" "$hours" "$set" > "$d/job.jdl"
    : > "$d/items.txt"
    echo "$B" >> "$LB/.last"
  elif [[ "${MEM[$B]}" != "$mem/$itr" ]]; then
    echo "ERROR: two lines for $set L=$L theta=$th with different memory/realisations per job; give them one" >&2; exit 1
  fi
  echo "$L $th $itr $i ${eps//,/:} $set $tag" >> "$d/items.txt"
}

each_job add_job
for B in $(cat "$LB/.last"); do echo "cat_batch: $B  $(wc -l < "$LB/$B/items.txt") jobs"; done
echo "cat_batch: $(count_jobs) jobs in $(wc -l < "$LB/.last") batches in ~/local/batches"
