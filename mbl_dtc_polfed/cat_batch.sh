#!/bin/bash
# Writes one HTCondor batch per (set, L, theta) of the job list, instead of one run file and one jdl per job:
#   ~/local/batches/<B>/run.sh      the one run script of the batch (takes the job's numbers as arguments)
#   ~/local/batches/<B>/job.jdl     one submit file, "queue ... from items.txt"
#   ~/local/batches/<B>/items.txt   one line per job: L theta r part eps(: separated) set tag
#   ~/local/batches/<B>/log/        err/out of the jobs and one cluster log
#   $BUDDY/batches/<B>/polfed_bands.jl   the julia script of this batch (a snapshot, so a later git pull never changes running jobs)
# <B> = mbldtc_<set>_L<L>_theta<theta>_<stamp> (set polfed, test or v21); mopthefloor_BAF.sh removes a batch once none of its jobs is in condor_q.
# Split (2026-11-15): a (set, L, theta) with more than MAXITEMS jobs (default 10000) is written as parts
#   <B>_p1, <B>_p2, ... of at most MAXITEMS jobs each, one condor_submit each. A submit file is all or nothing, and the
#   schedd refuses one above MAX_JOBS_PER_SUBMISSION (HTCondor default 20000): the V2.1 L=16 batches (30000 jobs) never got in.
#   The parts keep the prefix mbldtc_<set>_L<L>_theta<theta>_ that rm_jobs.sh, check.sh and restore_batches.sh match.
# The list of the batches written is in ~/local/batches/.last (read by condor_submit.sh).

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
STAMP="${STAMP:-$(date +%Y%m%d-%H%M%S)}"
MAXITEMS="${MAXITEMS:-10000}"
LB="${LB:-$HOME/local/batches}"; RB="$BUDDY/batches"
mkdir -p "$LB" "$RB"
: > "$LB/.last"
declare -A MEM CNT NDONE

# first pass: the number of jobs of every (set, L, theta), to know whether it has to be split
count_base() { local B="mbldtc_${9}_L${1}_theta${2}_${STAMP}"; CNT[$B]=$(( ${CNT[$B]:-0} + 1 )); }
each_job count_base

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
# no log copy to the cluster drive any more (Tony 2026-11-23: its file-number limit); the tail goes to the job's out file in ~/local/batches/<B>/log
tail -3 julia_\$tag.log
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
  local base="mbldtc_${set}_L${L}_theta${th}_${STAMP}" B d p np
  B=$base
  if (( ${CNT[$base]} > MAXITEMS )); then
    p=$(( ${NDONE[$base]:-0} / MAXITEMS + 1 )); np=$(( (${CNT[$base]} + MAXITEMS - 1) / MAXITEMS )); B="${base}_p${p}"
  fi
  NDONE[$base]=$(( ${NDONE[$base]:-0} + 1 ))
  [[ -n "${ONLY_PART:-}" && "${p:-1}" != "$ONLY_PART" ]] && return 0   # ONLY_PART=n: write only part n (L16/p1_theta0.0.sh)
  d="$LB/$B"
  if [[ -z "${MEM[$B]:-}" ]]; then
    MEM[$B]="$mem/$hours"
    mkdir -p "$d/log" "$RB/$B"
    cp "$POL_DIR/polfed_bands.jl" "$RB/$B/polfed_bands.jl"
    make_run "$B" > "$d/run.sh"; chmod +x "$d/run.sh"
    make_jdl "$B" "$L" "$th${p:+ (part $p/$np)}" "$mem" "$hours" > "$d/job.jdl"
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
