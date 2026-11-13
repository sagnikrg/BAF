#!/bin/bash
# Writes one HTCondor batch per (kind, L, theta) of the job list (BAF/HOWTO_cluster.md, rule 1), after mbl_dtc_polfed/cat_batch.sh:
#   ~/local/batches/<B>/run.sh      the one run script of the batch (takes the job's numbers as arguments)
#   ~/local/batches/<B>/job.jdl     one submit file, "queue ... from items.txt"
#   ~/local/batches/<B>/items.txt   one line per job: kind L theta g r0 r1 tag
#   ~/local/batches/<B>/log/        err/out of the jobs and one cluster log
#   $BUDDY/batches/<B>/code/        the julia scripts of this batch (a snapshot, so a later git pull never changes running jobs)
# <B> = sb_<kind>_L<L>_theta<theta>_<stamp>; mopthefloor_BAF.sh removes a batch once none of its jobs is in condor_q.
# The list of the batches written is in ~/local/batches/.last (read by condor_submit.sh).
#
# Job list line:  kind L theta g r0 r1 block cpus memGB hours   (realisations r0..r1 in jobs of <block>)
# All lines of one (kind, L, theta) must have the same cpus, memory and hours.
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
JOBS="${JOBS:-$here/jobs.txt}"
STAMP="${STAMP:-$(date +%Y%m%d-%H%M%S)}"
LB=~/local/batches; RB="$BUDDY/batches"
mkdir -p "$LB" "$RB"
: > "$LB/.last"
declare -A RES

make_run() {   # $1 = batch, $2 = cpus
cat <<EOF
#!/bin/bash
# one job of batch $1, arguments: kind L theta g r0 r1 tag
kind=\$1; L=\$2; th=\$3; g=\$4; a=\$5; b=\$6; tag=\$7

# Goes to the job work directory (never go on anywhere else, the end of this script empties it)
cd /jwd || exit 1
source /etc/profile

# the julia depot of the mbl_dtc jobs (julia 1.9.4, HDF5 0.17.2)
cp \$BUDDY/BAF/mbl_dtc/.julia/julia-1.9.4-08-08-24.tar.gz ./
tar -xf julia-1.9.4-08-08-24.tar.gz
rm -f julia-1.9.4-08-08-24.tar.gz
module load julia/1.9.4

# the scripts of this batch and the header
mkdir -p code
cp \$BUDDY/batches/$1/code/*.jl ./code/
cp -rf \$BUDDY/BAF/symmetry_breaking/.header ./

export OPENBLAS_NUM_THREADS=$2
export OMP_NUM_THREADS=$2
export JULIA_NUM_THREADS=1
export JULIA_PKG_OFFLINE=true

julia code/run_sb.jl ./out \$kind \$L \$th \$g \$a \$b > julia_\$tag.log 2>&1
echo "julia exit code \$?" >> julia_\$tag.log

mkdir -p /cephfs/user/sghosh/data/symmetry_breaking/logs
cp julia_\$tag.log /cephfs/user/sghosh/data/symmetry_breaking/logs/.
cp out/*.hdf5 /cephfs/user/sghosh/data/symmetry_breaking/.
rm -rf *
EOF
}

make_jdl() {   # $1 = batch, kind L theta cpus mem hours
cat <<EOF

#Job Script to be submitted using HTCondor (one batch: all jobs of items.txt)

Executable              = run.sh
JobBatchName            = sb $2: L=$3, theta=$4
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=(SUBHOST);

Arguments               = "\$(kind) \$(L) \$(th) \$(g) \$(a) \$(b) \$(tag)"

Universe                = vanilla

Transfer_executable     =
Transfer_input_files    =
Transfer_output_files   =


Error                   = log/err.\$(tag)
Output                  = log/out.\$(tag)
Log                     = log/cluster.log

Request_memory          = $6 GB
Request_cpus            = $5

Request_disk            = 2 GB

+CephFS_IO              = "low"

+MaxRuntimeHours        = $7

+ContainerOS            = "Debian12"

# Submit the jobs, one per line of items.txt
Queue kind,L,th,g,a,b,tag from items.txt

EOF
}

n=0
while read -r kind L th g r0 r1 block cpus mem hours _; do
  [[ -z "${kind:-}" || "$kind" == \#* ]] && continue
  B="sb_${kind}_L${L}_theta${th}_${STAMP}"; d="$LB/$B"
  if [[ -z "${RES[$B]:-}" ]]; then
    RES[$B]="$cpus/$mem/$hours"
    mkdir -p "$d/log" "$RB/$B/code"
    cp "$here"/code/*.jl "$RB/$B/code/"
    make_run "$B" "$cpus" > "$d/run.sh"; chmod +x "$d/run.sh"
    make_jdl "$B" "$kind" "$L" "$th" "$cpus" "$mem" "$hours" > "$d/job.jdl"
    : > "$d/items.txt"
    echo "$B" >> "$LB/.last"
  elif [[ "${RES[$B]}" != "$cpus/$mem/$hours" ]]; then
    echo "ERROR: two lines for $kind L=$L theta=$th with different cpus/memory/hours; give them one" >&2; exit 1
  fi
  for (( a=r0; a<=r1; a+=block )); do
    b=$(( a+block-1 )); (( b>r1 )) && b=$r1
    echo "$kind $L $th $g $a $b sb_${kind}_L${L}_th${th}_g${g}_r${a}-${b}" >> "$d/items.txt"
    n=$((n+1))
  done
done < <(tr -d '\r' < "$JOBS")

for B in $(cat "$LB/.last"); do echo "cat_batch: $B  $(wc -l < "$LB/$B/items.txt") jobs"; done
echo "cat_batch: $n jobs in $(wc -l < "$LB/.last") batches in ~/local/batches"
