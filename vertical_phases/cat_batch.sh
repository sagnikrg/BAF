#!/bin/bash
# Writes one HTCondor batch per line (L, g) of the job list (HOWTO_cluster.md rule 1), template mbl_dtc_polfed/cat_batch.sh:
#   ~/local/batches/<B>/run.sh      the one run script of the batch (takes the job's numbers as arguments)
#   ~/local/batches/<B>/job.jdl     one submit file, "queue ... from items.txt"
#   ~/local/batches/<B>/items.txt   one line per job: L g r0 r1 thetas(: separated or -) tag
#   ~/local/batches/<B>/log/        err/out of the jobs and one cluster log
#   $BUDDY/batches/<B>/             snapshot of code/*.jl and .header (a later git pull never changes running jobs)
# <B> = vertical_L<L>_g<g>_<stamp>; mopthefloor_BAF.sh removes a batch once none of its jobs is in condor_q.
#
# Job list line:   L  g  r0  r1  block  cpus  memory_GB  hours  [thetas]
#   one job per block of realisations r0..r1; thetas comma separated (default: the grid of run_vertical.jl)

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
JOBS="${JOBS:-$here/jobs.txt}"
STAMP="${STAMP:-$(date +%Y%m%d-%H%M%S)}"
LB=~/local/batches; RB="$BUDDY/batches"
TARBALL=julia-1.9.4-08-08-24.tar.gz
mkdir -p "$LB" "$RB"
: > "$LB/.last"

make_run() {   # $1 = batch, $2 = cpus
cat <<EOF
#!/bin/bash
# one job of batch $1, arguments: L g r0 r1 thetas(: separated or -) tag
L=\$1; g=\$2; a=\$3; b=\$4; th=\${5//:/,}; tag=\$6
[[ "\$th" == "-" ]] && th=""

# Goes to the job work directory (never go on anywhere else, the end of this script empties it)
cd /jwd || exit 1
source /etc/profile

# the julia depot of the mbl_dtc jobs (julia 1.9.4, HDF5 0.17.2)
cp \$BUDDY/BAF/mbl_dtc/.julia/${TARBALL} ./
tar -xf ${TARBALL}
rm -f ${TARBALL}
module load julia/1.9.4

# the code snapshot of this batch
cp \$BUDDY/batches/$1/*.jl ./
cp -rf \$BUDDY/batches/$1/.header ./

export OPENBLAS_NUM_THREADS=$2
export OMP_NUM_THREADS=$2
export JULIA_NUM_THREADS=1
export JULIA_PKG_OFFLINE=true

julia run_vertical.jl ./out \$L \$g \$a \$b \$th > julia_\$tag.log 2>&1
echo "julia exit code \$?" >> julia_\$tag.log

mkdir -p /cephfs/user/sghosh/data/vertical_phases/logs
cp julia_\$tag.log /cephfs/user/sghosh/data/vertical_phases/logs/.
cp out/*.hdf5 /cephfs/user/sghosh/data/vertical_phases/.
rm -rf *
EOF
}

make_jdl() {   # $1 = batch, L g cpus mem hours
cat <<EOF

#Job Script to be submitted using HTCondor (one batch: all jobs of items.txt)

Executable              = run.sh
JobBatchName            = vertical_phases: L=$2, g=$3
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=(SUBHOST);

Arguments               = "\$(L) \$(g) \$(a) \$(b) \$(th) \$(tag)"

Universe                = vanilla

Transfer_executable     =
Transfer_input_files    =
Transfer_output_files   =

Error                   = log/err.\$(tag)
Output                  = log/out.\$(tag)
Log                     = log/cluster.log

Request_memory          = $5  GB
Request_cpus            = $4

Request_disk            = 2 GB

+CephFS_IO  = "low"

+MaxRuntimeHours	= $6

+ContainerOS            = "Debian12"

# Submit the jobs, one per line of items.txt
Queue L,g,a,b,th,tag from items.txt

EOF
}

declare -A RES
while read -r L g r0 r1 block cpus mem hours thetas; do
  [[ -z "${L:-}" || "$L" == \#* ]] && continue
  B="vertical_L${L}_g${g}_${STAMP}"
  d="$LB/$B"
  if [[ -z "${RES[$B]:-}" ]]; then          # first line of this (L, g): write the batch
    RES[$B]="$cpus/$mem/$hours"
    mkdir -p "$d/log" "$RB/$B"
    cp "$here"/code/*.jl "$RB/$B/"
    cp -rf "$here/.header" "$RB/$B/"
    make_run "$B" "$cpus" > "$d/run.sh"; chmod +x "$d/run.sh"
    make_jdl "$B" "$L" "$g" "$cpus" "$mem" "$hours" > "$d/job.jdl"
    : > "$d/items.txt"
    echo "$B" >> "$LB/.last"
  elif [[ "${RES[$B]}" != "$cpus/$mem/$hours" ]]; then
    echo "ERROR: two lines for L=$L g=$g with different cpus/memory/hours; give them one" >&2; exit 1
  fi
  for (( a=r0; a<=r1; a+=block )); do
    b=$(( a+block-1 )); (( b>r1 )) && b=$r1
    th=${thetas//,/:}; echo "$L $g $a $b ${th:--} vertical_L${L}_g${g}_r${a}-${b}" >> "$d/items.txt"
  done
done < <(tr -d '\r' < "$JOBS")
for B in $(cat "$LB/.last"); do echo "cat_batch: $B  $(wc -l < "$LB/$B/items.txt") jobs"; done
