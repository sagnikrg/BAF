#!/bin/bash
# creates one run_file_<tag>.sh per job in ~/local (jobs from the job list, see common.sh)

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
cd ~/local/

make_run() {
local L=$1 th=$2 r=$3 j=$4 eps=$5 mem=$6 hours=$7 tag=$8 set=$9
cat <<EOF >run_file_${tag}.sh
#!/bin/bash

# Goes to the job work directory
cd /jwd
source /etc/profile

# Copies, extracts and removes the Julia tarball
cp \$BUDDY/BAF/mbl_dtc_polfed/.julia/${TARBALL} ./
tar -xf ${TARBALL}
rm -f ${TARBALL}

lscpu --json | grep "Model name" | awk -F '"' '{print \$8}' > model_name.txt

# Loads Julia
module load ${JULIA_MODULE}

# the script and the header (itensor and polfed only)
cp \$BUDDY/polfed_bands.jl ./
cp -rf \$BUDDY/BAF/mbl_dtc_polfed/.header ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Do the real thing here
julia polfed_bands.jl ${L} ${th} ${r} ${j} ${eps} 64 ${set} > julia_${tag}.log 2>&1
echo "julia exit code \$?" >> julia_${tag}.log
rm polfed_bands.jl

# Clean up the working directory
mkdir -p /cephfs/user/sghosh/data/mbl_dtc_polfed/logs
cp julia_${tag}.log /cephfs/user/sghosh/data/mbl_dtc_polfed/logs/.
cp *hdf5 /cephfs/user/sghosh/data/mbl_dtc_polfed/.
rm -rf *
EOF
}

each_job make_run
echo "cat_run: $(count_jobs) run files in ~/local"
