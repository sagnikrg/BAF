#!/bin/bash
# creates one run_file_<tag>.sh per job in ~/local (jobs from jobs.txt, see common.sh)

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"
cd ~/local/

make_run() {
local L=$1 th=$2 itr=$3 i=$4 eps=$5 mem=$6 tag=$7
cat <<EOF >run_file_${tag}.sh
#!/bin/bash

# Goes to the job work directory
cd /jwd
source /etc/profile

# Copies, extracts and removes the Julia tarball
cp \$BUDDY/BAF/mbl_dtc/.julia/julia-1.9.4-08-08-24.tar.gz ./
tar -xf julia-1.9.4-08-08-24.tar.gz
rm -f julia-1.9.4-08-08-24.tar.gz

lscpu --json | grep "Model name" | awk -F '"' '{print \$8}' > model_name.txt

# Loads Julia
module load julia/1.9.4

# the script and the header
cp \$BUDDY/eigen_win.jl ./
cp -rf \$BUDDY/BAF/mbl_dtc/.header ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Do the real thing here
julia eigen_win.jl ${L} ${th} ${itr} ${i} ${eps} > julia_${tag}.log 2>&1
echo "julia exit code \$?" >> julia_${tag}.log
rm eigen_win.jl

# Clean up the working directory
mkdir -p /cephfs/user/sghosh/data/mbl_dtc/logs
cp julia_${tag}.log /cephfs/user/sghosh/data/mbl_dtc/logs/.
cp *hdf5 /cephfs/user/sghosh/data/mbl_dtc/.
rm -rf *
EOF
}

each_job make_run
echo "cat_run: $(count_jobs) run files in ~/local"
