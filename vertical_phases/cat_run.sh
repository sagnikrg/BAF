#!/bin/bash
# Writes one run_file per job of jobs.txt into ~/local (see BAF.sh, repo="vertical_phases")
here=$(cd "$(dirname "$0")" && pwd)
cd ~/local/
grep -v '^\s*#' "${JOBS:-$here/jobs.txt}" | grep -v '^\s*$' | while read -r L g r0 r1 block cpus mem hours thetas; do
for (( a=r0; a<=r1; a+=block )); do
b=$(( a+block-1 )); (( b>r1 )) && b=$r1
name=vertical_L${L}_g${g}_r${a}-${b}
cat <<EOS >run_file_${name}.sh
#!/bin/bash
cd /jwd
source /etc/profile

# the julia depot of the mbl_dtc jobs (julia 1.9.4, HDF5 0.17.2), same as there
cp \$BUDDY/BAF/mbl_dtc/.julia/julia-1.9.4-08-08-24.tar.gz ./
tar -xf julia-1.9.4-08-08-24.tar.gz
rm -f julia-1.9.4-08-08-24.tar.gz
module load julia/1.9.4

cp \$BUDDY/BAF/vertical_phases/code/*.jl ./
cp -rf \$BUDDY/BAF/vertical_phases/.header ./

export OPENBLAS_NUM_THREADS=${cpus}
export OMP_NUM_THREADS=${cpus}
export JULIA_NUM_THREADS=1

julia run_vertical.jl ./out ${L} ${g} ${a} ${b} ${thetas}

mkdir -p /cephfs/user/sghosh/data/vertical_phases
cp out/*.hdf5 /cephfs/user/sghosh/data/vertical_phases/.
rm -rf *
EOS
done
done
