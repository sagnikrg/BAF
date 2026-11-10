#!/bin/bash
# Writes one run_file per job of the job list into ~/local (see BAF.sh)
here=$(cd "$(dirname "$0")" && pwd)
cd ~/local/
grep -v '^\s*#' "${JOBS:-$here/jobs.txt}" | grep -v '^\s*$' | while read -r kind L th g r0 r1 block cpus mem hours; do
for (( a=r0; a<=r1; a+=block )); do
b=$(( a+block-1 )); (( b>r1 )) && b=$r1
name=sb_${kind}_L${L}_th${th}_g${g}_r${a}-${b}
cat <<EOS >run_file_${name}.sh
#!/bin/bash
cd /jwd
source /etc/profile

# the julia depot of the mbl_dtc jobs (julia 1.9.4, HDF5 0.17.2), same as there
cp \$BUDDY/BAF/mbl_dtc/.julia/julia-1.9.4-08-08-24.tar.gz ./
tar -xf julia-1.9.4-08-08-24.tar.gz
rm -f julia-1.9.4-08-08-24.tar.gz
module load julia/1.9.4

mkdir -p code
cp \$BUDDY/BAF/symmetry_breaking/code/*.jl ./code/
cp -rf \$BUDDY/BAF/symmetry_breaking/.header ./

export OPENBLAS_NUM_THREADS=${cpus}
export OMP_NUM_THREADS=${cpus}
export JULIA_NUM_THREADS=1

julia code/run_sb.jl ./out ${kind} ${L} ${th} ${g} ${a} ${b}

mkdir -p /cephfs/user/sghosh/data/symmetry_breaking
cp out/*.hdf5 /cephfs/user/sghosh/data/symmetry_breaking/.
rm -rf *
EOS
done
done
