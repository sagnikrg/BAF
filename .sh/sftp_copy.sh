



#!/usr/bin/env bash
#set -euo pipefail


# global automated parameters:

repo="anderson_metastability"   
quan="eigendata"

#paths

source_path="cephfs/user/sghosh/data/"
dest_path="~/physik_mbl_dtc/${repo}/"




#for (( L=10; itr<=10; itr++ ))
#do

L=10

mkdir "${dest_path}/EIGENDATA/L${L}"
rsync -avz --progress "${source_path}/L${L}/${repo}_${quan}_${L}_*.hdf5" "${dest_path}/EIGENDATA/L${L}/"






#done

