



#!/usr/bin/env bash
#set -euo pipefail


# global automated parameters:

repo=anderson_metastability   
quan=eigendata

#paths


cd ~
mkdir -p physik_mbl_dtc/${repo}/EIGENDATA

for (( L=20; L<=24; L++ ))
do


mkdir -p physik_mbl_dtc/${repo}/EIGENDATA/L${L}
rsync -avz --progress $BUDDY/data/${repo}_${quan}_${L}_*.hdf5 physik_mbl_dtc/${repo}/EIGENDATA/L${L}/.

rm -r  $BUDDY/data/${repo}_${quan}_${L}_*.hdf5  



done

