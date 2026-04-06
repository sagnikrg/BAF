



#!/usr/bin/env bash
#set -euo pipefail


# global automated parameters:

repo=anderson_metastability   
quan=eigendata

#paths


cd ~

for (( L=11; L<=12; L++ ))
do

#L=10

mkdir  physik_mbl_dtc/${repo}/EIGENDATA/L${L}
rsync -avz --progress $BUDDY/data/${repo}_${quan}_${L}_*.hdf5 physik_mbl_dtc/${repo}/EIGENDATA/L${L}/.




done

