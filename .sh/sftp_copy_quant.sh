



#!/usr/bin/env bash
#set -euo pipefail


# global automated parameters:

repo=anderson_metastability   
quan=gamma

#paths


cd ~
mkdir physik_mbl_dtc/${repo}/${quan}

for (( L=16; L<=16; L++ ))
do


mkdir  physik_mbl_dtc/${repo}/${quan}/L${L}
rsync -avz --progress $BUDDY/data/${repo}_${quan}_${L}_*.hdf5 physik_mbl_dtc/${repo}/${quan}/L${L}/.

rm -r  $BUDDY/data/${repo}_${quan}_${L}_*.hdf5  



done

