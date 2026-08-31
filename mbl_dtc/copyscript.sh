#!/bin/bash

#Copies file Iteratively to physik_mbl_dtc


for i1 in {1..90};
do

src_file="/cephfs/user/sghosh/Data/$((i1+3372))_0/spin1_2dtc_L8_Eig.hdf5"
src_file_2="/cephfs/user/sghosh/Data/$((i1+3372))_0/spin1_2dtc_L8_eigvals_${i1}.hdf5"
dst_file="."

#mv "${src_file}" "${src_file_2}"
scp "${src_file_2}" "${dst_file}"

done

