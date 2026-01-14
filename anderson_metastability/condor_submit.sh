#!/bin/bash

cd ../..
L=8
Itrnumb=1

for (( itr=1; itr<=Itrnumb; itr++ ))
do

 condor_submit job_file_eigenscatterL${L}_${itr}.jdl

done
