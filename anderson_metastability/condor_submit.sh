#!/bin/bash

cd ~/local/


L=11
Itrnumb=10

for (( itr=1; itr<=Itrnumb; itr++ ))
do

 condor_submit job_file_anderson_metastability_L${L}_${itr}.jdl

done
