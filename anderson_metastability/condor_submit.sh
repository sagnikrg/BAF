#!/bin/bash

cd ~/local/


L=23
Itrnumb=60

for (( itr=1; itr<=Itrnumb; itr++ ))
do

 condor_submit job_file_anderson_metastability_L${L}_${itr}.jdl

done
