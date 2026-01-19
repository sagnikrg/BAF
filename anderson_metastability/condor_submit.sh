#!/bin/bash

cd ~/local/


L=9
Itrnumb=11

for (( itr=1; itr<=Itrnumb; itr++ ))
do

 condor_submit job_file_anderson_metastability_L${L}_${itr}.jdl

done
