#!/bin/bash

cd ~/local/


L=12
Itrnumb=9

for (( itr=1; itr<=Itrnumb; itr++ ))
do

 condor_submit job_file_anderson_metastability_L${L}_${itr}.jdl

done
