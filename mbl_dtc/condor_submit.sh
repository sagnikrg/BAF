#!/bin/bash

cd ~/local/


L=12
thetarun=0.2
Itrnumb=500


for (( i1=1; i1<=Itrnumb; i1++ ))
do

 condor_submit job_file_mbldtcL${L}_theta${thetarun}_${i1}.jdl

done
