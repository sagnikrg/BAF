#!/bin/bash

cd ~/local/

L=10
thetarun=0.1
Itrnumb=200


for (( i1=1; i1<=Itrnumb; i1++ ))
do

 condor_submit job_file_mbldtcL${L}_theta${thetarun}_${i1}.jdl

done
