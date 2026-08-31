#!/bin/bash

cd ~/local/

L=8
thetarun=0.0
Itrnumb=16

for (( i1=1; i1<=Itrnumb; i1++ ))
do

 condor_submit job_file_mbldtcL${L}_${i1}.jdl

done
