#!/bin/bash

cd ~/local/

L=8
thetarun=0.25
Itrnumb=4





for (( i1=1; i1<=Itrnumb; i1++ ))
do

 condor_submit job_file_mbldtc_L${L}_theta${thetarun}_${i1}.jdl

done
