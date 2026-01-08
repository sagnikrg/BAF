#!/bin/bash

cd ../

L=12

for i1 in {1..2000}
do

 condor_submit job_file_mbldtcL${L}_${i1}.jdl

done
