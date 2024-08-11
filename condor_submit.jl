#!/bin/bash

cd ../

L=8

for i1 in {1..10}
do

 condor_submit job_file_mbldtcL${L}_${i1}.jdl

done