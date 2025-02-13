#!/bin/bash

cd ../..

fp=0.1

for i1 in {1..1}
do

 condor_submit job_file_OPO_fp${fp}_${i1}.jdl

done
