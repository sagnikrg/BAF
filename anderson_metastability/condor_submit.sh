#!/bin/bash

cd ../..

Itrnumb=100

for itr in {1..1}
do

 condor_submit job_file_eigenscatterL${L}_${itr}.jdl

done
