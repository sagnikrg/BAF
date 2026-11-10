#!/bin/bash

cd ../..

L=8

for i1 in {1..1}
do

 condor_submit job_file_Ilattice_L${L}_${i1}.jdl

done
