#!/bin/bash

#creates job iteratively 

cd ~/local/


L=8
thetarun=0.0
Itrnumb=1


for (( i1=1; i1<=Itrnumb; i1++ ))
do




cat <<EOF >job_file_mbldtcL${L}_theta${thetarun}_${i1}.jdl

#Job Script to be submitted using HTCondor
        
Executable              = run_file_mbldtcL${L}_theta${thetarun}_${i1}.sh
JobBatchName            = mbldtc: L=${L}, theta=${thetarun}, Itr=${i1}                
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=$ENV(SUBHOST);
        
Arguments		= 8
        
Universe                = vanilla

Transfer_executable     = 
Transfer_input_files    = 
Transfer_output_files   =


Error                   = log/err.mbldtcL${L}_theta${thetarun}_${i1}
Output                  = log/out.mbldtcL${L}_theta${thetarun}_${i1}
Log                     = log/log.mbldtcL${L}_theta${thetarun}_${i1}

Request_memory          = 12  GB
Request_cpus            = 1
        
Request_disk            = 2 GB

+CephFS_IO  = "low"
        
+MaxRuntimeHours	= 48
        
+ContainerOS            = "Debian11"

# Submit job
Queue 1        
            
EOF




done



