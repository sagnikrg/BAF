#!/bin/bash

#creates job iteratively 

cd ~/local/

L=8
thetarun=0.0
Itrnumb=1

for (( i1=1; i1<=Itrnumb; i1++ ))
do




cat <<EOF >job_file_mbldtcL${L}_${i1}.jdl

#Job Script to be submitted using HTCondor
        
Executable              = run_file_mbldtcL${L}_${i1}.sh
JobBatchName            = mbldtc: L=${L}, theta=${thetarun} Itr=${i1}                
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=$ENV(SUBHOST);
        
Arguments		= 8
        
Universe                = vanilla

Transfer_executable     = True
Transfer_input_files    = 
Transfer_output_files   =


Error                   = log/err.\$(ClusterId).\$(Process)
Output                  = log/out.\$(ClusterId).\$(Process)
Log                     = log/log.\$(ClusterId).\$(Process)

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



