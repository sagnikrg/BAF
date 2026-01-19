#!/bin/bash

#creates job iteratively 

cd ~/local/

L=10
Itrnumb=10


for (( itr=1; itr<=Itrnumb; itr++ ))
do




cat <<EOF >job_file_anderson_metastability_L${L}_${itr}.jdl

#Job Script to be submitted using HTCondor


Executable              = run_file_anderson_metastability_L${L}_${itr}.sh
JobBatchName            = anderson_metastability: L=${L}, Itr=${itr}        
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=$ENV(SUBHOST);
        
Arguments		= 8
        
Universe                = vanilla

Transfer_executable     = 
Transfer_input_files    = 
Transfer_output_files   =


Error                   = log/err.anderson_metastability_L${L}_${itr}.log
Output                  = log/out.anderson_metastability_L${L}_${itr}.log
Log                     = 

Request_memory          = 12  GB
Request_cpus            = 1
        
Request_disk            = 4 GB

+CephFS_IO  = "low"
        
+MaxRuntimeHours	=   24
        
+ContainerOS        = "Debian12"

# Submit job
Queue 1        
            
EOF


done

