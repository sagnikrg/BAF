#!/bin/bash

#creates job iteratively 

cd ../..

L=8
Itrnumb=1


for itr in {1..Itrnumb}
do








cat <<EOF >job_file_eigenscatterL${L}_${itr}.jdl

#Job Script to be submitted using HTCondor


Executable              = run_file_eigenscatterL${L}_${itr}.sh
JobBatchName            = OPO: L=${L}, Itr=${itr}        
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=$ENV(SUBHOST);
        
Arguments		= 8
        
Universe                = vanilla

Transfer_executable     = 
Transfer_input_files    = 
Transfer_output_files   =


Error                   = log/err.eigenscatterL${L}_${itr}.log
Output                  = log/out.eigenscatterL${L}_${itr}.log
Log                     = 

Request_memory          = 12  GB
Request_cpus            = 1
        
Request_disk            = 4 GB

+CephFS_IO  = "low"
        
+MaxRuntimeHours	=   12
        
+ContainerOS        = "Debian12"

# Submit job
Queue 1        
            
EOF


done

