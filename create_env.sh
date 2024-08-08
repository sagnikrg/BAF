#!/bin/bash

#creates job iteratively 

cd ../data

L=8

for i1 in {1..1}
do




cat <<EOF >run_file_mbldtcL${L}_${i1}.sh
#!/bin/bash
     
  
     
# Goes to the job work direcrtory     
            
cd /jwd			
source /etc/profile



#Copies, Extracts and removes the Julia tarball 
            
cp \$BUDDY/BAF/julia/julia-1.9.4-08-08-24.tar.gz ./

tar -xf julia-1.9.4-08-08-24.tar.gz
rm -f julia-1.9.4-08-08-24.tar.gz



#Loads Julia
module load julia/1.9.4

#Creates folder for the Job

mkdir \${ClusterId}_\$Process
cd \${ClusterId}_\$Process


#Job submission 

cp \$BUDDY/mbldtc_L${L}_${i1}.jl ./
cp -rf \$BUDDY/BAF/header ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Do the real thing here
julia mbldtc_L${L}_${i1}.jl
rm mbldtc_L${L}_${i1}.jl


# copy results
#cp  /jwd/\${ClusterId}_\$Process/*.hdf5 /home/sghosh/physik/local/data/.
EOF







cat <<EOF >job_file_mbldtcL${L}_${i1}.jdl

#Job Script to be submitted using HTCondor
        
Executable              = run_file_mbldtcL${L}_${i1}.sh
        
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=$ENV(SUBHOST);
        
Arguments		= 8
        
Universe                = vanilla

Transfer_executable     = True
Transfer_input_files    = 
#Transfer_output_files   =


Error                   = log/err.\$(ClusterId).\$(Process)
Output                  = log/out.\$(ClusterId).\$(Process)
Log                     = log/log.\$(ClusterId).\$(Process)

Request_memory          = 8  GB
Request_cpus            = 1
        
Request_disk            = 2 GB

+CephFS_IO  = "low"
        
+MaxRuntimeHours	= 1
        
+ContainerOS            = "Debian11"

# Submit job
Queue 1        
            
EOF


done



