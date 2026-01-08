#!/bin/bash

#creates job iteratively 

cd ../..
L=0.1



for i1 in {1..1}
do







cat <<EOF >run_file_Ilattice_L${L}_${i1}.sh
#!/bin/bash
     
  
     
# Goes to the job work direcrtory     
            
cd /jwd			
source /etc/profile



#Copies, Extracts and removes the Julia tarball 
            
cp \$BUDDY/julia/julia-1.9.4-08-08-24.tar.gz ./

tar -xf julia-1.9.4-08-08-24.tar.gz
rm -f julia-1.9.4-08-08-24.tar.gz

 lscpu | grep 'Model name: '> model_name.txt



#Loads Julia
module load julia/1.9.4


#Job submission 

cp \$BUDDY/Ilattice_L_${L}_${i1}.jl ./
cp -rf \$BUDDY/BAF/InfoLattice/header ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Do the real thing here
julia Ilattice_L_${L}_${i1}.jl
rm Ilattice_L_${L}_${i1}.jl


# copy results
cp  *.hdf5 /cephfs/user/sghosh/data/.
rm -rf *
EOF







cat <<EOF >job_file_Ilattice_L${L}_${i1}.jdl

#Job Script to be submitted using HTCondor


Executable              = run_file_Ilattice_L${L}_${i1}.sh
JobBatchName            = Info Lattice: L=${L}, Itr=${i1}        
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=$ENV(SUBHOST);
        
Arguments		= 8
        
Universe                = vanilla

Transfer_executable     = 
Transfer_input_files    = 
Transfer_output_files   =


Error                   = log/err.InfoLattice_L_${L}_${i1}.log
Output                  = 
Log                     = 

Request_memory          = 10  GB
Request_cpus            = 1
        
Request_disk            = 2 GB

+CephFS_IO  = "low"
        
+MaxRuntimeHours	=   24
        
+ContainerOS        = "Debian12"

# Submit job
Queue 1        
            
EOF


done



