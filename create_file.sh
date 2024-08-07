#!/bin/bash

#creates job iteratively 

cd ..

L=8

for i1 in {1..1}
do




cat <<EOF >run_file_mbldtcL${L}_${i1}.sh
#!/bin/bash
     
  
     
# Goes to the job work direcrtory     
            
cd /jwd			
source /etc/profile

/home/sghosh/physik/

#Copies, Extracts and removes the Julia tarball 
            
cp /home/sghosh/physik/local/julia/julia-07-08-24.tar.gz ./

tar -xf julia-07-08-24.tar.gz
rm -f julia-07-08-24.tar.gz



#Loads Julia
module load julia/1.9.4

#Creates folder for the Job

mkdir \${ClusterId}_\$Process
cd \${ClusterId}_\$Process


#Job submission 

cp /home/sghosh/physik/local/mbldtc_L${L}_${i1}.jl ./
cp -rf /home/sghosh/physik/local/BAF/header ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Do the real thing here
julia mbldtcL${L}_${i1}.jl
rm mbldtcL{L}_${i1}.jl


# copy results
cp  /jwd/\${ClusterId}_\$Process/*.hdf5 /home/sghosh/physik/local/data/.
EOF







cat <<EOF >job_file_mbldtcL${L}_${i1}.jdl

#Job Script to be submitted using HTCondor
        
Executable              = run_file_mbldtcL${L}_${i1}.sh
        
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=$ENV(SUBHOST);
        
Arguments		= 8
        
Universe                = vanilla

Transfer_executable     = True
Transfer_input_files    = 
Transfer_output_files   =


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



#cd $BUDDY

for i1 in {1..1}
do



cat <<EOF >mbldtcL${L}_${i1}.jl


#########################################################################
# Computing Eigenstatistics of the MBL-DTC Unitary
#########################################################################

# Importing Headers

using Pkg					#   For Package Management
using Dates                 #   For Date and Time
using LinearAlgebra         #   For Linear Algebra
using Random                #   For Random Number Generation
using Kronecker             #   For Kronecker Product
using Arpack                #   For Eigenvalue Decomposition
using StatsBase             #   For Statistics
using HDF5                  #   For Saving Data


include("header/gates.jl")
include("header/brickwall.jl")
include("header/TransferMat.jl")
include("header/functions.jl")

L=${L};
theta=0.0;
epsilon=0.0; 



file=h5open("mbldtc_L8_theta_\$(theta)_${i1}.hdf5","cw")
attrs=attributes(file)


######################################
# Attributes:
######################################

	# Extracting Date and Time

	Dates.now()
	attrs["Date/Time"]=string(Dates.now())


	# Extracting Processor Type

    processor_type = Sys.CPU_NAME
    attrs["[ENV] Processor Type"] = string(processor_type)



	# Extracting Julia Version

	julia_version = VERSION
	attrs["[ENV] Julia Version"] = string(julia_version)

	# Container OS

	attrs["[ENV] Container OS"] = "Debian11"

	# Extracting Number of Threads

	num_threads = Threads.nthreads()
	attrs["[ENV] Number of Threads"] = string(num_threads)	

	# Code

	script_content = read("mbldtcL${L}_${i1}.jl", String)
	attrs["[ENV] Code"] = script_content

	# Modules

	module_list = Pkg.installed()
	attrs["[ENV] Modules"] = string(module_list)


	# Meta Data

#########################################################################
# Brickwall
#########################################################################


U=brickwall(L,theta,epsilon);


#########################################################################
# Eigenstatistics
#########################################################################

eigA,eigvecA=eigen(U);

#########################################################################
# Eigenstate Entanglement Entropy
#########################################################################

#########################################################################
# Lazarides-Luitz Staistics 
#########################################################################


file["L\$(L)/theta\$(theta)/epsilon"*first("\$(epsilon)",4)*"/Itr"]=eigA;
	
close(file)	
	






EOF

done 









