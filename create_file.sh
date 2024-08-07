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


#Copies, Extracts and removes the Julia tarball 
            
cp ~/local/julia/julia-07-08-24.tar.gz ./

tar -xf julia-07-08-24.tar.gz
rm -f julia-07-08-24.tar.gz



#Loads Julia
module load julia/1.9.4

#Creates folder for the Job

mkdir \${ClusterId}_\$Process
cd \${ClusterId}_\$Process


#Job submission 

cp ~/local/mbldtc_L${L}_${i1}.jl ./
cp -rf ~/local/BAF/header ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Do the real thing here
julia mbldtcL${L}_${i1}.jl
rm mbldtcL{L}_${i1}.jl


# copy results
cp  /jwd/\${ClusterId}_\$Process/*.hdf5 ~/local/data/.
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

using LinearAlgebra
using Random
using Kronecker
using Arpack


L=${L};
theta=0.0; 

#########################################################################
# Brickwall
#########################################################################

#########################################################################
# Eigenstatistics
#########################################################################

#########################################################################
# Eigenstate Entanglement Entropy
#########################################################################

#########################################################################
# Lazarides-Luitz Staistics 
#########################################################################


# getting the phase diagram

#for itr in 1:20
itr=1;
	for j in 1:14

file=h5open("EntanglementEntropyL10_theta_\$(thetalist[j])_${i1}.hdf5","w")

		for i in 1:1:51

			epsilon=(i-1)*0.02;
			thetamean=thetalist[j];
			Un=brickwall(L,thetamean,epsilon)
		
			for t in 1:16
				Un.=Un*Un
			end
		
			S,V,D=svd(reshape(Un[683,:],(32,32)));
			V=V.^2;


		file["SV/thetamea\$(thetamean)/epsilon"*first("\$(epsilon)",4)*"/Itr\$(itr)"]=V;
	end
	
close(file)	
	end
#end





EOF

done 









