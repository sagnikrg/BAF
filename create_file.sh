#!/bin/bash

#creates job iteratively 


for i1 in {1..1000}
do




cat <<EOF >run_file_hcee10_${i1}.sh
#!/bin/bash
     
  
     
# Goes to the job work direcrtory     
            
cd /jwd			
source /etc/profile


#Copies, Extracts and removes the Julia tarball 
            
cp \$BUDDY/julia-1.8.3.tar.gz ./

tar -xf julia-1.8.3.tar.gz
rm -f julia-1.8.3.tar.gz



#Loads Julia
module load julia/1.8.3

#Creates folder for the Job

mkdir \${ClusterId}_\$Process
cd \${ClusterId}_\$Process


#Job submission 

cp \$BUDDY/halfchainL10_${i1}.jl ./
cp -rf \$BUDDY/header ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Do the real thing here
julia halfchainL10_${i1}.jl
rm halfchainL10_${i1}.jl

# Make sure the folder exists
#mkdir -p \$BUDDY/log


# copy results
#cp -r /jwd/\${ClusterId}_\$Process \$BUDDY/Data
cp  /jwd/\${ClusterId}_\$Process/*.hdf5 \$BUDDY/Data
EOF







cat <<EOF >job_file_hcee8_${i1}.jdl

#Job Script to be submitted using HTCondor
        
Executable              = run_file_hcee10_${i1}.sh
        
Environment             = ClusterId=\$(ClusterId);Process=\$(Process);SubHost=$ENV(SUBHOST);
        
Arguments		= 8
        
Universe                = vanilla

Transfer_executable     = True
Transfer_input_files    = 
Transfer_output_files   =


Error                   = log/err.\$(ClusterId).\$(Process)
Output                  = log/out.\$(ClusterId).\$(Process)
Log                     = log/log.\$(ClusterId).\$(Process)

Request_memory          = 40  GB
Request_cpus            = 1
        
Request_disk            = 2 GB

+CephFS_IO  = "low"
        
+MaxRuntimeHours	= 48
        
+ContainerOS            = "Debian11"

# Submit job
Queue 1        
            
EOF


done



cd $BUDDY

for i1 in {1..1000}
do



cat <<EOF >halfchainL10_${i1}.jl


#########################################################################
# Computing Half Chain Entanglement Entropy for the Spin-1/2 DTC model 
#########################################################################

# Importing Headers

include("header/header.jl")



L=10;
thetalist=[0.0,0.01,0.02,0.03,0.04,0.06,0.07,0.08,0.09,0.11,0.12,0.13,0.14,0.16,0.17,0.18,0.19];


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









