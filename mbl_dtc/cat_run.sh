#!/bin/bash

#creates job iteratively 

cd ~/local/

L=10
thetarun=0.2
Itrnumb=230





for (( i1=1; i1<=Itrnumb; i1++ ))
do





cat <<EOF >run_file_mbldtc_h0_L${L}_theta${thetarun}_${i1}.sh
#!/bin/bash
     
  
     
# Goes to the job work direcrtory     
            
cd /jwd			  
source /etc/profile



#Copies, Extracts and removes the Julia tarball 
            
cp \$BUDDY/BAF/mbl_dtc/.julia/julia-1.9.4-08-08-24.tar.gz ./



tar -xf julia-1.9.4-08-08-24.tar.gz
rm -f julia-1.9.4-08-08-24.tar.gz

lscpu --json | grep "Model name" | awk -F '"' '{print $8}' > model_name.txt


#Loads Julia
module load julia/1.9.4




#Job submission 

cp \$BUDDY/mbldtc_h0_L${L}_theta${thetarun}_${i1}.jl ./
cp -rf \$BUDDY/BAF/mbl_dtc/.header ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Do the real thing here
julia mbldtc_h0_L${L}_theta${thetarun}_${i1}.jl
rm mbldtc_h0_L${L}_theta${thetarun}_${i1}.jl



# Clean up the working directory

cp  *hdf5 /cephfs/user/sghosh/data/.
rm -rf *
EOF





done



