#!/bin/bash

#creates job iteratively 

cd ~/local/

L=8
Itrnumb=1


for (( itr=1; itr<=Itrnumb; itr++ ))
do



cat <<EOF >run_file_anderson_metastability_L${L}_${itr}.sh
#!/bin/bash
     
  
     
# Goes to the job work direcrtory     
            
cd /jwd			
source /etc/profile



#Copies, Extracts and removes the Julia tarball 
            
cp \$BUDDY/BAF/anderson_metastability/.julia/julia_1.11.1_15.01.26.tar.gz ./

cp \$BUDDY/BAF/anderson_metastability/METADATA_rawdata.txt ./
cp \$BUDDY/BAF/anderson_metastability/METADATA_gamma.txt ./
cp \$BUDDY/BAF/anderson_metastability/METADATA_midspectra_gap.txt ./
cp \$BUDDY/BAF/anderson_metastability/METADATA_IPR_realbasis.txt ./
cp \$BUDDY/BAF/anderson_metastability/METADATA_IPR_eigenbasis.txt ./


tar -xf julia_1.11.1_15.01.26.tar.gz
rm -f julia_1.11.1_15.01.26.tar.gz

 lscpu | grep 'Model name: '> model_name.txt



#Loads Julia
module load julia/1.11.1


#Job submission 

cp \$BUDDY/anderson_metastability_L${L}_${itr}.jl ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Do the real thing here

julia anderson_metastability_L${L}_${itr}.jl
rm anderson_metastability_L${L}_${itr}.jl



# Clean up the working directory

cp  * /cephfs/user/sghosh/data/.
rm -rf *
EOF






done

