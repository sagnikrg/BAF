#!/bin/bash

#creates job iteratively 

cd ../..

L=8
Itrnumb=1


for itr in {1..Itrnumb}
do



cat <<EOF >run_file_eigenscatterL${L}_${itr}.sh
#!/bin/bash
     
  
     
# Goes to the job work direcrtory     
            
cd /jwd			
source /etc/profile



#Copies, Extracts and removes the Julia tarball 
            
cp \$BUDDY/julia/julia-1.11.1-anderson-Oct-30-25.tar.gz ./

tar -xf julia-1.11.1-anderson-Oct-30-25.tar.gz
rm -f julia-1.11.1-anderson-Oct-30-25.tar.gz

 lscpu | grep 'Model name: '> model_name.txt



#Loads Julia
module load julia/1.9.4


#Job submission 

cp \$BUDDY/eigenscatterL${L}_${itr}.jl ./
cp -rf \$BUDDY/BAF/anderson-metastability/.header ./

# set number of threads
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

# Do the real thing here
julia eigenscatterL${L}_${itr}.jl
rm eigenscatterL${L}_${itr}.jl


# copy results
cp  *.hdf5 /cephfs/user/sghosh/data/.
rm -rf *
EOF






done

