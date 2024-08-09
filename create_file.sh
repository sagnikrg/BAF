

cd ../

L=8

#cd $BUDDY

for i1 in {1..1}
do



cat <<EOF >mbldtc_L${L}_${i1}.jl


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









