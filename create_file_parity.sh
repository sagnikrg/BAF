

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
Itrnumb=100;

epsilonlist=[0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0]

global Parity_full=zeros(51)
global Parity_normal=zeros(51)



file=h5open("mbldtc_parity_L${L}_theta_\$(theta)_${i1}.hdf5","cw")
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

	script_content = read("mbldtc_L${L}_${i1}.jl", String)
	attrs["[ENV] Code"] = script_content

	# Modules

	module_list = Pkg.installed()
	attrs["[ENV] Modules"] = string(module_list)

	# Julia Environment

	attrs["[ENV] Julia Environment"] = "julia-1.9.4-08-08-24.tar.gz" 

	# Meta Data

	attrs["METADATA"] = "This Data file contains a full parameter scan of epsilon in epsilonlint of the MBL-DTC Unitary for L=\$(L) and theta=\$(theta). The data contains the overlap of catstates and preservation of the catstatedness under kick."

	# Author

	attrs["Author"] = "Sagnik Ghosh"

	# cluster

	attrs["Cluster"] = "BAF"

	# Parameters

	attrs["[Parameters] h"] = "2pi"
	attrs["[Parameters] L"] = L
	attrs["[Parameters] theta"] = theta
	attrs["[Parameters] epsilon"] = string(epsilonlist)

	attrs["[Parameters] Itrnumb"] = Itrnumb


######################################
# Code:
######################################	


for itr in 1:Itrnumb


for (i,epsilon) in enumerate(epsilonlist)
   
    U=brickwall(L,theta,epsilon);

    eigvals,eigvecs=eigen(brickwall(L,theta,epsilon));
    eigvals,eigvecs=phase_ordered_eigvecs(eigvals,eigvecs);

    eigvecsnew=copy(eigvecs)

    Lhalve=L-1

    for i in 1:2^Lhalve


        eigvecsnew[:,i]=(eigvecs[:,i]+eigvecs[:,i+2^Lhalve])/sqrt(2)
        eigvecsnew[:,i+2^Lhalve]=(eigvecs[:,i]-eigvecs[:,i+2^Lhalve])/sqrt(2)

    end
    
    parityeig=0.0
	parityeig_normal=0.0

    kicks=kick(L,epsilon)
	kicks_normal=kick(L,0.0)

    for i in 1:2^Lhalve
        parityeig+=abs(transpose(conj.(eigvecsnew[:,i+2^Lhalve])) *kicks*eigvecsnew[:,i])
        parityeig+=abs(transpose(conj.(eigvecsnew[:,i])) *kicks*eigvecsnew[:,i+2^Lhalve])
    end


	Parity_full[i]+=parityeig

	for i in 1:2^Lhalve
		parityeig_normal+=abs(transpose(conj.(eigvecsnew[:,i+2^Lhalve])) *kicks_normal*eigvecsnew[:,i])
		parityeig_normal+=abs(transpose(conj.(eigvecsnew[:,i])) *kicks_normal*eigvecsnew[:,i+2^Lhalve])
	end

	Parity_normal[i]+=parityeig_normal


	## epsilon loop ends here
end

## Itrnumb loop ends here
end

Parity_full=Parity_full/Itrnumb
Parity_normal=Parity_normal/(Itrnumb)

file["L\$(L)/theta\$(theta)/Parity_full"]=Parity_full;
file["L\$(L)/theta\$(theta)/Parity_normal"]=Parity_normal;
	
close(file)	
	






EOF

done 









