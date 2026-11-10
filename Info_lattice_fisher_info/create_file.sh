

cd $BUDDY

fp=0.1

#cd $BUDDY

#Itrnumb=1

for i1 in {1..1}
do



cat <<EOF >Ilattice_L_${L}_${i1}.jl




#########################################################################
# Computation of Information Lattice for the  DTC eigenstate
# for various epsilon
#########################################################################

# Importing Headers

using Pkg					#   For Package Management
using Dates                 #   For Date and Time
using HDF5                  #   For Saving Data


include("header/Headers.jl")
include("header/info_lattice.jl");

L=8
theta=0.0
epsilonlist=[0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0]





file=h5open("Infolattice_L_${L}_${i1}.hdf5","cw")
attrs=attributes(file)


######################################
# Attributes:
######################################

	# Extracting Date and Time

	start_time=Dates.now()
	attrs["[Benchmark] Date/Time"]=Dates.format(Dates.now(), " dd/MM/yyyy at HH:mm:ss")

	

	# Extracting Processor Type

    
	model_name = read_model_name("model_name.txt")
    attrs["[Benchmark] Processor Type"] = string(model_name)



	# Extracting Julia Version

	julia_version = VERSION
	attrs["[ENV] Julia Version"] = string(julia_version)

	# Container OS

	attrs["[ENV] Container OS"] = "Debian12"

	# Extracting Number of Threads

	num_threads = Threads.nthreads()
	attrs["[ENV] Number of Threads"] = string(num_threads)	

	# Code

	script_content = read("Ilattice_L_${L}_${i1}.jl", String)
	attrs["[ENV] Code"] = script_content

	# Header

	header_content = read("header/Brickwall.jl", String)
	attrs["[ENV] Header"] = header_content


	mod_content = read("header/info_lattice.jl", String)
	attrs["[ENV] Info Lattice"] = mod_content

	# Modules

	module_list = Pkg.installed()
	attrs["[ENV] Modules"] = string(module_list)

	# Julia Environment

	attrs["[ENV] Julia Environment"] = "julia-1.9.4-08-08-24.tar.gz" 

	# Meta Data

	attrs["METADATA"] = "This Data file contains scan over epsilon for the computation of Information Lattice for the DTC eigenstate. The information lattice is computed for all eigenstate seperately for a fixed disorder realization. The info lattice is saved for all eigenstates seperately for each epsilon value."

	# Author

	attrs["Author"] = "Sagnik Ghosh"

	# cluster

	attrs["Cluster"] = "BAF"

	# Parameters

	attrs["[Parameters] L"] = L
	
	attrs["[Parameters] theta"] = theta

	attrs["[Parameters] epsilon"] = string(epsilonlist)

	attrs["[Method] Eigen"] = "ED with Arpack"

######################################
# Code:
######################################	


for epsilon in epsilonlist

U=brickwall(L,theta, epsilon)
eigvals,eigvecs=eigen(U);


file["L$L/theta$theta/epsilon$epsilon/eigvals"]=eigvals

for i in 1:2^L
	
Psi=eigvecs[:,i]
info_lattice_t=info_lattice(Psi,L)

file["L$L/theta=$theta/epsilon$epsilon/eigvec/$i"]=info_lattice_t

end


end

end_time=Dates.now()

elapsed_time = end_time - start_time
attrs["[Benchmark] Elapsed Time"] = string(format_duration(elapsed_time))

close(file)


EOF
done 









