

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
tnum=800;
tnum2=40000;
Ntot=2^L;


epsilonlist=[0.0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.44,0.5,0.56,0.64,0.66,0.68,0.7,0.72,0.74,0.75,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1.0]






file=h5open("mbldtc-pi-tev-randombitstring_L\$(L)_theta_\$(theta)_${i1}.hdf5","cw")
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

	attrs["METADATA"] = "This Data file contains time evolution of half chain entanglement entropy and ZZ correlation for the combination of epsilon in epsilonlint of the MBL-DTC Unitary for L=\$(L) and theta=\$(theta)."
	
	
	# Author

	attrs["Author"] = "Sagnik Ghosh"

	# cluster

	attrs["Cluster"] = "BAF"

	# Parameters

	attrs["[Parameters] h"] = "2pi"
	attrs["[Parameters] L"] = L
	attrs["[Parameters] theta"] = theta
	attrs["[Parameters] epsilon"] = string(epsilonlist)
	attrs["[Parameters] Initial State"] = "randomMPS"

	attrs["[Parameters] Itrnumb"] = Itrnumb


######################################
# Code:
######################################	




for i in 1:length(epsilonlist)

	epsilon=epsilonlist[i]


	global halfchainee_t=fill(0.0, tnum)
	global zz_t=fill(0.0, tnum2,L)
	global zz_t_mean=fill(0.0, tnum2)

	for itr in 1:Itrnumb


		#########################################################################
		# Initial State
		#########################################################################


		sites = siteinds("S=1/2", L);
		dummysites = siteinds("S=1/2", L);

		# Initial State

		#psi = productMPS(sites,"↑"  )
		#psi = Neel_state(sites  )
		#psi = randomMPS(sites  )
		psi = random_BitString(sites  )

		Psi=psi[1]*psi[2]
		for i in 3:L
    		Psi=Psi*psi[i]
		end

		Psi=Psi/norm(Psi)
		norm(Psi)

		#########################################################################
		# Z gates
		#########################################################################
		
		gateZ=Array{ITensor}(undef, L)

		for k in 1:L

			gateZ[k]=ITensor(Z, dummysites[k],sites[k])

		end

		PsiZ=Array{ITensor}(undef, L)

		for k in 1:L

		PsiZ[k]=copy(Psi)
		PsiZ[k]=gateZ[k]*PsiZ[k]
		PsiZ[k]=PsiZ[k]*delta(sites[k],dummysites[k])

		end
    

		#########################################################################
		# Brickwall
		#########################################################################



		brick1=brickwall_tensor(L, theta, epsilon)
		brick=itensorise(brick1, sites, dummysites)


		#########################################################################
		# Time Evolution
		#########################################################################

		Psitemp=Array{ITensor}(undef, L)

		for t in 1:tnum
   
			halfchainee_t[t]=halfchainee(Psi,sites)

			for ttemp in 1:100

				Psi=brickwall_tev(Psi, brick, sites, dummysites)
				t_temp=(t-1)*100+ttemp
			

				if t_temp<=tnum2

				for k in 1:L


					PsiZ[k]=brickwall_tev(PsiZ[k], brick, sites, dummysites)

					Psitemp[k]=gateZ[k]*Psi
					Psitemp[k]=Psitemp[k]*delta(sites[k],dummysites[k])

					zz=inner(PsiZ[k],Psitemp[k])
					zz_t[t_temp,k]=real(zz)

				end

				end
			end
		
		end

		for t in 1:tnum2
   		 	zz_t_mean[t]=zz_t_mean[t]+mean(zz_t[t,:])
		end


		
		## Itr loop ends here
	end

	#########################################################################
	# Normalising the ZZ correlation, half chain entanglement entropy
	#########################################################################

	zz_t_mean=zz_t_mean/Itrnumb
	halfchainee_t=halfchainee_t/Itrnumb


	#########################################################################
	# Saving Data
	#########################################################################

	file["L\$(L)/theta\$(theta)/epsilon\$(epsilon)/HalfChainEE"]=halfchainee_t;

	file["L\$(L)/theta\$(theta)/epsilon\$(epsilon)/ZZ"]=zz_t_mean;

	## epsilon loop ends here
end

	
close(file)	
	






EOF

done 









