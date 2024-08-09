

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
Itrnumb=2;


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

	script_content = read("mbldtc_L${L}_${i1}.jl", String)
	attrs["[ENV] Code"] = script_content

	# Modules

	module_list = Pkg.installed()
	attrs["[ENV] Modules"] = string(module_list)


	# Meta Data

	attrs["[Parameters] L"] = L
	attrs["[Parameters] theta"] = theta
	attrs["[Parameters] epsilon"] = string(epsilonlist)

######################################
# Code:
######################################	


epsilonlist=[0.0,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.0]

levelspacing=fill(0.0,length(epsilonlist))
entanglement_ee=fill(0.0,length(epsilonlist))





for i in 1:length(epsilonlist)

	epsilon=epsilonlist[i]


	for itr in 1:\$(Itrnumb)
	
	#########################################################################
	# Brickwall
	#########################################################################


	U=brickwall(L,theta,epsilon);

	#h=rand(L)*pi/2;
	#Ind=collect(1:L)
	#ZRow=copy(kronlist(RZ.(h),Ind));




	A=brickwall(L,theta,epsilon)
	#A=A*ZRow





#########################################################################
# Eigenstatistics
#########################################################################
	eigA,eigvecA=eigen(A)

	#saving the eigenvalues
	file["L\$(L)/theta\$(theta)/epsilon"*first("\$(epsilon)",5)*"/Itr\$(itr)"]=eigA;
	
	#computing the level spacing
	levelspacing[i]=levelspacing[i]+LevelSpacingRatio(eigA)
    


#########################################################################
# Eigenstate Entanglement Entropy
#########################################################################

for j in 1:length(eigvecA)
	
	#computing the half chain entanglement entropy
	entanglement_ee[i]=entanglement_ee[i]+EntanglementEntropy(eigvecA[j],L/2)
	
end


#########################################################################
# Lazarides-Luitz Staistics 
#########################################################################

	end

end

levelspacing=levelspacing/\$(Itrnumb)


file["L\$(L)/theta\$(theta)/epsilon"*first("\$(epsilon)",4)*"/Itr"]=eigA;
	
close(file)	
	






EOF

done 









