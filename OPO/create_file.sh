

cd $BUDDY

fp=0.1

#cd $BUDDY

Itrnumb=1

for i1 in {1..Itrnumb}
do



cat <<EOF >opo_fp_${fp}_${i1}.jl


#########################################################################
# Computing Eigenstatistics of the MBL-DTC Unitary
#########################################################################

# Importing Headers

using Pkg					#   For Package Management
using Dates                 #   For Date and Time
using HDF5                  #   For Saving Data


include("header/Headers.jl")



# Parameters

g = 0.00118         		# Nonlinear interaction strength
kappa = 0.045    			# Decay coefficient
N = 512         			# Number of spatial points
dx = 0.6799         		# Spatial step size
L = dx*N        			# Length of the domain
k=1;

dt = 0.001       			# Time step size
T = 10.0         			# Total time

noise_strength = sqrt(kappa/dx)*sqrt(dt)  
							# Strength of the noise

f_0=\$(fp) 					# Pumping strength
k_p=1.4						# Pumping wave number
omega_p=-0.42;				# Pumping frequency




file=h5open("OPO_fp_\$(fp)_${i1}.hdf5","cw")
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

	attrs["METADATA"] = "This Data file contains a full parameter scan of epsilon in epsilonlint of the MBL-DTC Unitary for L=\$(L) and theta=\$(theta). The data contains the eigenvalues for each realisation, disorder averaged level spacing ratio, disorder averaged eigestate entanglement entropy and disorder averaged histograms diagonal, pi-diagonal and various offdiagonal elements of the Lazadires Matrix."

	# Author

	attrs["Author"] = "Sagnik Ghosh"

	# cluster

	attrs["Cluster"] = "BAF"

	# Parameters

	attrs["[Parameters] h"] = "0"
	attrs["[Parameters] L"] = L
	attrs["[Parameters] theta"] = theta
	attrs["[Parameters] epsilon"] = string(epsilonlist)

	attrs["[Parameters] Itrnumb"] = Itrnumb


######################################
# Code:
######################################	




for i in 1:length(epsilonlist)

	epsilon=epsilonlist[i]


	global histogram_diag=fill(0,2000)
	global histogram_pi_diag=fill(0,2000)

	global histogram_offdiag=fill(0,16,2000)
	global histogram_pi_offdiag=fill(0,16,2000)


	for itr in 1:Itrnumb


		#########################################################################
		# Brickwall
		#########################################################################



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

		for j in 1:Ntot
	
			#computing the half chain entanglement entropy
	
			entanglement_ee[i]=entanglement_ee[i]+EntanglementEntropy(eigvecA[:,j], L)
	
		end


		#########################################################################
		# Lazarides-Luitz Staistics 
		#########################################################################

		Corr=LazadiresDiagram(eigA,eigvecA);

		histogram_diag .+=histgram(diag(Corr), -1:0.001:1)
		histogram_pi_diag .+=histgram(pi_diag(Corr), -1:0.001:1)

		for j in 1:16
				histogram_offdiag[j,:] .+=histgram(offdiag(Corr,j), -1:0.001:1)
				histogram_pi_offdiag[j,:] .+=histgram(pi_offdiag(Corr,j), -1:0.001:1)
		end


		## Itr loop ends here
	end

	# normalising the Lazarides-Luitz statistics

	histogram_diag=histogram_diag/Itrnumb
	histogram_pi_diag=histogram_pi_diag/Itrnumb

	histogram_offdiag=histogram_offdiag/Itrnumb
	histogram_pi_offdiag=histogram_pi_offdiag/Itrnumb

	#saving the Lazarides-Luitz statistics

	file["L\$(L)/theta\$(theta)/epsilon"*first("\$(epsilon)",5)*"/Histogram/Diag"]=histogram_diag
	file["L\$(L)/theta\$(theta)/epsilon"*first("\$(epsilon)",5)*"/Histogram/PiDiag"]=histogram_pi_diag

	file["L\$(L)/theta\$(theta)/epsilon"*first("\$(epsilon)",5)*"/Histogram/OffDiag"]=histogram_offdiag
	file["L\$(L)/theta\$(theta)/epsilon"*first("\$(epsilon)",5)*"/Histogram/PiOffDiag"]=histogram_pi_offdiag


	## epsilon loop ends here
end

levelspacing=levelspacing/Itrnumb
entanglement_ee=entanglement_ee/(Itrnumb*Ntot)

file["L\$(L)/theta\$(theta)/Levelspacing"]=levelspacing;
file["L\$(L)/theta\$(theta)/EntanglementEE"]=entanglement_ee;
	
close(file)	
	






EOF

done 









