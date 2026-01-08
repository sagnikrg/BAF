

cd $BUDDY

L=8
Itrnumb=1

for itr in {1..Itrnumb}
do



cat <<EOF >eigenscatterL${L}_${itr}.jl

###################
# Headers
###################

include(".header/Headers.jl");

###################
# Parameters
###################
L=$(L)

WList=[ 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 , 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19 , 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4,  1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.63, 1.64, 1.65, 1.66, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87, 1.88, 1.89, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0,31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0];

γ=1.0
Itrnumber=1000


file_destination= h5open("anderson_eigendata$(L)_$(itr).hdf5","cw");

attrs=HDF5.attributes(file_destination)

######################################
# Attributes:
######################################

	# Extracting Date and Time

	#Dates.now()
	attrs["Date/Time"]=string(Dates.now())


	
	# Extracting Processor Type

    
	model_name = read_model_name("model_name.txt")
    attrs["[Benchmark] Processor Type"] = string(model_name)


	# Extracting Julia Version

	julia_version = VERSION
	attrs["[ENV] Julia Version"] = string(julia_version)

	# Extracting Number of Threads

	num_threads = Threads.nthreads()
	attrs["[ENV] Number of Threads"] = string(num_threads)	

	# Code

	script_content = read("eigenscatterL$(L).jl", String)
	attrs["[ENV] Code"] = script_content

    #Lindblad
    script_content = read(".header/.function/lindblad.jl", String)
    attrs["[ENV] Lindbladian Function"] = script_content
	# Modules

	module_list = Pkg.installed()
	attrs["[ENV] Modules"] = string(module_list)

    # Julia Environment

	attrs["[ENV] Julia Environment"] = "julia-1.9.4-08-08-24.tar.gz" 

	# Meta Data

	attrs["METADATA"] = "This file contains the raw eigenvalues of the Lindbladian for an 1D Anderson model, with open boundary conditions, and a Z channel dissipation on site 1. The data is organized by system size (L), disorder strength (W), and iteration number (itr). Each entry includes the eigenvalues, norms of the left and right eigenvectors, the specific disorder realization used, and the time taken to compute these values.
    
    
    We also compute the IPR on the fly of each eigenstates, and save the following eigenstates: ground state, state with maximal IPR, a state in the lower and upper band at L and 2L distance from the argmax IPR respectively"
	
	# Author

	attrs["Author"] = "Sagnik Ghosh"


	# cluster

	attrs["Cluster"] = "BAF"


	attrs["[Parameters] Itrnumb"] = \$(Itrnumber)
    attrs["[Parameters] L"] = \$(L)
    attrs["[Parameters] WList"] = string(WList)

    attrs["[Parameters] γ"] = string(γ)

close(file_destination)
###################
# Main Loop
###################



    for itr in 1:Itrnumber

       # println("L=",L," itr=",itr)
          
        
        file=h5open("anderson_eigendata$(L)_$(itr).hdf5","cw");
        for W in WList

            Time_begin_eigen=Dates.now()
	        
                #The Lindbladian:
                
                    Lind, μ = Lindbladian(W,L,γ);


                # Compute eigenvalues and right eigenvectors
                    eigenvalues, right_eigenvectors = eigen(Matrix(Lind));          
            
            Time_end_eigen=Dates.now()
            total_time_eigen=Time_end_eigen-Time_begin_eigen

                # Writing Benchmark time
                    file["L\$L/W\$(W)/itr\$(itr)/eigen_benchmark_time"] = string(format_duration(total_time_eigen)); 
              


            file["L\$L/W\$(W)/itr\$(itr)/eigenvalues"] = eigenvalues;
            file["L\$L/W\$(W)/itr\$(itr)/disorder_realisation"] = μ;

            #compute the left eigenvectors by looking at Lindblad complex conjugate transpose

            left_eigenvectors = right_eigenvectors';
            
            # Save results to HDF5 file

            file["L\$L/W\$(W)/itr\$(itr)/left_eigenvectors_norm"] = norm.(eachcol(left_eigenvectors)); 
            file["L\$L/W\$(W)/itr\$(itr)/right_eigenvectors_norm"] = norm.(eachcol(right_eigenvectors));


        #IPR in computational basis:

            global IPR_rightvec=fill(0.0, L^2)
            
            for j in 1:L^2

                IPR_rightvec[j]=IPR(right_eigenvectorsd[:,j])

            end


   

        # Writing IPR and Densities of some states

         file["L\$L/W\$(W)/itr\$(itr)/IPR"] = IPR_rightvec
         file["L\$L/W\$(W)/itr\$(itr)/Eigstate_ground"] = right_eigenvectors[:,end]
         file["L\$L/W\$(W)/itr\$(itr)/Eigstate_IPR_argmax"] = right_eigenvectors[:,argmax(IPR_rightvec)]


        end
        close(file)
        GC.gc()
    end


EOF
done 

