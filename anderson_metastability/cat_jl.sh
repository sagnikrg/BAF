

cd $BUDDY

L=8
Itrnumb=1

for (( itr=1; itr<=Itrnumb; itr++ ))
do



cat <<EOF >anderson_metastability_L${L}_${itr}.jl

###################
# Headers
###################



using LinearAlgebra     # for eigenvalues and eigenvectors
using StatsBase         # for statistical functions
using HDF5              # for saving data
using Dates             # for benchmark_time
using Pkg               # for package management
using SparseArrays     # for sparse matrix operations
using Random            # for random number generation




###########################---------------------------------
## Functions
###########################----------------------



# Lindbladian and Hamiltonian:


#--------------------------------------
# Hamiltonian:
#--------------------------------------

function H(W,L)
      μ = 2W * rand(L) .- W  # Random onsite potential uniformly distributed between -W and W
        
            sumμ = sum(μ)/2
            # Create the Hamiltonian matrix
            H = spzeros(L, L)
        
            # Fill the Hamiltonian matrix with disorder term (1 - 2 * c†c) and hopping term
            for i in 1:L
                H[i, i] = μ[i]-sumμ  # On-site potential term (-μ[i]/2 from (1 - 2c†c))
                if i < L
                    H[i, i+1] = 0.5  # Hopping term to the next site (scaled by 1/2 due to the transformation)
                    H[i+1, i] = 0.5  # Hopping term from the next site
                end
            end
        
            

            # Convert to dense matrix and create the Liouvillian for the Lindblad equation
            H_dense = Matrix(Hermitian(H))

            return H_dense, μ 
end



##----------------------------
# channels
##----------------------------

function Z_channel(L)
            identity_matrix = Matrix{Float64}(I, L, L)  # Identity matrix of size LxL
            n1 = spzeros(L, L)
            n1[1, 1] = 1.0                      # Corresponds to the occupation number operator c†c on site 1
            decay_channel_z =   2.0 * n1 -identity_matrix
     return decay_channel_z
end



##----------------------------
# Lindbladian
##----------------------------


function Lindbladian(W,L,γ)

        H_dense, μ =H(W,L);
        decay_channel_1=Z_channel(L)

        Hdim = size(H_dense)[1]
        identity_matrix = Matrix{Float64}(I, Hdim, Hdim)  # Corrected identity matrix creation
        
        H_transpose = transpose(H_dense)
        
        
        L_H = kron(identity_matrix, H_dense) - kron(H_transpose, identity_matrix)   # Commutator part

        
        conj_decay_channel_1 = conj(decay_channel_1)
        
        D1 = decay_channel_1' * decay_channel_1
        transpose_D1 = transpose(D1)
            
            
        L1 = kron(conj_decay_channel_1, decay_channel_1) - 0.5 * (kron(identity_matrix, D1) + kron(transpose_D1, identity_matrix))

            
        
        
        Lind = -1im * L_H + γ * (L1) # + L_X1 + L_Y1)  # Sum of all Lindblad terms

    return Lind, μ 
end



# Overloaded function with default γ=1.0


function Lindbladian(W,L)

        H_dense, μ =H(W,L);
        decay_channel_1=Z_channel(L)

        Hdim = size(H_dense)[1]
        identity_matrix = Matrix{Float64}(I, Hdim, Hdim)  # Corrected identity matrix creation
        
        H_transpose = transpose(H_dense)
        
        
        L_H = kron(identity_matrix, H_dense) - kron(H_transpose, identity_matrix)   # Commutator part

        
        conj_decay_channel_1 = conj(decay_channel_1)
        
        D1 = decay_channel_1' * decay_channel_1
        transpose_D1 = transpose(D1)
            
            
        L1 = kron(conj_decay_channel_1, decay_channel_1) - 0.5 * (kron(identity_matrix, D1) + kron(transpose_D1, identity_matrix))

        
        
        
        Lind = -1im * L_H + (L1) # + L_X1 + L_Y1)  # Sum of all Lindblad terms

    return Lind , μ 
end


# ------------------------------------------
# HDF5 Mods
# ------------------------------------------



# Function to read a single line from a file
function read_model_name(filename)
    open(filename, "r") do file
        return readline(file)
    end
end



function format_duration(duration::Millisecond)
    milliseconds=duration.value % 1000
    total_seconds = div(duration.value, 1000)  # Convert milliseconds to seconds
    hrs = div(total_seconds, 3600)             # Get hours
    mins = div(total_seconds % 3600, 60)       # Get minutes
    secs = total_seconds % 60                  # Get seconds
    return "\$(hrs) hour(s), \$(mins) minute(s), \$(secs) second(s), and \$(milliseconds) millisecond(s)"
end




#------------------------------------------------------------------------------------------

###################
# Parameters
###################
L=${L}

WList=[ 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 , 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19 , 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4,  1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.63, 1.64, 1.65, 1.66, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87, 1.88, 1.89, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0,31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0, 100.0];

γ=1.0
Itrnumber=10

mid_band_ind=2*L-1

#------------------------------------
# Output files:
#------------------------------------

file_destination_rawdata= h5open("anderson_metastability_eigendata_${L}_${itr}.hdf5","cw");
file_destination_gamma= h5open("anderson_metastability_gamma_${L}_${itr}_bootstrap.hdf5","cw");
file_destination_midscpectragap= h5open("anderson_metastability_midgspectragap_${L}_${itr}_bootstrap.hdf5","cw");



attrs=HDF5.attributes(file_destination_rawdata)

######################################
# Attributes: Raw Data File
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

	script_content = read("anderson_eigenscatter_L${L}_${itr}.jl", String)
	attrs["[ENV] Code"] = script_content

	# Modules

	module_list = Pkg.installed()
	attrs["[ENV] Modules"] = string(module_list)

    # Julia Environment

	attrs["[ENV] Julia Environment"] = "julia_1.11.1_30.10.25.tar.gz" 

	# Meta Data

	#attrs["METADATA"] = read("METADATA.txt", String)
	
	# Author

	attrs["Author"] = "Sagnik Ghosh"


	# cluster

	attrs["Cluster"] = "BAF"


	attrs["[Parameters] Itrnumb"] = Itrnumber
    attrs["[Parameters] L"] = L
    attrs["[Parameters] WList"] = string(WList)

    attrs["[Parameters] γ"] = string(γ)

close(file_destination_rawdata)
close(file_destination_gamma)
close(file_destination_midscpectragap)




###################
# Main Loop
###################



    for itr in 1:Itrnumber

       # println("L=",L," itr=",itr)

       global file_destination_rawdata=h5open("anderson_metastability_eigendata_${L}_${itr}.hdf5","cw");
       global file_destination_gamma=h5open("anderson_metastability_gamma_${L}_${itr}_bootstrap.hdf5","cw");
       global file_destination_midscpectragap=h5open("anderson_metastability_midgspectragap_${L}_${itr}_bootstrap.hdf5","cw");

        for W in WList

                Time_begin_eigen=Dates.now()
	        
                #The Lindbladian:
                
                    Lind, μ = Lindbladian(W,L,γ);


                # Compute eigenvalues and right eigenvectors
                    eigenvalues, right_eigenvectors = eigen(Matrix(Lind));          
            


                Time_end_eigen=Dates.now()
                total_time_eigen=Time_end_eigen-Time_begin_eigen

                # Writing Benchmark time
                    file_destination_rawdata["L\$L/W\$(W)/itr\$(itr)/eigen_benchmark_time"] = string(format_duration(total_time_eigen)); 
              


            file_destination_rawdata["L\$L/W\$(W)/itr\$(itr)/disorder_realisation"] = μ;
            file_destination_gamma["L\$L/W\$(W)/itr\$(itr)/disorder_realisation"] = μ;
            file_destination_midscpectragap["L\$L/W\$(W)/itr\$(itr)/disorder_realisation"] = μ;


            file_destination_rawdata["L\$L/W\$(W)/itr\$(itr)/eigenvalues"] = eigenvalues;
            


            #compute the left eigenvectors by looking at Lindblad complex conjugate transpose

            left_eigenvectors = right_eigenvectors';
            


            # Save results to HDF5 file

            file_destination_rawdata["L\$L/W\$(W)/itr\$(itr)/left_eigenvectors_norm"] = norm.(eachcol(left_eigenvectors)); 
            file_destination_rawdata["L\$L/W\$(W)/itr\$(itr)/right_eigenvectors_norm"] = norm.(eachcol(right_eigenvectors));



            # Slowest Decay Modes

                    eig_real=real.(eigenvalues)
                    eigsort=sort(eig_real)
                    file_destination_gamma["L\$L/W\$(W)/itr\$(itr)/slowest_decay_modes"] = eigsort[L^2-1];

            # midspectral gap

                    file_destination_midscpectragap["L\$L/W\$(W)/itr\$(itr)/midspectral_gap"] = (eigsort[mid_band_ind]-eigsort[mid_band_ind+1]);
                    
        end
        close(file_destination_rawdata)
        close(file_destination_gamma)
        close(file_destination_midscpectragap)

        GC.gc()
    end


EOF
done 

