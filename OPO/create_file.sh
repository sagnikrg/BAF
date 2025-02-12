

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

f_0=${fp} 					# Pumping strength
k_p=1.4						# Pumping wave number
omega_p=-0.42;				# Pumping frequency




file=h5open("OPO_fp_${fp}_${i1}.hdf5","cw")
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

	attrs["[ENV] Container OS"] = "Debian12"

	# Extracting Number of Threads

	num_threads = Threads.nthreads()
	attrs["[ENV] Number of Threads"] = string(num_threads)	

	# Code

	script_content = read("opo_fp_${fp}_${i1}.jl", String)
	attrs["[ENV] Code"] = script_content

	# Header

	header_content = read("header/Headers.jl", String)
	attrs["[ENV] Header"] = header_content

	# Modules

	module_list = Pkg.installed()
	attrs["[ENV] Modules"] = string(module_list)

	# Julia Environment

	attrs["[ENV] Julia Environment"] = "julia-1.9.4-08-08-24.tar.gz" 

	# Meta Data

	attrs["METADATA"] = "This Data file contains a RK4 simulation of the OPO model with a fixed pumping strength and a fixed noise strength for one iteration. We save the full dynamics for certain sites and the full wavefunction for certain time steps."

	# Author

	attrs["Author"] = "Sagnik Ghosh"

	# cluster

	attrs["Cluster"] = "BAF"

	# Parameters

	attrs["[Parameters] h"] = "0"
	attrs["[Parameters] fp"] = f_0
	
	attrs["[Parameters] Itrnumb"] = Itrnumb


######################################
# Code:
######################################	




	# Spatial grid
	x = range(0, L-dx, step=dx)

	# Generate complex Gaussian noise

	global u_final = fill(0.0*im, N)
	global tseries = fill(0.0*im, 10001)



	time_series_1 = Complex{Float64}[]
	time_series_2 = Complex{Float64}[]
	time_series_3 = Complex{Float64}[]
	time_series_4 = Complex{Float64}[]
	time_series_5 = Complex{Float64}[]
	time_series_6 = Complex{Float64}[]

	# Initial condition
	u = zeros(Complex{Float64}, N)

	# Time-stepping
	
	t = 0.0
	while t < T
    u = rk4_step_full(u, dt, t, g, f_0, kappa, dx).+noise_strength*randn(Complex{Float64}, N)
    t += dt
    push!(time_series_1, u[47])
    push!(time_series_2, u[126])
	push!(time_series_3, u[225])
	push!(time_series_4, u[324])
	push!(time_series_5, u[423])
	push!(time_series_6, u[512])

	# Saving the full wavefunction for certain time steps
	if t in [3.0, 3.5 ,4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
		dataset = string("psi_", Int(t))
		file[dataset] = u
    end


	file["psi_77/$(itr)"] = time_series_1
	file["psi_126/$(itr)"] = time_series_2
	file["psi_225/$(itr)"] = time_series_3
	file["psi_324/$(itr)"] = time_series_4
	file["psi_423/$(itr)"] = time_series_5
	file["psi_512/$(itr)"] = time_series_6



close(file)


done 









