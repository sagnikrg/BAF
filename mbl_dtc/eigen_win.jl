#########################################################################
# Computing Eigenstatistics of the MBL-DTC Unitary (random h)
# window points of the epsilon scan (thesis wishlist V1, option B)
#
#   julia eigen_win.jl L theta itr job eps
#   julia eigen_win.jl L theta itr job eps [set]
#     eps: comma separated epsilon list (epsilon = 1-g), from eps/theta_<theta>.txt or eps/grid.txt
#     set: win (default, the window points), grid (the old 0.025 grid) or scr (scrambling reference), only names the file
#
# Same physics and file layout as the production random h runs of 2026-09
# (mbldtc_L8_theta_*.hdf5 on E:), with two changes:
#   - epsilon groups keep the full value ("epsilon0.1125"), not the first 5 characters,
#     which would cut 0.1125 to "0.112"
#   - the stored gap wraps around the circle (min(d, 2pi-d)); the old one did not
#########################################################################

using Pkg					#   For Package Management
using Dates                 #   For Date and Time
using LinearAlgebra         #   For Linear Algebra
using Random                #   For Random Number Generation
using Kronecker             #   For Kronecker Product
using Arpack                #   For Eigenvalue Decomposition
using StatsBase             #   For Statistics
using HDF5                  #   For Saving Data


# The header is loaded without its plotting dependency: functions/histgram.jl has `using CairoMakie`
# (back since the header syncs of 2026-10-29/11-05), and CairoMakie is not in the julia 1.9.4 depot of
# the cluster. So `using CairoMakie` is dropped, and the pure plotting files (*_plots.jl, which use the
# L"..." strings of Makie) are not included; the other plot functions are defined but never called here.
function nomakie(ex)
    if ex isa Expr && ex.head == :using && any(a -> a isa Expr && a.args == [:CairoMakie], ex.args)
        return nothing
    elseif ex isa Expr && ex.head == :call && ex.args[1] == :include && length(ex.args) == 2 &&
           ex.args[2] isa String && endswith(ex.args[2], "_plots.jl")
        return nothing
    elseif ex isa Expr && ex.head == :call && ex.args[1] == :include && length(ex.args) == 2
        return Expr(:call, :include, :nomakie, ex.args[2])
    end
    return ex
end
include(nomakie, ".header/Header.jl")

# The two argument LazadiresDiagram of the production runs (header of 2026-08-31, commit 2ee4e53), from the
# eigenpairs already computed; the header since the sync of 2026-10-29 only has LazadiresDiagram(U), which
# would diagonalise U a second time. phase_ordered_eigvecs(EigA,Eigvec) is still in the header.
function LazadiresDiagram(EigA, Eigvec)
    EigvecNew=phase_ordered_eigvecs(EigA,Eigvec)[2];
    localdim=length(eigvals(Z));
    dim=convert(Int64,floor(log(size(EigvecNew)[1])/log(localdim)))
    Xi=copy(kron(Z,kron_power(I(localdim),(dim-1))));
    return real.(conj(transpose(EigvecNew))*Xi*EigvecNew)
end

L=parse(Int, ARGS[1]);
thetastr=ARGS[2];
theta=parse(Float64, thetastr);
Itrnumb=parse(Int, ARGS[3]);
job=parse(Int, ARGS[4]);
epsilonstr=split(ARGS[5], ",");
epsilonlist=parse.(Float64, epsilonstr);
set= length(ARGS) >= 6 ? ARGS[6] : "win";
Ntot=2^L;

global levelspacing=fill(0.0,length(epsilonlist))
global entanglement_ee=fill(0.0,length(epsilonlist))


file=h5open("mbldtc_$(set)_L$(L)_theta_$(thetastr)_$(job).hdf5","cw")
attrs=attributes(file)


######################################
# Attributes:
######################################

	attrs["Date/Time"]=string(Dates.now())
	attrs["[ENV] Processor Type"] = string(Sys.CPU_NAME)
	attrs["[ENV] Julia Version"] = string(VERSION)
	attrs["[ENV] Container OS"] = "Debian12"
	attrs["[ENV] Number of Threads"] = string(Threads.nthreads())
	attrs["[ENV] Code"] = read(@__FILE__, String)
	attrs["[ENV] Modules"] = string(Pkg.installed())
	attrs["[ENV] Julia Environment"] = "julia-1.9.4-08-08-24.tar.gz"
	attrs["METADATA"] = (occursin("scr",set) ? "Scrambling reference point theta=1, g=0.5 (thesis wishlist V1.3)" : occursin("grid",set) ? "Old 0.025 grid of the epsilon scan, more realisations (thesis wishlist V1.2)" : "Window points of the epsilon scan (thesis wishlist V1)")*" of the MBL-DTC Unitary for L=$(L) and theta=$(theta), random h. The data contains the eigenvalues and antipodal gaps for each realisation, disorder averaged level spacing ratio, disorder averaged eigenstate entanglement entropy and disorder averaged histograms of the diagonal, pi-diagonal and various offdiagonal elements of the Lazarides diagram."
	attrs["Author"] = "Sagnik Ghosh"
	attrs["Cluster"] = "BAF"
	attrs["[Parameters] h"] = "2pi"
	attrs["[Parameters] J"] = "pi"
	attrs["[Parameters] L"] = L
	attrs["[Parameters] theta"] = theta
	attrs["[Parameters] epsilon"] = string(epsilonlist)
	attrs["[Parameters] Itrnumb"] = Itrnumb
	attrs["[Parameters] job"] = job


######################################
# Code:
######################################

time_taken= @elapsed begin

for i in 1:length(epsilonlist)

	epsilon=epsilonlist[i]
	eg="L$(L)/theta$(theta)/epsilon"*epsilonstr[i]

	global histogram_diag=fill(0,2000)
	global histogram_pi_diag=fill(0,2000)

	global histogram_offdiag=fill(0,16,2000)
	global histogram_pi_offdiag=fill(0,16,2000)


	for itr in 1:Itrnumb

		#########################################################################
		# Brickwall
		#########################################################################

        h=rand(L)*2*pi ;
        J=rand(L-1)*pi;

        A=circuit_dtc(L, theta, epsilon, h, J);


		#########################################################################
		# Eigenstatistics
		#########################################################################

		eigA,eigvecA=eigen(A)

		# antipodal gap, distance taken around the circle
		omega=sort(angle.(eigA).+pi)
		gap=Float64[]
		for omega_i in omega
			pair=mod(pi+omega_i, 2pi)
			d=abs.(omega.-pair)
			push!(gap, minimum(min.(d, 2pi.-d)))
		end

		file[eg*"/Itr$(itr)/eig"]=eigA;
		file[eg*"/Itr$(itr)/gap"]=gap;

        attribs = HDF5.attributes(file[eg*"/Itr$(itr)"])
        attribs["J"]=J
        attribs["h"]=h

		levelspacing[i]=levelspacing[i]+MeanLevelSpacingRatio(eigA)


		#########################################################################
		# Eigenstate Entanglement Entropy
		#########################################################################

		for j in 1:Ntot
			entanglement_ee[i]=entanglement_ee[i]+EntanglementEntropy(eigvecA[:,j], L)
		end


		#########################################################################
		# Lazarides-Luitz Statistics
		#########################################################################

		Corr=LazadiresDiagram(eigA,eigvecA);

		histogram_diag .+=histgram_interval(diag(Corr), -1:0.001:1)
		histogram_pi_diag .+=histgram_interval(pi_diag(Corr), -1:0.001:1)

		for j in 1:16
			histogram_offdiag[j,:] .+=histgram_interval(offdiag(Corr,j), -1:0.001:1)
			histogram_pi_offdiag[j,:] .+=histgram_interval(pi_offdiag(Corr,j), -1:0.001:1)
		end

		## Itr loop ends here
	end

	histogram_diag=histogram_diag/Itrnumb
	histogram_pi_diag=histogram_pi_diag/Itrnumb
	histogram_offdiag=histogram_offdiag/Itrnumb
	histogram_pi_offdiag=histogram_pi_offdiag/Itrnumb

	file[eg*"/Histogram/Diag"]=histogram_diag
	file[eg*"/Histogram/PiDiag"]=histogram_pi_diag
	file[eg*"/Histogram/OffDiag"]=histogram_offdiag
	file[eg*"/Histogram/PiOffDiag"]=histogram_pi_offdiag

	## epsilon loop ends here
end
end

	levelspacing=levelspacing/Itrnumb
	entanglement_ee=entanglement_ee/(Itrnumb*Ntot)

	file["L$(L)/theta$(theta)/Levelspacing"]=levelspacing;
	file["L$(L)/theta$(theta)/EntanglementEE"]=entanglement_ee;

	attrs["Time Elapsed"]=time_taken
close(file)
