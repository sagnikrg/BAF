#########################################################################
# Computing Eigenstatistics of the MBL-DTC Unitary (random h), symmetry resolved
# theta scans on the two symmetric lines g = 0 (epsilon = 1) and g = 1 (epsilon = 0)
#
#   julia eigen_sym.jl L epsilon itr job thetas [outdir]
#     epsilon: 1.0 (g = 0, no kick) or 0.0 (g = 1, perfect kick)
#     thetas:  comma separated theta list (Vision's fine grid, Tony 2026-11-13)
#
# Same physics and file layout as eigen_win.jl (circuit_dtc with the noise on, full eigenvalues, antipodal gaps,
# level spacing ratio, eigenstate entanglement entropy, Lazarides histograms; no info lattice), with the scan over
# theta at one epsilon instead of over epsilon at one theta, and in addition the symmetry resolved spectrum:
#   epsilon = 1 (g = 0): "ztot"          [S^z, U] = 0, the L+1 sectors of m           (Sec. sec:symmetries)
#   epsilon = 0 (g = 1): "perfect_kick"  U maps m -> -m: U^2 inside m > 0 and U inside m = 0, each spectrum once
# The same h, J are used for all theta of one realisation (the noise of circuit_dtc is drawn at every call).
# Full spectrum by the Hermitian trick (residual checked, falls back to eigen(U)).
#########################################################################

using Pkg					#   For Package Management
using Dates                 #   For Date and Time
using LinearAlgebra         #   For Linear Algebra
using Random                #   For Random Number Generation
using StatsBase             #   For Statistics
using HDF5                  #   For Saving Data

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
include(nomakie, joinpath(@__DIR__, ".header/Header.jl"))
BLAS.set_num_threads(1)

function LazadiresDiagram(EigA, Eigvec)
    EigvecNew=phase_ordered_eigvecs(EigA,Eigvec)[2];
    localdim=length(eigvals(Z));
    dim=convert(Int64,floor(log(size(EigvecNew)[1])/log(localdim)))
    Xi=copy(kron(Z,kron_power(I(localdim),(dim-1))));
    return real.(conj(transpose(EigvecNew))*Xi*EigvecNew)
end

# eigenpairs of a unitary from a Hermitian combination; residual check, eigen(U) as the fallback
function unitary_eigen(U; alpha=0.3711, tol=1e-8)
    Q=eigen(Hermitian((cis(-alpha)*U+cis(alpha)*U')/2)).vectors
    d=vec(sum(conj(Q).*(U*Q), dims=1))
    res=maximum(norm(U*Q[:,k]-d[k]*Q[:,k]) for k in axes(Q,2))
    res>tol && (F=eigen(U); return F.values, F.vectors, true)
    return d, Q, false
end

L=parse(Int, ARGS[1]);
epsilonstr=ARGS[2];
epsilon=parse(Float64, epsilonstr);
Itrnumb=parse(Int, ARGS[3]);
job=parse(Int, ARGS[4]);
thetastr=split(ARGS[5], ",");
thetalist=parse.(Float64, thetastr);
outdir= length(ARGS) >= 6 ? ARGS[6] : ".";
Ntot=2^L;

symmetry = epsilon == 1.0 ? "ztot" : epsilon == 0.0 ? "perfect_kick" : error("epsilon must be 0.0 or 1.0")
eig_resolved(U) = symmetry == "ztot" ? eig_ztot(U; vectors=false)[1] : eig_perfect_kick(U)

Random.seed!(hash((L, epsilon, job, :symres)))

global levelspacing=fill(0.0,length(thetalist))
global levelspacing_sector=fill(0.0,length(thetalist))
global entanglement_ee=fill(0.0,length(thetalist))
global fell_back=0

mkpath(outdir)
file=h5open(joinpath(outdir, "mbldtc_symres_L$(L)_epsilon_$(epsilonstr)_$(job).hdf5"),"cw")
attrs=attributes(file)


######################################
# Attributes:
######################################

	attrs["Date/Time"]=string(Dates.now())
	attrs["[ENV] Processor Type"] = string(Sys.CPU_NAME)
	attrs["[ENV] Julia Version"] = string(VERSION)
	attrs["[ENV] Container OS"] = Sys.iswindows() ? "Windows (laptop)" : "Debian12"
	attrs["[ENV] Number of Threads"] = string(Threads.nthreads())
	attrs["[ENV] Code"] = read(@__FILE__, String)
	attrs["[ENV] Modules"] = string(Pkg.installed())
	attrs["[ENV] Julia Environment"] = string(Base.active_project())
	attrs["Symmetry"] = symmetry
	attrs["METADATA"] = "Symmetry resolved theta scan on the line g=$(1-epsilon) (epsilon=$(epsilon)) of the MBL-DTC Unitary for L=$(L), random h (thesis wishlist E10, Tony 2026-11-13, Vision's fine theta grid). Resolved symmetry: $(symmetry) ("*(symmetry=="ztot" ? "[S^z,U]=0, the L+1 sectors of m; time reversal makes the reference the COE" : "perfect kick, U maps m to -m; U^2 inside m>0 and U inside m=0, each spectrum once; reference CUE")*"). The data contains, for each realisation, the full eigenvalues, the antipodal gaps and the symmetry resolved eigenvalues (eig_sector, concatenated, with sector_sizes; for perfect_kick the m>0 entries are eigenvalues of U^2), and per theta the disorder averaged level spacing ratio of the full spectrum and of the resolved spectrum (pooled over sectors), the disorder averaged eigenstate entanglement entropy and the disorder averaged histograms of the diagonal, pi-diagonal and various offdiagonal elements of the Lazarides diagram."
	attrs["Author"] = "Sagnik Ghosh"
	attrs["Cluster"] = Sys.iswindows() ? "laptop" : "BAF"
	attrs["[Parameters] h"] = "2pi"
	attrs["[Parameters] J"] = "pi"
	attrs["[Parameters] L"] = L
	attrs["[Parameters] epsilon"] = epsilon
	attrs["[Parameters] theta"] = string(thetalist)
	attrs["[Parameters] Itrnumb"] = Itrnumb
	attrs["[Parameters] job"] = job


######################################
# Code:
######################################

# half chain entropy; gesdd (divide and conquer) fails to converge on some g = 1 eigenvectors (LAPACKException(1)),
# then the slower QR iteration
function ee_safe(v, L)
    M = reshape(v, 2^(L÷2), 2^(L-L÷2))
    S = try svdvals(M) catch; svd(M; alg = LinearAlgebra.QRIteration()).S end
    p = S[S .> 1e-10].^2
    -sum(p .* log.(p))
end

time_taken= @elapsed begin

hist0()=(fill(0,2000), fill(0,2000), fill(0,16,2000), fill(0,16,2000))
H=[hist0() for _ in thetalist]

for itr in 1:Itrnumb

    h=rand(L)*2*pi ;
    J=rand(L-1)*pi;

    for i in 1:length(thetalist)

        theta=thetalist[i]
        eg="L$(L)/epsilon$(epsilonstr)/theta"*thetastr[i]

        A=circuit_dtc(L, theta, epsilon, h, J);

        #########################################################################
        # Eigenstatistics
        #########################################################################

        eigA,eigvecA,fb=unitary_eigen(A)
        global fell_back+=fb

        omega=sort(angle.(eigA).+pi)
        gap=Float64[]
        for omega_i in omega
            pair=mod(pi+omega_i, 2pi)
            d=abs.(omega.-pair)
            push!(gap, minimum(min.(d, 2pi.-d)))
        end

        Es=eig_resolved(A)

        file[eg*"/Itr$(itr)/eig"]=eigA;
        file[eg*"/Itr$(itr)/gap"]=gap;
        file[eg*"/Itr$(itr)/eig_sector"]=reduce(vcat, Es);
        file[eg*"/Itr$(itr)/sector_sizes"]=length.(Es);

        attribs = HDF5.attributes(file[eg*"/Itr$(itr)"])
        attribs["J"]=J
        attribs["h"]=h

        levelspacing[i]+=MeanLevelSpacingRatio(eigA)
        levelspacing_sector[i]+=mean(LevelSpacingRatio(Es))
        file[eg*"/Itr$(itr)/r"]=MeanLevelSpacingRatio(eigA)
        file[eg*"/Itr$(itr)/r_sector"]=mean(LevelSpacingRatio(Es))

        #########################################################################
        # Eigenstate Entanglement Entropy
        #########################################################################

        ee=0.0
        for j in 1:Ntot
            ee+=ee_safe(eigvecA[:,j], L)
        end
        entanglement_ee[i]+=ee
        file[eg*"/Itr$(itr)/ee"]=ee/Ntot

        #########################################################################
        # Lazarides-Luitz Statistics
        #########################################################################

        Corr=LazadiresDiagram(eigA,eigvecA);
        hd,hpd,hod,hpod=H[i]
        hd .+=histgram_interval(diag(Corr), -1:0.001:1)
        hpd .+=histgram_interval(pi_diag(Corr), -1:0.001:1)
        for j in 1:16
            hod[j,:] .+=histgram_interval(offdiag(Corr,j), -1:0.001:1)
            hpod[j,:] .+=histgram_interval(pi_offdiag(Corr,j), -1:0.001:1)
        end
    end
    println(Dates.format(now(), "HH:MM:SS"), " L=$L eps=$epsilonstr job $job itr $itr/$Itrnumb"); flush(stdout)
end

for i in 1:length(thetalist)
    eg="L$(L)/epsilon$(epsilonstr)/theta"*thetastr[i]
    hd,hpd,hod,hpod=H[i]
    file[eg*"/Histogram/Diag"]=hd/Itrnumb
    file[eg*"/Histogram/PiDiag"]=hpd/Itrnumb
    file[eg*"/Histogram/OffDiag"]=hod/Itrnumb
    file[eg*"/Histogram/PiOffDiag"]=hpod/Itrnumb
end
end

	levelspacing=levelspacing/Itrnumb
	levelspacing_sector=levelspacing_sector/Itrnumb
	entanglement_ee=entanglement_ee/(Itrnumb*Ntot)

	file["L$(L)/epsilon$(epsilonstr)/Levelspacing"]=levelspacing;
	file["L$(L)/epsilon$(epsilonstr)/LevelspacingSector"]=levelspacing_sector;
	file["L$(L)/epsilon$(epsilonstr)/EntanglementEE"]=entanglement_ee;

	attrs["Time Elapsed"]=time_taken
	attrs["Fell back to eigen(U)"]=fell_back
close(file)
println("ALL DONE")
