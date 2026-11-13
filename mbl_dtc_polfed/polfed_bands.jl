#########################################################################
# V2: two antipodal bands of eigenpairs of the MBL-DTC unitary at large L, from normal POLFED
#
#   julia polfed_bands.jl L theta r job eps [nev] [set]
#     r   : realisation number; the disorder (h, J and the noise) is drawn from the seed 100000 L + r, so a
#           realisation is the same circuit at every epsilon and every theta (h, J are the first draws)
#     job : number of the job within the realisation (an L=16 realisation is split over epsilon)
#     eps : comma separated epsilon list (epsilon = 1-g)
#     nev : eigenpairs requested per band (default 64)
#     set : polfed (default) or test, only names the file
#
# For every epsilon: POLFED (ARPACK on the geometric sum g_k(U), Luitz 2021, parameters luitz_params) around the
# targets phi=0 and phi=pi, from the same random start. Stored per band: the converged quasienergies omega (complex),
# their residuals, and the half chain entanglement entropy of each eigenvector. The antipodal gap and <r> are taken
# from these in the analysis (only the levels whose antipode lies inside the other band count).
#
# Writes mbldtc_<set>_L<L>_theta_<theta>_r<r>_<job>.hdf5
#########################################################################

using Pkg, Dates, LinearAlgebra, Random, HDF5
BLAS.set_num_threads(1)
include(".header/.src/itensor/ITensor.jl")      # circuit_dtc(sites, theta, epsilon), brickwall_tev, random_vector
include(".header/.src/polfed/polfed.jl")        # FloquetOperator, polfed, floquet_matrix

# the model as in polfed/code/polfed.bench.jl
function model(L, theta, epsilon, seed)
    Random.seed!(seed)
    sites=siteinds("S=1/2", L)
    gE, gO, J, h = circuit_dtc(sites, theta, epsilon)
    return sites, gE, gO, J, h
end

# normal POLFED parameters of Luitz 2021, as in polfed.bench.jl
luitz_params(L)=(ncv=2^(L÷2+1), k=max(1, round(Int, 0.8*2^(L+1)/2^(L÷2+1))))

# half chain entanglement entropy of a dense eigenvector (an L legged ITensor)
function half_chain_ee(psi::ITensor, sites)
    L=length(sites)
    v=vec(Array(psi, sites...))
    p=svdvals(reshape(v, 2^(L÷2), 2^(L-L÷2))).^2
    p=p[p .> 1e-14]
    return -sum(p .* log.(p))
end

L=parse(Int, ARGS[1])
thetastr=ARGS[2]; theta=parse(Float64, thetastr)
r=parse(Int, ARGS[3])
job=parse(Int, ARGS[4])
epsilonstr=split(ARGS[5], ",")
epsilonlist=parse.(Float64, epsilonstr)
nev= length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 64
set= length(ARGS) >= 7 ? ARGS[7] : "polfed"
seed=100000*L + r
p=luitz_params(L)
nev=min(nev, p.ncv÷2)          # ARPACK needs ncv > nev; only small test sizes are capped (L=14: ncv=256, L=16: 512)

file=h5open("mbldtc_$(set)_L$(L)_theta_$(thetastr)_r$(r)_$(job).hdf5", "cw")
attrs=attributes(file)
attrs["Date/Time"]=string(Dates.now())
attrs["[ENV] Processor Type"]=string(Sys.CPU_NAME)
attrs["[ENV] Julia Version"]=string(VERSION)
attrs["[ENV] Container OS"]="Debian12"
attrs["[ENV] Code"]=read(@__FILE__, String)
attrs["[ENV] Modules"]=string(Pkg.dependencies() |> values |> collect |> x -> [(d.name, string(d.version)) for d in x if d.is_direct_dep])
attrs["METADATA"]="Two antipodal bands (targets phi=0 and phi=pi) of eigenpairs of the MBL-DTC unitary for L=$(L), theta=$(theta), random h, realisation $(r) (seed $(seed), the same circuit at every epsilon), from normal POLFED (ARPACK, Luitz 2021 parameters ncv=$(p.ncv), k=$(p.k)), nev=$(nev) per band. Per epsilon and band: omega (re, im) of the converged pairs, residual, half chain entanglement entropy S of each eigenvector, the Floquet steps and the time taken."
attrs["Author"]="Sagnik Ghosh"
attrs["Cluster"]="BAF"
attrs["[Parameters] h"]="2pi"
attrs["[Parameters] J"]="pi"
attrs["[Parameters] L"]=L
attrs["[Parameters] theta"]=theta
attrs["[Parameters] epsilon"]=string(epsilonlist)
attrs["[Parameters] realisation"]=r
attrs["[Parameters] seed"]=seed
attrs["[Parameters] job"]=job
attrs["[Parameters] nev"]=nev
attrs["[Parameters] ncv"]=p.ncv
attrs["[Parameters] k"]=p.k

time_taken=@elapsed for (i, epsilon) in enumerate(epsilonlist)
    eg="L$(L)/theta$(theta)/epsilon"*epsilonstr[i]
    haskey(file, eg) && continue
    sites, gE, gO, J, h = model(L, theta, epsilon, seed)
    U=FloquetOperator(gE, gO)
    for (band, phi) in (("band0", 0.0), ("bandpi", float(pi)))
        Random.seed!(seed + 1)
        t=@elapsed res=polfed(U, random_vector(sites); k=p.k, ncv=p.ncv, nev=nev, phi=phi, method=:arpack)
        S=[half_chain_ee(x, sites) for x in res.vectors]
        file[eg*"/$(band)/omega"]=hcat(real.(res.omega), imag.(res.omega))
        file[eg*"/$(band)/residual"]=collect(res.residuals)
        file[eg*"/$(band)/S"]=S
        a=attributes(file[eg*"/$(band)"])
        a["phi"]=phi; a["nsteps"]=res.nsteps; a["time"]=t
    end
    a=attributes(file[eg]); a["J"]=J; a["h"]=h
    flush(file)
    println("epsilon=$(epsilon) done"); flush(stdout)
end

attrs["Time Elapsed"]=time_taken
close(file)
println("maxrss_MB ", round(Sys.maxrss()/2^20, digits=1))   # peak memory of this job (condor samples MemoryUsage too rarely for short jobs)
