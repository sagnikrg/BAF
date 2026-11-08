#################################################################
# Production run: theta scans at fixed g, directly with the unitaries
#
#   julia -t 1 run_vertical.jl <outdir> <L> <g> <r0> <r1> [thetas]
#
#   realisations r0..r1, thetas comma separated (default THETAS below)
#   e.g. julia -t 1 run_vertical.jl data 12 0.5 1 1
#
# One HDF5 file per block: vertical_L<L>_g<g>_r<r0>-<r1>.hdf5, group theta<thetamean>/real<r>.
# The disorder (h, J, leg noise, theta noise) of realisation r is the same for all theta and
# is seeded by (L, g, r), and U_0 is the same circuit with theta = 0 exactly.
# Finished realisations are skipped, so a block can be restarted.
# The night run of 2026-11-07 used the older one file per (L, g) layout, vertical_L<L>_g<g>.hdf5;
# plot_vertical.jl reads both, but never mix the two in one data folder (same seeds).
#################################################################

include(joinpath(@__DIR__, "vertical.jl"))
include(joinpath(@__DIR__, "unitary_analysis.jl"))
using HDF5, Printf, Dates

# fine at small theta, where the peak at delta = 0 dies, up to the random drive at theta = 1
const THETAS = [0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.5, 0.75, 1.0]

function run_block(outdir, L, g, r0, r1, thetas)
    fname = joinpath(outdir, @sprintf("vertical_L%d_g%.2f_r%d-%d.hdf5", L, g, r0, r1))
    h5open(fname, "cw") do f
        if !haskey(HDF5.attributes(f), "L")
            a = HDF5.attributes(f)
            a["L"] = L; a["g"] = g; a["thetas"] = thetas
            a["logedges"] = collect(logedges())
            a["Date/Time"] = string(now()); a["[ENV] Processor Type"] = string(Sys.CPU_NAME)
            a["[ENV] Julia Version"] = string(VERSION); a["[ENV] BLAS threads"] = BLAS.get_num_threads()
        end
    end
    for r in r0:r1
        done = h5open(f -> all(haskey(f, "theta$(th)/real$(r)") for th in thetas), fname, "r")
        done && continue
        t0 = time()
        rng = MersenneTwister(hash((L, g, r)))
        p0 = draw_params(rng, L, 0.0, g)          # theta = noise only
        dth = p0.theta
        U0 = floquet_unitary(p0; theta = 0.0)
        phi, Q0, fell0 = unitary_eigen(U0)
        lsr0 = lsr(phi)
        U0 = nothing; GC.gc()
        for th in thetas
            key = "theta$(th)/real$(r)"
            U = floquet_unitary(p0; theta = th + dth)
            a = analyse_pair(phi, Q0, U)
            U = nothing; GC.gc()
            h5open(fname, "cw") do f
                haskey(f, key) && delete_object(f, key)
                grp = create_group(f, key)
                for k in (:s1, :s2, :sl, :c1, :m2, :l2, :c2, :hd, :hlog)
                    grp[String(k)] = getfield(a, k)
                end
                at = HDF5.attributes(grp)
                for k in (:diagw, :Shalf, :PRn, :Pself, :lsr)
                    at[String(k)] = getfield(a, k)
                end
                at["fell_back"] = Int(a.fell || fell0)
                at["theta"] = th + dth
                at["lsr_U0"] = lsr0
            end
        end
        @printf("%s L=%d g=%.2f real %d (%d-%d)  %.1f s\n", Dates.format(now(), "HH:MM:SS"), L, g, r, r0, r1, time() - t0)
        flush(stdout)
    end
end

outdir = ARGS[1]
mkpath(outdir)
L, g, r0, r1 = parse(Int, ARGS[2]), parse(Float64, ARGS[3]), parse(Int, ARGS[4]), parse(Int, ARGS[5])
thetas = length(ARGS) >= 6 ? parse.(Float64, split(ARGS[6], ",")) : THETAS
run_block(outdir, L, g, r0, r1, thetas)
println("ALL DONE")
