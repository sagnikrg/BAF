#####################################################################
# Chapter 7 (J.A.R.V.I.S.) on BAF: one job = realisations r0..r1 of one (kind, L, theta, g)
#
#   julia run_sb.jl <outdir> <kind> <L> <theta> <g> <r0> <r1>
#
# kind (same content as the laptop run.sb.jl, fresh seeds: the tag "baf" is in the hash):
#   full : open chain eigen -> K_obc(t), O2 numerators at the ends for nu = 1 and t/2,
#          O2 profiles from site 1 at t = 16, 64; periodic chain eigvals -> K_pbc(t)
#   scan : K_obc(t) and K_pbc(t) only
#   kb   : K_obc(t) and K_b(t), b = 2,3,4, field at site L/2 renewed every step within b
# Output: <outdir>/sb_L<L>_<kind>_theta<theta>_g<g>_r<r0>-<r1>.hdf5, same dataset names
# as the laptop files (group <kind>/theta<theta>/g<g>), plus "realisation_index".
#####################################################################

include(joinpath(@__DIR__, "sb.jl"))
using HDF5, Dates
BLAS.set_num_threads(1)

const TMAX = 256
const PROFILE_T = (16, 64)
outdir, kind = ARGS[1], Symbol(ARGS[2])
L = parse(Int, ARGS[3]); th = parse(Float64, ARGS[4]); g = parse(Float64, ARGS[5])
r0 = parse(Int, ARGS[6]); r1 = parse(Int, ARGS[7])
kind in (:full, :scan, :kb) || error("kind must be full, scan or kb")
mkpath(outdir)

# O2 profile numerators from site 1 to every y at the times PROFILE_T (as in run.sb.jl)
function o2_profile(lam, V, L)
    N = length(lam)
    zdiag(site) = [((s >> (L - site)) & 1) == 0 ? 1.0 : -1.0 for s in 0:N-1]
    Zx = V' * (zdiag(1) .* V)
    out1 = zeros(L, length(PROFILE_T)); outh = zeros(L, length(PROFILE_T))
    for y in 2:L
        Zy = V' * (zdiag(y) .* V)
        B = Zy .* transpose(Zx)
        for (k, t) in enumerate(PROFILE_T)
            lt = lam .^ t
            for (nu, out) in ((1, out1), (t ÷ 2, outh))
                c = cis(2pi * nu / t)
                s = zero(ComplexF64)
                @inbounds for m in 1:N, n in 1:N
                    q = c * lam[m] / lam[n]
                    gsum = abs(1 - q) < 1e-12 ? complex(float(t)) : (1 - lt[m] / lt[n]) / (1 - q)
                    s += B[n, m] * lt[n] * gsum
                end
                out[y, k] = abs2(s)
            end
        end
    end
    out1, outh
end

gs = string(round(g, digits=2))
name = "$(kind)/theta$(th)/g$(gs)"
file = joinpath(outdir, "sb_L$(L)_$(kind)_theta$(th)_g$(gs)_r$(r0)-$(r1).hdf5")
nr = r1 - r0 + 1
t0 = time()
Kobc = zeros(TMAX, nr); Kpbc = zeros(TMAX, nr)
n1 = zeros(TMAX, nr); nh = zeros(TMAX, nr)
p1 = zeros(L, length(PROFILE_T), nr); ph = similar(p1)
Kb = Dict(b => zeros(TMAX, nr) for b in (2, 3, 4))
for (i, ri) in enumerate(r0:r1)
    rng = MersenneTwister(hash((L, th, g, ri, kind, "baf")))
    r = Realisation(rng, L, th)
    W = floquet(r, g)
    if kind == :full
        F = eigen(W)
        a, b = o2_numerators(F.values, F.vectors, 1, L, TMAX)
        n1[:, i] = a; nh[:, i] = b
        Kobc[:, i] = sffs(F.values, TMAX)
        p1[:, :, i], ph[:, :, i] = o2_profile(F.values, F.vectors, L)
    else
        Kobc[:, i] = sffs(eigvals(W), TMAX)
    end
    if kind in (:full, :scan)
        Kpbc[:, i] = sffs(eigvals(floquet(Realisation(rng, L, th; pbc=true), g)), TMAX)
    else
        for b in (2, 3, 4); Kb[b][:, i] = sff_b(r, g, b, L ÷ 2, TMAX, rng); end
    end
    println("$name L=$L realisation $ri ($i/$nr) $(round(time() - t0, digits=1)) s"); flush(stdout)
end
h5open(file, "w") do f
    f["$name/K_obc"] = Kobc
    f["$name/realisation_index"] = collect(r0:r1)
    if kind in (:full, :scan); f["$name/K_pbc"] = Kpbc; end
    if kind == :full
        f["$name/O2num_nu1"] = n1; f["$name/O2num_half"] = nh
        f["$name/profile_nu1"] = p1; f["$name/profile_half"] = ph
        f["$name/profile_t"] = collect(PROFILE_T)
    end
    if kind == :kb
        for b in (2, 3, 4); f["$name/K_b$b"] = Kb[b]; end
        f["$name/field_site"] = L ÷ 2
    end
    attributes(f[name])["realisations"] = nr
    attributes(f[name])["seconds"] = time() - t0
end
println("done $file in $(round(time() - t0, digits=1)) s  $(now())")
