#####################################################################
# Symmetry breaking in the SFF transfer matrix (Thesis Chapter 7)
#
# Exact diagonalisation side of the Garratt-Chalker programme for the
# Khemani-Ippoliti circuit: averaged SFF with open and periodic
# boundaries, the clock order parameter correlator O2 in the sectors
# nu = 1 and nu = t/2, and the response to period-b symmetry breaking
# fields K_b.
#
#   W = Kick * W_even * W_odd * Fields
#
# exactly as circuit_dtc of the header for even L and open boundaries
# (tested in test.sb.jl), extended to odd L (every site kicked once)
# and to periodic boundaries (bond (L,1) in the even layer).
#####################################################################

using LinearAlgebra, Random, Statistics

const HEADER = joinpath(@__DIR__, "..", ".header", ".src")
include(joinpath(HEADER, "gates", "gates.jl"))
include(joinpath(HEADER, "functions", "kron.jl"))

const ZZ = kron(Z, Z)
const XXYY = kron(X, X) + kron(Y, Y)
const I2 = Matrix{ComplexF64}(I, 2, 2)

#####################################################################
# Disorder realisation
#####################################################################

struct Realisation
    L::Int
    pbc::Bool
    h::Vector{Float64}       # fields, uniform [0,2pi)
    J::Vector{Float64}       # Ising even couplings, uniform [0,pi)
    theta::Float64           # thetamean + N(0,pi/50), once per circuit
    dh::Matrix{Float64}      # 4 x nbonds leg noise, N(0,pi/50)
end

nbonds(L, pbc) = pbc ? L : L - 1

# same order of random draws as circuit_dtc(L, thetamean, epsilon) of the header
function Realisation(rng::AbstractRNG, L::Int, thetamean::Float64; pbc::Bool=false)
    h = rand(rng, L) * 2pi
    J = rand(rng, L) * pi
    theta = thetamean + randn(rng, 1)[] * pi / 50
    nb = nbonds(L, pbc)
    dh = zeros(4, nb)
    for j in 1:nb
        dh[:, j] = randn(rng, 4) * pi / 50
    end
    Realisation(L, pbc, h, J, theta, dh)
end

function bond_gate(r::Realisation, j::Int)
    d = r.dh[:, j]
    int1 = kron(RZ(d[1]), RZ(d[2]))
    int2 = exp(-im * r.J[j] * ZZ - im * r.theta / 2 * XXYY)
    int3 = kron(RZ(d[3]), RZ(d[4]))
    return int3 * int2 * int1
end

# basis permutation so that kron in the site order (2,3,...,L,1) can be
# used for the periodic even layer; site 1 is the most significant bit
function cyclic_perm(L)
    N = 2^L
    p = zeros(Int, N)
    for s in 0:N-1
        b1 = (s >> (L - 1)) & 1
        sp = ((s << 1) & (N - 1)) | b1      # bits (b2,...,bL,b1)
        p[s+1] = sp + 1
    end
    p
end

"""
    floquet(r, g; hx=nothing)

Floquet operator W = Kick Even Odd Fields. `g` is the thesis g (kick angle pi g).
`hx = (x, value)` replaces the field on site x (for the period-b fields).
"""
function floquet(r::Realisation, g::Float64; hx=nothing)
    L = r.L
    h = copy(r.h)
    hx === nothing || (h[hx[1]] = hx[2])
    Fs = [RZ(h[i]) for i in 1:L]
    Ks = [RX(pi * g) for i in 1:L]
    G = [bond_gate(r, j) for j in 1:nbonds(L, r.pbc)]

    # odd layer with the fields: blocks (1,2),(3,4),...  [+ (L) for odd L]
    oddb = Matrix{ComplexF64}[G[j] * kron(Fs[j], Fs[j+1]) for j in 1:2:L-1]
    isodd(L) && push!(oddb, Fs[L])

    if !r.pbc
        # kick with the even layer: blocks (1),(2,3),...,[(L)]
        evenb = Matrix{ComplexF64}[Ks[1]]
        append!(evenb, [kron(Ks[j], Ks[j+1]) * G[j] for j in 2:2:L-1])
        iseven(L) && push!(evenb, Ks[L])
        A = kron_list(evenb, 1:length(evenb))
    else
        @assert iseven(L) "periodic boundaries need even L"
        evenb = [kron(Ks[j], Ks[mod1(j + 1, L)]) * G[j] for j in 2:2:L]
        p = cyclic_perm(L)
        A = kron_list(evenb, 1:length(evenb))[p, p]   # built in site order (2,...,L,1)
    end
    return right_kron(A, oddb)
end

# A * kron(B_1,...,B_k) without forming the product densely
function right_kron(A::Matrix{ComplexF64}, blocks)
    # (A * K) = (K^T * A^T)^T ; apply K^T = kron(B_i^T) to the columns of A^T
    M = Matrix(transpose(A))
    N = size(M, 1)
    dims = [size(b, 1) for b in blocks]
    for (i, b) in enumerate(blocks)
        d = dims[i]
        inner = prod(dims[i+1:end]; init=1)       # less significant factors
        outer = N ÷ (d * inner)
        T3 = reshape(M, inner, d, outer * N)
        P = permutedims(T3, (2, 1, 3))            # d x inner x rest
        Q = transpose(b) * reshape(P, d, :)
        M = reshape(permutedims(reshape(Q, d, inner, outer * N), (2, 1, 3)), N, N)
    end
    return Matrix(transpose(M))
end

#####################################################################
# Spectral quantities
#####################################################################

# Tr W^t for t = 1..T from the eigenvalues
function traces(lam::AbstractVector{<:Complex}, T::Int)
    out = zeros(ComplexF64, T)
    z = copy(lam)
    for t in 1:T
        out[t] = sum(z)
        z .*= lam
    end
    out
end

sffs(lam, T) = abs2.(traces(lam, T))

# K_b(t) = |Tr (W^(b)...W^(1))^(t/b)|^2 at t = b, 2b, ...; NaN elsewhere
function sff_b(r::Realisation, g::Float64, b::Int, x::Int, T::Int, rng)
    V = Matrix{ComplexF64}(I, 2^r.L, 2^r.L)
    for k in 1:b
        V = floquet(r, g; hx=(x, rand(rng) * 2pi)) * V
    end
    lam = eigvals(V)
    tr = traces(lam, T ÷ b)
    out = fill(NaN, T)
    for m in 1:T÷b
        out[m*b] = abs2(tr[m])
    end
    out
end

"""
    o2_numerators(W, x, y, T)

For t = 1..T the numerators |sum_d e^{2 pi i nu d/t} Tr[W^{t-d} Z_y W^d Z_x]|^2
of the clock correlator O2 (to be divided by the averaged SFF), for nu = 1 and
nu = t/2 (even t only, NaN for odd t). Also returns the eigenvalues.
"""
function o2_numerators(W::AbstractMatrix, x::Int, y::Int, T::Int)
    F = eigen(W)
    n1, nh = o2_numerators(F.values, F.vectors, x, y, T)
    return n1, nh, F.values
end

function o2_numerators(lam::AbstractVector, V::AbstractMatrix, x::Int, y::Int, T::Int)
    L = round(Int, log2(length(lam)))
    zdiag(site) = [((s >> (L - site)) & 1) == 0 ? 1.0 : -1.0 for s in 0:2^L-1]
    Zx = V' * (zdiag(x) .* V)
    Zy = V' * (zdiag(y) .* V)
    B = Zy .* transpose(Zx)                     # B[n,m] = <n|Zy|m><m|Zx|n>
    E = [lam[n]^s for s in 0:T, n in eachindex(lam)]   # E[s+1,n] = lam_n^s
    Fm = E * B * transpose(E)                   # Fm[s+1,d+1] = Tr[W^s Zy W^d Zx]
    n1 = zeros(T)
    nh = fill(NaN, T)
    for t in 1:T
        s1 = zero(ComplexF64); sh = zero(ComplexF64)
        for d in 0:t-1
            f = Fm[t-d+1, d+1]
            s1 += cis(2pi * d / t) * f
            sh += (iseven(d) ? 1 : -1) * f
        end
        n1[t] = abs2(s1)
        iseven(t) && (nh[t] = abs2(sh))
    end
    return n1, nh
end
