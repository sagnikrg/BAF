#################################################################
# Vertical phases: Floquet heating in the vanilla DTC circuit
#
# The circuit is the one of the header (circuit_dtc), rebuilt here with
# all of its random parameters exposed, so that we can also write down
# the generators of every gate. The "energy" is the sum of the generators
# (the leading Magnus term, the analogue of the time averaged H_0 of
# Morningstar, Huse, Khemani, arXiv 2210.13444).
#################################################################

using LinearAlgebra, Random, Statistics

const HEADER = joinpath(@__DIR__, ".header", ".src")
include(joinpath(HEADER, "circuit", "time_crystals.jl"))   # brickwall_dtc, circuit_dtc, gates, kron

#################################################################
# Parameters of one realisation
#################################################################

struct DTCParams
    L::Int
    g::Float64                 # thesis convention g = 1 - epsilon
    theta::Float64             # theta = thetamean + noise (drawn once per circuit)
    h::Vector{Float64}         # Z field, uniform [0, 2pi)
    J::Vector{Float64}         # Ising, uniform [0, pi) (length L, the last one unused as in the header)
    delh::Matrix{Float64}      # 4 x (L-1) leg noise, legs 1,2 before and 3,4 after the gate
end

# Same order of random draws as circuit_dtc(L, thetamean, epsilon) -> brickwall_dtc
function draw_params(rng::AbstractRNG, L::Int, thetamean::Float64, g::Float64)
    h = rand(rng, L) * 2 * pi
    J = rand(rng, L) * pi
    theta = thetamean + randn(rng, 1)[] * pi / 50
    delh = zeros(4, L - 1)
    for j in 1:L-1
        delh[:, j] = randn(rng, 4) * pi / 50
    end
    DTCParams(L, g, theta, h, J, delh)
end

#################################################################
# The Floquet unitary (mirrors brickwall_dtc exactly)
#################################################################

function brickwall_gates(p::DTCParams; theta = p.theta, g = p.g)
    L = p.L
    ZZ = kron(Z, Z); XX = kron(X, X); YY = kron(Y, Y)
    FU = fill(fill(0.0 * im, 4, 4), L - 1)
    for j in 1:L-1
        int1 = kron(RZ(p.delh[1, j]), RZ(p.delh[2, j]))
        int2 = exp(-im * p.J[j] * ZZ - im * theta / 2 * (XX + YY))
        int3 = kron(RZ(p.delh[3, j]), RZ(p.delh[4, j]))
        FU[j] = int3 * int2 * int1
    end
    for i in 1:2:L-1
        FU[i] = FU[i] * kron(RZ(p.h[i]), RZ(p.h[i+1]))
    end
    gk = pi * g                               # = pi*(1-epsilon) of the header
    for i in 2:2:L-1
        FU[i] = kron(RX(gk), RX(gk)) * FU[i]
    end
    FU[1] = kron(RX(gk), I(2)) * FU[1]
    FU[end] = kron(I(2), RX(gk)) * FU[end]
    FU
end

floquet_unitary(p::DTCParams; kw...) = brickwall(brickwall_gates(p; kw...))

#################################################################
# Generators: Hbar = sum of the exponents of all gates
#   kick       e^{-i pi g/2 X}          -> (pi g/2) X_i
#   field      e^{-i h/2 Z}             -> h_i/2 Z_i
#   leg noise  e^{-i dh/2 Z}            -> dh/2 Z
#   Ising      e^{-i J ZZ}              -> J_j Z_j Z_{j+1}
#   XY         e^{-i theta/2 (XX+YY)}   -> theta/2 (X_j X_{j+1} + Y_j Y_{j+1})
#################################################################

function site_op(op, i, L)
    kron_list_ops = [k == i ? op : Matrix{ComplexF64}(I, 2, 2) for k in 1:L]
    reduce(kron, kron_list_ops)
end

function bond_op(op1, op2, j, L)
    ops = [k == j ? op1 : (k == j + 1 ? op2 : Matrix{ComplexF64}(I, 2, 2)) for k in 1:L]
    reduce(kron, ops)
end

# Diagonal part (field, noise, Ising) as a vector, X part and XY part as matrices
function generator_parts(p::DTCParams)
    L = p.L; D = 2^L
    # diagonal: z-values, site 1 is the leftmost kron factor (most significant bit)
    zval(s, i) = 1 - 2 * ((s >> (L - i)) & 1)
    hz = zeros(L)
    for i in 1:L
        hz[i] += p.h[i] / 2
    end
    for j in 1:L-1
        hz[j]   += (p.delh[1, j] + p.delh[3, j]) / 2
        hz[j+1] += (p.delh[2, j] + p.delh[4, j]) / 2
    end
    diag = zeros(D)
    for s in 0:D-1
        e = 0.0
        for i in 1:L
            e += hz[i] * zval(s, i)
        end
        for j in 1:L-1
            e += p.J[j] * zval(s, j) * zval(s, j + 1)
        end
        diag[s+1] = e
    end
    SX = zeros(ComplexF64, D, D)
    for i in 1:L
        SX .+= site_op(ComplexF64.(X), i, L)
    end
    SXY = zeros(ComplexF64, D, D)
    for j in 1:L-1
        SXY .+= bond_op(ComplexF64.(X), ComplexF64.(X), j, L) .+ bond_op(ComplexF64.(Y), ComplexF64.(Y), j, L)
    end
    diag, SX, SXY
end

# Hbar_0 (the zero drive U_0: Ising + field + kick) and Hbar (with the XY term)
function hbar(p::DTCParams; with_theta::Bool = true)
    diag, SX, SXY = generator_parts(p)
    H = Diagonal(diag) + (pi * p.g / 2) * SX
    if with_theta
        H = H + (p.theta / 2) * SXY
    end
    Hermitian(Matrix(H))
end

#################################################################
# Ladder analysis of the Floquet eigenstates in the eigenbasis of an energy
#################################################################

"""
ladder_analysis(U, H; nbins)

U = e^{-i lambda} on its eigenstates, H = sum_n E_n |n><n|.
For each Floquet eigenstate k, P_nk = |<n|psi_k>|^2 is its distribution in energy.
Returns
  lambda           quasienergies in (-pi, pi]
  E                eigenvalues of H
  R                ladder visibility |sum_n P_nk e^{i(E_n - lambda_k)}| per k
  Srung, Sens      entropy of the rung occupation and of the ladder ensemble (P_n ~ DOS) per k
  P                the full weight matrix (only if keepP)
"""
function ladder_analysis(U::AbstractMatrix, H::Hermitian; keepP::Bool = false, dos_sigma = nothing)
    F = eigen(U)
    lambda = -angle.(F.values)
    V = F.vectors
    EH = eigen(H)
    E = EH.values
    W = EH.vectors
    P = abs2.(W' * V)                         # P[n, k]
    D = length(E)
    # smoothed DOS for the ladder ensemble, Gaussian kernel of width ~ local mean spacing * 20
    sig = dos_sigma === nothing ? 0.25 : dos_sigma
    dos(x) = sum(exp.(-(x .- E) .^ 2 ./ (2 * sig^2))) / (sqrt(2pi) * sig)
    Emin, Emax = extrema(E)
    R = zeros(D); Srung = zeros(D); Sens = zeros(D); Pmax = zeros(D)
    for k in 1:D
        ph = E .- lambda[k]
        R[k] = abs(sum(P[:, k] .* cis.(ph)))
        nr = round.(Int, ph ./ (2pi))
        nmin, nmax = extrema(nr)
        occ = zeros(nmax - nmin + 1)
        for n in 1:D
            occ[nr[n]-nmin+1] += P[n, k]
        end
        Srung[k] = -sum(x -> x > 0 ? x * log(x) : 0.0, occ)
        Pmax[k] = maximum(occ)
        # ladder ensemble: rungs lambda + 2 pi n inside the spectrum, weight ~ DOS
        ens = Float64[]
        for n in floor(Int, (Emin - lambda[k]) / 2pi):ceil(Int, (Emax - lambda[k]) / 2pi)
            push!(ens, dos(lambda[k] + 2pi * n))
        end
        ens ./= sum(ens)
        Sens[k] = -sum(x -> x > 0 ? x * log(x) : 0.0, ens)
    end
    out = (lambda = lambda, E = E, R = R, Srung = Srung, Sens = Sens, Pmax = Pmax)
    keepP ? merge(out, (P = P,)) : out
end
