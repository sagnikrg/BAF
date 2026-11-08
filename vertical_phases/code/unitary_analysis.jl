#################################################################
# Analysis directly with the unitaries (no generators)
#
#   U_0 = U_F(theta = 0), the zero drive: Ising + field + kick, same disorder
#   V   = U_0^dagger U_F, the XX+YY part dressed by the disorder
#
# Analogues of Morningstar, Huse, Khemani (arXiv 2210.13444), Fig. 3 onwards:
#   V(phi_1, phi_2)   matrix elements of V between eigenstates of U_0 (Fig. 3)
#   V(dphi), G        G_{1<->2} = V^2 D_1 D_2 with D = 2^L / 2pi (Fig. 4)
#   P(delta)          weight of the U_F eigenstates on U_0 eigenstates at
#                     eigenphase distance delta = phi - lambda (Fig. 5a)
#   S_half            entropy of the weight in |delta| < pi/2 vs outside (Fig. 5b)
#################################################################

using LinearAlgebra, Statistics

wrap(x) = mod(x + pi, 2pi) - pi

"""
unitary_eigen(U; alpha)

Eigen decomposition of a unitary via the Hermitian matrix (e^{-i alpha} U + h.c.)/2, which is
several times cheaper than the general solver. Every vector is checked by its residual; if any is
off (accidental near degeneracy of cos(lambda + alpha)), falls back to eigen(U).
Returns (lambda, Q) with U Q = Q diag(e^{-i lambda}), lambda in (-pi, pi].
"""
function unitary_eigen(U::AbstractMatrix; alpha = 0.3711, tol = 1e-8)
    A = Hermitian((cis(-alpha) * U + cis(alpha) * U') / 2)
    Q = eigen(A).vectors
    UQ = U * Q
    d = vec(sum(conj(Q) .* UQ, dims = 1))
    res = maximum(norm(UQ[:, k] - d[k] * Q[:, k]) for k in axes(Q, 2))
    if res > tol
        F = eigen(U)
        return -angle.(F.values), F.vectors, true
    end
    return -angle.(d), Q, false
end

function lsr(phases)
    s = sort(mod.(phases, 2pi))
    d = diff(vcat(s, s[1] + 2pi))
    r = [min(d[i], d[i+1]) / max(d[i], d[i+1]) for i in 1:length(d)-1]
    mean(r)
end

struct Bins
    n1::Int    # bins of dphi in (-pi, pi]
    n2::Int    # bins of phi_1, phi_2 for the 2D map
    nd::Int    # bins of delta for the eigenstate distribution
end
const BINS = Bins(96, 48, 128)

binof(x, n) = clamp(floor(Int, (x + pi) / (2pi) * n) + 1, 1, n)

# log bins of |delta| for the peak at delta = 0: NLOG bins in log10|delta| from LOGMIN to log10(pi),
# everything below LOGMIN goes into the first one
const NLOG = 60
const LOGMIN = -5.0
logbinof(x) = clamp(floor(Int, (log10(max(abs(x), 1e-300)) - LOGMIN) / (log10(pi) - LOGMIN) * NLOG) + 1, 1, NLOG)
logedges() = 10 .^ range(LOGMIN, log10(pi), length = NLOG + 1)

"""
analyse_pair(phi, Q0, U)

phi, Q0: eigenphases and eigenvectors of U_0. U: the Floquet unitary at theta.
Returns a NamedTuple of per realisation reductions.
"""
function analyse_pair(phi, Q0, U; bins = BINS)
    D = length(phi)
    lam, Q, fell = unitary_eigen(U)
    M = Q0' * (U * Q0)                 # |M_mn| = |<m| U_0^dag U_F |n>|
    absM = abs.(M)
    # 1D in dphi = phi_n - phi_m, and 2D in (phi_m, phi_n), off diagonal only
    s1 = zeros(bins.n1); s2 = zeros(bins.n1); sl = zeros(bins.n1); c1 = zeros(Int, bins.n1)
    m2 = zeros(bins.n2, bins.n2); l2 = zeros(bins.n2, bins.n2); c2 = zeros(Int, bins.n2, bins.n2)
    b2 = [binof(x, bins.n2) for x in phi]
    @inbounds for n in 1:D, m in 1:D
        m == n && continue
        a = absM[m, n]
        la = log(max(a, 1e-300))
        b = binof(wrap(phi[n] - phi[m]), bins.n1)
        s1[b] += a; s2[b] += a^2; sl[b] += la; c1[b] += 1
        m2[b2[m], b2[n]] += a; l2[b2[m], b2[n]] += la; c2[b2[m], b2[n]] += 1
    end
    diagw = mean(abs2.(diag(M)))       # weight V keeps on the U_0 eigenstate itself
    # eigenstates of U_F in the U_0 basis
    P = abs2.(Q0' * Q)                 # P[m, k]
    hd = zeros(bins.nd); hlog = zeros(NLOG)
    Shalf = zeros(D); PR = zeros(D); Pself = zeros(D)
    @inbounds for k in 1:D
        pin = 0.0; ipr = 0.0
        for m in 1:D
            dlt = wrap(phi[m] - lam[k])
            w = P[m, k]
            hd[binof(dlt, bins.nd)] += w
            hlog[logbinof(dlt)] += w
            cos(dlt) > 0 && (pin += w)
            ipr += w^2
        end
        pin = clamp(pin, 0.0, 1.0)
        Shalf[k] = (pin <= 0 || pin >= 1) ? 0.0 : -(pin * log2(pin) + (1 - pin) * log2(1 - pin))
        PR[k] = 1 / ipr
        Pself[k] = maximum(@view P[:, k])
    end
    hd ./= D; hlog ./= D
    (s1 = s1, s2 = s2, sl = sl, c1 = c1, m2 = m2, l2 = l2, c2 = c2, diagw = diagw,
     hd = hd, hlog = hlog, Shalf = mean(Shalf), PRn = mean(PR) / D, Pself = mean(Pself),
     lsr = lsr(lam), fell = fell)
end
