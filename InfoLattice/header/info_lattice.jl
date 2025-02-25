###############################################
# This is a library with various functions to compute the information lattice of a quantum state.
# The plot function draw_infolattice can be found on plotmods
###############################################



# Von Neumann Entropy of a reduced density matrix with base 2

function von_neumann_entropy(eigenvalues; base::Real = 2.0) 
    eigenvalues= real(eigenvalues)
    # Filter out zero eigenvalues to avoid log(0) errors
    valid_eigenvalues = filter(x -> x > 0, eigenvalues)
    # Compute entropy
    entropy = -sum(λ -> λ * log(λ) / log(base), valid_eigenvalues)
    return entropy
end



# compute the correct reduced density matrix given psi, level and lattice size

function site_entropy(Psi, i, l, L)
    
psi=reshape(Psi,(2^(i-1),2^(l),2^(L-i-l+1)))
psi_permuted = permutedims(psi, (1, 3, 2))
psi=reshape(psi_permuted,(2^(L-l),2^l))

rho=psi'*psi
eigvals, eigvecs=eigen(rho)
#von_neumann_entropy(eigvals)

return von_neumann_entropy(eigvals)
end



# computing the bit information of a site

function site_bit(Psi, i, l, L)
    
    psi=reshape(Psi,(2^(i-1),2^(l),2^(L-i-l+1)))
    psi_permuted = permutedims(psi, (1, 3, 2))
    psi=reshape(psi_permuted,(2^(L-l),2^l))
    
    rho=psi'*psi
    eigvals, eigvecs=eigen(rho)
    #von_neumann_entropy(eigvals)
    
    return l-von_neumann_entropy(eigvals)
end

# implementation of site_bit with ITensors

using ITensors
using ITensorMPS


function site_bit(Psi, i, l, L; method=:ITensor)
    
    # computing the full density matrix and loading it in ITensor

    sites=siteinds("S=1/2",L)
    psi=reshape(Psi, 2,2,2,2,2,2,2,2)
    #psi=permutedims(psi, [8,7,6,5,4,3,2,1])
    psitensor= ITensor(psi, sites)
    rho = psitensor * dag(prime(psitensor))


    for j in 1:(i-1)
    
            rho=rho*delta(sites[j],sites[j]')
    end

    for j in (i+l):L
            rho=rho*delta(sites[j],sites[j]')
    end


    rho_array=array(rho)
    rho_array=reshape(rho_array, 2^(l), 2^(l))
    eigvals, eigvecs=eigen(rho_array')
    #von_neumann_entropy(eigvals)
    
    return l-von_neumann_entropy(real(eigvals))
    
end


# computing the bit information of a site

function site_bit(Psi, i, l, L)
    
    psi=reshape(Psi,(2^(i-1),2^(l),2^(L-i-l+1)))
    psi_permuted = permutedims(psi, (1, 3, 2))
    psi=reshape(psi_permuted,(2^(L-l),2^l))
    
    rho=psi'*psi
    eigvals, eigvecs=eigen(rho)
    #von_neumann_entropy(eigvals)
    
    return l-von_neumann_entropy(eigvals)
end


# the bare information with the connected part

function info_lattice_raw(Psi, L)
    # Predefine the type as a vector of Float64 vectors
    info_lattice_raw_t = Vector{Vector{Float64}}(undef, L)

    for l in 1:L
        bit_info = Float64[]  # Store bit information for sub-lattices of size `l`
        for i in 1:(L - l + 1)
            push!(bit_info, site_bit(Psi, i, l, L))#; method=:ITensor))  # Compute and store bit information
        end
        info_lattice_raw_t[l] = bit_info  # Assign explicitly typed row
    end

    return info_lattice_raw_t
end


#########################
# Info lattice of mutual information at various levels
#########################

function info_lattice(Psi, L)


    info_lattice_raw_t=info_lattice_raw(Psi,L)
    info_lattice_t = [zeros(Float64, size(subarray)) for subarray in info_lattice_raw_t]

    info_lattice_t[1].=info_lattice_raw_t[1]

    bit_info=Float64[]
    
    for i in 1:(L-1)
        push!(bit_info, info_lattice_raw_t[2][i]-info_lattice_raw_t[1][i]-info_lattice_raw_t[1][i+1])
    end

    info_lattice_t[2].=bit_info

    for l in 3:L
      
        bit_info=Float64[]
        for i in 1:(L-l+1)
            push!(bit_info, info_lattice_raw_t[l][i]-info_lattice_raw_t[l-1][i]-info_lattice_raw_t[l-1][i+1]+info_lattice_raw_t[l-2][i+1])
        end
        info_lattice_t[l].=bit_info
    end

    return info_lattice_t
end

