##############
# The functions for the Hamiltonian and the Lindbladian 
# in superoperator representation
# with a single Z decay channel at site 1
##############


#---- Dependencies------:

using SparseArrays

#--------------------



#--------------------------------------
# Hamiltonian:
#--------------------------------------

function H(W,L)
      μ = 2W * rand(L) .- W  # Random onsite potential uniformly distributed between -W and W
        
            sumμ = sum(μ)/2
            # Create the Hamiltonian matrix
            H = spzeros(L, L)
        
            # Fill the Hamiltonian matrix with disorder term (1 - 2 * c†c) and hopping term
            for i in 1:L
                H[i, i] = μ[i]-sumμ  # On-site potential term (-μ[i]/2 from (1 - 2c†c))
                if i < L
                    H[i, i+1] = 0.5  # Hopping term to the next site (scaled by 1/2 due to the transformation)
                    H[i+1, i] = 0.5  # Hopping term from the next site
                end
            end
        
            

            # Convert to dense matrix and create the Liouvillian for the Lindblad equation
            H_dense = Matrix(Hermitian(H))

            return H_dense, μ 
end



##----------------------------
# channels
##----------------------------

function Z_channel(L)
            identity_matrix = Matrix{Float64}(I, L, L)  # Identity matrix of size LxL
            n1 = spzeros(L, L)
            n1[1, 1] = 1.0                      # Corresponds to the occupation number operator c†c on site 1
            decay_channel_z =   2.0 * n1 -identity_matrix
     return decay_channel_z
end



##----------------------------
# Lindbladian
##----------------------------


function Lindbladian(W,L,γ)

        H_dense, μ =H(W,L);
        decay_channel_1=Z_channel(L)

        Hdim = size(H_dense)[1]
        identity_matrix = Matrix{Float64}(I, Hdim, Hdim)  # Corrected identity matrix creation
        
        H_transpose = transpose(H_dense)
        
        
        L_H = kron(identity_matrix, H_dense) - kron(H_transpose, identity_matrix)   # Commutator part

        
        conj_decay_channel_1 = conj(decay_channel_1)
        
        D1 = decay_channel_1' * decay_channel_1
        transpose_D1 = transpose(D1)
            
            
        L1 = kron(conj_decay_channel_1, decay_channel_1) - 0.5 * (kron(identity_matrix, D1) + kron(transpose_D1, identity_matrix))

            
        
        
        Lind = -1im * L_H + γ * (L1) # + L_X1 + L_Y1)  # Sum of all Lindblad terms

    return Lind, μ 
end



# Overloaded function with default γ=1.0


function Lindbladian(W,L)

        H_dense, μ =H(W,L);
        decay_channel_1=Z_channel(L)

        Hdim = size(H_dense)[1]
        identity_matrix = Matrix{Float64}(I, Hdim, Hdim)  # Corrected identity matrix creation
        
        H_transpose = transpose(H_dense)
        
        
        L_H = kron(identity_matrix, H_dense) - kron(H_transpose, identity_matrix)   # Commutator part

        
        conj_decay_channel_1 = conj(decay_channel_1)
        
        D1 = decay_channel_1' * decay_channel_1
        transpose_D1 = transpose(D1)
            
            
        L1 = kron(conj_decay_channel_1, decay_channel_1) - 0.5 * (kron(identity_matrix, D1) + kron(transpose_D1, identity_matrix))

        
        
        
        Lind = -1im * L_H + (L1) # + L_X1 + L_Y1)  # Sum of all Lindblad terms

    return Lind , μ 
end

