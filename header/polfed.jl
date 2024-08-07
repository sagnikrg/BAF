######################################################
# Header file dedicated to polynomial filtering
######################################################

using IterativeSolvers
using LinearAlgebra

# Function to perform Arnoldi iteration and return eigenvalues and eigenvectors
function arnoldi_iteration(A, k, v0)
    n = length(v0)
    V = zeros(ComplexF64, n, k)
        H = zeros(ComplexF64, k, k)
    V[:, 1] = v0 / norm(v0)

    for j in 1:k
        w = A * V[:, j]
        for i in 1:j
            H[i, j] = dot(V[:, i], w)
            w -= H[i, j] * V[:, i]
        end
        if j < k
            H[j+1, j] = norm(w)
            V[:, j+1] = w / H[j+1, j]
        end
    end

    return H, V
end

# Function to compute eigenvalues and eigenvectors using the Arnoldi iteration
function arnoldi_eig(A, k, v0)
    H, V = arnoldi_iteration(A, k, v0)
    eigvals, eigvecs = eigen(H)
    return  V * eigvecs
end


function extract_eigenvalues(A, V)
    k = size(V, 2)
    eigenvalues = zeros(ComplexF64, k)
    for i in 1:k
        v = V[:, i]
        Av = A * v
        eigenvalue = dot(v, Av) / dot(v, v)  # Rayleigh quotient
        eigenvalues[i] = eigenvalue
    end


    #kick out eigenvalues that are not converged
    eigenvalues = eigenvalues[abs.(eigenvalues) .> 1]
    
    return eigenvalues
end

function extract_converged_eigenvectors(A, V)
    k = size(V, 2)
    eigenvalues = zeros(ComplexF64, k)
    for i in 1:k
        v = V[:, i]
        Av = A * v
        eigenvalue = dot(v, Av) / dot(v, v)  # Rayleigh quotient
        eigenvalues[i] = eigenvalue
    end

    # Select eigenvectors corresponding to eigenvalues with absolute value greater than 1
    converged_indices = findall(abs.(eigenvalues) .> 1)
    converged_eigenvectors = V[:, converged_indices]
    
    return converged_eigenvectors
end

function arnoldi(A, k, max_iter, v0)
 V= zeros(ComplexF64, length(v0), k)

    for itr in 1:max_iter
        H, V = arnoldi_iteration(A, k, v0)
        eigvals, eigvecs = eigen(H)
        V = V * eigvecs
        v0 = V[:, 1]
    end


    return V
end

############################################
# Benchmarking Arnoldi iteration
############################################


# Function to benchmark the Arnoldi iteration and track convergence
function benchmark_arnoldi(A, max_iter, v0)
    n = length(v0)
    convergence_data = []

    for k in 1:max_iter
        H, V = arnoldi_iteration(A, k, v0)
        eigvals = eigen(H).values
        push!(convergence_data, (k, eigvals))
    end

    return convergence_data
end

using CairoMakie

function plot_convergence(convergence_data)
    iterations = [data[1] for data in convergence_data]
    mod_eigvals_list = [data[2] for data in convergence_data]

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Iteration", ylabel = "Modulus of Eigenvalue", title = "Convergence of Eigenvalues")

    for i in 1:length(mod_eigvals_list[1])
        lines!(ax, iterations, [abs.(mod_eigvals[i]) for mod_eigvals in mod_eigvals_list], label = "Eigenvalue $i")
    end

    axislegend(ax)
    fig
end
