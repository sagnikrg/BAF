using Plots

# Function to plot Anderson orbitals
function plot_anderson_orbitals(H)
    H_dense = Matrix(Hermitian(H))  # Make the matrix Hermitian to ensure real eigenvalues
    eigenvalues, eigenvectors = eigen(H_dense)
    L = size(H_dense, 1)

    fig = Figure(resolution = (800, 600), title = "Anderson Localization")
    ax = Axis(fig[1, 1], title="Anderson Orbitals", xlabel="Site Index", ylabel="Probability Density")
    
    for i in 1:L
        lines!(ax, abs.(eigenvectors[:, i]).^2, label="Eigenvector $(i)", linewidth=2)
    end
    
    Legend(fig[1, 2], ax, "Eigenvectors")
    fig
end

# Define parameters
L = 12  # number of lattice sites
W = 8.0  # Disorder strength

# Create Hamiltonian and plot orbitals
H = create_hamiltonian(L, W)
plot_anderson_orbitals(H)