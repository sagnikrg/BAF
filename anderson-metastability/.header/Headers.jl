




using LinearAlgebra     # for eigenvalues and eigenvectors
using StatsBase         # for statistical functions
using CairoMakie        # for plotting
using HDF5              # for saving data
using Dates             # for benchmark_time
using Pkg               # for package management

include(".function/lindblad.jl")
include(".function/IPR.jl")
include(".function/truncation.jl")

include("hdf5/hdf5Mods.jl")
