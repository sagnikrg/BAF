




using LinearAlgebra     # for eigenvalues and eigenvectors
using StatsBase         # for statistical functions
using HDF5              # for saving data
using Dates             # for benchmark_time
using Pkg               # for package management
using SparseArrays      # for SparseArrays


include(".function/lindblad.jl")
include(".function/IPR.jl")
include("hdf5/hdf5Mods.jl")
