

#Global Header File

using LinearAlgebra             #   Linear Algebra 
#using MKL                       #   MKL
#using Random                     #   Random RandomMatrices
#using HDF5                      #   For HDF5
#using CairoMakie                #   For Plotting
#using FFTW                      #   For FFT 
#using StatsBase
#using Kronecker

#using Arpack
#using ITensors
#using ITensorMPS
#using ITensorsVisualization     #   Packages for ITensors
#using BenchmarkTools

#include("gates.jl")
#include("brickwall.jl")
#include("plotMods.jl")
#include("functions.jl")

# Dependencies
include(".src/gates/gates.jl")
include(".src/functions/kron.jl")
include(".src/hdf5/hdf5_mods.jl")
include(".src/circuit/brickwall.jl")
include(".src/circuit/time_crystals.jl")
include(".src/functions/histgram.jl")
include(".src/eigen_statistics/eigen_statistics.jl")
include(".src/functions/entanglementee.jl")
include(".src/lazadires_diagram/lazadires_diagram.jl")
include(".src/lazadires_diagram/off-and-diagonals.jl")
