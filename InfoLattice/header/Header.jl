

#Global Header File

using LinearAlgebra             #   Linear Algebra 
#using MKL                       #   MKL
using Random                     #   Random RandomMatrices
using HDF5                      #   For HDF5
#using CairoMakie                #   For Plotting
#using FFTW                      #   For FFT 
#using ITensors,
#using ITensorsVisualization     #   Packages for ITensors
using StatsBase
using Kronecker
#using Distributions

include("Gates.jl")
include("Brickwall.jl")
#include("Plotmods.jl")
include("TransferMat.jl")
include("functions.jl")



# Function to read a single line from a file
function read_model_name(filename)
    open(filename, "r") do file
        return readline(file)
    end
end



function format_duration(duration::Millisecond)
    milliseconds=duration.value % 1000
    total_seconds = div(duration.value, 1000)  # Convert milliseconds to seconds
    hrs = div(total_seconds, 3600)             # Get hours
    mins = div(total_seconds % 3600, 60)       # Get minutes
    secs = total_seconds % 60                  # Get seconds
    return "$(hrs) hour(s), $(mins) minute(s), $(secs) second(s), and $(milliseconds) millisecond(s)"
end

