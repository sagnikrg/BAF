# Homemade functions for automated handling of HDF5 attributes

using HDF5
using Pkg
using Dates

include("attrsMods.jl")



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


