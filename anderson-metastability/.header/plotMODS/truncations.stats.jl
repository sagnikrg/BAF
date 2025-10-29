### recall to copy back the histgram.jl

########################################
# histogram, truncation etc.
########################################


function extract_number(s::AbstractString)
    m = match(r"\d+", s)
    return m !== nothing ? parse(Int, m.match) : nothing
end                                



############################################
# truncation
############################################

function replace_small_values(arr::Vector{Float64})
    new_arr = similar(arr)  # Create a new array of the same size and type
    for i in eachindex(arr)
        if arr[i] < 1e-15
            new_arr[i] = 1e-15
        else
            new_arr[i] = arr[i]
        end
    end
    return new_arr
end

############################################
# kick out small values
# get a new smaller array
############################################

function kick_out_small_values(arr::Vector{Float64})
    new_arr = Float64[]
    for i in eachindex(arr)
        if arr[i] > 1e-15
            push!(new_arr, arr[i])
        end
    end
    return new_arr
end



############################################
# crossing indices
############################################


function find_crossing_indices(curve1, curve2)
    crossing_indices = []

    for i in 2:length(curve1)
    # Check if the sign of the difference changes
        if (curve1[i] - curve2[i]) * (curve1[i-1] - curve2[i-1]) < 0
            push!(crossing_indices, i-1) # Store the first index of the pair where crossing occurs
        end
    end

    return crossing_indices
end