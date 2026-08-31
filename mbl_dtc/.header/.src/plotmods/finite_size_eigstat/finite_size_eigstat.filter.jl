
include("../truncation_stats.jl")




function plotdata_filter(quantity::Array{Array{Float64,1},1})

    L_tot = size(quantity)[1]

    Quantity_new=Array{Array{Float64,1},1}(undef, L_tot)

    for i in 1:L_tot
        Quantity_new[i]=kick_out_small_values(quantity[i])
    end

    return Quantity_new

end

################ code above is updated############################
#-----------------------------------------------------------------
############### code below is old ################################


# truncating quantity and errors simultaneously:

function plotdata_filter(quantity, error_up, error_down)

    L_tot, W_tot = size(quantity)

    Quantity_new=Array{Array{Float64,1},1}(undef, L_tot)
    Error_up_new=Array{Array{Float64,1},1}(undef, L_tot)
    Error_down_new=Array{Array{Float64,1},1}(undef, L_tot)


    for i in 1:L_tot
        Quantity_new_temp=kick_out_small_values(quantity[i,:])
        Error_up_new_temp=kick_out_small_values(error_up[i,:])
        Error_down_new_temp=kick_out_small_values(error_down[i,:])

        n_temp=minimum([length(Quantity_new_temp), length(Error_up_new_temp), length(Error_down_new_temp)])

        Quantity_new[i]=Quantity_new_temp[1:n_temp]
        Error_up_new[i]=Error_up_new_temp[1:n_temp]
        Error_down_new[i]=Error_down_new_temp[1:n_temp]

    end

    return Quantity_new, Error_up_new, Error_down_new

end
