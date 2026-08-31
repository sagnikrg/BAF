


function plotdata_log(quant::Array{Array{Float64,1},1})

    L_tot=length(quant)

    Quant_scaled=Array{Array{Float64,1},1}(undef, L_tot)

    for i in 1:L_tot
        Quant_scaled[i]=log10.(quant[i])
    end

    return Quant_scaled

end


# Wrapper:

function plotdata_log(quant1::Array{Array{Float64,1},1}, quant2::Array{Array{Float64,1},1}, quant3::Array{Array{Float64,1},1})

    return plotdata_log(quant1), plotdata_log(quant2), plotdata_log(quant3);

end