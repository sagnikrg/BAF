

function plotdata_rescale_4L(LList, quant::Array{Array{Float64,1},1})

    L_tot=length(LList)

    Quant_scaled=Array{Array{Float64,1},1}(undef, L_tot)

    for i in 1:L_tot
        Quant_scaled[i]=quant[i].*4^(LList[i])
    end

    return Quant_scaled

end



function plotdata_rescale_2L(LList, quant::Array{Array{Float64,1},1})

    L_tot=length(LList)

    Quant_scaled=Array{Array{Float64,1},1}(undef, L_tot)

    for i in 1:L_tot
        Quant_scaled[i]=quant[i].*2^(LList[i])
    end

    return Quant_scaled

end

#### wrappers

function plotdata_rescale_4L(LList, quant1::Array{Array{Float64,1},1}, quant2::Array{Array{Float64,1},1}, quant3::Array{Array{Float64,1},1})

    return plotdata_rescale_4L(LList, quant1), plotdata_rescale_4L(LList, quant2), plotdata_rescale_4L(LList, quant3);

end

function plotdata_rescale_2L(LList, quant1::Array{Array{Float64,1},1}, quant2::Array{Array{Float64,1},1}, quant3::Array{Array{Float64,1},1})

    return plotdata_rescale_2L(LList, quant1), plotdata_rescale_2L(LList, quant2), plotdata_rescale_2L(LList, quant3);

end
