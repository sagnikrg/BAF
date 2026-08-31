include("finite_size_eigstat.loaddata.jl")
include("finite_size_eigstat.filter.jl")
include("finite_size_eigstat.rescale.jl")
include("finite_size_eigstat.log.jl")

function plotdata_frl_4L(file_plotdata::HDF5.File, quantity_name::String)

file_attrs=HDF5.attributes(file_plotdata)
LList=read(file_attrs["[Parameters] LList"])

quantity=plotdata_load(file_plotdata, quantity_name)
quantity=-1 .* quantity;
quantity=plotdata_filter(quantity)
quantity=plotdata_rescale_4L(LList, quantity)
quantity=plotdata_log(quantity)
return quantity

end



function plotdata_frl_2L(file_plotdata::HDF5.File, quantity_name::String)

file_attrs=HDF5.attributes(file_plotdata)
LList=read(file_attrs["[Parameters] LList"])

quantity=plotdata_load(file_plotdata, quantity_name)
quantity=-1 .* quantity;
quantity=plotdata_filter(quantity)
quantity=plotdata_rescale_2L(LList, quantity)
quantity=plotdata_log(quantity)
return quantity

end




function plotdata_frl_4L_transposed(file_plotdata::HDF5.File, quantity_name::String)

file_attrs=HDF5.attributes(file_plotdata)
WList=read(file_attrs["[Parameters] WList"])

quantity=plotdata_load(file_plotdata, quantity_name)
quantity=-1 .* quantity;
quantity=plotdata_rescale_4L(LList, quantity)
quantity_transposed=plotdata_transpose(quantity)
quantity_transposed=plotdata_filter(quantity_transposed)
quantity_log=plotdata_log(quantity_transposed)
return quantity_log

end


function plotdata_frl_2L_transposed(file_plotdata::HDF5.File, quantity_name::String)

file_attrs=HDF5.attributes(file_plotdata)
WList=read(file_attrs["[Parameters] WList"])

quantity=plotdata_load(file_plotdata, quantity_name)
quantity=-1 .* quantity;
quantity=plotdata_rescale_2L(LList, quantity)
quantity_transposed=plotdata_transpose(quantity)
quantity_transposed=plotdata_filter(quantity_transposed)
quantity_log=plotdata_log(quantity_transposed)
return quantity_log

end
