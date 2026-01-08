##########################
# function to compute inverse participation ratio
##########################

function IPR(eigenvector)
    return sum(abs.(eigenvector).^4)/sum(abs.(eigenvector).^2)^2
end


function IPR(eigenvector::Vector{ComplexF64}, U_basis::Matrix{ComplexF64})
    eigenvector=U_basis*eigenvector
    return IPR(eigenvector)
end



################################
# Phase ordering
################################    

function phase_ordered_eigvecs(EigA,Eigvec)


    #########################################
    # Extracting eigenvalues and eigenvectors of the input matrix A:
    #########################################

    
    N=real.(EigA);


    #########################################
    # Defining a matrix to store the phase information of the eigenvectors:
    #########################################


    Ph=fill(0.0*im, length(EigA)+1,length(EigA));
    

    
    ########################################################
    # Storing the phase information of the eigenvectors in the matrix Ph:
    ########################################################

    Ph[1,:]=N;
    for i in 1:length(EigA);
        for j in 1:length(EigA)
            Ph[i+1,j]=Eigvec[i,j];
        end
    end
    
    
    ########################################################
    # Ordering the eigenstates from -pi to pi:
    ########################################################

    Phnew=copy(Ph[:,sortperm(real(Ph[1, :]), rev=true)]) # Phase orders the eigenstates from -pi to pi
    
    EigvecNew=copy(Eigvec);
    
    for i in 1:length(EigA);
        for j in 1:length(EigA)
            EigvecNew[i,j]=Phnew[i+1,j];
        end
    end

    EigA=Phnew[1,:];
    

    return EigA, EigvecNew

end


function heatmap_plt(Square)

###############################
# Square for Heatmaps  
###############################
    
# Desired aspect ratio
aspect_ratio = 1.1

# Set the figure width (for example, 1000 pixels)
fig_width = 1000
fig_height = fig_width / aspect_ratio
fig = Figure(resolution = (fig_width, fig_height), fontsize=16)#, title = L"Eigenstructre: $<\omega'|\sigma^Z|\omega>$ for L=8")


ax=Axis(fig[1,1] ,xlabel = L"$|n>$",          #Label for X  
ylabel = L"$<m|$",         #Label for Y
#title = L"Eigenstructre: $<\omega'|\sigma^Z|\omega>$ for L=8",      #Plot Title
xticks = (1:L),                #Xticks
yticks = (1:L) )



    xs = range(1, L, length = L)              #Axis Range X
    ys = range(1, L, length = L)              #Axis Range Y    
    zs = Square
    
    #The Heatmap
    CairoMakie.heatmap!(ax, xs, ys, zs,              #Yticks
    colormap = Reverse(:deep))

                                                             
    Colorbar(fig[1,2], colormap=Reverse(:deep),                       #Colorbar
                 ticks = 0.0:0.1:1.0)            #Colorbar ticks
    
    fig
    
end


function heatmap_plt(Square, cmap)

###############################
# Square for Heatmaps (with custom colormap)
###############################
    
# Desired aspect ratio
aspect_ratio = 1.1

# Set the figure width (for example, 1000 pixels)
fig_width = 1000
fig_height = fig_width / aspect_ratio
fig = Figure(resolution = (fig_width, fig_height), fontsize=16)#, title = L"Eigenstructre: $<\omega'|\sigma^Z|\omega>$ for L=8")


ax=Axis(fig[1,1] ,xlabel = L"$|n>$",          #Label for X  
ylabel = L"$<m|$",         #Label for Y
#title = L"Eigenstructre: $<\omega'|\sigma^Z|\omega>$ for L=8",      #Plot Title
xticks = (1:L),                #Xticks
yticks = (1:L) )



    xs = range(1, L, length = L)              #Axis Range X
    ys = range(1, L, length = L)              #Axis Range Y    
    zs = Square
    
    #The Heatmap
    CairoMakie.heatmap!(ax, xs, ys, zs,              #Yticks
    colormap = cmap)

                                                             
    Colorbar(fig[1,2], colormap=cmap,                       #Colorbar
                 ticks = 0.0:0.1:1.0)            #Colorbar ticks
    
    fig
    
end


function heatmap_plt(Square, type ,cmap)

###############################
# Square for Heatmaps (with custom colormap)
###############################
    
# Desired aspect ratio
aspect_ratio = 1.1

# Set the figure width (for example, 1000 pixels)
fig_width = 1000
fig_height = fig_width / aspect_ratio
fig = Figure(resolution = (fig_width, fig_height), fontsize=22)


ax=Axis(fig[1,1] ,xlabel = L"$|n>$",          #Label for X  
ylabel = L"$<m|$",         #Label for Y
title = "Heatmap of "*string(type)*" for L=$(L)",
#title = L"Eigenstructre: $<\omega'|\sigma^Z|\omega>$ for L=8",      #Plot Title
xticks = (1:L),                #Xticks
yticks = (1:L) )



    xs = range(1, L, length = L)              #Axis Range X
    ys = range(1, L, length = L)              #Axis Range Y    
    zs = Square
    
    #The Heatmap
    CairoMakie.heatmap!(ax, xs, ys, zs,              #Yticks
    colormap = cmap)

                                                             
    Colorbar(fig[1,2], colormap=cmap,                       #Colorbar
                 ticks = 0.0:0.1:1.0)            #Colorbar ticks
    
    fig
    
end