
viridis=[colorant"#440154", colorant"#482878", colorant"#3e4989", colorant"#31688e", colorant"#26828e", colorant"#1f9e89", colorant"#35b779", colorant"#6ece58", colorant"#b5de2b", colorant"#fde725"]

function LazadiresDiagram(A)

    #########################################
    # Extracting local hilbert space dimension:
    #########################################

    localdim=length(eigvals(Z));
    
    #########################################
    # Extracting eigenvalues and eigenvectors of the input matrix A:
    #########################################

    EigA=eigvals(A);
    N=angle.(EigA);
    Eigvec=eigvecs(A);
    
    #########################################
    # Defining a matrix to store the phase information of the eigenvectors:
    #########################################


    Ph=fill(0.0*im, length(EigA)+1,length(EigA));
    
    #########################################
    # Extracting Number of qubits from the input matrix A:
    #########################################

    dim=convert(Int64,floor(log(length(EigA))/log(localdim)))
    

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

    Phnew=copy(Ph[:,sortperm(real(Ph[1, :]))]); # Phase orders the eigenstates from -pi to pi
    EigvecNew=Eigvec;
    
    for i in 1:length(EigA);
        for j in 1:length(EigA)
            EigvecNew[i,j]=Phnew[i+1,j];
        end
    end
    
    ########################################################
    # Constructing the correlation matrix:
    ########################################################

    Corr=fill(0.0, length(EigA),length(EigA));
    
    ########################################################
    # Defining the symmetry operator whose correlation we want to calculate:
    ########################################################
    
    Xi=copy(kronecker(Z,kronecker(I(localdim),(dim-1))));
    Corr=abs.(conj(transpose(EigvecNew))*Xi*EigvecNew)
    
###############################
#Phase ordered <n|Z|m>  
###############################


xs = range(-pi, pi, length = 257)              #Axis Range X
ys = range(-pi, pi, length = 257)              #Axis Range Y    
zs = Corr                                    #The Heatmap
fig,ax,hm=CairoMakie.heatmap(xs, ys, zs,
     axis=(; xlabel = L"$T\omega$",          #Label for X  
             ylabel = L"$T\omega'$",         #Label for Y
             title = L"Eigenstructre: $<\omega'|\sigma^Z|\omega>$ for L=8",      #Plot Title
             xticks = (-3:3),                #Xticks
             yticks = (-3:3) ),              #Yticks
             colormap = Reverse(:deep))      #Colormap
                                                         
Colorbar(fig[1,2],hm,                        #Colorbar
             ticks = 0.0:0.1:1.0)            #Colorbar ticks

fig


end



function CircPlot(EigA)


##############################################
#Plot of Arnoldi Eigenvalues 
############################################## 


thet=range( 0, 2*pi, length = 500);
xcirc=cos.(thet)
ycirc=sin.(thet)


scene,ax,ts=CairoMakie.lines(xcirc, ycirc,
     axis=( ; xlabel = L"Re{λ}",                                     #Label for X  
             ylabel = L"Im{λ}",             #Label for Y
             title = "Scatter Plot of Polfed Eigenvalues L=8",
             aspect = 1)
             , color = :black, linewidth = 1, linestyle = :dash,
             )
                                                              
xs = real.(EigA);              #Axis Range X (log scale)
zs = imag.(EigA);              #Axis Range X (log scale)

CairoMakie.scatter!(xs, zs, color = colorant"#00539a", markersize = 17, strokecolor = :black, strokewidth=0.1)

scene

end