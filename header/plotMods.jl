
viridis=[colorant"#440154", colorant"#482878", colorant"#3e4989", colorant"#31688e", colorant"#26828e", colorant"#1f9e89", colorant"#35b779", colorant"#6ece58", colorant"#b5de2b", colorant"#fde725"]

function LazadiresDiagramPlot(Corr)

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