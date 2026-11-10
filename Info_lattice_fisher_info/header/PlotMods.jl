#####################################
# Color Schemes
#####################################

using CairoMakie
using GeometryBasics
using Colors
using ColorSchemes



viridis=[colorant"#440154", colorant"#482878", colorant"#3e4989", colorant"#31688e", colorant"#26828e", colorant"#1f9e89", colorant"#35b779", colorant"#6ece58", colorant"#b5de2b", colorant"#fde725"]


# gradients
greens=[colorant"#cbfff8",colorant"#c4f9f2",colorant"#bdf3ec",colorant"#aee8e0",colorant"#A7E2DA",colorant"#addbd5",colorant"#a5d2cd",colorant"#9dc9c5",colorant"#8db8b5",colorant"#85afad",colorant"#7da6a5",colorant"#769e9d",colorant"#6e9695",colorant"#678d8d",colorant"#5f8585",colorant"#587d7e",colorant"#517576",colorant"#4a6d6e",colorant"#436567",colorant"#3c5d60",colorant"#355558",colorant"#2f4e51",colorant"#28464a",colorant"#223f43",colorant"#1b383c",colorant"#153135"]

blues=[colorant"#08203E", colorant"#0B2441", colorant"#0E2745", colorant"#112B48", colorant"#142E4B", colorant"#17324E", colorant"#1A3552", colorant"#1D3955", colorant"#203C58", colorant"#23405B", colorant"#26435F", colorant"#294762", colorant"#2C4A65", colorant"#2F4E69", colorant"#31526C", colorant"#34556F", colorant"#375972", colorant"#3A5C76", colorant"#3D6079", colorant"#40637C", colorant"#43677F", colorant"#466A83", colorant"#496E86", colorant"#4C7189", colorant"#4F758C", colorant"#527890", colorant"#557C93"]

reds=[colorant"#3D0808", colorant"#400B0B", colorant"#440E0E", colorant"#471111", colorant"#4A1414", colorant"#4E1717", colorant"#511A1A", colorant"#541D1D", colorant"#582020", colorant"#5B2323", colorant"#5E2626", colorant"#622929", colorant"#652C2C", colorant"#692F2F", colorant"#6C3232", colorant"#6F3535", colorant"#733838", colorant"#763B3B", colorant"#793E3E", colorant"#7D4141", colorant"#804444", colorant"#834747", colorant"#874A4A", colorant"#8A4D4D", colorant"#8D5050", colorant"#915353", colorant"#945656"]

yellows=[colorant"#3D2F08", colorant"#43340A", colorant"#48390B", colorant"#4E3D0D", colorant"#53420E", colorant"#594710", colorant"#5E4C11", colorant"#645113", colorant"#6A5515", colorant"#6F5A16", colorant"#755F18", colorant"#7A6419", colorant"#80691B", colorant"#866E1D", colorant"#8B721E", colorant"#917720", colorant"#967C21", colorant"#9C8123", colorant"#A18624", colorant"#A78A26", colorant"#AD8F28", colorant"#B29429", colorant"#B8992B", colorant"#BD9E2C", colorant"#C3A22E", colorant"#C8A72F", colorant"#CEAC31"]

ccafeblues=[colorant"#07151b",colorant"#fadcdc",colorant"#bd8d9e",colorant"#504459",colorant"#edf0c2",colorant"#584931",colorant"#c7d3aa",colorant"#ca3925"]



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


##############################################
# Information Lattice Plot
##############################################

function draw_infolattice(infolattice)
        cmap = ColorScheme(get(ColorSchemes.tempo, range(0.05, 0.95, length=1000)))
        # other possible schemes acton
        fig = Figure(size = (950, 670), fontsize = 20, font = "computer-modern")
        ax = Axis(fig[1, 1], aspect = DataAspect(), 
                xticks = (1.5:0.5:(length(infolattice[1])+0.5), string.(0:0.5:(length(infolattice[1])-1))), 
                yticks = ((1:(length(infolattice[1]))).*(sqrt(3)/2), string.(0:(length(infolattice[1])-1))),
                xlabel = L"\text{sites }(n)",
                ylabel = L"\text{level }(l)",
                #title = "Information Lattice"
                )
        
    
    
        a = Point2f(1, 0)
        b = Point2f(0.5, sqrt(3)/2)
    
        xlims!(ax, 0, length(infolattice[1])+2)
        ylims!(ax, 0, length(infolattice) )
    
        radius = 0.5
        hex_vertices = [Point2f(radius * cos(θ), radius * sin(θ)) for θ in π/6:π/3:2π]
    
    
        for (nb, line) in enumerate(infolattice)
            for (na, val) in enumerate(line)
                coord = na * a + nb * b
                poly = Polygon([v + coord for v in hex_vertices])
    
                poly!(ax, poly, color = val, colormap = cmap, colorrange = (0.0, 1.0), transparency = 0.8)
            end
        end
        
    
        #coloraxis = Colorbar(fig[1, 2], colormap = cmap, label = "Value", size = 20)
            
    # Creating a colorbar with fixed length and specified tick placement
    colorbar = Colorbar(fig[1, 2], colormap = cmap, label = L"I^l_n",
                       # labelpadding = 10,  # Padding between the colorbar and the label
                        height = Relative(0.7),  # Sets the colorbar to half the height of the axis
                        width = 25,  # Sets the colorbar to 10% of the width of the axis
                        ticks = LinRange(0, 1,11),  # Five equally spaced ticks
                        valign = :top)  # Vertically top aligned
                        
    
    
    fig[1, 2] = colorbar
    
        return fig
    end