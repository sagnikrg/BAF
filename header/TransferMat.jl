
#############################################################################################
#Computing the Right Translation Matrtix
#It is based on the observation that the bases are left shifted 
#in their binary form when they are multiplied with the corresponding level dimension.
#############################################################################################

function Tr(level,n)
    #n = 3;              ## number of q-dits
    #level = 2;          ## local dimension 
    dim = level^n;      ## Dimension of Hilbert Space



    A=fill(0+im*0, dim,dim);

    for i in 1:dim-1
        j=1+mod((i-1)*level,(dim-1)) ;
        A[i,j]=1;
    end

    A[dim,dim]=1;

    A
end
;



################################
#function for Histogram
################################

function histgram(data, intervals)                                                                                
    ranges = intervals[1:end-1] .=> intervals[2:end]                                                               
    bins = [similar(data, 0) for _ in 1:length(ranges)]  
    Nbins=fill(1,length(ranges) )                                                                         
    for x in data                                                                                                  
        for (i, (a, b)) in pairs(ranges)                                                                                       
            if a <= x < b                                                                                          
                push!(bins[i], x)                                                                                  
                break                                                                                              
            end                                                                                                    
        end                                                                                                        
    end
    for i in 1:length(ranges)
        Nbins[i]=length(bins[i]);
    end

    return Nbins                                                                                                    
end                                                                                                                


################################
# Level Spacing Statistics
################################ 


function LevelSpacing(EigA)

    N=fill(0.0,length(EigA))        # N is the array of eigenvalues
   
    N=angle.(EigA);                 # N is the array of eigenvalues in angle form
    N= sort(N,rev=true);            # N is the array of eigenvalues in angle form in descending order
    EigA1=copy(N);                  # EigA1 is the array of eigenvalues in angle form in descending order
   
    ls = deleteat!(EigA1,1);        # ls is the array of eigenvalues in angle form in descending order without the first element
    la = deleteat!(N,length(N));    # la is the array of eigenvalues in angle form in descending order without the last element            
    
    m=copy(la-ls)/mean(la-ls);      # m is the array of level spacing statistics

    m

end    


function LevelSpacingRatio(EigA)

    N=fill(0.1,length(EigA))
    N=angle.(EigA);

    N= sort(N,rev=true);


    EigA1=copy(N);
    ls = deleteat!(EigA1,1);
    la = deleteat!(N,length(N));
    m=copy(la-ls)/mean(la-ls)

    #histogram(m)

    #Computing Level Spacing Ratio

    n=fill(0.1,length(EigA)-2);

        for i in 1:length(n)
            #n[i]=m[i+1]/m[i];
            n[i]= minimum([m[i], m[i+1]])/maximum([m[i], m[i+1]]);
        end

    #histogram(n)
mean(n)
end    


function LevelSpacingExtended(EigA)

    En=angle.(EigA);
    Em=copy(En);
        
    Eigdiff=fill(0.0*im, 2^L, 2^L);
        
    for j in 1:2^L
        for k in 1:2^L
            Eigdiff[j,k]=En[j]-Em[k];
        end
    end
    
    m    

end


################################
# finding the crossing indices
################################

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




function pi_diag(Corr)
    L=size(Corr)[1];
    Lh=convert(Int64,L/2)
    
    A=[]
    for i in 1:Lh
        A=push!(A, Corr[i,i+Lh])
        A=push!(A, Corr[i+Lh,i])
    end

    A
end



function offdiag(Corr, k)
    L=size(Corr)[1];
    
    A=[]
    for i in 1:L-k
        A=push!(A, Corr[i,i+k])
        A=push!(A, Corr[i+k,i])
    end

    A
end


function pi_offdiag(Corr,k)
    L=size(Corr)[1];
    Lh=convert(Int64,L/2)
    
    A=[]
    for i in 1:Lh+k
        A=push!(A, Corr[i,i+Lh-k])
        A=push!(A, Corr[i+Lh-k,i])
    end

    for i in 1:Lh-k
    A=push!(A, Corr[i,i+Lh+k])
    A=push!(A, Corr[i+Lh+k,i])
    end

    A
end


function entanglement_entropy(eigenstate)
    # Compute the density matrix of the eigenstate
    density_matrix = eigenstate * eigenstate'

    # Compute the eigenvalues of the density matrix
    eigenvalues = eigvals(density_matrix)

    # Compute the entanglement entropy
    entropy = -sum(eigenvalues .* log.(eigenvalues))

    return entropy
end