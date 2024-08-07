
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

function histgram(data)                                                                                
    n = length(data)
    n_start=data[1]
    n_end=data[end]

#    intervals=n_start:0.000025:n_end
    intervals=LinRange(n_start,n_end,100)
 
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

    Nbins=Nbins/length(data)
    LRange=collect(intervals[1:end-1])
    return LRange, Nbins                                                                                                    
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

################################
# SFF 
################################

function sff(En,t)
    Em=copy(En);
		  
    Eigdiff=fill(0.0*im, length(En), length(En));

    for j in 1:length(En)
        for k in 1:length(En)
            Eigdiff[j,k]=(En[j]/Em[k])^t;
        end
    end

    return real(sum(Eigdiff))
end


################################
# Phase ordering
################################    

function phase_ordered_eigvecs(A)

  #########################################
    # Extracting local hilbert space dimension:
    #########################################

    localdim=length(eigvals(Z));
    
    #########################################
    # Extracting eigenvalues and eigenvectors of the input matrix A:
    #########################################

    EigA,Eigvec=eigen(A);
    N=angle.(EigA);
    
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

    Phnew=copy(Ph[:,sortperm(real(Ph[1, :]))]); # Phase orders the eigenstates from -pi to pi
    EigvecNew=Eigvec;
    
    for i in 1:length(EigA);
        for j in 1:length(EigA)
            EigvecNew[i,j]=Phnew[i+1,j];
        end
    end
    
    return EigA,EigvecNew

end


################################
# LazadiresDiagram 
################################



function LazadiresDiagram(U)

    EigvecNew=phase_ordered_eigvecs(U)[2];

    #########################################
    # Extracting Number of qubits from the input matrix A:
    #########################################

   
    localdim=length(eigvals(Z));
    dim=convert(Int64,floor(log(size(EigvecNew)[1])/log(localdim)))


    ########################################################
    # Constructing the correlation matrix:
    ########################################################

    Corr=fill(0.0, size(EigvecNew));
    
    ########################################################
    # Defining the symmetry operator whose correlation we want to calculate:
    ########################################################
    
    Xi=copy(kronecker(Z,kronecker(I(localdim),(dim-1))));
    Corr=real.(conj(transpose(EigvecNew))*Xi*EigvecNew)
    
Corr


end

################################
# Detuning of eigenvector overlap
################################



function Detuning(U_Dressed)

    EigvecNew=phase_ordered_eigvecs(U_Dressed);
    
    
    localdim=length(eigvals(Z));
    dim=10

    ########################################################
    # Constructing the correlation matrix:
    ########################################################

    Corr=fill(0.0, size(EigvecNew[2]));
    Corr_detu=fill(0.0, size(EigvecNew[1]));
    ########################################################
    # Defining the symmetry operator whose correlation we want to calculate:
    ########################################################
    Xi=copy(kronecker(Z,kronecker(I(localdim),(dim-1))));
    
    #Xi=copy(kronecker(X,(dim)));
    Corr=abs.(conj(transpose(EigvecNew[2]))*Xi*EigvecNew[2])
    Corr_detu=offdiag(Corr)
    
Corr_detu


end


################################
# Detuning with maximal overlap
################################

function Detuning_Max(EigvecNewA, U_Dressed)

    Corr=fill(0.0, size(U_Dressed));
    Corr_max=fill(0.0, size(EigvecNewA));
    Corr_ind=Int64[];

    EigvecNewB=phase_ordered_eigvecs(U_Dressed)[2];
    
    Corr=abs.(conj(transpose(EigvecNewB))*EigvecNewA)
    for i in 1:size(U_Dressed)[1]
        Corr_max[i]=maximum(Corr[i,:])
        push!(Corr_ind,argmax(Corr[i,:]))
    end
    
return Corr_ind,Corr_max


end

