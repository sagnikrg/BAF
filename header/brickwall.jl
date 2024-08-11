using ITensors
ITensors.disable_warn_order()
#######################################
# Computing Kronecker Products
#######################################

# Needs LinearAlgebra


# For Powers of Kronecker Product:

function kron_product(X, n)
    result = X
    for i in 2:n
        result = kron(result, X)
    end
    return result
end



# For Listed Kronecker Product:

function kronlist(List, ind)
    A=List[ind[1]];
    for i in 2:length(ind)
    A=kron(A,List[i]);
    end
    A 
end


########################################################
# The Brickwall as a Function 
########################################################

# For now we restrict to spin 1/2 systems, i.e. q-dits with level=2
 

function brickwall(L,thetamean,epsilon)
  
    #Constructing the background Z field
    
        #h=rand(L)*pi;
        #Ind=collect(1:L)
        #ZRow=copy(kronlist(RZ.(h),Ind));

 
        #Constructing the Random Brickwall 


            J=rand(L)*pi;                     #Ising Even Disorder on the two body gates

            fonez=copy(kron(Z,Z))               #
            ftwox=copy(kron(X,X));              #For the two body XX+YY gates
            ftwoy=copy(kron(Y,Y));              #

    
            thetadev=pi/50;
            theta=thetamean+randn(1)[]*thetadev;               #Interaction (Normal sampling)
                                      
            delh=randn(4)*pi/50;                                          #Imperfection in Z tuning (Normal sampling)

        FU=fill(fill(0.1+im, 4,4), L);


        for j in 1:length(FU)
        
            int1=kron(RZ(delh[1]),RZ(delh[2]));
                int2=exp(-im*J[j]*fonez-im*theta/2*(ftwox+ftwoy));
                int3=kron(RZ(delh[3]),RZ(delh[4]));

                FU[j]=int3*int2*int1;
            end

        Indodd=collect(1:2:L-1);    
        Indeven=collect(2:2:L-1);    
        
        # for even L

        if L%2==0
           

        UOdd=copy(kronlist(FU,Indodd));
        UEven=copy(kron(I(2),kronlist(FU,Indeven),I(2)));
 
        end

        # for odd L

        if L%2==1
        
        UOdd=copy(kron(kronlist(FU,Indodd)),I(2));
        UEven=copy(kron(I(2),kronlist(FU,Indeven)));
        end
            


        #Constructing the X Kicks 
 
            g=pi*(1-epsilon);
            XRow=copy(kron_product(RX(g),L));


          A=XRow*UEven*UOdd;
          

 A    
end
;
#check header

########################################################
# The Brickwall as a Function 
########################################################

# For now we restrict to spin 1/2 systems, i.e. q-dits with level=2
 


########################################################

# The Brickwall as a Tensor Function
# format

# |     |     |     |     |     |     |     |  
#---------------------------------------------
#       |  2  |     |  4  |     |  6  |         
#---------------------------------------------
# |  1  |     |  3  |     |  5  |     |  7  |
#---------------------------------------------
# |     |     |     |     |     |     |     |

########################################################


function brickwall_tensor(L,thetamean,epsilon)

 #########################################
 # Background Disorder for the Z field
 #########################################

   h=rand(L)*2pi;              
   ZRow=RZ.(h);


 ###########################################
 #Constructing the Random Brickwall 
 ###########################################

        J=rand(L)*pi;             #Ising Even Disorder on the two body gates

        fonez=copy(kron(Z,Z))
        ftwox=copy(kron(X,X));  #For the two body XX+YY gates
        ftwoy=copy(kron(Y,Y));


        thetadev=pi/50;
        theta=thetamean+randn(1)[]*thetadev;  #Interaction
                                  
    FU=fill(fill(0.0*im, 4,4), L-1);
    FUTensor=fill(fill(0.0*im, 2,2,2,2), L-1);

    for j in 1:length(FU)
    
            delh=randn(4)*pi/50;         #Imperfection in Z tuning
            int1=kron(RZ(delh[1]),RZ(delh[2]));
            int2=exp(-im*J[j]*fonez-im*theta/2*(ftwox+ftwoy));
            int3=kron(RZ(delh[3]),RZ(delh[4]));

            FU[j]=int3*int2*int1;
    end

 ##########################################
 #Constructing the X Kicks 
 ##########################################


    g=pi*(1-epsilon);



 #############################################
 # Defining two body gates as array of Tensors        
 ##############################################


    #for i in 1:2:L-1
        #FU[i]=FU[i]*kron(ZRow[i],ZRow[i+1]);
        #FUTensor[i]=reshape(Int,2,2,2,2)
    #end

    for i in 2:2:L-1
        FU[i]=kron(RX(g),RX(g))*FU[i];
        #FUTensor[i]=reshape(Int,2,2,2,2)
    end

 ############################################        
 # Open Boundary Condition:
 ############################################

 FU[1]=kron(RX(g),I(2))*FU[1];
 FU[L-1]=kron(I(2),RX(g))*FU[L-1];


 for i in 1:2:L-1
    #FU[i]=FU[i]*kron(ZRow[i],ZRow[i+1]);
    FUTensor[i]=reshape(FU[i],2,2,2,2)
end

for i in 2:2:L-1
    #FU[i]=kron(RX(g),RX(g))*FU[i];
    FUTensor[i]=reshape(FU[i],2,2,2,2)
end
 # Returning Brickwall:        
 FUTensor

end
;



function brickwall_bare(L,thetamean)
  
    #Constructing the background Z field
    
        h=rand(L)*pi/2;
        Ind=collect(1:L)
        ZRow=copy(kronlist(RZ.(h),Ind));

 
        #Constructing the Random Brickwall 


            J=rand(L)*pi/2;                     #Ising Even Disorder on the two body gates

            fonez=copy(kron(Z,Z))               #
            ftwox=copy(kron(X,X));              #For the two body XX+YY gates
            ftwoy=copy(kron(Y,Y));              #

    
            thetadev=pi/50;
            theta=thetamean+randn(1)[]*thetadev;               #Interaction
                                      
            delh=randn(4)*pi/50;                                          #Imperfection in Z tuning

        FU=fill(fill(0.1+im, 4,4), L);


        for j in 1:length(FU)
        
            int1=kron(RZ(delh[1]),RZ(delh[2]));
                int2=exp(-im*J[j]*fonez-im*theta/2*(ftwox+ftwoy));
                int3=kron(RZ(delh[3]),RZ(delh[4]));

                FU[j]=int3*int2*int1;
            end

        Indodd=collect(1:2:L-1);    
        Indeven=collect(2:2:L-1);    
        
        # for even L

        if L%2==0
           

        UOdd=copy(kronlist(FU,Indodd));
        UEven=copy(kron(I(2),kronlist(FU,Indeven),I(2)));
 
        end

        # for odd L

        if L%2==1
        
        UOdd=copy(kron(kronlist(FU,Indodd)),I(2));
        UEven=copy(kron(I(2),kronlist(FU,Indeven)));
        end
            



          A=UEven*UOdd*ZRow;
          

 A    
end
;

function kick(L,epsilon)
    g=pi*(1-epsilon);
    XRow=copy(kron_product(RX(g),L));
    XRow

end

    ########################################
    # time evolution with ITensor
    ########################################



    
function itensorise(FUTensor, sites, dummysites)
    L=length(FUTensor);
    gates = ITensor[]
    for i in 1:2:L
        push!(gates, ITensor(FUTensor[i],dummysites[i+1],dummysites[i],sites[i],sites[i+1]))
    end
  

    for i in 2:2:L
        push!(gates, ITensor(FUTensor[i],sites[i+1],sites[i],dummysites[i],dummysites[i+1]))
    end
  
    gates
end


##################################
# Note in this case the order of the gates in brick has changed.
# Unlike tev_1 case we first put all the odd gates and then the even gates.
# This has to be accounted appropriately in the tev function.
##################################

function brickwall_tev(Psi, brick, sites, dummysites)

    L=length(brick)
    Lhalf=Int((L+1)/2)
    for i in 1:Lhalf
        Psi=brick[i]*Psi
    end

    #for i in 1:L+1
    #    Psi=Psi*delta(sites[i],dummysites[i])
    #end


    for i in (Lhalf+1):L
        Psi=brick[i]*Psi
    end
    
    #for i in 2:L
        Psi=Psi*delta(sites[1],dummysites[1])
        Psi=Psi*delta(sites[L+1],dummysites[L+1])
    #end
    
   
    # psi=apply(brick[1:2:L-1], psi);  # Apply the odd gates
   # psi=apply(brick[2:2:L-1], psi);  # Apply the even gates
    Psi
end



###################################
# Néel state
###################################

function Neel_state(sites)
    # Create a site set for N spin-half sites
   

    # Initialize the MPS as a Néel state
    state = productMPS(sites, n -> isodd(n) ? "↑" : "↓")

    return state
end


###################################
# Random Bit string
###################################

function random_BitString(sites)
    # Create a site set for N spin-half sites

    N=length(sites)
  # Create an array to store the random states
random_states = Vector{String}(undef, N)

# Fill the array with random "↑" or "↓" states
for i in 1:N
    random_states[i] = rand(Bool) ? "↑" : "↓"
end

# Create the MPS representing the product state with the random bit string
psi = productMPS(sites, random_states)

    return psi
end


function EntanglementEntropy(eigvec, l)

    lhalf=Int(l/2)
    
    eigreshape=reshape(eigvec,2^lhalf,2^lhalf)
    
    U,S,V=svd(eigreshape)
    
    ee=0.0
    
    for i in 1:2^lhalf
        if S[i]>1e-10
        ee=ee+-S[i]^2*log(S[i]^2)
        end
    end
    
    return ee
    
end
    



function halfchainee(Psi, sites)

    L=length(sites)
    Lhalf=Int(L/2)
    N=2^Lhalf;


    lft=combiner(sites[1:Lhalf], tags="left half")
    rgt=combiner(sites[Lhalf+1:L], tags="right half")
    Psi_svd=Psi*lft
    Psi_svd=Psi_svd*rgt


    U,V,D=svd(Psi_svd, inds(Psi_svd)[1]);
    ee=0.0

for i in 1:N
    if V[i,i]>1e-10
    ee=ee+-V[i,i]^2*log(V[i,i]^2)
    end
end

ee
end

function halfchainee_compile()

L=6;
epsilon=0.9;
thetamean=0.0;


# Initialise the product state and the gates

sites = siteinds("S=1/2", L);
dummysites = siteinds("S=1/2", L);
psi = productMPS(sites, fill("↑", L));

Psi=psi[1]*psi[2]
for i in 3:L
    Psi=Psi*psi[i]
end

brick=brickwall_tensor(L,thetamean,epsilon);
brick=itensorise(brick,sites,dummysites);

e_tev=Array{Float64}(undef, 0)

for t in 1:100
    Psi=brickwall_tev(Psi, brick, sites, dummysites);
    push!(e_tev,halfchainee(Psi,sites))
end

end
   