# Packages for Random Unitaries
# This will be later expanded to include All the Standard Random Unitary Ensembles as well as Harr Random Matrices


#Hermite Ensembles (with one global symmetry)

#For these ensembles the SD of the diagonals has to be double of that of the offdiagoanal ones
#Following standard methods we choose the SD of diagonal to be 1


l = Normal(0, 0.50);   #off diagonal elements
d = Normal(0, 1.00);   #diagonal elements


function GOE(n)
 
    A=rand(l,n,n);

    for i in 1:n
        A[i,i]=rand(d);
    end

    RH=(A+transpose(A))/2;      #Symmetrisation
    RH
end

     


function GUE(n)
 
    A=rand(l,n,n)+im*rand(l,n,n);

    for i in 1:n
        A[i,i]=rand(d)+im*rand(d);
    end

    RH=(A+A')/2;                #Symmetrisation
    U=exp(-im*RH/2);

    U
end

     
    