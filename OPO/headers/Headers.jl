
using LinearAlgebra
using Random



function complex_noise(n, strength)
    return (randn(n) + 1im * randn(n)) * strength
end



function second_derivative_matrix(L, dx)
   # dx = 1.0  # Assuming unit spacing for simplicity; adjust as necessary
    dx2 = dx^2

    # Create a sparse matrix
    D2 = zeros(L, L)

    for i in 1:L
        D2[i, i] = -2 / dx2  # Main diagonal
        if i == 1
            D2[i, L] = 1 / dx2  # Periodic boundary condition for left end
        else
            D2[i, i-1] = 1 / dx2  # Lower diagonal
        end
        
        if i == L
            D2[i, 1] = 1 / dx2  # Periodic boundary condition for right end
        else
            D2[i, i+1] = 1 / dx2  # Upper diagonal
        end
    end

    return D2
end


function fourth_derivative_matrix(L, dx)
   # dx = 1.0  # Assuming unit spacing for simplicity; adjust as necessary
    dx4 = dx^4

    # Create a sparse matrix
    D4 = zeros(L, L)

    for i in 1:L
        D4[i, i] = 6 / dx4  # Main diagonal

        # Handling the first and last indices specially for periodic boundary conditions
        indices = [i-2, i-1, i+1, i+2]
        coefficients = [-1, -4, -4, -1] / dx4

        for offset in 1:length(indices)
            j = indices[offset]
            # Implement periodic boundary conditions
            if j < 1
                j += L
            elseif j > L
                j -= L
            end
            D4[i, j] = coefficients[offset]
        end
    end

    return D4
end


function pump(f_0, k_p, omega_p, t, x)
    return f_0 .* exp.(im.*(k_p .*x .- omega_p * t))
end

function Diff(n,dx)
    d2= second_derivative_matrix(n, dx)
    d4= fourth_derivative_matrix(n, dx)
    D= d2+sqrt(d4 +4*I(n))

    return D
end

# Diffusion + Nonlinear interaction + Fourth derivative + Decay operator with noise
function gpe_full(u, D, g, k, kappa, dx, pump)
 
    #d2u = zeros(Complex{Float64}, n)
    #d4u = zeros(Complex{Float64}, n)

    #second derivative

    #for i in 2:n-1
    #    d2u[i] =   (u[i+1] - 2*u[i] + u[i-1]) / dx^2
      #  d4u[i] = k * (u[i+2] - 4*u[i+1] + 6*u[i] - 4*u[i-1] + u[i-2]) / dx^4
    #end


    # Apply periodic boundary conditions
    #d2u[1] =  (u[2] - 2*u[1] + u[n]) / dx^2
    #d2u[n] =  (u[1] - 2*u[n] + u[n-1]) / dx^2
    

    # fourth derivative
    #for i in 2:n-1
    #    d4u[i] =   (d2u[i+1] - 2*d2u[i] + d2u[i-1]) / dx^2
      #  d4u[i] = k * (u[i+2] - 4*u[i+1] + 6*u[i] - 4*u[i-1] + u[i-2]) / dx^4
    #end


    # Apply periodic boundary conditions
    #d4u[1] =  (d2u[2] - 2*d2u[1] + d2u[n]) / dx^2
    #d4u[n] =  (d2u[1] - 2*d2u[n] + d2u[n-1]) / dx^2
    

    #d2u=D*d2u
    #q4u=sqrt.(d4u .+4)

    
    Du= D*u

    # Combine the derivatives and add the nonlinear, noise, and decay terms
    return im .* (Du)./2 -im* g .*( abs.(u).^2 .-(1/dx)).* u  -  kappa * u .+im.*pump  
end



# RK4 step for complex-valued functions with full dynamics
function rk4_step_full(u, dt, t, g, f_0, kappa, dx)
    k=0
    k_p=1.4
    omega_p=0.42
    n = length(u)
    D=Diff(n,dx)


    pump_1 = pump(f_0, k_p, omega_p, t, x)

    k1 = gpe_full(u, D, g, k, kappa, dx, pump_1)

    pump_2= pump(f_0, k_p, omega_p, t+0.5*dt, x)
    k2 = gpe_full(u + 0.5*dt*k1, D, g, k, kappa, dx,  pump_2)

    pump_3= pump(f_0, k_p, omega_p, t+0.5*dt, x)
    k3 = gpe_full(u + 0.5*dt*k2, D, g, k, kappa, dx,  pump_3)

    pump_4= pump(f_0, k_p, omega_p, t+dt, x)
    k4 = gpe_full(u + dt*k3, D, g, k, kappa, dx,  pump_4)
    return u + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
end