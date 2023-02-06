#Problem set 2
using Plots
import Pkg 
Pkg.add("Parameters")
Pkg.add("Roots")
using Parameters
using Roots
#Parameters
@with_kw struct Params
    α::Float64 = 1/3 ;
    z::Float64 = 1;
    σ::Float64 = 2;
    η::Float64 = 1;
    β::Float64 = 0.95;
    δ::Float64 = 1; #Assume full depreciation
    max_ier::Int64 = 2000; #Maximum number of iterations
    dis_tol::Float64 = 1E-8 ;
    h_tol::Float64 = 1E-8; #for Howard's policy iteration
end
p = Params();

#Calculate part 4: 
klss = ((1-p.β+p.β*p.δ)/(p.β*p.z*p.α))^(1/(p.α-1)); #k/l ratio in the steady state
l_ss = 0.4;
χ = (((p.z*klss^p.α)- p.δ*klss)^(-p.σ)*p.z*(1-p.α)*klss^p.α)/(l_ss^(p.η+p.σ));
println("χ = ", χ)

#Steady state capital:
k_ss = klss*l_ss;

#Part 5:
#Make the grid for k
function kgrid(k_ss,n_k)
    k_grid = range(1E-5,2*k_ss, length=n_k)
    return k_grid
end

function utility(k,kp,χ,p::Params)
    @unpack σ, η, α, δ, z = p
    #Solve for the static choice of labor supply by bisection
    #define the foc for l
    f(l_sol)= (z*k^α*l_sol^(1-α)+(1-δ)k-kp)^(1-σ)*(1-α)*z*k^α*l_sol^(-α) - χ*l_sol^η
    l_sol = find_zero(f,(0,10), Bisection())
    if 0 <= l_sol <= 1
        l = l_sol
    else
        l=1
    end 
    c = z*k^α*l^(1-α) + (1-δ)*k - kp
    if c>0
        return c^(1-σ)/(1-σ) - χ* l^(1+η)/(1+η)
    else
        return -Inf
    end
end

# Euler Error
function Euler_Error(k,kp,kpp,l,p::Params)
    # Return percentage error in Euler Equation
    @unpack z, α, δ, β, σ = p
    LHS = (z*k^α*l^(1-α) + (1-δ)*k - kp)^(-σ)
    RHS =  β*(z*kp^α^l^(1-α) +(1-δ)*kp - kpp)^(-σ)*(α*z*kp^(α-1)*l^(1-α)+1-σ)
    return (RHS/LHS-1)*100
end

##(a) Plain VFI
# Define function for value update and policy function
function T_operator(V_old, U_mat, p::Params)
    @unpack β = p
    n_k = length(V_old)
    V, G_kp = findmax(U_mat .+  β*repeat(V_old, n_k, 1) , dims = 2 )
    G_kp = [G_kp[i][2] for i in 1:n_k] 
    return V, G_kp
end

# Solve VFI with grid search and loops
function VFI_grid(T::Function, k_grid, p::Params)
    # VFI parameters
    @unpack max_iter, dist_tol = p
    # Initialize variables for loop
    n_k    = length(k_grid) ; 
    V_old  = zeros(n_k)     ; 
    V_dist = 1              ; 
    for iter = 1:max_iter
        # Update value function
        V_new, G_kp = T(V_old)
        # Update distance and iterations
        V_dist = maximum(abs.(V_new./V_old.-1))
        # Update old function
        V_old  = V_new
        if V_dist<= dist_tol 
            return V_new, G_kp
        end
    end
    # If loop ends there was no convergence -> Error!
    error("No convergence")
end




