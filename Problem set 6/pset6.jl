import Pkg
using Random, Distributions
using Parameters
using Statistics
using Plots
using LinearAlgebra
using Interpolations
using Dierckx
using ForwardDiff
using Optim
using Roots

#Set seed
Random.seed!(420)
#Parameters
@with_kw struct Par 
    #Model parameters
    β::Float64=0.98;
    α::Float64= 1/3;
    z0::Float64 = 1; #Reference level for productivity
    σ::Float64=2;
    η::Float64=1;
    δ::Float64=0.05;
    ρ::Float64= 0.9;
    σ_η::Float64= 0.01;
    #VFI parameters
    max_iter::Int64 = 2000; 
    dist_tol::Float64 = 1E-9;
    #Policy functions
    H_tol::Float64 = 1E-9;
    N_H::Int64 = 20;
    #Minimum consumption for numerical optimization
    c_min::Float64 = 1E-16;
end
p = Par();
# Steady state values
function SS_values(p::Par)
    @unpack z0,α,β,δ = p
    l_ss = 0.4 # Labor
    klss = ((1-β+β*δ)/(β*z0*α))^(1/(α-1)); #k/l ratio in the steady state
    k_ss = l_ss*klss #Capital
    c_ss = l_ss*(z0*klss^α-δ*klss) #Consumption
    y_ss = c_ss + δ*k_ss #Production
    r_ss = z0*α*klss^(α-1) #Rental rate
    w_ss = z0*(1-α)*klss^α #Wage
    return l_ss, klss, k_ss, c_ss, y_ss, r_ss, w_ss
end
# Test steady state function
l_ss, klss, k_ss, c_ss, y_ss, r_ss, w_ss = SS_values(p);
#Find χ so that l_ss = 0.4:
function find_chi(p::Par,klss,l_ss)
    @unpack z0, α, η, σ, δ = p
    χ = ((z0*klss^α- δ*klss)^(-σ)*z0*(1-α)*klss^α)/(l_ss^(η+σ))
    return χ
end
global χ = find_chi(p,klss,l_ss);

# Function to make grid for capital 
function Make_K_Grid(n_k,θ_k,p::Par)
    # Get SS
    l_ss, klss, k_ss, c_ss, y_ss, r_ss, w_ss = SS_values(p)
    # Lower and upper bounds
    lb = 1E-5
    ub = 2*k_ss
    # Get k_grid
    if θ_k≠1
        k_grid = PolyRange(lb,ub;θ=θ_k,N=n_k)
    else
    k_grid = range(lb,ub,length=n_k)
    end
    # Return
    return k_grid
end

# Function that returns the percentage error in the Euler equation
function Euler_Error(k,kp,kpp,l)
    @unpack z0, α, δ, β, σ = p
    LHS = (z0.*k.^α.*l^(1-α) .+ (1-δ).*k .- kp).^(-σ)
    RHS =  β.*(z0.*kp.^α.*l^.(1-α) +(1-δ).*kp .- kpp).^(-σ).*(α.*z0.*kp.^(α-1).*l.^(1-α).+(1-σ))
    return real((RHS./LHS.-1).*100)
end

#Define utility function
function utility(z,k,kp,l, p::Par)
    @unpack α, σ, δ, η, c_min = p
    c = max(z.*k.^α.*l.^(1-α) + (1-δ).*k .- kp,c_min)
    u_c = c.^(1-σ)./(1-σ)
    u_l = χ.*l.^(1+η)./(1+η)
    return u_c .- u_l
end 

#Derivative of utility function wrt labor l
function d_utility_l(k,z,kp,l,p::Par)
    @unpack α,δ,σ,η,c_min = p
    c = z*(k^α)*(l^(1-α))+(1-δ)*k-kp
    d_u = 0
    if c>c_min
        d_u = c^(-σ)
    else
        d_u = c_min^(-σ)
    end
    return d_u*z*(k^α)*(1-α)*(l^(-α))-χ*(l^η)
end

# Derivative of utility function wrt capital k'
function d_utility_kp(k,z,kp,l,p::Par)
    @unpack α,δ,σ,η,c_min = p
    c = z*(k^α)*(l^(1-α))+(1-δ)*k-kp
    d_u = 0
    if c>c_min
        d_u = c^(-σ)
    else
        d_u = c_min^(-σ)
    end
    return -d_u
end

#Part 1a:
# Function to distretize AR(1) markov process with Rouwenhorst (1995)
function Rouwenhorst95(N,p::Par)
    @unpack ρ,σ_η=p
    # INPUTS:
        # ρ: persistence of unerlying AR(1) process where log(z') = ρlog(z)+η
        # σ_z: Std dev of inovation η in AR(1) process where η∼N(0,σ^2)
        # N: Size of grid for discrete process
    # OUTPUT:
        # z: All possible values of discretized AR(1) process, equally spaced grid of size N
        # Π: Matrix of transition probabilities
        # PDF_z: Stationary PDF of z
    #---------------------------------------------------------------------------
    Π = zeros(N,N)
    Π_Nm = zeros(N-1,N-1)
    P = (1+ρ)/2
    ϕ = σ_η*(sqrt((N-1)/(1-ρ^2)))
    z = range(-ϕ,ϕ;length=N)
    if N==2
        Π = [P 1-P;1-P P]
    else
        Π_Nm = Rouwenhorst95(N-1,p)[2]
        o = zeros(N-1)
        Π = P*[Π_Nm o; o' 0] + (1-P)*[o Π_Nm; 0 o'] + (1-P)*[o' 0; Π_Nm o] + P*[0 o';o Π_Nm]
        Π = Π./repeat(sum(Π,dims=2),1,N)
    end
    PDF_z = pdf.(Binomial(N-1,0.5),(0:N-1))
    return (z,Π,PDF_z)
end


