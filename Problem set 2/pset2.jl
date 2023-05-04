using Plots
import Pkg 
Pkg.add("Parameters")
Pkg.add("Roots")
Pkg.add("Optim")
using Parameters
using Roots
using Optim
#Parameters
@with_kw struct Par
    α::Float64 = 1/3 ;
    z::Float64 = 1;
    σ::Float64 = 2;
    η::Float64 = 1;
    β::Float64 = 0.98;
    δ::Float64 = 1; #Assume full depreciation
    max_iter::Int64 = 2000; #Maximum number of iterations
    dis_tol::Float64 = 1E-4 ;
    # Howard's Policy Iterations
    H_tol::Float64 = 1E-4; #Tolerance for policy function iteration
    n_h = 20;
end
p = Par();
# Steady state values
function SS_values()
    @unpack z,α,β,δ = p
    l_ss = 0.4 # Labor
    klss = ((1-p.β+p.β*p.δ)/(p.β*p.z*p.α))^(1/(p.α-1)); #k/l ratio in the steady state
    k_ss = l_ss*klss #Capital
    c_ss = l_ss*(z*klss^α-δ*klss) #Consumption
    y_ss = c_ss + δ*k_ss #Production
    r_ss = z*α*klss^(α-1) #Rental rate
    w_ss = z*(1-α)*klss^α #Wage
    return l_ss, klss, k_ss, c_ss, y_ss, r_ss, w_ss
end
# Test steady state function
l_ss, klss, k_ss, c_ss, y_ss, r_ss, w_ss = SS_values();

#Find χ so that l_ss = 0.4:
function find_chi()
    @unpack z, α, η, σ, δ = p
    χ = ((z*klss^α- δ*klss)^(-σ)*z*(1-α)*klss^α)/(l_ss^(η+σ))
    return χ
end
println("χ = ", find_chi())

#Define utility function
function utility(k,kp,l)
    @unpack α, z, σ, δ, η = p
    χ = find_chi()
    c = z*k^α*l^(1-α) + (1-δ)*k - kp
    u_c = c^(1-σ)/(1-σ)
    u_l = χ*l^(1+η)/(1+η)
    if c > 0 
        return u_c - u_l
    else 
        return -Inf
    end 
end 

#Part 1a:
#Direct maximization over (k',l)
function Make_K_Grid(n_k)
    # Get SS
    l_ss, klss, k_ss, c_ss, y_ss, r_ss, w_ss = SS_values()
    # Get k_grid
    k_grid = range(1E-5,2*k_ss;length=n_k) ; # Equally spaced grid between 0 and 2*k_ss
    # Return
    return k_grid
end

#Search for the labor supply choice l given k and k': 
function l_optim(k,kp)
    n_l = 50
    l_aux = zeros(n_l)
    l_grid = range(1E-5, 1; length = n_l)
    for i = 1: n_l
        l_aux[i] = utility(k,kp,l_grid[i])
    end
    V_l, index = findmax(l_aux)
    return l_grid[index]
end

function VFI_grid_loop(n_k)
    @unpack max_iter, dis_tol = p
    println(" ")
    println("------------------------")
    println("VFI with grid search and loops - n_k=$n_k")
    # Get SS
    l_ss, klss, k_ss, c_ss, y_ss, r_ss, w_ss = SS_values()
    # Get k_grid
    k_grid   = Make_K_Grid(n_k)  # Equally spaced grid between 0 and 2*k_ss
    # Initialize variables for loop
    V_old  = zeros(n_k)  # Initial value, a vector of zeros
    iter   = 1           # Iteration index
    V_dist = 1
    converge = 1        # Initialize distance
    data_k = zeros(100,6);
    while iter<=max_iter && V_dist>dis_tol
        println("---------------iteration-------------")
        # Update value function
        V_new, G_kp, G_c = T_grid_loop(V_old,k_grid)
        # Update distance and iterations
        V_dist = maximum(abs.(V_new./V_old.-1))
        iter  += 1
        # Update old function
        V_old  = V_new
        # Report progress
        if mod(iter,50)==0
            data_k[:,converge] = V_old
            converge += 1
            println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
        end
    end
    # Check solution
    if (V_dist<=dis_tol)
        # Recover value and policy functions
        V, G_kp, G_c = T_grid_loop(V_old,k_grid)
        # Return
        println("VFI with grid search and loops - Completed - n_k=$n_k")
        println("Iterations = $iter and Distance = ",100*V_dist,"%")
        println("------------------------")
        println(" ")
        data_k[:,converge] = V
        return V, G_kp, G_c, k_grid
    else
        println("Error in VFI with loops - Solution not found ")
        return V_old, G_kp, G_c, k_grid
        #error("Error in VFI with loops - Solution not found")
    end
end

# Define function for Value update and policy functions
function T_grid_loop(V_old,k_grid)
    @unpack β, α, δ, z = p
    n_k  = length(k_grid)
    V    = zeros(n_k)
    G_kp = fill(0,n_k)
    G_c  = zeros(n_k) #positive
    for i = 1:n_k
        V_aux = zeros(n_k) ; # Empty vector for auxiliary value of V(i,j)
        for j = 1:n_k
            l_choice = l_optim(k_grid[i],k_grid[j])
            V_aux[j] = utility(k_grid[i],k_grid[j],l_choice) + β*V_old[j]
        end
        # Choose maximum value given current capital k_i
        V[i], G_kp[i] = findmax(V_aux)
        G_l = l_optim(k_grid[i],k_grid[G_kp[i]])
        G_c[i]   = z*k_grid[i]^α*G_l^(1-α) + (1-δ)*k_grid[i] - k_grid[G_kp[i]]
    end
    return V, G_kp, G_c
end

@time V_100, G_kp_100, G_c_100, k_grid_100 = VFI_grid_loop(100)
data = zeros(100, 1)

for i = 1:100
    data[i]= k_grid_100[G_kp_100[i]]
end

#Plot policy function
plot(k_grid_100, data,
label = "g(k)", legend = :bottomright, width = [10], title = "Policy Function")
plot!(k_grid_100,k_grid_100, legend = :bottomright,label="")
plot!([SS_values()[2]], seriestype="vline", label="Steady State")
xlabel!("Capital")
ylabel!("Capital")
#savefig("HW2_5a_policy")

#Value function graph
plot(k_grid_100, [data_k[:,1]],
label = "iter 50", legend = :bottomright, width = [10], title = "Value Function")
plot!(k_grid_100, [data_k[:,2]],
label = "iter 100", legend = :bottomright, width = [10])
plot!(k_grid_100, [data_k[:,3]],
label = "iter 150", legend = :bottomright, width = [10])
plot!(k_grid_100, [data_k[:,4]],
label = "iter 200", legend = :bottomright, width = [10])
plot!(k_grid_100, [data_k[:,5]],
label = "iter 250", legend = :bottomright, width = [10])
plot!(k_grid_100, [data_k[:,6]],
label = "limit", legend = :bottomright, width = [10])
xlabel!("Capital")
ylabel!("Value")
#savefig("HW2_5a_Value")

function Euler_Error(k,kp,kpp,l)
    @unpack z, α, δ, β, σ = p
    LHS = (z*k^α*l^(1-α) + (1-δ)*k - kp)^(-σ)
    RHS =  β*(z*kp^α^l^(1-α) +(1-δ)*kp - kpp)^(-σ)*(α*z*kp^(α-1)*l^(1-α)+1-σ)
    return (RHS/LHS-1)*100
end

function VFI_Analytical_Results(n_k,G_kp)
    # Get Grid
    k_grid = Make_K_Grid(n_k)
    # Analytical value function
    # Euler error of numerical policy function on grid
    Euler = zeros(n_k)
    for i=1:n_k
        k   = k_grid[i]
        kp  = k_grid[G_kp[i]]
        kpp = k_grid[G_kp[G_kp[i]]]
        l1 = l_optim(k,kp)
        l2 = l_optim(kp,kpp) #should be the same as l1
        Euler[i] = Euler_Error(k,kp,kpp,l1)
    end
    # Return
    return Euler
end

Euler_100 = VFI_Analytical_Results(100,G_kp_100)

#Euler graph
plot(k_grid_100, Euler_100,linetype=:scatter,marker=(:diamond,3),
label = "Error", title= "Euler Equation Error(%)", legend = :bottomright)
xlabel!("Capital")
ylabel!("Percentage Points")
#savefig("HW2_5a_Euler")


#Howard's policy iteration
function VFI_HPI_grid_mat(n_H,n_k)
    @unpack max_iter, dis_tol = p
    println(" ")
    println("------------------------")
    println("VFI with Howard's Policy Iteration - n_k=$n_k")
    # Get SS
    l_ss, klss, k_ss, c_ss, y_ss, r_ss, w_ss = SS_values()
    # Get k_grid
    k_grid   = range(1E-5,2*k_ss;length=n_k) ; # Equally spaced grid between 0 and 2*k_ss
    # Utility matrix

    U_mat = [utility(k_grid[i],k_grid[j],l_optim(k_grid[i],k_grid[j])) for i in 1:n_k, j in 1:n_k]
    # Initialize variables for loop
    V_old  = zeros(n_k) ; # Initial value, a vector of zeros
    iter   = 0          ; # Iteration index
    V_dist = 1          ; # Initialize distance
    while iter<=max_iter && V_dist>dis_tol
        # Update value function
        println("Loop iter= $iter")
        V_new, G_kp, G_c = HT_grid_mat(V_old,U_mat,k_grid,n_h)
        # Update distance and iterations
        V_dist = maximum(abs.(V_new./V_old.-1))
        iter  += 1
        # Update old function
        V_old  = V_new
        # Report progress
        if mod(iter,50)==0
            println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
        end
    end
    # Check solution
    if (V_dist<=dis_tol)
        # Recover value and policy functions
        V, G_kp, G_c = T_grid_loop(V_old,k_grid)
        # Return
        println("VFI with  Howard's Policy Iteration - Completed - n_k=$n_k")
        println("Iterations = $iter and Distance = ",100*V_dist,"%")
        println("------------------------")
        println(" ")
        return V, G_kp, G_c, k_grid
    else
        error("Error in VFI with Howard's Policy Iteration - Solution not found")
    end
end

function HT_grid_mat(V_old,U_mat,k_grid,n_h)
    @unpack β, H_tol, n_h, z, α, δ = p
    # Get Policy Function
    n_k    = length(V_old)
    V,G_kp = findmax( U_mat .+ β*repeat(V_old',n_k,1) , dims=2 )
    V_old  = V
    # "Optimal" U for Howard's iteration
        U_vec = U_mat[G_kp]
    # Howard's policy iteration
    # G_kp is a Cartesian Index
    for i=1:n_h
        V = U_vec .+ β*repeat(V_old',n_k,1)[G_kp]
        if maximum(abs.(V./V_old.-1))<=H_tol
            break
        end
        V_old = V
    end
    # Recover Policy Functions
    G_kp   = [G_kp[i][2] for i in 1:n_k] # G_kp is a Cartesian index
    l_choice2 = zeros(n_k)
    for i=1:n_k
        l_choice2[i] = l_optim(k_grid[i],k_grid[G_kp[i]])
   end

    G_c    = (z*(k_grid.^α).*l_choice2.^(1-α)) + (1-δ)*k_grid .- k_grid[G_kp]
    # Return output
    return V, G_kp, G_c
end

#@time V_50_h, G_kp_50_h, G_c_50_h, k_grid_50_h = VFI_HPI_grid_mat(50,50)
#@time V_200_h, G_kp_200_h, G_c_200_h, k_grid_200_h = VFI_HPI_grid_mat(200,200)
@time V_500_h, G_kp_500_h, G_c_500_h, k_grid_500_h = VFI_HPI_grid_mat(500,500)

Euler_500_h = VFI_Analytical_Results(500,G_kp_500_H)

#Euler graph HPI
plot(k_grid_500_h, Euler_500_h,linetype=:scatter,marker=(:diamond,3),
label = "HPI Euler Error", title= "HPI Euler Equation Error(%)", legend = :topright)
xlabel!("Capital")
ylabel!("Percentage Points")
#savefig("HW2_5b_Euler")

#Value function graph HPI
plot(k_grid_500_h, V_500_h,
label = "n_k = 500", legend = :bottomright, widths = [10], title = " HPI Value Function")
xlabel!("Capital")
ylabel!("Value")
#savefig("HW2_5b_value")

#Plot policy function HPI
data_h = zeros(500)
for i = 1:500
    data_h[i]= k_grid_500_h[G_kp_500_h[i]]
end

plot(k_grid_500_h, data_h,
label = "g(k)", legend = :bottomright, widths = [10], title = "HPI Policy Function")
plot!(k_grid_500_h,k_grid_500_h, legend = :bottomright,label="")
plot!([SS_values()[2]], seriestype="vline", label="Steady State")
xlabel!("Capital")
ylabel!("Capital")
#savefig("HW2_5b_policy")




#MacQueen-Porteus Bounds
# Solve VFI with MPB
function Solve_VFI_MPB(n_k)
    # Get Grid
    k_grid = Make_K_Grid(n_k)
    # Utility matrix
    U_mat = [utility(k_grid[i],k_grid[j], l_optim(k_grid[i],k_grid[j])) for i in 1:n_k, j in 1:n_k]
    # Solve VFI
    V, G_kp, G_c = VFI_grid_MPB(x->T_grid_mat(x,U_mat,k_grid),k_grid)
    # Return Solution
    return V,G_kp, G_c, k_grid
end

# Fixed point with MPB
function VFI_grid_MPB(T::Function,k_grid)
    @unpack dis_tol, β = p
    # Initialize variables for loop
    n_k    = length(k_grid) ; # Number of grid nodes
    V_old  = zeros(n_k)     ; # Initial value, a vector of zeros
    iter   = 0              ; # Iteration index
    V_dist = 1              ; # Initialize distance
    println(" ")
    println("------------------------")
    println("VFI - Grid Search - MPB - n_k=$n_k")
    for iter=1:max_iter
        # Update value function
        V_new, G_kp, G_c = T(V_old)
        # MPB and Distance
        MPB_l  = β/(1-β)*minimum(V_new-V_old)
        MPB_h  = β/(1-β)*maximum(V_new-V_old)
        V_dist = MPB_h - MPB_l
        # Update old function
        V_old  = V_new
        # Report progress
        if mod(iter,100)==0
            println("   VFI Loop: iter=$iter, dist=",100*V_dist,"%")
        end
        # Check Convergence
        if (V_dist<=dis_tol)
            # Recover value and policy functions
            V = V_old .+ (MPB_l+MPB_h)/2
            # Return
            println("VFI - Grid Search - MPB - n_k=$n_k")
            println("Iterations = $iter and Distance = ",100*V_dist)
            println("------------------------")
            println(" ")
            return V, G_kp, G_c
        end
    end
    # Report error for non-convergence
    error("Error in VFI - Grid Search - MPB - Solution not found")
end

# Define function for Value update and policy functions
function T_grid_mat(V_old,U_mat,k_grid)
    n_k    = length(V_old)
    V,G_kp = findmax( U_mat .+ β*repeat(V_old',n_k,1) , dims=2 )
    G_kp   = [G_kp[i][2] for i in 1:n_k] # G_kp is a Cartesian index
    opl = zeros(n_k)
    for i=1:n_k
        opl[i] = l_optim(k_grid[i],k_grid[G_kp[i]])
    end
    G_c    = (z*(k_grid.^α).*opl.^(1-α)) + (1-δ)*k_grid.- k_grid[G_kp]
    return V, G_kp, G_c
end

 @time V_600_MPB, G_kp_600_MPB, G_c_600_MPB, k_grid_600_MPB = Solve_VFI_MPB(600)

 plot(k_grid_600_MPB, V_600_MPB,
 label = "n_k = 600", legend = :bottomright, widths = [10], title = " MPB Value Function")
 xlabel!("Capital")
 ylabel!("Value")
 #savefig("HW2_5c_Value")

 Euler_600_MPB = VFI_Analytical_Results(600,G_kp_600_MPB)
 plot(k_grid_600_MPB, Euler_600_MPB,linetype=:scatter,marker=(:diamond,3),
 label = "MPB Euler Error", title= "MPB Euler Equation Error(%)", legend = :topright)
 xlabel!("Capital")
 ylabel!("Percentage Points")
 #savefig("HW2_5c_Euler")

 data_MPB = zeros(600)
 for i = 1:600
     data_MPB[i]= k_grid_600_MPB[G_kp_600_MPB[i]]
 end
 plot(k_grid_600_MPB, data_MPB,
 label = "g(k)", legend = :bottomright, widths = [10], title = "MPB Policy Function")
 plot!(k_grid_600_MPB,k_grid_600_MPB, legend = :bottomright,label="")
 plot!([SS_values()[2]], seriestype="vline", label="Steady State")
 xlabel!("Capital")
 ylabel!("Capital")
#savefig("HW2_5c_policy")
