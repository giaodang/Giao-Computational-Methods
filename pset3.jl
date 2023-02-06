using Plots
import Pkg
Pkg.add("Dierckx")
using Dierckx
Pkg.add("Interpolations")
using Interpolations
Pkg.add("ForwardDiff")
using ForwardDiff
Pkg.add("LaTeXStrings")
using LaTeXStrings
Pkg.add("Latexify")
using Latexify
Pkg.add("PyPlot")

##Part 1
#Define the utility functions
function u1(x)
    log(x)
end
function u2(x)
    sqrt(x)
end
function u(x,σ)
    x^(1-σ)/(1-σ) 
end
u3(x)=u(x,2);
u4(x)=u(x,5);
u5(x)=u(x,10);
#Define the grid
grid = range(0.05,2, length = 1000);
#Polynomial approximation - Newton basis: 
function diff(x::Array, y::Array)
    m = length(x)
    a = [y[i] for i in 1:m]
    
    for j in 2:m
        for i in reverse(collect(j:m))
            a[i]=(a[i]-a[i-1])/(x[i]-x[i-(j-1)])
        end
    end
    return(a)
end

function newton(x::Array,y::Array,z)
    m=length(x) 
    a=diff(x,y)
    sum=a[1]
    pr=1.0
    for j in 1:(m-1)
        pr=pr*(z-x[j])
        sum=sum+a[j+1]*pr
    end
    return sum
end

label = ["Log","Square root","CRRA" *latexinline("σ=2"), "CRRA"*latexinline("σ=5"),"CRRA"*latexinline("σ=10")]
funs = [u1,u2,u3,u4,u5]
for i in 1:5
    name = label[i]
    for n = (4,6,11,21)
        fn = funs[i].(grid)
        xi = collect(range(0.05,2;length=n)) ;
        yi = funs[i].(xi)
        interp=map(z->newton(xi,yi,z),grid)
        gr()
        plot(title="Interpolation $name n=$n - Newton Polynomial")
        plot!(grid,fn,linewidth=3,label = "Function: $name",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(grid,interp,linewidth=3,label="Interpolation")
        plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
        savefig("graphs/Newton $name _n_$n")
    end
end

#Use linear spline


#Interpolation error for Newton polynomials
function poly_inter_err(a,b,n,f::Function)
    a , b  = min(a,b), max(a,b)
    h = (b-a)/n
    x = [a+(i-1)*h for i in 1:(n+1)]
    Π = h^(n+1)*factorial(n)/4 
    ξ = abs(maximum(f.(x)))
    err = Π*ξ/factorial(n+1)
    return err
end


err=[poly_inter_err.(0.05,2,i,funs) for i in [4,6,11]] #I dont consider 21 since 21! is too large
Err=[err[1] err[2] err[3]] 
copy_to_clipboard(true)
mdtable(Err, side=label, head=["n=4","n=6","n=11"],fmt = FancyNumberFormatter(4))|> print
mdtable(Err, side=label, head=["n=4","n=6","n=11"],fmt = FancyNumberFormatter(4))|> display

#Plot the accuracy of interpolation as a function of grid size

##Part 2:
function polygrid(n::Integer,a::Float64,b::Float64,θ::Float64)
    grid = collect(range(0,1; length=n))
    xi = a .+ (b-a).*grid.^θ
    return xi
end

function interp_err_SN(F::Array,I::Array)
    error = maximum(abs.(I.-F))/maximum(abs.(F))
end
