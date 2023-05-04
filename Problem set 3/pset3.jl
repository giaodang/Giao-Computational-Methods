using Plots
import Pkg
Pkg.add("Interpolations")
using Interpolations
Pkg.add("ForwardDiff")
using ForwardDiff
Pkg.add("LaTeXStrings")
using LaTeXStrings
Pkg.add("Latexify")
using Latexify
Pkg.add("PyPlot")
using PyPlot
Pkg.add("LinearAlgebra")
using LinearAlgebra

##Part 1
#Define utility functions

u1(x) = log(x)
u2(x) = (x)^(1/2)
u3(x,σ) = x^(1-σ)/(1-σ)
u4(x) = u3(x,2)
u5(x) = u3(x,5)
u6(x) = u3(x,10)

#Define corresponding derivatives
u1_d(x) = 1/x
u2_d(x) = (1/2)*x^(-1/2)
u3_d(x,σ) = x^(-σ)
u4_d(x) = x^(-2)
u5_d(x) = x^(-5)
u6_d(x) = x^(-10)

grid = range(0.05,2;length=1000)


#Newton polynomials

    function diff(x::Array,y::Array)
        m = length(x)
        a = [y[i] for i in 1:m]

        for j in 2:m
            for i in reverse(collect(j:m))
                a[i]=(a[i]-a[i-1])/(x[i]-x[i-(j-1)])
            end
        end
        return a
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

nom = ["Log","Square",L"CRRA σ=2", L"CRRA σ=5", L"CRRA σ=10"]
funs = [u1,u2,u4,u5,u6]
for i in 1:5
    name = nom[i]
    for n = (4,6,11,21)
        fn = funs[i].(aaa)

        xi = collect(range(0.05,2;length=n)) ;
        yi = funs[i].(xi)

        interp=map(z->newton(xi,yi,z),aaa)

        gr()
        plot(title="Interpolation $name n=$n - Newton Polynomial")
        plot!(aaa,fn,linewidth=3,label = "Function: $name",legend=(0.75,0.75),foreground_color_legend = nothing,background_color_legend = nothing)
        plot!(aaa,interp,linewidth=3,label="Interpolation")
        plot!(xi,yi,linetype=:scatter,marker=(:diamond,9),markercolor=RGB(0.5,0.1,0.1),label = "Data")
        savefig("graphs/Newton $name _n_$n")

    end
end
                



                        









