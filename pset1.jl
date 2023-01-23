##Part a
using Pkg
Pkg.add("Plots")
alpha = 1/3
beta = 0.95
z=1
n=100
#Define vectors to save iteration values for each variable: 
kpath = zeros(Float64, n, 1)
ypath = zeros(Float64, n, 1)
cpath = zeros(Float64, n-1, 1)
rpath = zeros(Float64, n, 1)
wpath = zeros(Float64, n, 1)
#Denote the initial steady state values: 
kpath[1] = (beta*alpha*z)^(1/(1-alpha))
ypath[1] = z*kpath[1]^alpha
rpath[1] = alpha*z*kpath[1]^(alpha-1)
wpath[1] = (1-alpha)*z*kpath[1]^alpha
cpath[1] = z*kpath[1]^alpha - kpath[1]
#Capital shock: 
kpath[2] = 0.8kpath[1]
#Iterations:
for i = 3:100
    kpath[i] = beta*alpha*z*kpath[i-1]^alpha
end

for i = 2:100
    ypath[i] = z*kpath[i]^alpha
end

for i = 2:100
    rpath[i] = alpha*z*kpath[i]^(alpha-1)
end

for i = 2:100
    wpath[i] = (1-alpha)*z*kpath[i]^alpha
end

for i = 2:99
    cpath[i] = z*kpath[i]^alpha - kpath[i+1]
end

#Graphs:
using Plots
gr()
plot(kpath, label = "Capital k")

using Plots
gr()
plot(ypath,color = :green, label = "Production y")

using Plots
gr()
plot(cpath, color =:orange, label ="Consumption c")

using Plots
gr()
plot(wpath, color =:purple, label ="Wage w")

using Plots
gr()
plot(rpath, color = :red, label = "Rental rate r")

##Part b
using Pkg
Pkg.add("Plots")
z2=1.05
#Define vectors to save iteration values for each variable: 
kpath2 = zeros(Float64, n, 1)
ypath2 = zeros(Float64, n, 1)
cpath2 = zeros(Float64, n-1, 1)
rpath2 = zeros(Float64, n, 1)
wpath2 = zeros(Float64, n, 1)
#Denote the initial steady state values: 
kpath2[1] = kpath[1]
ypath2[1] = ypath[1]
rpath2[1] = rpath[1]
wpath2[1] = wpath[1]
cpath2[1] = cpath[1]
#Iterations: 
for i = 2:100 
    kpath2[i] = beta*alpha*z2*kpath[i-1]^alpha
end

for i = 2:100
    ypath2[i] = z2*kpath2[i]^alpha
end

for i = 2:100
    rpath2[i] = alpha*z2*kpath2[i]^(alpha-1)
end

for i = 2:100
    wpath2[i] = (1-alpha)*z2*kpath2[i]^alpha
end

for i = 2:99
    cpath2[i] = z2*kpath2[i]^alpha - kpath[i+1]
end

#Graphs:
using Plots
gr()
plot(kpath2, label = "Capital k")

using Plots
gr()
plot(ypath2,color = :green, label = "Production y")

using Plots
gr()
plot(cpath2, color =:orange, label ="Consumption c")

using Plots
gr()
plot(wpath2, color =:purple, label ="Wage w")

using Plots
gr()
plot(rpath2, color = :red, label = "Rental rate r")


