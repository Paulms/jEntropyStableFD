# Test problem 3
# Cauchy Parabolic Problem (Karlsen, Koley & Risebro)
# u_t+(u²/2)_x=K(u)_xx
# [-π/2, π]
# Initial Condition:
# uo =   sin(x)

# Parameters:
const CFL = 0.9
const Tend = 1.0
const M = 8000       #Cells in reference solution
const α = 0.3

#Functions:
K(u) = 0.5*max(u,0.0)^2
Flux(u) = u^2/2
FluxN(ul, ur) = (ur^2 + ul*ur + ul^2)/6.0
cdt(u, CFL, dx) = CFL/(1/dx*maximum(abs(u))+1/dx^2*2*maximum(abs(u)))

#Setup initial Conditions
function setup_initial(N)
  dx = (3/2*π)/(N-1)
  xx = [i*dx-dx-π/2 for i in 1:N]
  uinit = zeros(N)
  uinit = sin(xx)
  return dx,xx, uinit
end

include("numeric_schemes.jl")
## Save reference data
N = M
dx, xx, uinit = setup_initial(N)
uu3 = Entropy_conservative(uinit,dx,CFL,N,Tend)
writedlm("test_3_reference.txt", [xx uu3], '\t')

reference = readdlm("test_3_reference.txt")
steps = [200,400,800,1600,3200]
errors = zeros(2,5)
for (i,step) in enumerate(steps)
  println("Testing with ", step, " steps")
  N = step
  dx, xx, uinit = setup_initial(N)
  uu = Engquist_Osher(uinit,dx,CFL,N,Tend)  #MS
  error = estimate_error(reference[:,2], M, uu, N)
  println("Error: ", error)
  errors[1,i] = error
  uu2 = Entropy_conservative(uinit,dx,CFL,N,Tend) #ESC
  error = estimate_error(reference[:,2], M, uu2, N)
  println("Error: ", error)
  errors[2,i] = error
end

#Compute order
order = log2(errors[:,1:4]./errors[:,2:5])

## Compute errors
using DataFrames
df = DataFrame(errors);
names!(df, map(Symbol,steps));
df[:method] = ["MS","ESC"];
dfo = DataFrame(order);
names!(dfo, map(Symbol,steps[2:5]));
dfo[:method] = ["MS","ESC"];
println(df)
println(dfo)

N=100
dx, xx, uinit = setup_initial(N)
uu = Engquist_Osher(uinit,dx,CFL,N,Tend)
uu2 =  Entropy_conservative(uinit,dx,CFL,N,Tend)
uu3 = Entropy_conservative(uinit,dx,CFL,N,Tend,FORWARD_EULER,α*dx,true) #ESC-0.3
#Plot
using(Plots)
plot(xx, uinit, lab="u0")
plot!(xx, uu, lab="MS")
plot!(xx, uu2,lab="ESC")
plot!(reference[:,1], reference[:,2], lab="REF")
