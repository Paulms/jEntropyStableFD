# Test problem 1
# Burger's Equation
# u_t+(u²/2)_x=K(u)_xx
# Initial Condition:
# uo =   (1-x²)²   if -1<x<1
#         0         otherwise
# Parameters:
const μ = 0.01
const CFL = 0.9
const Tend = 0.5
const α = 0.1
const M = 16000       #Cells in reference solution

#Functions:
K(u) = μ*u^2
Flux(u) = u^2/2
FluxN(ul, ur) = (ur^2 + ul*ur + ul^2)/6.0
function kvisc(ul, ur)
  if (abs(ul)< eps() && abs(ur) <eps())
    0.0
  else
    μ*4/3*(ul^2+ul*ur+ur^2)/(ul+ur)
  end
end

#Setup initial Conditions
function setup_initial(N)
  dx = 4.0/(N-1)
  dt = CFL/(1.0*(1/dx+4*μ/dx^2))
  ntime = Int(floor(Tend/dt))
  xx = [i*dx-dx-2 for i in 1:N]
  uinit = zeros(N)
  for (i,x) in enumerate(xx)
    if (x >=-1 && x<=1)
      uinit[i] = (1-x^2)^2
    end
  end
  return dx, dt, ntime, xx, uinit
end

include("numeric_schemes.jl")
#Save reference data
# N = M
# dx, dt, ntime, xx, uinit = setup_initial(N)
# uu3 = Entropy_conservative(uinit,dx,dt,N,ntime)
# writedlm("burger_1_reference.txt", [xx uu3], '\t')

reference = readdlm("burger_1_reference.txt")
steps = [200,400,800,1600,3200]
errors = zeros(5,5)
for (i,step) in enumerate(steps)
  println("Testing with ", step, " steps")
  N = step
  dx, dt, ntime, xx, uinit = setup_initial(N)
  uu = Engquist_Osher(uinit,dx,dt,N,ntime)  #MS
  error = estimate_error(reference[:,2], M, uu, N)
  println("Error: ", error)
  errors[1,i] = error
  uu2 = Entropy_conservative(uinit,dx,dt,N,ntime) #ESC
  error = estimate_error(reference[:,2], M, uu2, N)
  println("Error: ", error)
  errors[2,i] = error
  uu3 = Entropy_nonconservative(uinit,dx,dt,N,ntime) #ESNC
  error = estimate_error(reference[:,2], M, uu3, N)
  println("Error: ", error)
  errors[3,i] = error
  uu4 = Entropy_conservative(uinit,dx,dt,N,ntime, TVD_RK2) #ESC2
  error = estimate_error(reference[:,2], M, uu4, N)
  println("Error: ", error)
  errors[4,i] = error
  uu5 = Entropy_nonconservative(uinit,dx,dt,N,ntime, TVD_RK2) #ESNC2
  error = estimate_error(reference[:,2], M, uu5, N)
  println("Error: ", error)
  errors[5,i] = error
end

#Compute order
order = log2(errors[:,1:4]./errors[:,2:5])

## Display Errors and Order:
# using DataFrames
# df = DataFrame(errors);
# names!(df, map(Symbol,steps));
# df[:method] = ["MS","ESC","ESNC","ESC2","ESNC2"];
# dfo = DataFrame(order);
# names!(dfo, map(Symbol,steps[2:5]));
# dfo[:method] = ["MS","ESC","ESNC","ESC2","ESNC2"];
# println(df)
# println(dfo)

N=400
dx, dt, ntime, xx, uinit = setup_initial(N)
uu = Engquist_Osher(uinit,dx,dt,N,ntime)
uu2 =  Entropy_conservative(uinit,dx,dt,N,ntime, TVD_RK2)
#Plot
using(Plots)
plot(xx, uinit, lab="u0")
plot!(xx, uu, lab="MS")
plot!(xx, uu2,lab="ESCN2")
plot!(reference[:,1], reference[:,2], lab="REF")
