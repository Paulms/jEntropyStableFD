# Test problem 1
#Burger's Equation
#Parameters
const μ = 0.01
const CFL = 0.9
const Tend = 0.5

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
  dx = 4.0/(N+1)
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
N = 16000
dx, dt, ntime, xx, uinit = setup_initial(N)
uu3 = Entropy_conservative(uinit,dx,dt,N,ntime)
writedlm("burger_1_reference.txt", [xx uu3], '\t')

reference = readdlm("burger_1_reference.txt")
steps = [200,400,800,1600,3200]
errors = zeros(2,5)
for (i,step) in enumerate(steps)
  println("Testing with ", step, " steps")
  N = step
  dx, dt, ntime, xx, uinit = setup_initial(N)
  uu = Engquist_Osher(uinit,dx,dt,N,ntime)
  error = estimate_error(reference, xx, uu, dx, N)
  println("Error: ", error)
  errors[1,i] = error
  uu2 = Entropy_conservative(uinit,dx,dt,N,ntime)
  error = estimate_error(reference, xx, uu2, dx, N)
  println("Error: ", error)
  errors[2,i] = error
end

#Compute order
order = log2(errors[:,1:4]./errors[:,2:5])

N=400
dx, dt, ntime, xx, uinit = setup_initial(N)
uu = Engquist_Osher(uinit,dx,dt,N,ntime)
uu2 = Entropy_conservative(uinit,dx,dt,N,ntime)
#Plot
using(Plots)
plot(xx, uinit)
plot!(xx, uu)
plot!(xx, uu2)
plot!(reference[:,1], reference[:,2])
