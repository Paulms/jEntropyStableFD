# Test problem 4
# Parabolic System

# Parameters:
const CFL = 0.9
const Tend = 1.0
#const M = 8000       #Cells in reference solution
const μ = 0.1

#Functions:
Flux(u) = u.^2/2
kk(u) = μ*sum(u.^2)*[1.0 0.0; 0.0 1.0]
FluxN(ul, ur) = (ur.^2 + ul.*ur + ul.^2)/6.0
function cdt(u, CFL, dx)
  uu = zeros(size(u,1))
  for j = 1:size(u,1)
    uu[j] = sqrt(sum(kk(u[j,:]).^2))
  end
  return CFL/(1/dx*maximum(sqrt(u[:,1].^2 + u[:,2].^2))+1/dx^2*2*maximum(uu))
end
kvisc(ul,ur) = μ*(sum(ul.^2 + ur.^2))/2.0*I#[1.0 0.0;0.0 1.0]

#Setup initial Conditions
function setup_initial(N)
  dx = 5.0/N
  # (-2.5, 2.5)
  xx = [i*dx+dx/2-2.5 for i in (0:N-1)]
  uinit = zeros(N, 2)
  for i = 1:N
    if -1.5 < xx[i] < -1.3
      uinit[i,1] = 1.0
      uinit[i,2] = 0.0
    elseif -0.5 < xx[i] < 0.5
      uinit[i,1] = 1.0
      uinit[i,2] = 1.0
    elseif 1.3 < xx[i] < 1.5
      uinit[i,1] = 0.0
      uinit[i,2] = 1.0
    else
     uinit[i,1] = 0.0
     uinit[i,2] = 0.0
   end
  end
  return dx,xx, uinit
end

include("numeric_schemes_nd.jl")
# Save reference data
# N = M
# dx, xx, uinit = setup_initial(N)
# uu3 = Entropy_conservative(uinit,dx,CFL,N,Tend)
# writedlm("test_4_reference.txt", [xx uu3], '\t')

# reference = readdlm("test_4_reference.txt")
# steps = [200,400,800,1600,3200]
# errors = zeros(2,5)
# for (i,step) in enumerate(steps)
#   println("Testing with ", step, " steps")
#   N = step
#   dx, xx, uinit = setup_initial(N)
#   uu = Engquist_Osher(uinit,dx,CFL,N,Tend)  #MS
#   error = estimate_error(reference[:,2], M, uu, N)
#   println("Error: ", error)
#   errors[1,i] = error
#   uu2 = Entropy_conservative(uinit,dx,CFL,N,Tend) #ESC
#   error = estimate_error(reference[:,2], M, uu2, N)
#   println("Error: ", error)
#   errors[2,i] = error
# end
#
# #Compute order
# order = log2(errors[:,1:4]./errors[:,2:5])
#
# ## Compute errors
# using DataFrames
# df = DataFrame(errors);
# names!(df, map(Symbol,steps));
# df[:method] = ["MS","ESC"];
# dfo = DataFrame(order);
# names!(dfo, map(Symbol,steps[2:5]));
# dfo[:method] = ["MS","ESC"];
# println(df)
# println(dfo)

N=1000
dx, xx, uinit = setup_initial(N)
@time uu3 = Entropy_nonconservative_nd(uinit,dx,CFL,N,Tend) #ESNC

#Plot
using(Plots)
plot(xx, uinit[:,1], lab="u1o",line=(:dot,2))
plot!(xx, uu3[:,1],lab="ESNC u1")
plot(xx, uinit[:,2], lab="u2o",line=(:dot,2))
plot!(xx, uu3[:,2],lab="ESNC u2")
