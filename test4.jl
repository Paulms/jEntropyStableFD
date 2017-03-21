# Test problem 4
# Parabolic System

# Parameters:
const CFL = 0.9
const Tend = 2.0
#const M = 8000       #Cells in reference solution
const μ = 0.1

#Functions:
Flux(u) = u.^2/2
kk(u) = μ*[norm(u)^2 0.0; 0.0 norm(u)^2]
FluxN(ul, ur) = (ur.^2 + ul.*ur + ul.^2)/6.0
cdt(u, CFL, dx) = CFL/(1/dx*maximum(norm(u))+1/dx^2*2*maximum(norm(kk(u))))
kvisc(ul,ur) = μ*(norm(ul)^2 + norm(ur)^2)/2.0*eye(2)

#Setup initial Conditions
function setup_initial(N)
  dx = 5/N
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
uu3 = Entropy_nonconservative_nd(uinit,dx,CFL,N,Tend) #ESNC
#Plot
using(Plots)
plot(xx, uinit, lab="u0",line=(:dot,2))
plot!(xx, uu3[:,1],lab="ESNC")
