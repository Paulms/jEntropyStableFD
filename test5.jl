# Test problem 5
# Shallow water system with flat bottom

# Parameters:
const CFL = 0.9
const Tend = 0.2
#const M = 8000       #Cells in reference solution
const μ = 1e-3
const gr = 9.8    #gravedad

#Functions:
kk(u) = μ*sum(vv(u).^2)*[1.0 0.0; 0.0 1.0]
vv(u) = [gr*u[1]-0.5*(u[2]/u[1])^2 (u[2]/u[1])]
FluxN(ul, ur) = [0.25*(ur[1]+ul[1])*(ur[2]/ur[1]+ul[2]/ul[1]); (0.5*gr*0.25*(ur[1]+ul[1])^2+0.25*0.5*(ur[1]+ul[1])*(ur[2]/ur[1]+ul[2]/ul[1])^2)]
function cdt(u, CFL, dx)
  uu = zeros(size(u,1))
  for j = 1:size(u,1)
    uu[j] = sqrt(sum(kk(u[j,:]).^2))
  end
  h = u[:,1]; q = u[:,2]
  fu = sqrt(1+(gr*h-q.^2./h.^2).^2+4*q.^2./h.^2)
  return CFL/(1/dx*maximum(fu)+1/dx^2*2*maximum(uu))
end
kvisc(ul,ur) = μ*(sum(vv(ul).^2 + vv(ur).^2))/2.0*[1.0 0.0;0.0 1.0]

#Setup initial Conditions
function setup_initial(N)
  dx = 10.0/N
  # (-5, 5)
  xx = [i*dx+dx/2-5.0 for i in (0:N-1)]
  uinit = zeros(N, 2)
  for i = 1:N
    if xx[i] < 0.0
      uinit[i,1] = 2.0
    else
     uinit[i,1] = 1.0
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

N=500
dx, xx, uinit = setup_initial(N)
@time uu3 = Entropy_nonconservative_nd(uinit,dx,CFL,N,Tend, TVD_RK2) #ESNC

#Plot
using(Plots)
plot(xx, uinit[:,1], lab="ho",line=(:dot,2))
plot!(xx, uu3[:,1],lab="ESNC h")
plot(xx, uinit[:,2]./uinit[:,1], lab="uo",line=(:dot,2))
plot!(xx, uu3[:,2]./uu3[:,1],lab="ESNC u")
