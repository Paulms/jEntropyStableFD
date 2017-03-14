# Test problem 2
# u_t+(u²)_x=(k(u)u_x)_x

# Parameters:
const CFL = 0.9
const Tend = 0.15
const α = 0.2
const M = 8000       #Cells in reference solution
global ϵ = 0.0

#Functions:
function K(u)
  if (u <= 0.5)
    0.0
  elseif (0.5 < u <0.6)
    1.25*u^2-1.25*u+5/16
  else
    0.25*u-11/80
  end
end
Flux(u) = u^2
FluxN(ul, ur) = (ur^2 + ul*ur + ul^2)/3.0
function r(u)
  if (u <= 0.5)
    0.0
  elseif (0.5<u<0.6)
    5.0/6*u^3-5.0/8*u^2+5/96
  else
    u^2/8-91/2400
  end
end

function kvisc(ul, ur)
  if (ul <= 0.5 && ur <= 0.5)
    0.0
  elseif (0.5 < ul < 0.6 && 0.5 < ur < 0.6)
    2.0/(ul+ur)*(5.0/6*(ul^2+ul*ur+ur^2)-1/8*(ul+ur))
  elseif (ul >= 0.6 && ur >= 0.6)
    2.0/(ul+ur)*(1.0/8*(ul+ur))
  else
    2/(ul+ur)*(r(ur)-r(ul))/(ur-ul)
  end
end

# function kvisc(ul, ur)
#   if (abs(ul-0.5) < ϵ && abs(ur-0.5) < ϵ && (ul-0.5)*(ur-0.5) < 0)
#     0.0
#   elseif (abs(ul-0.6) < ϵ && abs(ur-0.6) < ϵ && (ul-0.6)*(ur-0.6) < 0)
#     0.25/0.6
#   end
# end


#Setup initial Conditions
function setup_initial(N)
  dx = 2.0/(N-1)
  dt = CFL/(2.0/dx+0.5/dx^2)
  global ϵ
  ϵ = α*dx
  ntime = Int(floor(Tend/dt))
  xx = [i*dx-dx-1 for i in 1:N]
  uinit = zeros(N)
  for (i,x) in enumerate(xx)
    if (x <=-0.5)
      uinit[i] = 0.0
    elseif (-0.5<x<-0.3)
      uinit[i] = 5*(x+0.5)
    elseif (-0.3<x<0.3)
      uinit[i] = 1.0
    elseif (0.3<x<0.5)
      uinit[i] = 5*(0.5-x)
    else
      uinit[i] = 0.0
    end
  end
  return dx, dt, ntime, xx, uinit
end

include("numeric_schemes.jl")
#Save reference data
# N = M
# dx, dt, ntime, xx, uinit = setup_initial(N)
# @time uu3 = Entropy_conservative(uinit,dx,dt,N,ntime,FORWARD_EULER,α*dx,true)
# writedlm("test_2_reference.txt", [xx uu3], '\t')

reference = readdlm("test_2_reference.txt")
steps = [200,400,800,1600,3200]
errors = zeros(3,5)
for (i,step) in enumerate(steps)
  println("Testing with ", step, " steps")
  N = step
  dx, dt, ntime, xx, uinit = setup_initial(N)
  uu = Engquist_Osher(uinit,dx,dt,N,ntime)  #MS
  error = estimate_error(reference[:,2], M, uu, N)
  println("Error: ", error)
  errors[1,i] = error
  uu2 = Entropy_conservative(uinit,dx,dt,N,ntime,FORWARD_EULER,α*dx,true) #ESC-0.2
  error = estimate_error(reference[:,2], M, uu2, N)
  println("Error: ", error)
  errors[2,i] = error
  uu3 = Entropy_nonconservative(uinit,dx,dt,N,ntime,FORWARD_EULER,α*dx,true) #ESNC-0.2
  error = estimate_error(reference[:,2], M, uu3, N)
  println("Error: ", error)
  errors[3,i] = error
end

#Compute order
order = log2(errors[:,1:4]./errors[:,2:5])

## Display Errors and Order:
# using DataFrames
# df = DataFrame(errors);
# names!(df, map(Symbol,steps));
# df[:method] = ["MS","ESC-0.2","ESNC-02"];
# dfo = DataFrame(order);
# names!(dfo, map(Symbol,steps[2:5]));
# dfo[:method] = ["MS","ESC-0.2","ESNC-0.2"];
# println(df)
# println(dfo)

N=400
dx, dt, ntime, xx, uinit = setup_initial(N)
uu = Engquist_Osher(uinit,dx,dt,N,ntime)
uu2 = Entropy_conservative(uinit,dx,dt,N,ntime,FORWARD_EULER,α*dx,true) #ESC-0.2
uu3 = Entropy_nonconservative(uinit,dx,dt,N,ntime,FORWARD_EULER,α*dx,true)
#writedlm("test_2_400.txt", [xx uu uu2 uu3], '\t')
#Plot
using(Plots)
plot(xx, uinit, lab="u0")
plot!(xx, uu, lab="MS")
plot!(xx, uu2,lab="ESC-0.2")
plot!(xx, uu3,lab="ESNC-0.2")
plot!(reference[:,1], reference[:,2], lab="REF")
