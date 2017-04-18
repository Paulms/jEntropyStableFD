#Setup initial Conditions
function setup_initial(N)
  # We use ghost cells
  dx=2*π/N
  xx = [i*dx+dx/2 for i in 0:(N-1)]
  uinit = zeros(N)
  uinit = exp.(-4*(xx-π*1/2).^2)-exp.(-4*(xx-π*3/2).^2)
  return dx, xx, uinit
end
CFL = 0.25
N = 256
kk(u) = 0
cdt(u, CFL, dx) = CFL/(1/dx*maximum(abs(u)))
Flux(u) = u^2/2
Tend = 3
include("kt_scheme.jl")
dx, xx, uinit = setup_initial(N)
@time uu3 =  KT(uinit,dx,CFL/4,N,Tend, TVD_RK2)
@time uu4 =  KT2(uinit,dx,CFL/4,N,Tend, TVD_RK2)

using(Plots)
plot(xx, uinit, lab="u0",line=(:dot,2))
plot!(xx, uu3, lab="KT")
plot!(xx, uu4, lab="KT2")
