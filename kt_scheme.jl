using ForwardDiff

@enum StepMethod FORWARD_EULER TVD_RK2
@enum BoundaryCondition ZERO_FLUX PERIODIC

function KT(uinit,dx,dt,N,Tend,tempSteps = FORWARD_EULER,
  Θ = 1.0, boundary = ZERO_FLUX)
  uu = copy(uinit)
  #Print progress
  percentage = 0
  limit = Tend/5
  println("Starting KT scheme")
  t = 0.0
  utemp = zeros(N)
  utemp2 = zeros(N)
  while t <= Tend
    uold = copy(uu)
    dt = cdt(uold, CFL, dx)
    if (tempSteps == FORWARD_EULER)
      update_KT(uu, uold, N, dx, dt, Θ, boundary)
    elseif (tempSteps == TVD_RK2)
      #FIRST Step
      update_KT(utemp, uold, N, dx, dt, Θ, boundary)
      #Second Step
      update_KT(utemp2, utemp, N, dx, dt, Θ, boundary)
      uu = 0.5*(uold + utemp2)
    end
    # Print Progress
    if (t > limit)
      percentage = percentage + 20
      limit = limit + Tend/5
      println(percentage, "% completed")
    end
    t = t+dt
  end
  println("Completed...")
  return uu
end

function minmod(a,b,c)
  if (a > 0 && b > 0 && c > 0)
    min(a,b,c)
  elseif (a < 0 && b < 0 && c < 0)
    max(a,b,c)
  else
    zero(a)
  end
end

function update_KT(uu, uold, N, dx, dt, Θ, boundary)
  #Compute diffusion
  # TODO: Periodic condition is not working
  if (boundary == ZERO_FLUX)
    uleft = uold[1]; uright = uold[N]
    flag = 0.0
  else
    uleft = uold[N]; uright = uold[1]
    flag = 1.0
  end
  KK = map(kk, uold)
  Kleft = kk(uleft); Kright = kk(uright)
  λ = dt/dx
  #update vector
  # 1. slopes
  ∇u = zeros(N)
  for j = 2:(N-1)
    ∇u[j] = minmod(Θ*(uold[j]-uold[j-1]),(uold[j+1]-uold[j-1])/2,Θ*(uold[j+1]-uold[j]))
  end
  # Local speeds of propagation
  Jf = x -> ForwardDiff.derivative(Flux,x)
  aa=max(abs(Jf.(uold[1:N-1]+0.5*∇u[1:N-1])),abs(Jf.(uold[2:N]-0.5*∇u[2:N])))
  println(aa)
  #Flux slopes
  u_l = zeros(N-1)
  u_r = zeros(N-1)
  for j = 2:N
    u_l[j-1] = uold[j-1] + (0.5-λ*aa[j-1])*∇u[j-1]
    u_r[j-1] = uold[j] - (0.5-λ*aa[j-1])*∇u[j]
  end
  ∇f_l = zeros(N-1)
  ∇f_r = zeros(N-1)
  for j = 2:(N-2)
    ∇f_l = minmod(Θ*(Flux(u_l[j])-Flux(u_l[j-1])),(Flux(u_l[j+1])-Flux(u_l[j-1])),
    Θ*(Flux(u_l[j+1])-Flux(u_l[j])))
    ∇f_r = minmod(Θ*(Flux(u_r[j])-Flux(u_r[j-1])),(Flux(u_r[j+1])-Flux(u_r[j-1])),
    Θ*(Flux(u_r[j+1])-Flux(u_r[j])))
  end
  # Predictor solution values
  Ψ_l = u_l - λ/2*∇f_l
  Ψ_r = u_r - λ/2*∇f_r

  # Aproximate cell averages
  Ψ = zeros(N)
  Ψr = zeros(N)
  for i = 2:(N-1)
    Ψ[i] = uold[i] - λ/2*(aa[i]-aa[i-1])*∇u[i]-λ/(1-λ*(aa[i]+aa[i-1]))*
    (Flux(Ψ_l[i])-Flux(Ψ_r[i-1]))
  end
  for i = 1:(N-1)
    Ψr[i] = 0.5*(uold[i]+uold[i+1])+(1-λ*aa[i])/4*(∇u[i]-∇u[i+1])-1/(2*aa[i])*
    (Flux(Ψ_r[i])-Flux(Ψ_l[i]))
  end

  # Discrete derivatives
  ∇Ψ = zeros(N-1)
  for i = 2:(N-2)
    ∇Ψ=2\dx*minmod(Θ*(Ψ_r[i]-Ψ[i])/(1+λ*(aa[i]-aa[i-1])),
    (Ψ[i+1]-Ψ[i])/(2+λ*(aa[i]-aa[i-1]-aa[i+1])),
    Θ*(Ψ[i+1]-Ψ_r[i])/(1+λ*(aa[i]-aa[i+1])))
  end

  # Numerical Fluxes
  hh = zeros(N-1)
  hh = 0.5*(Flux.(Ψ_r)+Flux.(Ψ_l))-0.5*(uold[2:N]-uold[1:N-1])./aa+
  aa.*(1-λ*aa)/4.*(∇u[2:N]-∇u[1:N-1]) + λ*dx/2*(aa).^2.*∇Ψ
  ∇u_ap = (∇u[2:N]-∇u[1:N-1])/dx
  pp = 0.5*(KK[1:N-1].*∇u_ap + KK[2:N].*∇u_ap)

  j = 1
  uu[j] = uold[j] - dt/dx * (hh[j] - pp[j])

  for j = 2:(N-1)
    uu[j] = uold[j] - dt/dx * (hh[j]-hh[j-1]-(pp[j]-pp[j-1]))
  end
  j = N
  uu[j] = uold[j] - dt/dx*(-hh[j-1]+pp[j-1])
end
