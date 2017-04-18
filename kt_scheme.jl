using ForwardDiff

@enum StepMethod FORWARD_EULER TVD_RK2 RK4
@enum BoundaryCondition ZERO_FLUX PERIODIC

function KT(uinit,dx,CFL,N,Tend,tempSteps = FORWARD_EULER, boundary = ZERO_FLUX)
  uu = copy(uinit)
  #Print progress
  percentage = 0
  limit = Tend/5
  println("Starting KT scheme")
  t = 0.0
  rhs1 = zeros(N)
  rhs2 = zeros(N)
  rhs3 = zeros(N)
  rhs4 = zeros(N)
  while t <= Tend
    uold = copy(uu)
    dt = cdt(uold, CFL, dx)
    if (tempSteps == FORWARD_EULER)
      update_KT(rhs1, uold, N, dx, boundary)
      uu = uold + dt*rhs1
    elseif (tempSteps == TVD_RK2)
      #FIRST Step
      update_KT(rhs1, uold, N, dx, boundary)
      #Second Step
      update_KT(rhs2, uold+dt*rhs1, N, dx, boundary)
      uu = 0.5*(uold + uold + dt*rhs1 + dt*rhs2)
    elseif (tempSteps == RK4)
      #FIRST Step
      update_KT(rhs1, uold, N, dx, boundary)
      #Second Step
      update_KT(rhs2, uold+dt/2*rhs1, N, dx, boundary)
      #Third Step
      update_KT(rhs3, uold+dt/2*rhs2, N, dx, boundary)
      #Fourth Step
      update_KT(rhs4, uold+dt*rhs3, N, dx, dt, boundary)
      uu = uold + dt/6 * (rhs1+2*rhs2+2*rhs3+rhs4)
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

function KT2(uinit,dx,CFL,N,Tend,tempSteps = FORWARD_EULER, boundary = ZERO_FLUX,Θ=1)
  uu = copy(uinit)
  #Print progress
  percentage = 0
  limit = Tend/5
  println("Starting KT2 scheme")
  t = 0.0
  rhs1 = zeros(N)
  rhs2 = zeros(N)
  rhs3 = zeros(N)
  rhs4 = zeros(N)
  while t <= Tend
    uold = copy(uu)
    dt = cdt(uold, CFL, dx)
    if (tempSteps == FORWARD_EULER)
      update_KT2(rhs1, uold, N, dx, dt, Θ, boundary)
      uu = uold + dt*rhs1
    elseif (tempSteps == TVD_RK2)
      #FIRST Step
      update_KT2(rhs1, uold, N, dx, dt, Θ, boundary)
      #Second Step
      update_KT2(rhs2, uold+dt*rhs1, N, dx, dt, Θ, boundary)
      uu = 0.5*(uold + uold + dt*rhs1 + dt*rhs2)
    elseif (tempSteps == RK4)
      #FIRST Step
      update_KT2(rhs1, uold, N, dx, boundary)
      #Second Step
      update_KT2(rhs2, uold+dt/2*rhs1, N, dx, dt/2, Θ, boundary)
      #Third Step
      update_KT2(rhs3, uold+dt/2*rhs2, N, dx, dt/2, Θ, boundary)
      #Fourth Step
      update_KT2(rhs4, uold+dt*rhs3, N, dx, dt, dt, Θ, boundary)
      uu = uold + dt/6 * (rhs1+2*rhs2+2*rhs3+rhs4)
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
function minmod(a,b)
  0.5*(sign(a)+sign(b))*min(abs(a),abs(b))
end

function update_KT(rhs, uold, N, dx, boundary)
  ∇u = zeros(N)
  for j = 2:N-1
    ∇u[j] = minmod((uold[j] - uold[j-1])/dx,(uold[j+1] - uold[j])/dx)
  end
  uplus = uold[2:N] - dx/2*∇u[2:N]
  uminus = uold[1:N-1] + dx/2*∇u[1:N-1]
  Jf = x -> ForwardDiff.derivative(Flux,x)
  aa = max.(abs(Jf.(uplus)),abs(Jf.(uminus)))

  # Numerical Fluxes
  hh = zeros(N-1)
  hh = 0.5*(Flux.(uplus)+Flux.(uminus))-0.5*(aa.*(uplus-uminus))
  ∇u_ap = (uold[2:N]-uold[1:N-1])/dx
  KK = map(kk, uold)
  pp = 0.5*(KK[1:N-1].*∇u_ap + KK[2:N].*∇u_ap)
  j = 1
  rhs[j] = - 1/dx * (hh[j] - pp[j])
  for j = 2:(N-1)
    rhs[j] = - 1/dx * (hh[j]-hh[j-1]-(pp[j]-pp[j-1]))
  end
  j = N
  rhs[j] = - 1/dx *(-hh[j-1]+pp[j-1])
end


function update_KT2(rhs, uold, N, dx, dt,Θ, boundary)
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
  λ = dt/dx
  #update vector
  # 1. slopes
  ∇u = zeros(N)
  for j = 2:(N-1)
    ∇u[j] = minmod(Θ*(uold[j]-uold[j-1]),(uold[j+1]-uold[j-1])/2,Θ*(uold[j+1]-uold[j]))
  end
  # Local speeds of propagation
  Jf = x -> ForwardDiff.derivative(Flux,x)
  aa=max.(abs(Jf.(uold[1:N-1]+0.5*∇u[1:N-1])),abs(Jf.(uold[2:N]-0.5*∇u[2:N])))
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
    ∇f_l[j] = minmod(Θ*(Flux(u_l[j])-Flux(u_l[j-1])),(Flux(u_l[j+1])-Flux(u_l[j-1]))/2,
    Θ*(Flux(u_l[j+1])-Flux(u_l[j])))
    ∇f_r[j] = minmod(Θ*(Flux(u_r[j])-Flux(u_r[j-1])),(Flux(u_r[j+1])-Flux(u_r[j-1]))/2,
    Θ*(Flux(u_r[j+1])-Flux(u_r[j])))
  end
  # Predictor solution values
  Φ_l = u_l - λ/2*∇f_l
  Φ_r = u_r - λ/2*∇f_r

  # Aproximate cell averages
  Ψr = zeros(N-1)
  Ψ = zeros(N)
  for i = 1:(N-1)
    if (aa[i] != 0)
      Ψr[i] = 0.5*(uold[i]+uold[i+1])+(1-λ*aa[i])/4*(∇u[i]-∇u[i+1])-1/(2*aa[i])*
      (Flux(Φ_r[i])-Flux(Φ_l[i]))
    else
      Ψr[i] = 0.5*(uold[i]+uold[i+1])
    end
  end
  for i = 2:(N-1)
    Ψ[i] = uold[i] - λ/2*(aa[i]-aa[i-1])*∇u[i]-λ/(1-λ*(aa[i]+aa[i-1]))*
    (Flux(Φ_l[i])-Flux(Φ_r[i-1]))
  end

  # Discrete derivatives
  ∇Ψ = zeros(N-1)
  for i = 2:(N-2)
    ∇Ψ=2\dx*minmod(Θ*(Φ_r[i]-Ψ[i])/(1+λ*(aa[i]-aa[i-1])),
    (Ψ[i+1]-Ψ[i])/(2+λ*(2*aa[i]-aa[i-1]-aa[i+1])),
    Θ*(Ψ[i+1]-Φ_r[i])/(1+λ*(aa[i]-aa[i+1])))
  end

  # Numerical Fluxes
  hh = zeros(N-1)
  hh = 0.5*(Flux.(Φ_r)+Flux.(Φ_l))-0.5*(uold[2:N]-uold[1:N-1]).*aa+
  aa.*(1-λ*aa)/4.*(∇u[2:N]+∇u[1:N-1]) + λ*dx/2*(aa).^2.*∇Ψ
  ∇u_ap = (uold[2:N]-uold[1:N-1])/dx
  pp = 0.5*(KK[1:N-1].*∇u_ap + KK[2:N].*∇u_ap)
  j = 1
  rhs[j] = - 1/dx * (hh[j] - pp[j])
  for j = 2:(N-1)
    rhs[j] = - 1/dx * (hh[j]-hh[j-1]-(pp[j]-pp[j-1]))
  end
  j = N
  rhs[j] =  -1/dx*(-hh[j-1]+pp[j-1])
end
