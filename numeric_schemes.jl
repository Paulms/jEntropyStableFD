@enum StepMethod FORWARD_EULER TVD_RK2
function Engquist_Osher(uinit,dx,dt,N,ntime)
  uu = copy(uinit)
  uleft = uu[1]; uright = uu[N]
  fplusleft = Flux(uleft); fminusright = Flux(uright)
  Kleft = K(uleft); Kright = K(uright)
  #Print progress
  percentage = 0
  limit = ntime/5
  println("Starting monotone scheme")
  for t = 1:ntime
    uold = copy(uu)
    #Compute fluxes
    fplus = zeros(N); fminus = zeros(N)
    for j = 1:N
      if (uold[j] > 0.0)
        fplus[j] = Flux(uold[j])
      else
        fminus[j] = Flux(uold[j])
      end
    end
    KK = map(K, uold)
    #update vector
    j = 1
    uu[j] = uold[j] - dt/dx * (fplus[j] + fminus[j+1] - fplusleft-fminus[j]) +
    dt/dx^2*(KK[j+1] - 2*KK[j] + Kleft)
    for j = 2:(N-1)
      uu[j] = uold[j] - dt/dx * (fplus[j] + fminus[j+1] - fplus[j-1]-fminus[j]) +
      dt/dx^2*(KK[j+1] - 2*KK[j] + KK[j-1])
    end
    j = N
    uu[j] = uold[j] - dt/dx * (fplus[j] + fminusright - fplus[j-1]-fminus[j]) +
    dt/dx^2*(Kright - 2*KK[j] + KK[j-1])
    # Print Progress
    if (t > limit)
      percentage = percentage + 20
      limit = limit + ntime/5
      println(percentage, "% completed")
    end
  end
  println("Completed...")
  return uu
end

function Entropy_conservative(uinit,dx,dt,N,ntime, tempSteps = FORWARD_EULER, ϵ = 0.0, Extra_Viscosity = false)
  uu = copy(uinit)
  uleft = uu[1]; uright = uu[N]
  Kleft = K(uleft); Kright = K(uright)
  #Print progress
  percentage = 0
  limit = ntime/5
  println("Starting entropy conservative scheme")
  utemp = zeros(N)
  utemp2 = zeros(N)
  for t = 1:ntime
    uold = copy(uu)
    if (tempSteps == FORWARD_EULER)
      update_uu_EC(uu, uold, N, dx, dt, uleft, uright, Kleft, Kright, ϵ, Extra_Viscosity)
    elseif (tempSteps == TVD_RK2)
      #FIRST Step
      update_uu_EC(utemp, uold, N, dx, dt, uleft, uright, Kleft, Kright, ϵ, Extra_Viscosity)
      #Second Step
      update_uu_EC(utemp2, utemp, N, dx, dt, uleft, uright, Kleft, Kright, ϵ, Extra_Viscosity)
      uu = 0.5*(uold + utemp2)
    end
    # Print Progress
    if (t > limit)
      percentage = percentage + 20
      limit = limit + ntime/5
      println(percentage, "% completed")
    end
  end
  println("Completed...")
  return uu
end

function update_uu_EC(uu, uold, N, dx, dt, uleft, uright, Kleft, Kright, ϵ, Extra_Viscosity)
  #Compute diffusion
  KK = map(K, uold)
  #update vector
  j = 1
  uu[j] = uold[j] - dt/dx * (FluxN(uold[j], uold[j+1])-FluxN(uleft, uold[j])) +
  dt/dx^2*(KK[j+1] - 2*KK[j] + Kleft) +
  ϵ*dt/dx^2*(Extra_Viscosity ? uold[j+1]-2*uold[j]+uleft : 0.0)

  for j = 2:(N-1)
    uu[j] = uold[j] - dt/dx * (FluxN(uold[j], uold[j+1])-FluxN(uold[j-1], uold[j])) +
    dt/dx^2*(KK[j+1] - 2*KK[j] + KK[j-1])+
    ϵ*dt/dx^2*(Extra_Viscosity ? uold[j+1]-2*uold[j]+uold[j-1] : 0.0)
  end
  j = N
  uu[j] = uold[j] - dt/dx * (FluxN(uold[j], uright)-FluxN(uold[j-1], uold[j])) +
  dt/dx^2*(Kright - 2*KK[j] + KK[j-1])+
  ϵ*dt/dx^2*(Extra_Viscosity ? uright-2*uold[j]+uold[j-1]:0.0)
end

function Entropy_nonconservative(uinit,dx,dt,N,ntime, tempSteps = FORWARD_EULER, ϵ = 0.0, Extra_Viscosity = false)
  uu = copy(uinit)
  uleft = uu[1]; uright = uu[N]
  Kleft = K(uleft); Kright = K(uright)
  #Print progress
  percentage = 0
  limit = ntime/5
  println("Starting entropy non-conservative scheme")
  utemp = zeros(N)
  utemp2 = zeros(N)
  for t = 1:ntime
    uold = copy(uu)
    #update vector
    if (tempSteps == FORWARD_EULER)
      update_uu_NC(uu, uold, N, dx, dt, uleft, uright, ϵ, Extra_Viscosity)
    elseif (tempSteps == TVD_RK2)
      #FIRST Step
      update_uu_NC(utemp, uold, N, dx, dt, uleft, uright,ϵ, Extra_Viscosity)
      #Second Step
      update_uu_NC(utemp2, utemp, N, dx, dt, uleft, uright, ϵ, Extra_Viscosity)
      uu = 0.5*(uold + utemp2)
    end
    # Print Progress
    if (t > limit)
      percentage = percentage + 20
      limit = limit + ntime/5
      println(percentage, "% completed")
    end
  end
  println("Completed...")
  return uu
end

function update_uu_NC(uu, uold, N, dx, dt, uleft, uright, ϵ, Extra_Viscosity)
  j = 1
  uu[j] = uold[j] - dt/dx * (FluxN(uold[j], uold[j+1])-FluxN(uleft, uold[j])) +
  dt/dx^2*(kvisc(uold[j],uold[j+1])*(uold[j+1]-uold[j]) - kvisc(uleft,uold[j])*(uold[j]-uleft)) +
  ϵ*dt/dx^2*(Extra_Viscosity ? uold[j+1]-2*uold[j]+uleft : 0.0)
  for j = 2:(N-1)
    uu[j] = uold[j] - dt/dx * (FluxN(uold[j], uold[j+1])-FluxN(uold[j-1], uold[j])) +
    dt/dx^2*(kvisc(uold[j],uold[j+1])*(uold[j+1]-uold[j]) - kvisc(uold[j-1],uold[j])*(uold[j]-uold[j-1]))+
    ϵ*dt/dx^2*(Extra_Viscosity ? uold[j+1]-2*uold[j]+uold[j-1] : 0.0)
  end
  j = N
  uu[j] = uold[j] - dt/dx * (FluxN(uold[j], uright)-FluxN(uold[j-1], uold[j])) +
  dt/dx^2*(kvisc(uold[j],uright)*(uright-uold[j]) - kvisc(uold[j-1],uold[j])*(uold[j]-uold[j-1]))+
  ϵ*dt/dx^2*(Extra_Viscosity ? uright-2*uold[j]+uold[j-1]:0.0)
end

function estimate_error(reference,M, uu,N)
  uexact = zeros(N)
  R = Int(round(M/N))
  for i = 1:N
      uexact[i] = 1.0/R*sum(reference[R*(i-1)+1:R*i])
  end
  sum(1.0/N*abs(uu - uexact))
end
