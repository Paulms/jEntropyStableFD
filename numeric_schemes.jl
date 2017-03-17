@enum StepMethod FORWARD_EULER TVD_RK2
function Engquist_Osher(uinit,dx,CFL,N,Tend)
  uu = copy(uinit)
  #Print progress
  percentage = 0
  limit = Tend/5
  println("Starting monotone scheme")
  t = 0.0
  while t <= Tend
    uold = copy(uu)
    dt = cdt(uold, CFL, dx)
    #Compute fluxes
    uleft = uold[1]; uright = uold[N]
    fplusleft = (uleft > 0) ? Flux(uleft) : 0.0
    fminusright = (uright > 0) ? 0.0 : Flux(uright)
    Kleft = K(uleft); Kright = K(uright)
    
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
      limit = limit + Tend/5
      println(percentage, "% completed")
    end
    t = t + dt
  end
  println("Completed...")
  return uu
end

function Entropy_conservative(uinit,dx,dt,N,Tend, tempSteps = FORWARD_EULER, ϵ = 0.0, Extra_Viscosity = false)
  uu = copy(uinit)
  #Print progress
  percentage = 0
  limit = Tend/5
  println("Starting entropy conservative scheme")
  t = 0.0
  utemp = zeros(N)
  utemp2 = zeros(N)
  while t <= Tend
    uold = copy(uu)
    dt = cdt(uold, CFL, dx)
    if (tempSteps == FORWARD_EULER)
      update_uu_EC(uu, uold, N, dx, dt, ϵ, Extra_Viscosity)
    elseif (tempSteps == TVD_RK2)
      #FIRST Step
      update_uu_EC(utemp, uold, N, dx, dt, ϵ, Extra_Viscosity)
      #Second Step
      update_uu_EC(utemp2, utemp, N, dx, dt, ϵ, Extra_Viscosity)
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

function update_uu_EC(uu, uold, N, dx, dt, ϵ, Extra_Viscosity)
  #Compute diffusion
  uleft = uold[1]; uright = uold[N]
  KK = map(K, uold)
  Kleft = K(uleft); Kright = K(uright)
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

function Entropy_nonconservative(uinit,dx,dt,N,Tend, tempSteps = FORWARD_EULER, ϵ = 0.0, Extra_Viscosity = false)
  uu = copy(uinit)
  #Print progress
  percentage = 0
  limit = Tend/5
  println("Starting entropy non-conservative scheme")
  t = 0.0
  utemp = zeros(N)
  utemp2 = zeros(N)
  while  t<= Tend
    uold = copy(uu)
    dt = cdt(uold, CFL, dx)
    #update vector
    if (tempSteps == FORWARD_EULER)
      update_uu_NC(uu, uold, N, dx, dt, ϵ, Extra_Viscosity)
    elseif (tempSteps == TVD_RK2)
      #FIRST Step
      update_uu_NC(utemp, uold, N, dx, dt, ϵ, Extra_Viscosity)
      #Second Step
      update_uu_NC(utemp2, utemp, N, dx, dt, ϵ, Extra_Viscosity)
      uu = 0.5*(uold + utemp2)
    end
    # Print Progress
    if (t > limit)
      percentage = percentage + 20
      limit = limit +Tend/5
      println(percentage, "% completed")
    end
    t = t + dt
  end
  println("Completed...")
  return uu
end

function update_uu_NC(uu, uold, N, dx, dt, ϵ, Extra_Viscosity)
  uleft = uold[1]; uright = uold[N]
  Kleft = K(uleft); Kright = K(uright)
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
