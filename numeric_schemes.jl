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

function Entropy_conservative(uinit,dx,dt,N,ntime, ϵ = 0.0, Extra_Viscosity = false)
  uu = copy(uinit)
  uleft = uu[1]; uright = uu[N]
  Kleft = K(uleft); Kright = K(uright)
  #Print progress
  percentage = 0
  limit = ntime/5
  println("Starting entropy conservative scheme")
  for t = 1:ntime
    uold = copy(uu)
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

function estimate_error(reference, xx, uu, dx, N)
  dxr = reference[2,1] - reference[1,1]
  uexact = zeros(N)
  for i = 1:N
    j = Int(round((xx[i] - reference[1,1])/dxr))+1
    uexact[i] = reference[j,2]
  end
  sum(dx*abs(uu - uexact))
end
