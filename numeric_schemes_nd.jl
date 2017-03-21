@enum StepMethod FORWARD_EULER TVD_RK2
function Entropy_nonconservative_nd(uinit,dx,dt,N,Tend, tempSteps = FORWARD_EULER, ϵ = 0.0, Extra_Viscosity = false)
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
      update_uu_NCd(uu, uold, N, dx, dt, ϵ, Extra_Viscosity)
    elseif (tempSteps == TVD_RK2)
      #FIRST Step
      update_uu_NCd(utemp, uold, N, dx, dt, ϵ, Extra_Viscosity)
      #Second Step
      update_uu_NCd(utemp2, utemp, N, dx, dt, ϵ, Extra_Viscosity)
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

function update_uu_NCd(uu, uold, N, dx, dt, ϵ, Extra_Viscosity)
  uleft = uold[1,:]; uright = uold[N,:]
  j = 1
  uu[j,:] = uold[j,:] - dt/dx*(FluxN(uold[j,:], uold[j+1,:])-FluxN(uleft, uold[j,:])) +
  dt/dx^2*(kvisc(uold[j,:],uold[j+1,:])*(uold[j+1,:]-uold[j,:]) - kvisc(uleft,uold[j,:])*(uold[j,:]-uleft)) +
  ϵ*dt/dx^2*(Extra_Viscosity ? uold[j+1,:]-2*uold[j,:]+uleft : 0.0)
  for j = 2:(N-1)
    uu[j,:] = uold[j,:] - dt/dx*(FluxN(uold[j,:], uold[j+1,:])-FluxN(uold[j-1,:], uold[j,:])) +
    dt/dx^2*(kvisc(uold[j,:],uold[j+1,:])*(uold[j+1,:]-uold[j,:]) - kvisc(uold[j-1,:],uold[j,:])*(uold[j,:]-uold[j-1,:]))+
    ϵ*dt/dx^2*(Extra_Viscosity ? uold[j+1,:]-2*uold[j,:]+uold[j-1,:] : 0.0)
  end
  j = N
  uu[j,:] = uold[j,:] - dt/dx*(FluxN(uold[j,:], uright)-FluxN(uold[j-1,:], uold[j,:])) +
  dt/dx^2*(kvisc(uold[j,:],uright)*(uright-uold[j,:]) - kvisc(uold[j-1,:],uold[j,:])*(uold[j,:]-uold[j-1,:]))+
  ϵ*dt/dx^2*(Extra_Viscosity ? uright-2*uold[j,:]+uold[j-1,:]:0.0)
end
